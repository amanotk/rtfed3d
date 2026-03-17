// -*- C++ -*-

///
/// @brief Orszag-Tang Vortex Problem
///
/// $Id: vortex.cpp,v 267ff9a5277d 2016/03/16 09:14:48 amano $
///
#include "global.hpp"
#include "solver.hpp"
#include "integrator.hpp"
#include "io_hdf5.hpp"
#include "bc_periodic.hpp"
#include "cmdparser.hpp"
#include "configparser.hpp"

class Vortex;
typedef Vortex T_global;
typedef RK3 T_integrator;
typedef MC2 T_solver;

///
/// @brief Orszag-Tang Vortex Problem
///
class Vortex : public Global
{
public:
  /// constructor
  Vortex(const int shape[3], const int nb, const float64 cc,
         const float64 gamma, const float64 delh)
    : Global(shape, nb, cc, gamma)
  {
    physt = 0.0;
    delx = delh;
    dely = delh;
    delz = delh;
    mpiutil::getLocalRange(0, delz, zrange);
    mpiutil::getLocalRange(1, dely, yrange);
    mpiutil::getLocalRange(2, delx, xrange);

    for(int ix=Lbx-Nb; ix <= Ubx+Nb ;ix++) {
      xig(ix) = xrange[0] + delx*(ix-Lbx+0.5);
    }

    for(int iy=Lby-Nb; iy <= Uby+Nb ;iy++) {
      yig(iy) = yrange[0] + dely*(iy-Lby+0.5);
    }

    for(int iz=Lbz-Nb; iz <= Ubz+Nb ;iz++) {
      zig(iz) = zrange[0] + delz*(iz-Lbz+0.5);
    }
  }

  /// set up initial condition
  void init(ConfigParser &config)
  {
    float64 kx[3], ky[3];

    int dir  = config.getAs<int>("dir");

    boundary = new PeriodicBoundary(Nz, Ny, Nx, Nb);

    switch(dir) {
    case 0: // x-y plane
      kx[0] = 2*common::pi;
      kx[1] = 0.0;
      kx[2] = 0.0;
      ky[0] = 0.0;
      ky[1] = 2*common::pi;
      ky[2] = 0.0;
      break;
    case 1: // z-x plane
      kx[0] = 0.0;
      kx[1] = 0.0;
      kx[2] = 2*common::pi;
      ky[0] = 2*common::pi;
      ky[1] = 0.0;
      ky[2] = 0.0;
      break;
    case 2: // y-z plane
      kx[0] = 0.0;
      kx[1] = 2*common::pi;
      kx[2] = 0.0;
      ky[0] = 0.0;
      ky[1] = 0.0;
      ky[2] = 2*common::pi;
      break;
    default:
      tfm::format(std::cerr, "Error: input direction (%d) is invalid\n", dir);
      exit(-1);
      break;
    }

    int kn = 2-dir;
    int kt = (kn+1) % 3;
    int ks = (kn+2) % 3;

    // parameters
    float64 mpe = config.getAs<float64>("mpe");
    float64 lp  = config.getAs<float64>("lp");

    // proton
    qp  = c / (lp * sqrt(common::pi4*mpe));
    mp  = 1.0;
    qmp = qp/mp;

    // electron
    qe  =-qp;
    me  = mp / mpe;
    qme = qe/me;

    for(int iz=Lbz-1; iz <= Ubz+1 ;iz++) {
      for(int iy=Lby-1; iy <= Uby+1 ;iy++) {
        for(int ix=Lbx-1; ix <= Ubx+1 ;ix++) {
          float64 b0    = 1.0;
          float64 phix  = kx[0]*xig(ix) + kx[1]*yig(iy) + kx[2]*zig(iz);
          float64 phiy  = ky[0]*xig(ix) + ky[1]*yig(iy) + ky[2]*zig(iz);
          float64 sinx  = sin(phix);
          float64 siny  = sin(phiy);
          float64 sin2x = sin(2*phix);
          float64 rho   = gam*gam / common::pi4;
          float64 uz    = c*b0/(qmp*rho)*(cos(2*phix)+0.5*cos(phiy));
          float64 vx    =-0.5*siny;
          float64 vy    =+0.5*sinx;
          float64 rgm   = sqrt( (1+uz*uz)/(1-vx*vx-vy*vy) );
          float64 vz    = uz/rgm;
          float64 p     = gam / common::pi4;
          float64 bx    =-b0*siny;
          float64 by    = b0*sin2x;
          float64 bz    = 0.0;

          // proton
          uf(iz,iy,ix,0)     = rho / (mp + me);
          uf(iz,iy,ix,kn+1)  =+uz;
          uf(iz,iy,ix,kt+1)  = vx * rgm;
          uf(iz,iy,ix,ks+1)  = vy * rgm;
          uf(iz,iy,ix,4)     = p / 2;
          // electron
          uf(iz,iy,ix,5)     = rho / (mp + me);
          uf(iz,iy,ix,kn+6)  =-uz;
          uf(iz,iy,ix,kt+6)  = vx * rgm;
          uf(iz,iy,ix,ks+6)  = vy * rgm;
          uf(iz,iy,ix,9)     = p / 2;
          // e-field
          ebc(iz,iy,ix,kn)   =-(vx*by - vy*bx)*rc;
          ebc(iz,iy,ix,kt)   =-(vy*bz - vz*by)*rc;
          ebc(iz,iy,ix,ks)   =-(vz*bx - vx*bz)*rc;
          // b-field
          ebc(iz,iy,ix,kn+3) = bz;
          ebc(iz,iy,ix,kt+3) = bx;
          ebc(iz,iy,ix,ks+3) = by;
        }
      }
    }

    int_c2f<2>(ebc, ueb);
    boundary->set_field(*this, ueb);
    boundary->set_fluid(*this, uf, ebc);

    // initialize temporary
    vf  = uf;
    veb = ueb;

    // output message
    msg += "# Execute simulation with the following parameters\n";
    msg += tfm::format("%-20s = %20d\n", "direction", dir);
    msg += tfm::format("%-20s = %20.12e\n", "speed of light", c);
    msg += tfm::format("%-20s = %20.12e\n", "delx", delx);
    msg += tfm::format("%-20s = %20.12e\n", "dely", dely);
    msg += tfm::format("%-20s = %20.12e\n", "delz", delz);
    msg += tfm::format("%-20s = %20.12e\n", "gamma", gam);
    msg += tfm::format("%-20s = %20.12e\n", "mass ratio", mpe);
    msg += tfm::format("%-20s = %20.12e\n", "inertia length", lp);
  }
};

//
// *** main loop ***
//
int main(int argc, char** argv)
{
  const int Nb = 2;

  int shape[3], domain[3];
  int nw, Ng, Nt, dir;
  float64 cc, gamma, tmax, cfl, delh;
  std::string cfg, out;

  // parse command line
  CmdParser parser("vortex.cfg");
  parser.parse_check(argc, argv);

  cfg       = parser.get<std::string>("config");
  domain[0] = parser.get<int>("zdomain");
  domain[1] = parser.get<int>("ydomain");
  domain[2] = parser.get<int>("xdomain");

  // read configuration file
  ConfigParser config(cfg.c_str());
  Global::debugstream.redirect(config.getAs<std::string>("debug"));
  out   = config.getAs<std::string>("out") + ".h5";
  dir   = config.getAs<int>("dir");
  Nt    = config.getAs<int>("Nt");
  Ng    = config.getAs<int>("Ng");
  tmax  = config.getAs<float64>("tmax");
  cfl   = config.getAs<float64>("cfl");
  cc    = config.getAs<float64>("cc");
  gamma = config.getAs<float64>("gamma");
  delh  = 1.0/static_cast<float64>(Ng);

  shape[0] = dir == 0 ? 1 : Ng;
  shape[1] = dir == 1 ? 1 : Ng;
  shape[2] = dir == 2 ? 1 : Ng;
  mpiutil::initialize(&argc, &argv, shape, domain);

  // global
  T_global global(shape, Nb, cc, gamma, delh);
  global.init(config);

  // solver
  T_solver solver(global);

  // integrator
  T_integrator integrator(global);

  // show message
  global.print_message(std::cerr);

  // initial data
  HDF5IO io;
  io.write_field(global, out, true);

  // main loop
  nw = 1;
  global.delt = tmax / Nt;
  while(global.physt < tmax) {
    global.delt *= global.get_dt_factor(cfl);

    // push
    integrator.push(global.delt, global, solver);
    global.physt += global.delt;

    // diagnostics
    if( global.physt + 0.5*global.delt >= tmax*nw / Nt ) {
      float64 Ef, Ep;

      global.energy(Ef, Ep);
      tfm::format(std::cout, "%15.8e %15.8e %15.8e %15.8e\n",
                  global.physt, Ef, Ep, Ef+Ep);
      io.write_field(global, out);
      nw = nw + 1;
    }
  }

  mpiutil::finalize();

  return 0;
}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
