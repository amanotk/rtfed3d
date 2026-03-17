// -*- C++ -*-

///
/// @brief Tearing Instability
///
/// $Id: tearing.cpp,v e4aecbe07eb4 2016/02/02 06:28:09 amano $
///
#include "global.hpp"
#include "solver.hpp"
#include "integrator.hpp"
#include "io_hdf5.hpp"
#include "bc_periodic.hpp"
#include "cmdparser.hpp"
#include "configparser.hpp"

class Tearing;
typedef Tearing T_global;
typedef RK3 T_integrator;
typedef MC2 T_solver;

///
/// @brief Tearing Instabiity
///
class Tearing : public Global
{
public:
  /// constructor
  Tearing(const int shape[3], const int nb, const float64 cc,
          const float64 gamma)
    : Global(shape, nb, cc, gamma)
  {
    const float64 L = 2.0;
    int global_shape[3], global_offset[3];

    mpiutil::getGlobalShape(global_shape);
    mpiutil::getGlobalOffset(global_offset);

    physt = 0.0;
    delz  = 1*L / static_cast<float64>(global_shape[0]);
    dely  = 2*L / static_cast<float64>(global_shape[1]);
    delx  = 1*L / static_cast<float64>(global_shape[2]);

    zrange[0] = delz*global_offset[0] - 0.5*L;
    zrange[1] = delz*Nz + zrange[0];
    zrange[2] = zrange[1] - zrange[0];
    yrange[0] = dely*global_offset[1] - 1.0*L;
    yrange[1] = dely*Ny + yrange[0];
    yrange[2] = yrange[1] - yrange[0];
    xrange[0] = delx*global_offset[2] - 0.5*L;
    xrange[1] = delx*Nx + xrange[0];
    xrange[2] = xrange[1] - xrange[0];

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
    float64 sigma = config.getAs<float64>("sigma");
    float64 mpe   = config.getAs<float64>("mpe");
    float64 lp    = config.getAs<float64>("lp");
    float64 lcs   = config.getAs<float64>("lcs");
    float64 b0    = sqrt(common::pi4*sigma);
    float64 bb    = 1.0e-3 * b0;
    float64 rho   = 1.0;
    float64 pr    = 0.1;
    float64 S     = config.getAs<float64>("S");

    //
    // set resistivity
    //
    // NOTE:
    //     The magnetic Reynolds number is defined by
    //         S = 4 pi L / eta c
    //
    eta = common::pi4 * lcs/(c*S);

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
          float64 tanh1 = tanh((yig(iy)-1)/lcs);
          float64 tanh2 = tanh((yig(iy)+1)/lcs);
          float64 sech1 = 1/cosh((yig(iy)-1)/lcs);
          float64 sech2 = 1/cosh((yig(iy)+1)/lcs);
          float64 sinx  = sin(common::pi*xig(ix));
          float64 u0    = c*b0/(common::pi4*rho*qp*lcs);

          // proton
          uf(iz,iy,ix,0)     = rho / (mp + me);
          uf(iz,iy,ix,kn+1)  =-u0*(sech1*sech1 - sech2*sech2);
          uf(iz,iy,ix,kt+1)  =-u0*(sech1*tanh1 - sech2*tanh2);
          uf(iz,iy,ix,ks+1)  = 0.0;
          uf(iz,iy,ix,4)     = pr / 2;
          // electron
          uf(iz,iy,ix,5)     = rho / (mp + me);
          uf(iz,iy,ix,kn+6)  =+u0*(sech1*sech1 - sech2*sech2);
          uf(iz,iy,ix,kt+6)  =+u0*(sech1*tanh1 - sech2*tanh2);
          uf(iz,iy,ix,ks+6)  = 0.0;
          uf(iz,iy,ix,9)     = pr / 2;
          // e-field
          ebc(iz,iy,ix,kn)   = 0.0;
          ebc(iz,iy,ix,kt)   = 0.0;
          ebc(iz,iy,ix,ks)   = 0.0;
          // b-field
          ebc(iz,iy,ix,kn+3) = b0*(sech1 - sech2);
          ebc(iz,iy,ix,kt+3) = b0*(tanh1 - tanh2 + 1);
          ebc(iz,iy,ix,ks+3) = bb*sinx;
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
  CmdParser parser("tearing.cfg");
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

  shape[0] = 1;
  shape[1] = Ng*2;
  shape[2] = Ng;
  mpiutil::initialize(&argc, &argv, shape, domain);

  // global
  T_global global(shape, Nb, cc, gamma);
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
