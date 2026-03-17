// -*- C++ -*-

///
/// @brief Current Sheets
///
/// $Id: sheets.cpp,v f3d1b0f8cdaa 2015/09/07 13:33:49 amano $
///
#include "global.hpp"
#include "solver.hpp"
#include "integrator.hpp"
#include "io_hdf5.hpp"
#include "bc_periodic.hpp"
#include "cmdparser.hpp"
#include "configparser.hpp"

class Sheets;
typedef Sheets T_global;
typedef RK3 T_integrator;
typedef MC2 T_solver;

///
/// @brief Current Sheets
///
class Sheets : public Global
{
public:
  /// constructor
  Sheets(const int shape[3], const int nb, const float64 cc,
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
    boundary = new PeriodicBoundary(Nz, Ny, Nx, Nb);

    // parameters
    float64 mpe  = config.getAs<float64>("mpe");
    float64 lp   = config.getAs<float64>("lp");
    float64 beta = config.getAs<float64>("beta");
    float64 v0   = config.getAs<float64>("v0");

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
          float64 b0    = sqrt(common::pi4);
          float64 rho   = 1.0;
          float64 vx    = v0 * cos(common::pi * yig(iy));
          float64 vy    = 0.0;
          float64 vz    = 0.0;
          float64 rgm   = 1.0 / sqrt(1 - vx*vx - vy*vy - vz*vz);
          float64 p     = 0.5 * beta;
          float64 bx    = 0.0;
          float64 by    = fabs(xig(ix) - 1.0) < 0.5 ? b0 : -b0;
          float64 bz    = 0.0;

          // proton
          uf(iz,iy,ix,0)  = rho / (mp + me);
          uf(iz,iy,ix,1)  = vx * rgm;
          uf(iz,iy,ix,2)  = vy * rgm;
          uf(iz,iy,ix,3)  = vz * rgm;
          uf(iz,iy,ix,4)  = p / 2;
          // electron
          uf(iz,iy,ix,5)  = rho / (mp + me);
          uf(iz,iy,ix,6)  = vx * rgm;
          uf(iz,iy,ix,7)  = vy * rgm;
          uf(iz,iy,ix,8)  = vz * rgm;
          uf(iz,iy,ix,9)  = p / 2;
          // e-field
          ebc(iz,iy,ix,0) =-(vy*bz - vz*by)*rc;
          ebc(iz,iy,ix,1) =-(vz*bx - vx*bz)*rc;
          ebc(iz,iy,ix,2) =-(vx*by - vy*bx)*rc;
          // b-field
          ebc(iz,iy,ix,3) = bx;
          ebc(iz,iy,ix,4) = by;
          ebc(iz,iy,ix,5) = bz;
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
  CmdParser parser("sheets.cfg");
  parser.parse_check(argc, argv);

  cfg       = parser.get<std::string>("config");
  domain[0] = parser.get<int>("zdomain");
  domain[1] = parser.get<int>("ydomain");
  domain[2] = parser.get<int>("xdomain");

  // read configuration file
  ConfigParser config(cfg.c_str());
  Global::debugstream.redirect(config.getAs<std::string>("debug"));
  out   = config.getAs<std::string>("out") + ".h5";
  Nt    = config.getAs<int>("Nt");
  Ng    = config.getAs<int>("Ng");
  tmax  = config.getAs<float64>("tmax");
  cfl   = config.getAs<float64>("cfl");
  cc    = config.getAs<float64>("cc");
  gamma = config.getAs<float64>("gamma");
  delh  = 2.0/static_cast<float64>(Ng);

  shape[0] = 1;
  shape[1] = Ng;
  shape[2] = Ng;
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
