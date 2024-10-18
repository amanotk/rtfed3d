// -*- C++ -*-

///
/// @brief Blast Wave Problem
///
/// $Id$
///
#include "global.hpp"
#include "solver.hpp"
#include "integrator.hpp"
#include "io_hdf5.hpp"
#include "bc_periodic.hpp"
#include "cmdparser.hpp"
#include "configparser.hpp"

class Blastwave;
typedef Blastwave T_global;
typedef RK3 T_integrator;
typedef MC2 T_solver;

///
/// @brief Blast Wave Problem
///
class Blastwave : public Global
{
public:
  /// constructor
  Blastwave(const int shape[3], const int nb,
            const float64 cc, const float64 gamma)
    : Global(shape, nb, cc, gamma)
  {
    const float64 L = 12.0;
    int global_shape[3], global_offset[3];

    mpiutil::getGlobalShape(global_shape);
    mpiutil::getGlobalOffset(global_offset);

    physt = 0.0;
    delz  = L / static_cast<float64>(global_shape[0]);
    dely  = L / static_cast<float64>(global_shape[1]);
    delx  = L / static_cast<float64>(global_shape[2]);

    zrange[0] = delz*global_offset[0] - 0.5*L;
    zrange[1] = delz*Nz + zrange[0];
    zrange[2] = zrange[1] - zrange[0];
    yrange[0] = dely*global_offset[1] - 0.5*L;
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
    float64 R  = 0.8;
    float64 dR = 1.0 - R;

    boundary = new PeriodicBoundary(Nz, Ny, Nx, Nb);

    // parameters
    float64 mpe = config.getAs<float64>("mpe");
    float64 lp  = config.getAs<float64>("lp");
    float64 ri  = config.getAs<float64>("ri");
    float64 pi  = config.getAs<float64>("pi");
    float64 ro  = config.getAs<float64>("ro");
    float64 po  = config.getAs<float64>("po");
    float64 bx  = config.getAs<float64>("bx");
    float64 by  = config.getAs<float64>("by");
    float64 bz  = config.getAs<float64>("bz");

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
          float64 r = sqrt(xig(ix)*xig(ix) + yig(iy)*yig(iy) + zig(iz)*zig(iz));

          // density and pressure
          float64 rho = (ro - ri)*(r - R)/dR + ri;
          float64 p   = (po - pi)*(r - R)/dR + pi;
          rho = limiter::min(ri, limiter::max(ro, rho));
          p   = limiter::min(pi, limiter::max(po, p));

          // proton
          uf(iz,iy,ix,0)  = rho / (mp + me);
          uf(iz,iy,ix,1)  = 0.0;
          uf(iz,iy,ix,2)  = 0.0;
          uf(iz,iy,ix,3)  = 0.0;
          uf(iz,iy,ix,4)  = p / 2;
          // electron
          uf(iz,iy,ix,5)  = rho / (mp + me);
          uf(iz,iy,ix,6)  = 0.0;
          uf(iz,iy,ix,7)  = 0.0;
          uf(iz,iy,ix,8)  = 0.0;
          uf(iz,iy,ix,9)  = p / 2;
          // e-field
          ebc(iz,iy,ix,0) = 0.0;
          ebc(iz,iy,ix,1) = 0.0;
          ebc(iz,iy,ix,2) = 0.0;
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
    msg += tfm::format("%-20s = %20.12e\n", "density inside", ri);
    msg += tfm::format("%-20s = %20.12e\n", "pressure inside", pi);
    msg += tfm::format("%-20s = %20.12e\n", "density outside", ro);
    msg += tfm::format("%-20s = %20.12e\n", "pressure outside", po);
    msg += tfm::format("%-20s = %20.12e\n", "Bx" , bx);
    msg += tfm::format("%-20s = %20.12e\n", "By" , by);
    msg += tfm::format("%-20s = %20.12e\n", "Bz" , bz);
  }
};

//
// *** main loop ***
//
int main(int argc, char** argv)
{
  const int Nb = 2;

  int shape[3], domain[3];
  int nw, Nx, Ny, Nz, Nt;
  float64 cc, gamma, tmax, cfl;
  std::string cfg, out;

  // parse command line
  CmdParser parser("blastwave.cfg");
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
  Nx    = config.getAs<int>("Nx");
  Ny    = config.getAs<int>("Ny");
  Nz    = config.getAs<int>("Nz");
  tmax  = config.getAs<float64>("tmax");
  cfl   = config.getAs<float64>("cfl");
  cc    = config.getAs<float64>("cc");
  gamma = config.getAs<float64>("gamma");

  shape[0] = Nz;
  shape[1] = Ny;
  shape[2] = Nx;
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
