// -*- C++ -*-

///
/// @brief Resistive Diffusion Problem
///
/// $Id: diffusion.cpp,v 975a18d6e205 2016/02/12 14:01:29 amano $
///
#include "global.hpp"
#include "solver.hpp"
#include "integrator.hpp"
#include "io_hdf5.hpp"
#include "bc_periodic.hpp"
#include "bc_conducting.hpp"
#include "cmdparser.hpp"
#include "configparser.hpp"

class Diffusion;
typedef Diffusion T_global;
typedef RK3 T_integrator;
typedef MC2 T_solver;

///
/// @brief Boundary condition
///
class Boundary : public PeriodicBoundary
{
protected:
  ConductingBoundary wall;

public:
  Boundary(const int Nz, const int Ny, const int Nx, const int Nb)
    : PeriodicBoundary(Nz, Ny, Nx, Nb)
  {
  }

  virtual void set_field(Global &g, T_vector &eb, int nb=-1)
  {
    const int Nb = (nb < 0) ? g.Nb : nb;

    PeriodicBoundary::set_field(g, eb, Nb);

    wall.set_field_x(g, eb, Nb);
  }

  virtual void set_fluid(Global &g, T_vector &uf, T_vector &eb, int nb=-1)
  {
    const int Nb = (nb < 0) ? g.Nb : nb;

    PeriodicBoundary::set_fluid(g, uf, eb, Nb);

    wall.set_fluid_x(g, uf, eb, Nb);
  }
};

///
/// @brief Resistive Diffusion Problem
///
class Diffusion : public Global
{
public:
  /// constructor
  Diffusion(const int shape[3], const int nb, const float64 cc,
            const float64 gamma)
    : Global(shape, nb, cc, gamma)
  {
    const float64 L = 3.0;
    int global_shape[3], global_offset[3];

    mpiutil::getGlobalShape(global_shape);
    mpiutil::getGlobalOffset(global_offset);

    physt = 0.0;
    delz  = L / static_cast<float64>(global_shape[0]);
    dely  = L / static_cast<float64>(global_shape[1]);
    delx  = L / static_cast<float64>(global_shape[2]);

    zrange[0] = 0.0;
    zrange[1] = delz*Nz + zrange[0];
    zrange[2] = zrange[1] - zrange[0];
    yrange[0] = 0.0;
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
    boundary = new Boundary(Nz, Ny, Nx, Nb);

    // parameters
    float64 mpe  = 1.0;
    float64 D    = 0.01;
    float64 lp   = config.getAs<float64>("lp");
    float64 b0   = sqrt(common::pi4);

    // resistivity
    eta  = D * common::pi4 / (c*c);

    // proton
    qp  = c / (lp * sqrt(common::pi4*mpe));
    mp  = 1.0;
    qmp = qp/mp;

    // electron
    qe  =-qp;
    me  = mp / mpe;
    qme = qe/me;

    float64 u0   = c*b0/(common::pi4*qmp * sqrt(common::pi*D));

    for(int iz=Lbz-1; iz <= Ubz+1 ;iz++) {
      for(int iy=Lby-1; iy <= Uby+1 ;iy++) {
        for(int ix=Lbx-1; ix <= Ubx+1 ;ix++) {
          float64 xx    = 0.5 * xig(ix) / sqrt(D);
          float64 rho   = 1.0;
          float64 p     = 50.0;
          float64 uz    = u0 * exp(-xx*xx);
          float64 by    = b0 * erf(xx);

          // proton
          uf(iz,iy,ix,0)  = rho / 2;
          uf(iz,iy,ix,1)  = 0.0;
          uf(iz,iy,ix,2)  = 0.0;
          uf(iz,iy,ix,3)  =+uz;
          uf(iz,iy,ix,4)  = p / 2;
          // electron
          uf(iz,iy,ix,5)  = rho / 2;
          uf(iz,iy,ix,6)  = 0.0;
          uf(iz,iy,ix,7)  = 0.0;
          uf(iz,iy,ix,8)  =-uz;
          uf(iz,iy,ix,9)  = p / 2;
          // e-field
          ebc(iz,iy,ix,0) = 0.0;
          ebc(iz,iy,ix,1) = 0.0;
          ebc(iz,iy,ix,2) = 0.0;
          // b-field
          ebc(iz,iy,ix,3) = 0.0;
          ebc(iz,iy,ix,4) = by;
          ebc(iz,iy,ix,5) = 0.0;
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

  int shape[3], domain[3], period[3];
  int nw, Ng, Nt, dir;
  float64 cc, gamma, tmax, cfl;
  std::string cfg, out;

  // parse command line
  CmdParser parser("diffusion.cfg");
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

  shape[0] = 1;
  shape[1] = 1;
  shape[2] = Ng;
  period[0] = 1;
  period[1] = 1;
  period[2] = 0;
  mpiutil::initialize(&argc, &argv, shape, domain, period);

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
