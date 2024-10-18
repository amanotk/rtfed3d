// -*- C++ -*-

///
/// @brief Magnetic Reconnection Problem
///
/// $Id: reconnection.cpp,v 267ff9a5277d 2016/03/16 09:14:48 amano $
///
#include "global.hpp"
#include "solver.hpp"
#include "integrator.hpp"
#include "io_hdf5.hpp"
#include "bc_periodic.hpp"
#include "bc_conducting.hpp"
#include "cmdparser.hpp"
#include "configparser.hpp"

class Reconnection;
typedef Reconnection T_global;
typedef RK3 T_integrator;
typedef MC2 T_solver;

///
/// @brief Boundary condition
///
class Boundary : public PeriodicBoundary
{
protected:
  int walldir;
  ConductingBoundary wall;

public:
  Boundary(const int Nz, const int Ny, const int Nx, const int Nb,
           const int dir=0)
  : PeriodicBoundary(Nz, Ny, Nx, Nb), walldir(dir)
  {
  }

  virtual void set_field(Global &g, T_vector &eb, int nb=-1)
  {
    const int Nb = (nb < 0) ? g.Nb : nb;

    PeriodicBoundary::set_field(g, eb, Nb);

    switch(walldir) {
    case 0:
      wall.set_field_y(g, eb, Nb);
      break;
    case 1:
      wall.set_field_x(g, eb, Nb);
      break;
    case 2:
      wall.set_field_z(g, eb, Nb);
      break;
    }
  }

  virtual void set_fluid(Global &g, T_vector &uf, T_vector &eb, int nb=-1)
  {
    const int Nb = (nb < 0) ? g.Nb : nb;

    PeriodicBoundary::set_fluid(g, uf, eb, Nb);

    switch(walldir) {
    case 0:
      wall.set_fluid_y(g, uf, eb, Nb);
      break;
    case 1:
      wall.set_fluid_x(g, uf, eb, Nb);
      break;
    case 2:
      wall.set_fluid_z(g, uf, eb, Nb);
      break;
    }
  }
};

///
/// @brief Magnetic Reconnection Problem
///
class Reconnection : public Global
{
public:
  /// constructor
  Reconnection(const int shape[3], const int nb,
               const float64 cc, const float64 gamma)
    : Global(shape, nb, cc, gamma)
  {
    physt = 0.0;
  }

  /// set up initial condition
  void init(ConfigParser &config)
  {
    int global_offset[3];
    int dir   = config.getAs<int>("dir");
    int Ng    = config.getAs<int>("Ng");
    float64 L = config.getAs<float64>("L");
    blitz::TinyVector<float64,3> kx, ky, kz;

    mpiutil::getGlobalOffset(global_offset);

    boundary = new Boundary(Nz, Ny, Nx, Nb, dir);

    // set grid
    {
      float64 dh = L / static_cast<float64>(Ng);

      delx = dh;
      dely = dh;
      delz = dh;

      switch(dir) {
      case 0: // x-y plane
        kx = 1.0, 0.0, 0.0;
        ky = 0.0, 1.0, 0.0;
        kz = 0.0, 0.0, 1.0;
        xrange[0] = delx * global_offset[2] - L;
        xrange[1] = xrange[0] + Nx * delx;
        xrange[2] = xrange[1] - xrange[0];
        yrange[0] = dely * global_offset[1] - L * 0.5;
        yrange[1] = yrange[0] + Ny * dely;
        yrange[2] = yrange[1] - yrange[0];
        zrange[0] = delz * global_offset[0];
        zrange[1] = zrange[0] + Nz * delz;
        zrange[2] = zrange[1] - zrange[0];
        break;
      case 1: // z-x plane
        ky = 1.0, 0.0, 0.0;
        kz = 0.0, 1.0, 0.0;
        kx = 0.0, 0.0, 1.0;
        zrange[0] = delz * global_offset[0] - L;
        zrange[1] = zrange[0] + Nz * delz;
        zrange[2] = zrange[1] - zrange[0];
        xrange[0] = delx * global_offset[2] - L * 0.5;
        xrange[1] = xrange[0] + Nx * delx;
        xrange[2] = xrange[1] - xrange[0];
        yrange[0] = dely * global_offset[1];
        yrange[1] = yrange[0] + Ny * dely;
        yrange[2] = yrange[1] - yrange[0];
        break;
      case 2: // y-z plane
        kz = 1.0, 0.0, 0.0;
        kx = 0.0, 1.0, 0.0;
        ky = 0.0, 0.0, 1.0;
        yrange[0] = dely * global_offset[1] - L;
        yrange[1] = yrange[0] + Ny * dely;
        yrange[2] = yrange[1] - yrange[0];
        zrange[0] = delz * global_offset[0] - L * 0.5;
        zrange[1] = zrange[0] + Nz * delz;
        zrange[2] = zrange[1] - zrange[0];
        xrange[0] = delx * global_offset[2];
        xrange[1] = xrange[1] + Nx * delx;
        xrange[2] = xrange[1] - xrange[0];
        break;
      default:
        tfm::format(std::cerr, "Error: input direction (%d) is invalid\n", dir);
        break;
      }

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

    int kn = 2-dir;
    int kt = (kn+1) % 3;
    int ks = (kn+2) % 3;

    // parameters
    float64 mpe   = config.getAs<float64>("mpe");
    float64 sigma = config.getAs<float64>("sigma");
    float64 dA    = config.getAs<float64>("dA");
    float64 lcs   = config.getAs<float64>("lcs");
    float64 nbg   = config.getAs<float64>("nbg");
    float64 S     = config.getAs<float64>("S");
    float64 wpp   = 1.0/sqrt(sigma);

    //
    // set resistivity
    //
    // NOTE:
    //     The magnetic Reynolds number is defined by
    //         S = 4 pi L sigma^(1/2) / eta c
    //
    eta = common::pi4 * lcs/(c*S) * sqrt(sigma);

    // mass
    mp = 1.0;
    me = 1.0 / mpe;

    // magnetic field
    float64 b0 = sqrt(common::pi4 * sigma);

    // charge
    qp  = 1.0 / sqrt(common::pi4 * sigma);
    qe  =-qp;
    qmp = qp/mp;
    qme = qe/me;

    // fluid quantities
    float64 n0   = 1.0;
    float64 pr0  = 0.25 * sigma;
    float64 prbg = pr0 * nbg;
    float64 uz0  = (c*b0)/(2*common::pi4*lcs * n0*qp);

    // mode number for perturbation
    float64 dk = common::pi/L;
    float64 db = 0.25 * dA * b0;

    for(int iz=Lbz-Nb; iz <= Ubz+Nb ;iz++) {
      for(int iy=Lby-Nb; iy <= Uby+Nb ;iy++) {
        for(int ix=Lbx-Nb; ix <= Ubx+Nb ;ix++) {
          float64 xx    = kx[0]*xig(ix) + kx[1]*yig(iy) + kx[2]*zig(iz);
          float64 yy    = ky[0]*xig(ix) + ky[1]*yig(iy) + ky[2]*zig(iz);
          float64 zz    = kz[0]*xig(ix) + kz[1]*yig(iy) + kz[2]*zig(iz);
          float64 coshy = cosh(yy/lcs);
          float64 tanhy = tanh(yy/lcs);
          float64 sechy = 1.0/(coshy*coshy);
          float64 bx    =-db * cos(dk*xx) * sin(dk*yy);
          float64 by    =+db * sin(dk*xx) * cos(dk*yy);
          float64 bz    = 0.0;
          float64 uz    =-uz0 * sechy/(sechy + nbg/n0);

          // proton
          uf(iz,iy,ix,0)     = n0*sechy + nbg;
          uf(iz,iy,ix,kn+1)  =+uz;
          uf(iz,iy,ix,kt+1)  = 0.0;
          uf(iz,iy,ix,ks+1)  = 0.0;
          uf(iz,iy,ix,4)     = pr0 * sechy + prbg;
          // electron
          uf(iz,iy,ix,5)     = n0*sechy + nbg;
          uf(iz,iy,ix,kn+6)  =-uz;
          uf(iz,iy,ix,kt+6)  = 0.0;
          uf(iz,iy,ix,ks+6)  = 0.0;
          uf(iz,iy,ix,9)     = pr0 * sechy + prbg;
          // e-field
          ebc(iz,iy,ix,kn)   = 0.0;
          ebc(iz,iy,ix,kt)   = 0.0;
          ebc(iz,iy,ix,ks)   = 0.0;
          // b-field
          ebc(iz,iy,ix,kn+3) = bz;
          ebc(iz,iy,ix,kt+3) = bx + b0*tanhy;
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
  }
};

//
// *** main loop ***
//
int main(int argc, char** argv)
{
  const int Nb = 2;

  int shape[3], domain[3], period[3];
  int nw, Nx, Ny, Nz, Ng, Nt, dir;
  float64 cc, gamma, tmax, delt;
  std::string cfg, out;

  // parse command line
  CmdParser parser("reconnection.cfg");
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
  delt  = config.getAs<float64>("delt");
  cc    = config.getAs<float64>("cc");
  gamma = config.getAs<float64>("gamma");

  period[0] = 1;
  period[1] = 1;
  period[2] = 1;

  switch(dir) {
  case 0:
    shape[0]  = 1;
    shape[1]  = Ng;
    shape[2]  = 2*Ng;
    period[1] = 0;
    break;
  case 1:
    shape[0]  = 2*Ng;
    shape[1]  = 1;
    shape[2]  = Ng;
    period[2] = 0;
    break;
  case 2:
    shape[0]  = Ng;
    shape[1]  = 2*Ng;
    shape[2]  = 1;
    period[0] = 0;
    break;
  }

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
  global.delt = delt;
  while(global.physt < tmax) {
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
