// -*- C++ -*-

///
/// @brief 1D Riemann Problem
///
/// Initial left and right states must be in MHD equilibria.
///
/// $Id: riemann1d.cpp,v 641a96a86755 2015/09/14 19:35:09 amano $
///
#include "global.hpp"
#include "solver.hpp"
#include "integrator.hpp"
#include "io_hdf5.hpp"
#include "bc_periodic.hpp"
#include "bc_symmetric.hpp"
#include "cmdparser.hpp"
#include "configparser.hpp"

class Riemann;
typedef Riemann T_global;
typedef RK3 T_integrator;
typedef MC2 T_solver;

///
/// @brief Boundary condition
///
class Boundary : public PeriodicBoundary
{
protected:
  int symdir;
  SymmetricBoundary symmetric;

public:
  Boundary(const int Nz, const int Ny, const int Nx, const int Nb,
           const int dir=0)
  : PeriodicBoundary(Nz, Ny, Nx, Nb), symdir(dir)
  {
  }

  virtual void set_field(Global &g, T_vector &eb, int nb=-1)
  {
    const int Nb = (nb < 0) ? g.Nb : nb;

    PeriodicBoundary::set_field(g, eb, Nb);

    switch(symdir) {
    case 0:
      symmetric.set_field_z(g, eb, Nb);
      break;
    case 1:
      symmetric.set_field_y(g, eb, Nb);
      break;
    case 2:
      symmetric.set_field_x(g, eb, Nb);
      break;
    }
  }

  virtual void set_fluid(Global &g, T_vector &uf, T_vector &eb, int nb=-1)
  {
    const int Nb = (nb < 0) ? g.Nb : nb;

    PeriodicBoundary::set_fluid(g, uf, eb, Nb);

    switch(symdir) {
    case 0:
      symmetric.set_fluid_z(g, uf, eb, Nb);
      break;
    case 1:
      symmetric.set_fluid_y(g, uf, eb, Nb);
      break;
    case 2:
      symmetric.set_fluid_x(g, uf, eb, Nb);
      break;
    }
  }
};

///
/// @brief 1D Riemann Problem
///
class Riemann : public Global
{
public:
  /// constructor
  Riemann(const int shape[3], const int nb, const float64 cc,
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
    float64 ul[8], ur[8], *uu[2][2][2];

    int dir = config.getAs<int>("dir");

    boundary = new Boundary(Nz, Ny, Nx, Nb, dir);

    switch(dir) {
    case 0: // along z dir
      uu[0][0][0] = ul;
      uu[1][0][0] = ur;
      break;
    case 1: // along y dir
      uu[0][0][0] = ul;
      uu[0][1][0] = ur;
      break;
    case 2: // along x dir
      uu[0][0][0] = ul;
      uu[0][0][1] = ur;
      break;
    default:
      tfm::format(std::cerr, "Error: input direction (%d) is invalid\n", dir);
      exit(-1);
      break;
    }

    // index for vector quantities
    int kn = 2 - dir;
    int kt = (kn+1) % 3;
    int ks = (kn+2) % 3;

    float64 bn  = config.getAs<float64>("bn");

    // left state
    ul[0]    = config.getAs<float64>("n_l");
    ul[kn+1] = config.getAs<float64>("un_l");
    ul[kt+1] = config.getAs<float64>("ut_l");
    ul[ks+1] = 0.0;
    ul[4]    = config.getAs<float64>("p_l");
    ul[kn+5] = bn;
    ul[kt+5] = config.getAs<float64>("bt_l");
    ul[ks+5] = 0.0;
    // right state
    ur[0]    = config.getAs<float64>("n_r");
    ur[kn+1] = config.getAs<float64>("un_r");
    ur[kt+1] = config.getAs<float64>("ut_r");
    ur[ks+1] = 0.0;
    ur[4]    = config.getAs<float64>("p_r");
    ur[kn+5] = bn;
    ur[kt+5] = config.getAs<float64>("bt_r");
    ur[ks+5] = 0.0;

    // parameters
    float64 mpe = config.getAs<float64>("mpe");
    float64 lp  = config.getAs<float64>("lp");
    float64 tau = config.getAs<float64>("tau");

    eta = lp*lp / tau;

    // proton
    qp  = c / (lp * sqrt(common::pi4*mpe));
    mp  = 1.0;
    qmp = qp/mp;

    // electron
    qe  =-qp;
    me  = mp / mpe;
    qme = qe/me;

    // set initial condition
    for(int iz=Lbz-Nb; iz <= Ubz+Nb ;iz++) {
      int jz = ((dir == 0) && (zig(iz) > 0.5)) ? 1 : 0;
      for(int iy=Lby-Nb; iy <= Uby+Nb ;iy++) {
        int jy = ((dir == 1) && (yig(iy) > 0.5)) ? 1 : 0;
        for(int ix=Lbx-Nb; ix <= Ubx+Nb ;ix++) {
          int jx = ((dir == 2) && (xig(ix) > 0.5)) ? 1 : 0;
          float64 *v = uu[jz][jy][jx];
          // proton
          uf(iz,iy,ix,0) = v[0] / (mp + me);
          uf(iz,iy,ix,1) = v[1];
          uf(iz,iy,ix,2) = v[2];
          uf(iz,iy,ix,3) = v[3];
          uf(iz,iy,ix,4) = v[4] / 2;
          // electron
          uf(iz,iy,ix,5) = v[0] / (mp + me);
          uf(iz,iy,ix,6) = v[1];
          uf(iz,iy,ix,7) = v[2];
          uf(iz,iy,ix,8) = v[3];
          uf(iz,iy,ix,9) = v[4] / 2;
          // e-field
          float64 ux = v[1];
          float64 uy = v[2];
          float64 uz = v[3];
          float64 bx = v[5];
          float64 by = v[6];
          float64 bz = v[7];
          float64 gm = sqrt(1 + ux*ux + uy*uy + uz*uz);
          ebc(iz,iy,ix,0) =-(uy*bz - uz*by)/gm;
          ebc(iz,iy,ix,1) =-(uz*bx - ux*bz)/gm;
          ebc(iz,iy,ix,2) =-(ux*by - uy*bx)/gm;
          // b-field
          ebc(iz,iy,ix,3) = v[5];
          ebc(iz,iy,ix,4) = v[6];
          ebc(iz,iy,ix,5) = v[7];
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
    msg += "--- left state ---\n";
    msg += tfm::format("%-20s = %20.12e\n", "rho", ul[0]);
    msg += tfm::format("%-20s = %20.12e\n", "Ux" , ul[1]);
    msg += tfm::format("%-20s = %20.12e\n", "Uy" , ul[2]);
    msg += tfm::format("%-20s = %20.12e\n", "Uz" , ul[3]);
    msg += tfm::format("%-20s = %20.12e\n", "P"  , ul[4]);
    msg += tfm::format("%-20s = %20.12e\n", "Bx" , ul[5]);
    msg += tfm::format("%-20s = %20.12e\n", "By" , ul[6]);
    msg += tfm::format("%-20s = %20.12e\n", "Bz" , ul[7]);
    msg += "--- right state ---\n";
    msg += tfm::format("%-20s = %20.12e\n", "rho", ur[0]);
    msg += tfm::format("%-20s = %20.12e\n", "Ux" , ur[1]);
    msg += tfm::format("%-20s = %20.12e\n", "Uy" , ur[2]);
    msg += tfm::format("%-20s = %20.12e\n", "Uz" , ur[3]);
    msg += tfm::format("%-20s = %20.12e\n", "P"  , ur[4]);
    msg += tfm::format("%-20s = %20.12e\n", "Bx" , ur[5]);
    msg += tfm::format("%-20s = %20.12e\n", "By" , ur[6]);
    msg += tfm::format("%-20s = %20.12e\n", "Bz" , ur[7]);
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
  float64 cc, gamma, tmax, cfl, delh;
  std::string cfg, out;

  // parse command line
  CmdParser parser("riemann1d.cfg");
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

  shape[0]  = dir == 0 ? Ng : 1;
  shape[1]  = dir == 1 ? Ng : 1;
  shape[2]  = dir == 2 ? Ng : 1;
  period[0] = dir == 0 ? 0 : 1;
  period[1] = dir == 1 ? 0 : 1;
  period[2] = dir == 2 ? 0 : 1;
  mpiutil::initialize(&argc, &argv, shape, domain, period);

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
