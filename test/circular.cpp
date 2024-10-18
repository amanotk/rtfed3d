// -*- C++ -*-

///
/// @brief Propagation of Circularly Polarized Electromagnetic Wave
///
/// $Id: circular.cpp,v 267ff9a5277d 2016/03/16 09:14:48 amano $
///
#include "global.hpp"
#include "solver.hpp"
#include "integrator.hpp"
#include "io_hdf5.hpp"
#include "bc_periodic.hpp"
#include "cmdparser.hpp"
#include "configparser.hpp"

class Circular;
typedef Circular T_global;
typedef RK3 T_integrator;
typedef MC2 T_solver;

namespace {

///
/// dispersion solver
///
class DispersionSolver
{
private:
  int maxiter;    ///< # maximum iteration
  float64 tol;    ///< tolerance
  float64 wpp;    ///< lab frame proton plasma frequency
  float64 wpe;    ///< lab frame electron plasma frequency
  float64 wcp;    ///< non-relativistic proton cyclotron frequency
  float64 wce;    ///< non-relativistic electron cyclotron frequency

public:
  // constructor
  DispersionSolver(float64 wp, float64 wc, float64 mpe)
  {
    maxiter = 200;
    tol     = 1.0e-12;

    wpp = wp * wp;
    wcp = wc;
    wpe = wpp * mpe;
    wce =-wcp * mpe;
  }

  void f(float64 w, float64 k, float64 gp, float64 ge,
         float64 &ff, float64 &df)
  {
    float64 wwcp = gp*w + wcp;
    float64 wwce = ge*w + wce;

    ff = w*w - k*k - wpp * w/wwcp - wpe * w/wwce;
    df = 2*w - wpp * wcp/(wwcp*wwcp) - wpe * wce/(wwce*wwce);
  }

  void g(float64 x, float64 w, float64 k, float64 wc, float64 eta,
         float64 &gg, float64 &dg)
  {
    float64 ww = wc/w;
    float64 vp = eta * w/k;
    float64 a = ww * ww;
    float64 b = vp * vp * a;
    float64 y = sqrt(1 + x);

    gg = x*x + (1 + a - b) * x + 2 * ww * x*y         - b;
    dg = 2*x + (1 + a - b)     +     ww * (2 + 3*x)/y;
  }

  // calculate Lorentz factor
  float64 gamma(float64 w, float64 k, float64 eta, float64 wc, float64 gm)
  {
    bool status = false;
    float64 x = gm*gm - 1;

    for(int i=0; i < maxiter ;i++) {
      float64 dx, gg, dg;
      g(x, w, k, wc, eta, gg, dg);

      dx = -gg/dg;
      x  = x + dx;

      if( fabs(dx) < tol ) {
        status = true;
        break;
      }
    }

    if( status == false ) {
      std::cerr << "Failed to find gamma factor\n";
      exit(-1);
    }

    return sqrt(1 + x);
  }

  // solve dispersion relation
  void solve(float64 k, float64 eta, float64 &w, float64 &gp, float64 &ge)
  {
    bool status = false;

    for(int i=0; i < maxiter ;i++) {
      float64 dw, ff, df;
      gp = gamma(w, k, eta, wcp, gp);
      ge = gamma(w, k, eta, wce, ge);

      f(w, k, gp, ge, ff, df);
      dw = -ff/df;
      w  = w + dw;

      if( fabs(dw) < tol ) {
        status = true;
        break;
      }
    }

    if( status == false ) {
      std::cerr << "Failed to find frequency\n";
      exit(-1);
    }
  }
};

}

///
/// @brief Propagation of Circularly Polarized Electromagnetic Wave
///
class Circular : public Global
{
protected:
  int global_shape[3];
  float64 period;
  float64 xlen;
  float64 ylen;
  float64 zlen;

public:
  /// constructor
  Circular(const int shape[3], const int nb,
           const float64 cc, const float64 gamma, const float64 delh)
    : Global(shape, nb, cc, gamma)
  {
    mpiutil::getGlobalShape(global_shape);

    physt = 0.0;
    zlen = delh * global_shape[0];
    ylen = delh * global_shape[1];
    xlen = delh * global_shape[2];
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
    int N = limiter::max(Nx, Ny, Nz);
    float64 kk, kx, ky, kz;
    float64 ev0[3], ev1[3], ev2[3], ev3[3];

    int mode = config.getAs<int>("mode");

    boundary = new PeriodicBoundary(Nz, Ny, Nx, Nb);

    if( N == Nx ) {
      ev0[0] = 0.0;
      ev0[1] = 0.0;
      ev0[2] = 1.0;
    } else if( N == Ny ) {
      ev0[0] = 1.0;
      ev0[1] = 0.0;
      ev0[2] = 0.0;
    } else if( N == Nz ) {
      ev0[0] = 0.0;
      ev0[1] = 1.0;
      ev0[2] = 0.0;
    }

    kx = global_shape[2] > 1 ? mode*common::pi2/xlen : 0.0;
    ky = global_shape[1] > 1 ? mode*common::pi2/ylen : 0.0;
    kz = global_shape[0] > 1 ? mode*common::pi2/zlen : 0.0;
    kk = sqrt(kx*kx + ky*ky + kz*kz);

    // parallel to the wavenumber vector
    ev1[0] = kx / kk;
    ev1[1] = ky / kk;
    ev1[2] = kz / kk;

    // perp 1
    ev2[0] = ev0[1]*ev1[2] - ev0[2]*ev1[1];
    ev2[1] = ev0[2]*ev1[0] - ev0[0]*ev1[2];
    ev2[2] = ev0[0]*ev1[1] - ev0[1]*ev1[0];

    // perp 2
    ev3[0] = ev1[1]*ev2[2] - ev1[2]*ev2[1];
    ev3[1] = ev1[2]*ev2[0] - ev1[0]*ev2[2];
    ev3[2] = ev1[0]*ev2[1] - ev1[1]*ev2[0];

    // normalize unit vectors
    {
      float64 en2 = 1/sqrt(ev2[0]*ev2[0] + ev2[1]*ev2[1] + ev2[2]*ev2[2]);
      float64 en3 = 1/sqrt(ev3[0]*ev3[0] + ev3[1]*ev3[1] + ev3[2]*ev3[2]);

      ev2[0] *= en2;
      ev2[1] *= en2;
      ev2[2] *= en2;
      ev3[0] *= en3;
      ev3[1] *= en3;
      ev3[2] *= en3;
    }

    // parameters
    float64 mpe   = config.getAs<float64>("mpe");
    float64 sigma = config.getAs<float64>("sigma");
    float64 h0    = config.getAs<float64>("h0");
    float64 T0    = (gam-1)/gam * (h0 - 1);
    float64 b0    = sqrt(common::pi4 * sigma * h0) * c;
    float64 eta   = config.getAs<float64>("eta");
    float64 w0    = config.getAs<float64>("w0");
    float64 g0    = config.getAs<float64>("g0");
    float64 gp    = g0;
    float64 ge    = g0;

    // proton
    qp  = c * sqrt(h0 / common::pi4);
    mp  = 1.0;
    qmp = qp/mp;

    // electron
    qe  =-qp;
    me  = mp / mpe;
    qme = qe/me;

    // solve dispersion relation
    float64 k0    = kk;
    float64 wpp   = c;
    float64 wcp   = sqrt(sigma) * wpp;
    float64 wce   =-wcp * mpe;

    DispersionSolver disp(wpp, wcp, mpe);
    disp.solve(k0, eta, w0, gp, ge);

    // proper density
    float64 np0 = 1.0 / gp;
    float64 ne0 = 1.0 / ge;

    // transverse velocity
    float64 upt =-wcp/(w0 + wcp/gp) * w0/k0 * eta;
    float64 uet =-wce/(w0 + wce/ge) * w0/k0 * eta;

    // pressure
    float64 pp0 = np0 * T0 * mp*cc;
    float64 pe0 = ne0 * T0 * me*cc;

    // E field
    float64 bp   = eta * b0;
    float64 ep   = bp * w0/(k0*c);

    for(int iz=Lbz-1; iz <= Ubz+1 ;iz++) {
      for(int iy=Lby-1; iy <= Uby+1 ;iy++) {
        for(int ix=Lbx-1; ix <= Ubx+1 ;ix++) {
          float64 phi  = kx*xig(ix) + ky*yig(iy) + kz*zig(iz);
          float64 sinp = sin(phi);
          float64 cosp = cos(phi);
          float64 vp[3], ve[3], ee[3], bb[3];

          vp[0] = 0.0;
          vp[1] =+upt*cosp;
          vp[2] =-upt*sinp;
          ve[0] = 0.0;
          ve[1] =+uet*cosp;
          ve[2] =-uet*sinp;
          ee[0] = 0.0;
          ee[1] =-ep*sinp;
          ee[2] =-ep*cosp;
          bb[0] = b0;
          bb[1] =+bp*cosp;
          bb[2] =-bp*sinp;

          // proton
          uf(iz,iy,ix,0)  = np0;
          uf(iz,iy,ix,1)  = ev1[0]*vp[0] + ev2[0]*vp[1] + ev3[0]*vp[2];
          uf(iz,iy,ix,2)  = ev1[1]*vp[0] + ev2[1]*vp[1] + ev3[1]*vp[2];
          uf(iz,iy,ix,3)  = ev1[2]*vp[0] + ev2[2]*vp[1] + ev3[2]*vp[2];
          uf(iz,iy,ix,4)  = pp0;
          // electron
          uf(iz,iy,ix,5)  = ne0;
          uf(iz,iy,ix,6)  = ev1[0]*ve[0] + ev2[0]*ve[1] + ev3[0]*ve[2];
          uf(iz,iy,ix,7)  = ev1[1]*ve[0] + ev2[1]*ve[1] + ev3[1]*ve[2];
          uf(iz,iy,ix,8)  = ev1[2]*ve[0] + ev2[2]*ve[1] + ev3[2]*ve[2];
          uf(iz,iy,ix,9)  = pe0;
          // e-field
          ebc(iz,iy,ix,0) = ev1[0]*ee[0] + ev2[0]*ee[1] + ev3[0]*ee[2];
          ebc(iz,iy,ix,1) = ev1[1]*ee[0] + ev2[1]*ee[1] + ev3[1]*ee[2];
          ebc(iz,iy,ix,2) = ev1[2]*ee[0] + ev2[2]*ee[1] + ev3[2]*ee[2];
          // b-field
          ebc(iz,iy,ix,3) = ev1[0]*bb[0] + ev2[0]*bb[1] + ev3[0]*bb[2];
          ebc(iz,iy,ix,4) = ev1[1]*bb[0] + ev2[1]*bb[1] + ev3[1]*bb[2];
          ebc(iz,iy,ix,5) = ev1[2]*bb[0] + ev2[2]*bb[1] + ev3[2]*bb[2];
        }
      }
    }

    int_c2f<2>(ebc, ueb);
    boundary->set_field(*this, ueb);
    boundary->set_fluid(*this, uf, ebc);

    // initialize temporary
    vf  = uf;
    veb = ueb;

    // execution time
    period = 2*common::pi/w0;

    // output message
    float64 wp = sqrt(common::pi4 * (qp*qp/(h0*mp) + qe*qe/(h0*me)));
    msg += "# Execute simulation with the following parameters\n";
    msg += tfm::format("%-20s = %20.12e\n", "speed of light", c);
    msg += tfm::format("%-20s = %20.12e\n", "delx", delx);
    msg += tfm::format("%-20s = %20.12e\n", "dely", dely);
    msg += tfm::format("%-20s = %20.12e\n", "delz", delz);
    msg += tfm::format("%-20s = %20.12e\n", "gamma", gam);
    msg += tfm::format("%-20s = %20.12e\n", "mass ratio", mpe);
    msg += tfm::format("%-20s = %20.12e\n", "inertia length", c/wp);
    msg += tfm::format("%-20s = %20.12e\n", "amplitude", eta);
    msg += tfm::format("%-20s = %20.12e\n", "frequency", w0);
    msg += tfm::format("%-20s = %20.12e\n", "wavenumber", k0);
    msg += tfm::format("%-20s = %20.12e\n", "up", sqrt(gp*gp-1));
    msg += tfm::format("%-20s = %20.12e\n", "ue", sqrt(ge*ge-1));

    msg += tfm::format("%-20s = %24.15e\n", "omega", w0);
    msg += tfm::format("%-20s = %24.15e\n", "gamp-1", gp-1);
    msg += tfm::format("%-20s = %24.15e\n", "game-1", ge-1);
  }

  // return maximum time
  float64 getPeriod()
  {
    return period;
  }
};

//
// *** main loop ***
//
int main(int argc, char** argv)
{
  const int Nb = 2;

  int shape[3], domain[3];
  int nw, Nx, Ny, Nz, Nt, cycle;
  float64 cc, gamma, cfl, tmax, len, delh;
  std::string cfg, out;

  // parse command line
  CmdParser parser("circular.cfg");
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
  cycle = config.getAs<int>("cycle");
  cfl   = config.getAs<float64>("cfl");
  len   = config.getAs<float64>("len");
  cc    = config.getAs<float64>("cc");
  gamma = config.getAs<float64>("gamma");

  delh = len / limiter::max(Nx, Ny, Nz);

  shape[0] = Nz;
  shape[1] = Ny;
  shape[2] = Nx;
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
  tmax = cycle * global.getPeriod();
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
