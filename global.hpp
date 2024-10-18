// -*- C++ -*-
#ifndef _GLOBAL_HPP_
#define _GLOBAL_HPP_

///
/// @file global.hpp
/// @brief Definition of global variables
///
/// $Id: global.hpp,v 641a96a86755 2015/09/14 19:35:09 amano $
///
#include <blitz/array.h>
#include "common.hpp"
#include "limiter.hpp"
#include "mpiutil.hpp"
#include "primitive.hpp"
#include "MersenneTwister.hpp"
#include "debug.hpp"

#include "decls.hpp"

/// global variables
class Global
{
protected:
  /// constructor
  Global(const int shape[3], const int nb,
         const float64 lightspeed, const float64 gamma)
    : boundary(0)
  {
    Nb  = nb;
    Nz  = shape[0];
    Ny  = shape[1];
    Nx  = shape[2];
    calc_grid_bounds(Nz, Nb, Mz, Lbz, Ubz);
    calc_grid_bounds(Ny, Nb, My, Lby, Uby);
    calc_grid_bounds(Nx, Nb, Mx, Lbx, Ubx);

    xig.resize(Mx);
    yig.resize(My);
    zig.resize(Mz);

    uf.resize(Mz, My, Mx, 10);
    vf.resize(Mz, My, Mx, 10);
    ff.resize(Mz, My, Mx, 3, 10);
    ebc.resize(Mz, My, Mx, 6);
    ueb.resize(Mz, My, Mx, 6);
    veb.resize(Mz, My, Mx, 6);
    feb.resize(Mz, My, Mx, 6);

    c   = lightspeed;
    cc  = c*c;
    rc  = 1/c;
    rcc = 1/cc;
    gam = gamma;
    G   = gam/(gam-1);
    eta = 0.0;
    msg = "";
    diverr.resize(Mz, My, Mx, 2);
  }

  // for internal use
  void write_parameter(std::ofstream &ofs);
  template <class T_array>
  void write_debug_array(std::ofstream &ofs, T_array &x);

public:
  // utility
  static common::DebugStream debugstream;
  std::string msg;

  //
  // boundary condition
  //
  BaseBoundary *boundary;

  //
  // array sizes
  //
  int Nb;                        ///< # grid points for boundary margin
  int Nx;                        ///< # grid points in x
  int Ny;                        ///< # grid points in y
  int Nz;                        ///< # grid points in z
  int Mx;                        ///< # array length in x
  int My;                        ///< # array length in y
  int Mz;                        ///< # array length in z
  int Lbx;                       ///< lower bound of array in x
  int Ubx;                       ///< upper bound of array in x
  int Lby;                       ///< lower bound of array in y
  int Uby;                       ///< upper bound of array in y
  int Lbz;                       ///< lower bound of array in z
  int Ubz;                       ///< upper bound of array in z

  //
  // simulation parameters
  //
  float64  c;                    ///< speed of light
  float64  cc;                   ///< c^2
  float64  rc;                   ///< 1 / c
  float64  rcc;                  ///< 1 / c^2
  float64  gam;                  ///< polytropic index
  float64  G;                    ///< polytropic index
  float64  eta;                  ///< resistivity
  float64  physt;                ///< physical time
  float64  delt;                 ///< time step
  float64  delx;                 ///< grid size in x
  float64  dely;                 ///< grid size in y
  float64  delz;                 ///< grid size in z
  float64  xrange[3];            ///< spatial domain in x dir.
  float64  yrange[3];            ///< spatial domain in y dir.
  float64  zrange[3];            ///< spatial domain in z dir.

  //
  // proton and electron quantities
  //
  float64  qe;                   ///< electron charge
  float64  qp;                   ///< proton charge
  float64  me;                   ///< electron mass
  float64  mp;                   ///< proton mass
  float64  qme;                  ///< charge / mass of electron
  float64  qmp;                  ///< charge / mass of proton

  //
  // array definition
  //
  T_coord xig;                   ///< integer grid points in x
  T_coord yig;                   ///< integer grid points in x
  T_coord zig;                   ///< integer grid points in x
  T_vector uf;                   ///< fluid and EM-field quantities at cell center
  T_vector vf;                   ///< temporary
  T_tensor ff;                   ///< numerical flux for fluid
  T_vector ebc;                  ///< EM-field at cell center
  T_vector ueb;                  ///< EM-field at face center
  T_vector veb;                  ///< temporary
  T_vector feb;                  ///< numerical flux for EM-field
  T_vector diverr;               ///< divergence error (for debugging purpose)

  //
  // member functions
  //
  /// show message to the stream
  void print_message(std::ostream &ofs);

  /// calculate energies
  void energy(float64 &Ef, float64 &Ep);

  /// interpolate EM-field fromc cell center to face center
  template <int R> void int_c2f(T_vector &ebc, T_vector &ueb);

  /// interpolate EM-field fromc face center to cell center
  template <int R> void int_f2c(T_vector &ebc, T_vector &ueb);

  /// convert conservative to primitive variables
  void primitive(float64 uc[10], float64 eb[6], float64 uf[10]);

  /// convert primitive to conservative variables
  void conservative(float64 uf[10], float64 eb[6], float64 uc[10]);

  /// calculate right-hand-side for the fluid equations
  void rhs(float64 dt, float64 uf[10], float64 eb[6], float64 rhs[10]);

  /// calculate flux in x direction for fluid equations
  void xflux(float64 uf[10], float64 eb[6], float64 fx[10]);

  /// calculate flux in y direction for fluid equations
  void yflux(float64 uf[10], float64 eb[6], float64 fy[10]);

  /// calculate flux in z direction for fluid equations
  void zflux(float64 uf[10], float64 eb[6], float64 fz[10]);

  // TODO: to be reorganized
  float64 get_dt_factor(float64 cfl);
  void write_data(std::ofstream &ofs, bool init=false);
  void write_debug_data(std::ofstream &ofs, bool init=false);
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
