// -*- C++ -*-
#ifndef _SOLVER_HPP_
#define _SOLVER_HPP_

///
/// @breif Implementation of Solvers
///
/// $Id: solver.hpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
#include "global.hpp"

///
/// Base Solver
///
class BaseSolver
{
protected:
  static const int Nc = 16;
  T_vector ur;
  T_vector ul;
  T_vector ebx;
  T_vector eby;
  T_vector ebz;

public:
  // constructor
  BaseSolver(const int Mz, const int My, const int Mx)
  {
    ur.resize(Mz, My, Mx, Nc);
    ul.resize(Mz, My, Mx, Nc);
    ebx.resize(Mz, My, Mx, 6);
    eby.resize(Mz, My, Mx, 6);
    ebz.resize(Mz, My, Mx, 6);
  }

  // destructor
  virtual ~BaseSolver()
  {
  }

  // calculate numerical flux for fluid
  virtual void flux_fluid(float64 dt, Global &g,
                          T_vector &uf, T_vector &ebc, T_vector &ueb,
                          T_tensor &ff) = 0;

  // calculate numerical flux for field
  virtual void flux_field(float64 dt, Global &g, T_vector &ueb,
                          T_vector &feb) = 0;


  void xflux_fluid_hll(float64 dt, Global &g, T_tensor &ff, T_vector &ebx,
                       const int Nb);

  void yflux_fluid_hll(float64 dt, Global &g, T_tensor &ff, T_vector &eby,
                       const int Nb);

  void zflux_fluid_hll(float64 dt, Global &g, T_tensor &ff, T_vector &ebz,
                       const int Nb);

  void xflux_field_hll(float64 dt, Global &g, T_vector &feb);

  void yflux_field_hll(float64 dt, Global &g, T_vector &feb);

  void zflux_field_hll(float64 dt, Global &g, T_vector &feb);
};

//
// include definition of specific solvers
//
#include "mc2.hpp"

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
