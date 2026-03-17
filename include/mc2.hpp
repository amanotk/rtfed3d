// -*- C++ -*-
#ifndef _MC2_HPP_
#define _MC2_HPP_

///
/// @brief MC2 Solver
///
/// $Id: mc2.hpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
#include "global.hpp"
#include "solver.hpp"
#include "bc.hpp"

///
/// MC2 Solver
///
class MC2 : public BaseSolver
{
protected:
  static const int Nb = 2;

public:
  // constructor
  MC2(Global &g)
  : BaseSolver(g.Mz, g.My, g.Mx)
  {
  }

  // interpolation with MC2 limiter
  void mc2int(float64 um, float64 uc, float64 up,
              float64 &ul, float64 &ur)
  {
    float64 du = limiter::mc2(uc - um, up - uc);

    ul = uc + 0.5*du;
    ur = uc - 0.5*du;
  }

  // calculate numerical flux for fluid
  virtual void flux_fluid(float64 dt, Global &g,
                          T_vector &uf, T_vector &ebc, T_vector &ueb,
                          T_tensor &ff);

  // calculate numerical flux for field
  virtual void flux_field(float64 dt, Global &g, T_vector &ueb,
                          T_vector &feb);

  // interpolate from cell center to face center in x dir
  void xint_c2f(Global &g, T_vector &uf, T_vector &ebc, T_vector &ueb);

  // interpolate from cell center to face center in y dir
  void yint_c2f(Global &g, T_vector &uf, T_vector &ebc, T_vector &ueb);

  // interpolate from cell center to face center in z dir
  void zint_c2f(Global &g, T_vector &uf, T_vector &ebc, T_vector &ueb);

  // interpolate from face center to edge center in x dir
  void xint_f2e(Global &g, T_vector &eby, T_vector &ebz);

  // interpolate from face center to edge center in y dir
  void yint_f2e(Global &g, T_vector &ebz, T_vector &ebx);

  // interpolate from face center to edge center in z dir
  void zint_f2e(Global &g, T_vector &ebx, T_vector &eby);
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
