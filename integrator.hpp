// -*- C++ -*-
#ifndef _INTEGRATE_HPP_
#define _INTEGRATE_HPP_

///
/// @brief Time Integration Schemes
///
/// $Id$
///
#include "global.hpp"
#include "solver.hpp"

///
/// Base Integrator
///
class BaseIntegrator
{
public:
  virtual void push(float64 dt, Global &g, BaseSolver &solver) = 0;
};

///
/// 3rd order TVD Runge-Kutta integrator
///
class RK3 : public BaseIntegrator
{
private:

  // substep
  void push_substep(float64    dt,
                    Global     &g,
                    float64    coeff[4],
                    T_vector   &uf,
                    T_vector   &vf,
                    T_tensor   &ff,
                    T_vector   &ebc,
                    T_vector   &ueb,
                    T_vector   &veb,
                    T_vector   &feb,
                    BaseSolver &solver);

public:
  RK3(Global &g)
  {
  }

  virtual void push(float64 dt, Global &g, BaseSolver &solver);
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
