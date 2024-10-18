// -*- C++ -*-
#ifndef _DECLS_HPP_
#define _DECLS_HPP_

///
/// @brief Forward declarations
///
/// $Id: decls.hpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
typedef blitz::Array<float64,1> T_coord;
typedef blitz::Array<float64,3> T_scalar;
typedef blitz::Array<float64,4> T_vector;
typedef blitz::Array<float64,5> T_tensor;

// forward declaration
class Global;
class BaseBoundary;

namespace {

// constant
const float64 rpi4 = 1/common::pi4;

/// Calculate Lorentz factor
float64 lorentz(float64 ux, float64 uy, float64 uz)
{
  return sqrt(1 + ux*ux + uy*uy + uz*uz);
}

void calc_grid_bounds(int N, int Nb, int &M, int &Lb, int &Ub)
{
  M  = N + 2*Nb;
  Lb = Nb;
  Ub = N + Nb - 1;
}

}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
