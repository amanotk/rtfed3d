// -*- C++ -*-
#ifndef _PRIMITIVE_HPP_
#define _PRIMITIVE_HPP_

///
/// @brief Primitive Solver for Relativistic Hydrodynamics
///
/// $Id: primitive.hpp,v 6f2924be4891 2014/09/01 02:53:29 amano $
///

namespace
{

///
/// @brief Convert conservative to primitive variables
///
/// The Brown method with specialization to relativistic hydrodynamic
/// equations is used to obtain primitive variables.
///
/// * Known Problems
/// - cancellation of significant digits occurs at gamma >> 1
///
/// * References
/// - Zenitani et al., 2009, Astrophys. J., 696, 1385
/// - Nunohiro & Hirano, 2003, Trans. Japan Soc. Ind. Appl. Math., 13, 159
/// - Nunohiro et al., 1996, Trans. Japan Soc. Ind. Appl. Math., 6, 173
///
/// * Parameters
/// - uc : conservative variables
/// - up : primitive variables to be recovered
/// - cc : speed of light
/// - G  : Gamma/(Gamma-1)
///
void rhd_primitive(float64 uc[], float64 up[],
                   float64 mm, float64 cc, float64 G)
{
  static const float64 cubic_re = 0.5;
  static const float64 cubic_im = 0.5*sqrt(3.0);

  float64 U = sqrt(uc[1]*uc[1] + uc[2]*uc[2] + uc[3]*uc[3]);
  float64 Y = U*cc/uc[4];
  float64 Z = uc[0]*cc*cc/uc[4];
  float64 denom = G*G*(1+Y)*(1-Y);
  float64 a =-(2*G*Y*Z) / denom;
  float64 b =+(G*G - 2*G*(G-1)*Y*Y - Z*Z) / denom;
  float64 c =-(2*(G-1)*Y*Z) / denom;
  float64 d =-((G-1)*(G-1)*Y*Y) / denom;
  float64 root_cubic, root_quartic;

  //
  // the equation to be solved is the following quartic equation
  //   x^4 + a*x^3 + b*x^2 + c*x + d = 0
  // where a < 0, c < 0, d < 0
  //

  //
  // solve a reduced cubic equation (according to the Brown method)
  // careful examination of sign of coefficients shows that this equation has
  // one real and two complex conjudate roots.
  //
  {
    float64 xreal[3], ximag[3], xabs[3];
    float64 aa =-b;
    float64 bb =+a*c - 4*d;
    float64 cc =+d*(4*b-a*a) - c*c;
    float64 dd = aa/3.0;
    float64 ee = dd*dd;
    float64 pp = bb/3.0 - ee;
    float64 qq = ee*dd - 0.5*(dd*bb - cc);
    float64 DD = common::maximum(pp*pp*pp + qq*qq, 0.0);
    float64 ww = cbrt(-qq - copysign(1.0, qq)*sqrt(DD)); // DD>0
    float64 uu = ww - pp/ww;
    float64 vv = ww + pp/ww;

    xreal[0] =+uu - dd;
    ximag[0] =+0.0;
    xreal[1] =-cubic_re*uu - dd;
    ximag[1] =+cubic_im*vv;
    xreal[2] =-cubic_re*uu - dd;
    ximag[2] =-cubic_im*vv;

    //
    // correct the obtained real root if needed
    // (no correction needed to complex roots for our purpose)
    //
    xabs[0] = xreal[0]*xreal[0] + ximag[0]*ximag[0];
    xabs[1] = xreal[1]*xreal[1] + ximag[1]*ximag[1];
    if( xabs[0] < xabs[1] ) {
      xreal[0] =-cc/xabs[1];
    }
    root_cubic = xreal[0];
  }

  //
  // obtain roots by using a root of the above cubic equation.
  //
  {
    //
    // v   = a root of the reduced cubic equation
    // w^2 = v^2 - 4*d
    // u^2 = a^2 - 4*(v-b)
    // u*w = a*v - 2*c
    //
    // Note: v ~ b for gamma << 1 (non-relativistic velocities)
    //
    float64 vv = root_cubic;        // real root of cubic equation
    float64 ww = sqrt(vv*vv - 4*d); // d < 0
    float64 uu = (a*vv - 2*c)/ww;   // u = sqrt(a*a + 4*(v-b)) is dangerous

    //
    // p = (a + u)/2
    // q = (v + w)/2
    // r = (a - u)/2
    // s = (v - w)/2
    //
    // a = p + r
    // b = p*r + q + s
    // c = p*s + q*r
    // d = q*s
    //
    // it can be shown that p > 0, q > 0, r < 0, s < 0
    //
    float64 p, q, r, s;

    if( vv > 0 ) {
      q = 0.5*(vv + ww);
      s = d/q;
    } else {
      s = 0.5*(vv - ww);
      q = d/s;
    }

    r = 0.5*(a - uu);
    p = (c - q*r)/s;

    //
    // roots of the the quartic equation are given by the roots of the
    // following  two quadratic equations:
    //   x^2 + p*x + q = 0
    //   x^2 + r*x + s = 0
    // but the real roots are obtained by the latter. The physical solution is
    // a positive one.
    //
    root_quartic =-0.5*(r - sqrt(r*r -4*s));
  }

  // recover primitive variables
  float64 u = root_quartic;
  float64 g = sqrt(1 + u*u);
  U = 1.0/(U + 1.0e-20);
  up[0] = uc[0] / g;
  up[1] = u * uc[1]*U;
  up[2] = u * uc[2]*U;
  up[3] = u * uc[3]*U;
  up[4] = up[0] * (g*(1 - g*Z))/((g*g*G-1)*Z);

  // mass density to number density
  up[0] /= mm;
}

}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
