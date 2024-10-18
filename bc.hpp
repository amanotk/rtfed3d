// -*- C++ -*-
#ifndef _BC_HPP_
#define _BC_HPP_

///
/// @brief Boundary Condition
///
/// $Id: bc.hpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
#include "decls.hpp"

///
/// Base class for handling boundary condition
///
class BaseBoundary
{
protected:
  BaseBoundary() {};
  virtual ~BaseBoundary() {}

public:
  virtual void set_field(Global &g, T_vector &eb, int nb=-1) {}
  virtual void set_fluid(Global &g, T_vector &uf, T_vector &eb, int nb=-1) {}
  virtual void set_field_flux(Global &g, T_vector &flux, int nb=-1) {}
  virtual void set_fluid_flux(Global &g, T_tensor &flux, int nb=-1) {}
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
