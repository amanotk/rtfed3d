// -*- C++ -*-
#ifndef _BC_PERIODIC_HPP_
#define _BC_PERIODIC_HPP_

///
/// @brief Periodic Boundary Condition
///
/// $Id: bc_periodic.hpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
#include "global.hpp"
#include "bc.hpp"
#include "mpiutil.hpp"

///
/// Implementation of periodic boundary condition via MPI
///
class PeriodicBoundary : public BaseBoundary
{
protected:
  blitz::TinyVector<int,3> m_lb;
  blitz::TinyVector<int,3> m_ub;
  int      m_bufnum;
  float64 *m_bufsnd[3][2];
  float64 *m_bufrcv[3][2];

  template <class T_array, class T_shape>
  void pack(float64 *buf, int &pos, T_array &array,
            T_shape &Lb, T_shape &Ub, bool is_boundary=false);

  template <class T_array, class T_shape>
  void unpack(float64 *buf, int &pos, T_array &array,
              T_shape &Lb, T_shape &Ub, bool is_boundary=false);

public:
  PeriodicBoundary(const int Nz, const int Ny, const int Nx, const int Nb);

  virtual ~PeriodicBoundary();

  // set boundary condition for EM-field at face center
  virtual void set_field(Global &g, T_vector &eb, int nb=-1);

  // set boundary condition for fluid and EM-field at cell center
  virtual void set_fluid(Global &g, T_vector &uf, T_vector &eb, int nb=-1);
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
