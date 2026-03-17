// -*- C++ -*-

///
/// @brief Implementation of time integration schemes
///
/// $Id: integrator.cpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
#include "integrator.hpp"

void RK3::push(float64 dt, Global &g, BaseSolver &solver)
{
  const float64 tvd_rk3[3][4] =
    { 1.0,         0.0,    1.0,      0.0,
      3.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0,
      1.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0};

  // convert to conservative variables
  for(int iz=g.Lbz; iz <= g.Ubz ;iz++) {
    for(int iy=g.Lby; iy <= g.Uby ;iy++) {
      for(int ix=g.Lbx; ix <= g.Ubx ;ix++) {
        float64 uc[10];
        // copy
        for(int k=0; k < 10 ;k++) {
          uc[k] = g.uf(iz,iy,ix,k);
        }
        // convert to conservative
        g.conservative(uc, &g.ebc(iz,iy,ix,0), &g.uf(iz,iy,ix,0));
      }
    }
  }

  // *** 1st step ***
  push_substep(dt, g, const_cast<float64*>(tvd_rk3[0]),
               g.uf, g.vf, g.ff, g.ebc, g.ueb, g.veb, g.feb, solver);

  // *** 2nd step ***
  push_substep(dt, g, const_cast<float64*>(tvd_rk3[1]),
               g.uf, g.vf, g.ff, g.ebc, g.ueb, g.veb, g.feb, solver);

  // *** 3rd step ***
  push_substep(dt, g, const_cast<float64*>(tvd_rk3[2]),
               g.uf, g.vf, g.ff, g.ebc, g.ueb, g.veb, g.feb, solver);


  // update solution
  g.uf  = g.vf;
  g.ueb = g.veb;
}

void RK3::push_substep(float64    dt,
                       Global     &g,
                       float64    coeff[4],
                       T_vector   &uf,
                       T_vector   &vf,
                       T_tensor   &ff,
                       T_vector   &ebc,
                       T_vector   &ueb,
                       T_vector   &veb,
                       T_vector   &feb,
                       BaseSolver &solver)
{
  const float64 dtx = 1 /g.delx;
  const float64 dty = 1 /g.dely;
  const float64 dtz = 1 /g.delz;
  const float64 cdtx = g.c * dtx;
  const float64 cdty = g.c * dty;
  const float64 cdtz = g.c * dtz;
  const float64 pidt = common::pi4;

  // calculate fluxes
  solver.flux_fluid(dt, g, vf, ebc, veb, ff);
  solver.flux_field(dt, g, veb, feb);

  //
  // update fluid
  //
  for(int iz=g.Lbz; iz <= g.Ubz ;iz++) {
    for(int iy=g.Lby; iy <= g.Uby ;iy++) {
      for(int ix=g.Lbx; ix <= g.Ubx ;ix++) {
        float64 vc[10], rhs[10];
        // conservative variables
        g.conservative(&vf(iz,iy,ix,0), &ebc(iz,iy,ix,0), vc);
        // RHS
        g.rhs(dt, &vf(iz,iy,ix,0), &ebc(iz,iy,ix,0), rhs);
        // update conservative variables
        for(int k=0; k < 10 ;k++) {
          vc[k] =
            + coeff[0] * uf(iz,iy,ix,k)
            + coeff[1] * vc[k]
            + coeff[2] *
            ( ( dtx*(ff(iz,iy,ix-1,0,k) - ff(iz,iy,ix,0,k)) +
                dty*(ff(iz,iy-1,ix,1,k) - ff(iz,iy,ix,1,k)) +
                dtz*(ff(iz-1,iy,ix,2,k) - ff(iz,iy,ix,2,k)) )
              + rhs[k] );
        }
        // copy
        for(int k=0; k < 10 ;k++) {
          vf(iz,iy,ix,k) = vc[k];
        }
      }
    }
  }

  //
  // update EM field
  //
  for(int iz=g.Lbz-1; iz <= g.Ubz ;iz++) {
    for(int iy=g.Lby-1; iy <= g.Uby ;iy++) {
      for(int ix=g.Lbx-1; ix <= g.Ubx ;ix++) {
        // Ex
        veb(iz,iy,ix,0) =
          + coeff[0] * ueb(iz,iy,ix,0)
          + coeff[1] * veb(iz,iy,ix,0)
          + coeff[2] *
          ( + cdty*(feb(iz,iy,ix,5) - feb(iz,iy-1,ix,5))
            - cdtz*(feb(iz,iy,ix,4) - feb(iz-1,iy,ix,4))
            - pidt * ff(iz,iy,ix,0,5) );
        // Ey
        veb(iz,iy,ix,1) =
          + coeff[0] * ueb(iz,iy,ix,1)
          + coeff[1] * veb(iz,iy,ix,1)
          + coeff[2] *
          ( + cdtz*(feb(iz,iy,ix,3) - feb(iz-1,iy,ix,3))
            - cdtx*(feb(iz,iy,ix,5) - feb(iz,iy,ix-1,5))
            - pidt * ff(iz,iy,ix,1,5) );
        // Ez
        veb(iz,iy,ix,2) =
          + coeff[0] * ueb(iz,iy,ix,2)
          + coeff[1] * veb(iz,iy,ix,2)
          + coeff[2] *
          ( + cdtx*(feb(iz,iy,ix,4) - feb(iz,iy,ix-1,4))
            - cdty*(feb(iz,iy,ix,3) - feb(iz,iy-1,ix,3))
            - pidt * ff(iz,iy,ix,2,5) );
        // Bx
        veb(iz,iy,ix,3) =
          + coeff[0] * ueb(iz,iy,ix,3)
          + coeff[1] * veb(iz,iy,ix,3)
          + coeff[2] *
          ( - cdty*(feb(iz,iy,ix,2) - feb(iz,iy-1,ix,2))
            + cdtz*(feb(iz,iy,ix,1) - feb(iz-1,iy,ix,1)) );
        // By
        veb(iz,iy,ix,4) =
          + coeff[0] * ueb(iz,iy,ix,4)
          + coeff[1] * veb(iz,iy,ix,4)
          + coeff[2] *
          ( - cdtz*(feb(iz,iy,ix,0) - feb(iz-1,iy,ix,0))
            + cdtx*(feb(iz,iy,ix,2) - feb(iz,iy,ix-1,2)) );
        // Bz
        veb(iz,iy,ix,5) =
          + coeff[0] * ueb(iz,iy,ix,5)
          + coeff[1] * veb(iz,iy,ix,5)
          + coeff[2] *
          ( - cdtx*(feb(iz,iy,ix,1) - feb(iz,iy,ix-1,1))
            + cdty*(feb(iz,iy,ix,0) - feb(iz,iy-1,ix,0)) );
      }
    }
  }

  g.boundary->set_field(g, veb);

  //
  // recover primitive variables
  //
  g.int_f2c<2>(ebc, veb);

  for(int iz=g.Lbz; iz <= g.Ubz ;iz++) {
    for(int iy=g.Lby; iy <= g.Uby ;iy++) {
      for(int ix=g.Lbx; ix <= g.Ubx ;ix++) {
        float64 vc[10];
        // copy
        for(int k=0; k < 10 ;k++) {
          vc[k] = vf(iz,iy,ix,k);
        }
        // convert to primitive
        g.primitive(vc, &ebc(iz,iy,ix,0), &vf(iz,iy,ix,0));
      }
    }
  }

  g.boundary->set_fluid(g, vf, ebc);
}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
