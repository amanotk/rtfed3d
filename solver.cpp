// -*- C++ -*-

///
/// @brief Implementation for Base Solver
///
/// $Id: solver.cpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
#include "solver.hpp"

void BaseSolver::xflux_fluid_hll(float64 dt, Global &g, T_tensor &ff, T_vector &ebx,
                                 const int Nb)
{
  const float64 c  = g.c;

  for(int iz=g.Lbz-Nb; iz <= g.Ubz+Nb ;iz++) {
    for(int iy=g.Lby-Nb; iy <= g.Uby+Nb ;iy++) {
      for(int ix=g.Lbx-1; ix <= g.Ubx ;ix++) {
        float64 fxr[10], fxl[10], ucr[10], ucl[10];
        float64 *uur = &ur(iz,iy,ix,0);
        float64 *uul = &ul(iz,iy,ix,0);
        float64 vmax = c;
        float64 vmin = c;
        float64 rv   = 1.0 / (vmax + vmin);

        // HLL flux
        g.xflux(&uur[0], &uur[10], fxr);
        g.xflux(&uul[0], &uul[10], fxl);
        g.conservative(&uur[0], &uur[10], ucr);
        g.conservative(&uul[0], &uul[10], ucl);

        for(int k=0; k < 10 ;k++) {
          ff(iz,iy,ix,0,k) = (vmax*fxl[k] + vmin*fxr[k] -
                              vmax*vmin*(ucr[k] - ucl[k])) * rv * dt;
        }

        // HLL average
        for(int k=0; k < 6 ;k++) {
          ebx(iz,iy,ix,k) = (vmax*uul[10+k] + vmin*uur[10+k]) * rv;
        }
      }
    }
  }
}

void BaseSolver::yflux_fluid_hll(float64 dt, Global &g, T_tensor &ff, T_vector &eby,
                                 const int Nb)
{
  const float64 c  = g.c;

  for(int iz=g.Lbz-Nb; iz <= g.Ubz+Nb ;iz++) {
    for(int iy=g.Lby-1; iy <= g.Uby ;iy++) {
      for(int ix=g.Lbx-Nb; ix <= g.Ubx+Nb ;ix++) {
        float64 fyr[10], fyl[10], ucr[10], ucl[10];
        float64 *uur = &ur(iz,iy,ix,0);
        float64 *uul = &ul(iz,iy,ix,0);
        float64 vmax = c;
        float64 vmin = c;
        float64 rv   = 1.0 / (vmax + vmin);

        // HLL flux
        g.yflux(&uur[0], &uur[10], fyr);
        g.yflux(&uul[0], &uul[10], fyl);
        g.conservative(&uur[0], &uur[10], ucr);
        g.conservative(&uul[0], &uul[10], ucl);

        for(int k=0; k < 10 ;k++) {
          ff(iz,iy,ix,1,k) = (vmax*fyl[k] + vmin*fyr[k] -
                              vmax*vmin*(ucr[k] - ucl[k])) * rv * dt;
        }

        // HLL average
        for(int k=0; k < 6 ;k++) {
          eby(iz,iy,ix,k) = (vmax*uul[10+k] + vmin*uur[10+k]) * rv;
        }
      }
    }
  }
}

void BaseSolver::zflux_fluid_hll(float64 dt, Global &g, T_tensor &ff, T_vector &ebz,
                                 const int Nb)
{
  const float64 c  = g.c;

  for(int iz=g.Lbz-1; iz <= g.Ubz ;iz++) {
    for(int iy=g.Lby-Nb; iy <= g.Uby+Nb ;iy++) {
      for(int ix=g.Lbx-Nb; ix <= g.Ubx+Nb ;ix++) {
        float64 fzr[10], fzl[10], ucr[10], ucl[10];
        float64 *uur = &ur(iz,iy,ix,0);
        float64 *uul = &ul(iz,iy,ix,0);
        float64 vmax = c;
        float64 vmin = c;
        float64 rv   = 1.0 / (vmax + vmin);

        // HLL flux
        g.zflux(&uur[0], &uur[10], fzr);
        g.zflux(&uul[0], &uul[10], fzl);
        g.conservative(&uur[0], &uur[10], ucr);
        g.conservative(&uul[0], &uul[10], ucl);

        for(int k=0; k < 10 ;k++) {
          ff(iz,iy,ix,2,k) = (vmax*fzl[k] + vmin*fzr[k] -
                              vmax*vmin*(ucr[k] - ucl[k])) * rv * dt;
        }

        // HLL average
        for(int k=0; k < 6 ;k++) {
          ebz(iz,iy,ix,k) = (vmax*uul[10+k] + vmin*uur[10+k]) * rv;
        }
      }
    }
  }
}

void BaseSolver::xflux_field_hll(float64 dt, Global &g, T_vector &feb)
{
  const float64 c  = g.c;
  const float64 rc = g.rc;

  for(int iz=g.Lbz-1; iz <= g.Ubz ;iz++) {
    for(int iy=g.Lby-1; iy <= g.Uby ;iy++) {
      for(int ix=g.Lbx-1; ix <= g.Ubx ;ix++) {
        float64 *ebr = &ur(iz,iy,ix,0);
        float64 *ebl = &ul(iz,iy,ix,0);

        // y-face : Ez, Bz
        {
          float64 vmax = c;
          float64 vmin = c;
          float64 rv   = 1.0 / (vmax + vmin);

          feb(iz,iy,ix,2) += (0.5*(vmax*ebl[1] + vmin*ebr[1]) +
                              rc * vmax*vmin*(ebr[2] - ebl[2])) * rv * dt;
          feb(iz,iy,ix,5) += (0.5*(vmax*ebl[3] + vmin*ebr[3]) -
                              rc * vmax*vmin*(ebr[0] - ebl[0])) * rv * dt;
        }

        // z-face : Ey, By
        {
          float64 vmax = c;
          float64 vmin = c;
          float64 rv   = 1.0 / (vmax + vmin);

          feb(iz,iy,ix,1) += (0.5*(vmax*ebl[4] + vmin*ebr[4]) -
                              rc * vmax*vmin*(ebr[7] - ebl[7])) * rv * dt;
          feb(iz,iy,ix,4) += (0.5*(vmax*ebl[6] + vmin*ebr[6]) +
                              rc * vmax*vmin*(ebr[5] - ebl[5])) * rv * dt;
        }
      }
    }
  }
}

void BaseSolver::yflux_field_hll(float64 dt, Global &g, T_vector &feb)
{
  const float64 c  = g.c;
  const float64 rc = g.rc;

  for(int iz=g.Lbz-1; iz <= g.Ubz ;iz++) {
    for(int iy=g.Lby-1; iy <= g.Uby ;iy++) {
      for(int ix=g.Lbx-1; ix <= g.Ubx ;ix++) {
        float64 *ebr = &ur(iz,iy,ix,0);
        float64 *ebl = &ul(iz,iy,ix,0);

        // z-face : Ex, Bx
        {
          float64 vmax = c;
          float64 vmin = c;
          float64 rv   = 1.0 / (vmax + vmin);

          feb(iz,iy,ix,0) += (0.5*(vmax*ebl[0] + vmin*ebr[0]) +
                              rc * vmax*vmin*(ebr[3] - ebl[3])) * rv * dt;
          feb(iz,iy,ix,3) += (0.5*(vmax*ebl[2] + vmin*ebr[2]) -
                              rc * vmax*vmin*(ebr[1] - ebl[1])) * rv * dt;
        }

        // x-face : Ez, Bz
        {
          float64 vmax = c;
          float64 vmin = c;
          float64 rv   = 1.0 / (vmax + vmin);

          feb(iz,iy,ix,2) += (0.5*(vmax*ebl[5] + vmin*ebr[5]) -
                              rc * vmax*vmin*(ebr[6] - ebl[6])) * rv * dt;
          feb(iz,iy,ix,5) += (0.5*(vmax*ebl[7] + vmin*ebr[7]) +
                              rc * vmax*vmin*(ebr[4] - ebl[4])) * rv * dt;
        }
      }
    }
  }
}

void BaseSolver::zflux_field_hll(float64 dt, Global &g, T_vector &feb)
{
  const float64 c  = g.c;
  const float64 rc = g.rc;

  for(int iz=g.Lbz-1; iz <= g.Ubz ;iz++) {
    for(int iy=g.Lby-1; iy <= g.Uby ;iy++) {
      for(int ix=g.Lbx-1; ix <= g.Ubx ;ix++) {
        float64 *ebr = &ur(iz,iy,ix,0);
        float64 *ebl = &ul(iz,iy,ix,0);

        // x-face : Ey, By
        {
          float64 vmax = c;
          float64 vmin = c;
          float64 rv   = 1.0 / (vmax + vmin);

          feb(iz,iy,ix,1) += (0.5*(vmax*ebl[1] + vmin*ebr[1]) +
                              rc * vmax*vmin*(ebr[2] - ebl[2])) * rv * dt;
          feb(iz,iy,ix,4) += (0.5*(vmax*ebl[3] + vmin*ebr[3]) -
                              rc * vmax*vmin*(ebr[0] - ebl[0])) * rv * dt;
        }

        // y-face : Ex, Bx
        {
          float64 vmax = c;
          float64 vmin = c;
          float64 rv   = 1.0 / (vmax + vmin);

          feb(iz,iy,ix,0) += (0.5*(vmax*ebl[4] + vmin*ebr[4]) -
                              rc * vmax*vmin*(ebr[7] - ebl[7])) * rv * dt;
          feb(iz,iy,ix,3) += (0.5*(vmax*ebl[6] + vmin*ebr[6]) +
                              rc * vmax*vmin*(ebr[5] - ebl[5])) * rv * dt;
        }
      }
    }
  }
}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
