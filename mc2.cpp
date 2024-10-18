// -*- C++ -*-

///
/// @brief Implementation of MC2 solver
///
/// $Id: mc2.cpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
#include "mc2.hpp"

void MC2::flux_fluid(float64 dt, Global &g,
                     T_vector &uf, T_vector &ebc, T_vector &ueb, T_tensor &ff)
{
  xint_c2f(g, uf, ebc, ueb);
  xflux_fluid_hll(dt, g, ff, ebx, Nb);

  yint_c2f(g, uf, ebc, ueb);
  yflux_fluid_hll(dt, g, ff, eby, Nb);

  zint_c2f(g, uf, ebc, ueb);
  zflux_fluid_hll(dt, g, ff, ebz, Nb);

  g.boundary->set_fluid_flux(g, ff);
}

void MC2::flux_field(float64 dt, Global &g, T_vector &ueb, T_vector &feb)
{
  // clear flux
  feb = 0.0;

  xint_f2e(g, eby, ebz);
  xflux_field_hll(dt, g, feb);

  yint_f2e(g, ebz, ebx);
  yflux_field_hll(dt, g, feb);

  zint_f2e(g, ebx, eby);
  zflux_field_hll(dt, g, feb);

  g.boundary->set_field_flux(g, ueb);
}

void MC2::xint_c2f(Global &g, T_vector &uf, T_vector &ebc, T_vector &ueb)
{
  for(int iz=g.Lbz-Nb; iz <= g.Ubz+Nb ;iz++) {
    for(int iy=g.Lby-Nb; iy <= g.Uby+Nb ;iy++) {
      for(int ix=g.Lbx-1; ix <= g.Ubx+1 ;ix++) {
        float64 *um = &uf(iz,iy,ix-1,0);
        float64 *uc = &uf(iz,iy,ix  ,0);
        float64 *up = &uf(iz,iy,ix+1,0);
        float64 *vm = &ebc(iz,iy,ix-1,0);
        float64 *vc = &ebc(iz,iy,ix  ,0);
        float64 *vp = &ebc(iz,iy,ix+1,0);
        float64 *wr = &ueb(iz,iy,ix-1,0);
        float64 *wl = &ueb(iz,iy,ix  ,0);
        float64 *fr = &ur(iz,iy,ix-1,0);
        float64 *fl = &ul(iz,iy,ix  ,0);

        // fluid
        for(int k=0; k < 10 ;k++) {
          mc2int(um[k], uc[k], up[k], fl[k], fr[k]);
        }

        // EM-field
        {
          // Ex
          fl[10] = wl[0];
          fr[10] = wr[0];

          // Ey
          mc2int(vm[1], vc[1], vp[1], fl[11], fr[11]);

          // Ez
          mc2int(vm[2], vc[2], vp[2], fl[12], fr[12]);

          // Bx
          fl[13] = wl[3];
          fr[13] = wr[3];

          // By
          mc2int(vm[4], vc[4], vp[4], fl[14], fr[14]);

          // Bz
          mc2int(vm[5], vc[5], vp[5], fl[15], fr[15]);
        }
      }
    }
  }
}

void MC2::yint_c2f(Global &g, T_vector &uf, T_vector &ebc, T_vector &ueb)
{
  for(int iz=g.Lbz-Nb; iz <= g.Ubz+Nb ;iz++) {
    for(int iy=g.Lby-1; iy <= g.Uby+1 ;iy++) {
      for(int ix=g.Lbx-Nb; ix <= g.Ubx+Nb ;ix++) {
        float64 *um = &uf(iz,iy-1,ix,0);
        float64 *uc = &uf(iz,iy  ,ix,0);
        float64 *up = &uf(iz,iy+1,ix,0);
        float64 *vm = &ebc(iz,iy-1,ix,0);
        float64 *vc = &ebc(iz,iy  ,ix,0);
        float64 *vp = &ebc(iz,iy+1,ix,0);
        float64 *wr = &ueb(iz,iy-1,ix,0);
        float64 *wl = &ueb(iz,iy  ,ix,0);
        float64 *fr = &ur(iz,iy-1,ix,0);
        float64 *fl = &ul(iz,iy  ,ix,0);

        // fluid
        for(int k=0; k < 10 ;k++) {
          mc2int(um[k], uc[k], up[k], fl[k], fr[k]);
        }

        // EM field
        {
          // Ex
          mc2int(vm[0], vc[0], vp[0], fl[10], fr[10]);

          // Ey
          fl[11] = wl[1];
          fr[11] = wr[1];

          // Ez
          mc2int(vm[2], vc[2], vp[2], fl[12], fr[12]);

          // Bx
          mc2int(vm[3], vc[3], vp[3], fl[13], fr[13]);

          // By
          fl[14] = wl[4];
          fr[14] = wr[4];

          // Bz
          mc2int(vm[5], vc[5], vp[5], fl[15], fr[15]);
        }
      }
    }
  }
}

void MC2::zint_c2f(Global &g, T_vector &uf, T_vector &ebc, T_vector &ueb)
{
  for(int iz=g.Lbz-1; iz <= g.Ubz+1 ;iz++) {
    for(int iy=g.Lby-Nb; iy <= g.Uby+Nb ;iy++) {
      for(int ix=g.Lbx-Nb; ix <= g.Ubx+Nb ;ix++) {
        float64 *um = &uf(iz-1,iy,ix,0);
        float64 *uc = &uf(iz  ,iy,ix,0);
        float64 *up = &uf(iz+1,iy,ix,0);
        float64 *vm = &ebc(iz-1,iy,ix,0);
        float64 *vc = &ebc(iz  ,iy,ix,0);
        float64 *vp = &ebc(iz+1,iy,ix,0);
        float64 *wr = &ueb(iz-1,iy,ix,0);
        float64 *wl = &ueb(iz  ,iy,ix,0);
        float64 *fr = &ur(iz-1,iy,ix,0);
        float64 *fl = &ul(iz  ,iy,ix,0);

        // fluid
        for(int k=0; k < 10 ;k++) {
          mc2int(um[k], uc[k], up[k], fl[k], fr[k]);
        }

        // EM field
        {
          // Ex
          mc2int(vm[0], vc[0], vp[0], fl[10], fr[10]);

          // Ey
          mc2int(vm[1], vc[1], vp[1], fl[11], fr[11]);

          // Ez
          fl[12] = wl[2];
          fr[12] = wr[2];

          // Bx
          mc2int(vm[3], vc[3], vp[3], fl[13], fr[13]);

          // By
          mc2int(vm[4], vc[4], vp[4], fl[14], fr[14]);

          // Bz
          fl[15] = wl[5];
          fr[15] = wr[5];
        }
      }
    }
  }
}

void MC2::xint_f2e(Global &g, T_vector &eby, T_vector &ebz)
{
  for(int iz=g.Lbz-1; iz <= g.Ubz+1 ;iz++) {
    for(int iy=g.Lby-1; iy <= g.Uby+1 ;iy++) {
      for(int ix=g.Lbx-1; ix <= g.Ubx+1 ;ix++) {
        float64 *fr = &ur(iz,iy,ix-1,0);
        float64 *fl = &ul(iz,iy,ix  ,0);

        // y-face
        {
          float64 *vm = &eby(iz,iy,ix-1,0);
          float64 *vc = &eby(iz,iy,ix  ,0);
          float64 *vp = &eby(iz,iy,ix+1,0);

          // Ey
          mc2int(vm[1], vc[1], vp[1], fl[0], fr[0]);

          // Ez
          mc2int(vm[2], vc[2], vp[2], fl[1], fr[1]);

          // By
          mc2int(vm[4], vc[4], vp[4], fl[2], fr[2]);

          // Bz
          mc2int(vm[5], vc[5], vp[5], fl[3], fr[3]);
        }

        // z-face
        {
          float64 *vm = &ebz(iz,iy,ix-1,0);
          float64 *vc = &ebz(iz,iy,ix  ,0);
          float64 *vp = &ebz(iz,iy,ix+1,0);

          // Ey
          mc2int(vm[1], vc[1], vp[1], fl[4], fr[4]);

          // Ez
          mc2int(vm[2], vc[2], vp[2], fl[5], fr[5]);

          // By
          mc2int(vm[4], vc[4], vp[4], fl[6], fr[6]);

          // Bz
          mc2int(vm[5], vc[5], vp[5], fl[7], fr[7]);
        }
      }
    }
  }
}

void MC2::yint_f2e(Global &g, T_vector &ebz, T_vector &ebx)
{
  for(int iz=g.Lbz-1; iz <= g.Ubz+1 ;iz++) {
    for(int iy=g.Lby-1; iy <= g.Uby+1 ;iy++) {
      for(int ix=g.Lbx-1; ix <= g.Ubx+1 ;ix++) {
        float64 *fr = &ur(iz,iy-1,ix,0);
        float64 *fl = &ul(iz,iy  ,ix,0);

        // z-face
        {
          float64 *vm = &ebz(iz,iy-1,ix,0);
          float64 *vc = &ebz(iz,iy  ,ix,0);
          float64 *vp = &ebz(iz,iy+1,ix,0);

          // Ex
          mc2int(vm[0], vc[0], vp[0], fl[0], fr[0]);

          // Ez
          mc2int(vm[2], vc[2], vp[2], fl[1], fr[1]);

          // Bx
          mc2int(vm[3], vc[3], vp[3], fl[2], fr[2]);

          // Bz
          mc2int(vm[5], vc[5], vp[5], fl[3], fr[3]);
        }

        // x-face
        {
          float64 *vm = &ebx(iz,iy-1,ix,0);
          float64 *vc = &ebx(iz,iy  ,ix,0);
          float64 *vp = &ebx(iz,iy+1,ix,0);

          // Ex
          mc2int(vm[0], vc[0], vp[0], fl[4], fr[4]);

          // Ez
          mc2int(vm[2], vc[2], vp[2], fl[5], fr[5]);

          // Bx
          mc2int(vm[3], vc[3], vp[3], fl[6], fr[6]);

          // Bz
          mc2int(vm[5], vc[5], vp[5], fl[7], fr[7]);
        }
      }
    }
  }
}

void MC2::zint_f2e(Global &g, T_vector &ebx, T_vector &eby)
{
  for(int iz=g.Lbz-1; iz <= g.Ubz+1 ;iz++) {
    for(int iy=g.Lby-1; iy <= g.Uby+1 ;iy++) {
      for(int ix=g.Lbx-1; ix <= g.Ubx+1 ;ix++) {
        float64 *fr = &ur(iz-1,iy,ix,0);
        float64 *fl = &ul(iz  ,iy,ix,0);

        // x-face
        {
          float64 *vm = &ebx(iz-1,iy,ix,0);
          float64 *vc = &ebx(iz  ,iy,ix,0);
          float64 *vp = &ebx(iz+1,iy,ix,0);

          // Ex
          mc2int(vm[0], vc[0], vp[0], fl[0], fr[0]);

          // Ey
          mc2int(vm[1], vc[1], vp[1], fl[1], fr[1]);

          // Bx
          mc2int(vm[3], vc[3], vp[3], fl[2], fr[2]);

          // By
          mc2int(vm[4], vc[4], vp[4], fl[3], fr[3]);
        }

        // y-face
        {
          float64 *vm = &eby(iz-1,iy,ix,0);
          float64 *vc = &eby(iz  ,iy,ix,0);
          float64 *vp = &eby(iz+1,iy,ix,0);

          // Ex
          mc2int(vm[0], vc[0], vp[0], fl[4], fr[4]);

          // Ey
          mc2int(vm[1], vc[1], vp[1], fl[5], fr[5]);

          // Bx
          mc2int(vm[3], vc[3], vp[3], fl[6], fr[6]);

          // By
          mc2int(vm[4], vc[4], vp[4], fl[7], fr[7]);
        }
      }
    }
  }
}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
