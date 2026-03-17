// -*- C++ -*-

///
/// @file global.cpp
/// @brief Implementation of global class
///
/// $Id: global.cpp,v 641a96a86755 2015/09/14 19:35:09 amano $
///
#include "global.hpp"

common::DebugStream Global::debugstream;

/// return time step
float64 Global::get_dt_factor(float64 cfl)
{
  float64 dtmax, dtexp, delta;

  delta = limiter::min(delx, dely, delz);
  dtmax = cfl * delta * rc;
  dtexp = ceil(log2(delt/dtmax));
  return pow(2.0, -dtexp);
}

/// show message
void Global::print_message(std::ostream &ofs)
{
  ofs << msg;
}

/// write simulation parameter
void Global::write_parameter(std::ofstream &ofs)
{
  size_t isize = sizeof(int32);
  size_t dsize = sizeof(float64);

  ofs.write(reinterpret_cast<const char*>(&Nb), isize);
  ofs.write(reinterpret_cast<const char*>(&Nx), isize);
  ofs.write(reinterpret_cast<const char*>(&Ny), isize);
  ofs.write(reinterpret_cast<const char*>(&Nz), isize);
  ofs.write(reinterpret_cast<const char*>(&delx), dsize);
  ofs.write(reinterpret_cast<const char*>(&dely), dsize);
  ofs.write(reinterpret_cast<const char*>(&delz), dsize);
  for(int ix=Lbx; ix <= Ubx ;ix++) {
    ofs.write(reinterpret_cast<const char*>(&xig(ix)), dsize);
  }
  for(int iy=Lby; iy <= Uby ;iy++) {
    ofs.write(reinterpret_cast<const char*>(&yig(iy)), dsize);
  }
  for(int iz=Lbz; iz <= Ubz ;iz++) {
    ofs.write(reinterpret_cast<const char*>(&zig(iz)), dsize);
  }
  // physical parameters
  ofs.write(reinterpret_cast<const char*>(&c), dsize);
  ofs.write(reinterpret_cast<const char*>(&gam), dsize);
  // proton
  ofs.write(reinterpret_cast<const char*>(&qp), dsize);
  ofs.write(reinterpret_cast<const char*>(&mp), dsize);
  // electron
  ofs.write(reinterpret_cast<const char*>(&qe), dsize);
  ofs.write(reinterpret_cast<const char*>(&me), dsize);
}

/// save data to the disk
void Global::write_data(std::ofstream &ofs, bool init)
{
  static int Nwc;
  size_t isize = sizeof(int32);
  size_t dsize = sizeof(float64);

  if( init ) {
    Nwc = 1;
    write_parameter(ofs);
  }

  // time
  ofs.write(reinterpret_cast<const char*>(&physt), dsize);
  ofs.write(reinterpret_cast<const char*>(&delt), dsize);

  // fluid
  for(int iz=Lbz; iz <= Ubz ;iz++) {
    for(int iy=Lby; iy <= Uby ;iy++) {
      for(int ix=Lbx; ix <= Ubx ;ix++) {
        ofs.write(reinterpret_cast<const char*>(&uf(iz,iy,ix,0)), 10*dsize);
      }
    }
  }

  // electromagentic field at cell center
  for(int iz=Lbz; iz <= Ubz ;iz++) {
    for(int iy=Lby; iy <= Uby ;iy++) {
      for(int ix=Lbx; ix <= Ubx ;ix++) {
        ofs.write(reinterpret_cast<const char*>(&ebc(iz,iy,ix,0)), 6*dsize);
      }
    }
  }

  // electromagnetic field at face center
  for(int iz=Lbz-1; iz <= Ubz ;iz++) {
    for(int iy=Lby-1; iy <= Uby ;iy++) {
      for(int ix=Lbx-1; ix <= Ubx ;ix++) {
        ofs.write(reinterpret_cast<const char*>(&ueb(iz,iy,ix,0)), 6*dsize);
      }
    }
  }

  // end-of-record marker
  ofs.write(reinterpret_cast<const char*>(&Nwc), isize);
  ofs.flush();
  Nwc++;
}

/// dump whole array to the disk
template <class T_array>
void Global::write_debug_array(std::ofstream &ofs, T_array &x)
{
  // write rank
  int rank = x.rank();
  ofs.write(reinterpret_cast<const char*>(&rank), sizeof(int));

  // write array shape
  for(int r=0; r < rank ;r++) {
    int ext = x.extent(r);
    ofs.write(reinterpret_cast<const char*>(&ext), sizeof(int));
  }

  // dump data content
  {
    typename T_array::T_numtype *ptr = x.data();

    ofs.write(reinterpret_cast<const char*>(ptr), sizeof(*ptr)*x.size());
  }

  ofs.flush();
}

/// save debugging data to the disk
void Global::write_debug_data(std::ofstream &ofs, bool init)
{
  static int Nwc;
  size_t isize = sizeof(int32);
  size_t dsize = sizeof(float64);

  if( init ) {
    Nwc = 1;
    write_parameter(ofs);
  }

  // time
  ofs.write(reinterpret_cast<const char*>(&physt), dsize);
  ofs.write(reinterpret_cast<const char*>(&delt), dsize);

  // dump data
  write_debug_array(ofs, uf);
  write_debug_array(ofs, ff);
  write_debug_array(ofs, ebc);
  write_debug_array(ofs, ueb);
  write_debug_array(ofs, feb);

  // calculate div(B) and div(E) error
  {
    float64 rdx = 1/delx;
    float64 rdy = 1/dely;
    float64 rdz = 1/delz;

    for(int iz=Lbz; iz <= Ubz ;iz++) {
      for(int iy=Lby; iy <= Uby ;iy++) {
        for(int ix=Lbx; ix <= Ubx ;ix++) {
          float64 *up = &uf(iz,iy,ix,0);
          float64 *ue = &uf(iz,iy,ix,5);
          float64 gp = lorentz(up[1], up[2], up[3]);
          float64 ge = lorentz(ue[1], ue[2], ue[3]);
          diverr(iz,iy,ix,0) =
            - common::pi4 * (qp*gp*up[0] + qe*ge*ue[0])
            + rdx*(ueb(iz,iy,ix,0) - ueb(iz,iy,ix-1,0))
            + rdy*(ueb(iz,iy,ix,1) - ueb(iz,iy-1,ix,1))
            + rdz*(ueb(iz,iy,ix,2) - ueb(iz-1,iy,ix,2));
          diverr(iz,iy,ix,1) =
            + rdx*(ueb(iz,iy,ix,3) - ueb(iz,iy,ix-1,3))
            + rdy*(ueb(iz,iy,ix,4) - ueb(iz,iy-1,ix,4))
            + rdz*(ueb(iz,iy,ix,5) - ueb(iz-1,iy,ix,5));
        }
      }
    }
  }
  write_debug_array(ofs, diverr);

  // end-of-record marker
  ofs.write(reinterpret_cast<const char*>(&Nwc), isize);
  ofs.flush();
  Nwc++;
}

/// calculate energy
void Global::energy(float64 &Ef, float64 &Ep)
{
  // electromagnetic field
  Ef = 0.0;
  for(int iz=Lbz; iz <= Ubz ;iz++) {
    for(int iy=Lby; iy <= Uby ;iy++) {
      for(int ix=Lbx; ix <= Ubx ;ix++) {
        float64 *eb = &ebc(iz,iy,ix,0);
        Ef += 0.5 *(eb[0]*eb[0] + eb[1]*eb[1] + eb[2]*eb[2] +
                    eb[3]*eb[3] + eb[4]*eb[4] + eb[5]*eb[5]) * rpi4;
      }
    }
  }

  // plasma
  Ep = 0.0;
  for(int iz=Lbz; iz <= Ubz ;iz++) {
    for(int iy=Lby; iy <= Uby ;iy++) {
      for(int ix=Lbx; ix <= Ubx ;ix++) {
        float64 *up = &uf(iz,iy,ix,0);
        float64 *ue = &uf(iz,iy,ix,5);
        float64 gp = lorentz(up[1], up[2], up[3]);
        float64 ge = lorentz(ue[1], ue[2], ue[3]);
        float64 wp = up[0]*mp*cc + G*up[4];
        float64 we = ue[0]*me*cc + G*ue[4];
        Ep += gp*gp*wp + ge*ge*we - (up[4] + ue[4]);
      }
    }
  }

  // take sum of all PEs
  {
    float64 le[2], ge[2];

    le[0] = Ef;
    le[1] = Ep;

    MPI_Allreduce(le, ge, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD);

    Ef = ge[0];
    Ep = ge[1];
  }
}

template <>
void Global::int_c2f<2>(T_vector &ebc, T_vector &ueb)
{
  // Ex, Bx
  for(int iz=Lbz-Nb; iz <= Ubz+Nb ;iz++) {
    for(int iy=Lby-Nb; iy <= Uby+Nb ;iy++) {
      for(int ix=Lbx-Nb; ix <= Ubx+Nb-1 ;ix++) {
        ueb(iz,iy,ix,0) = 0.5*(ebc(iz,iy,ix,0) + ebc(iz,iy,ix+1,0));
        ueb(iz,iy,ix,3) = 0.5*(ebc(iz,iy,ix,3) + ebc(iz,iy,ix+1,3));
      }
    }
  }

  // Ey, By
  for(int iz=Lbz-Nb; iz <= Ubz+Nb ;iz++) {
    for(int iy=Lby-Nb; iy <= Uby+Nb-1 ;iy++) {
      for(int ix=Lbx-Nb; ix <= Ubx+Nb ;ix++) {
        ueb(iz,iy,ix,1) = 0.5*(ebc(iz,iy,ix,1) + ebc(iz,iy+1,ix,1));
        ueb(iz,iy,ix,4) = 0.5*(ebc(iz,iy,ix,4) + ebc(iz,iy+1,ix,4));
      }
    }
  }

  // Ez, Bz
  for(int iz=Lbz-Nb; iz <= Ubz+Nb-1 ;iz++) {
    for(int iy=Lby-Nb; iy <= Uby+Nb ;iy++) {
      for(int ix=Lbx-Nb; ix <= Ubx+Nb ;ix++) {
        ueb(iz,iy,ix,2) = 0.5*(ebc(iz,iy,ix,2) + ebc(iz+1,iy,ix,2));
        ueb(iz,iy,ix,5) = 0.5*(ebc(iz,iy,ix,5) + ebc(iz+1,iy,ix,5));
      }
    }
  }
}

template <>
void Global::int_f2c<2>(T_vector &ebc, T_vector &ueb)
{
  // Ex, Bx
  for(int iz=Lbz-Nb; iz <= Ubz+Nb ;iz++) {
    for(int iy=Lby-Nb; iy <= Uby+Nb ;iy++) {
      for(int ix=Lbx-Nb+1; ix <= Ubx+Nb ;ix++) {
        ebc(iz,iy,ix,0) = 0.5*(ueb(iz,iy,ix,0) + ueb(iz,iy,ix-1,0));
        ebc(iz,iy,ix,3) = 0.5*(ueb(iz,iy,ix,3) + ueb(iz,iy,ix-1,3));
      }
    }
  }

  // Ey, By
  for(int iz=Lbz-Nb; iz <= Ubz+Nb ;iz++) {
    for(int iy=Lby-Nb+1; iy <= Uby+Nb ;iy++) {
      for(int ix=Lbx-Nb; ix <= Ubx+Nb ;ix++) {
        ebc(iz,iy,ix,1) = 0.5*(ueb(iz,iy,ix,1) + ueb(iz,iy-1,ix,1));
        ebc(iz,iy,ix,4) = 0.5*(ueb(iz,iy,ix,4) + ueb(iz,iy-1,ix,4));
      }
    }
  }

  // Ez, Bz
  for(int iz=Lbz-Nb+1; iz <= Ubz+Nb ;iz++) {
    for(int iy=Lby-Nb; iy <= Uby+Nb ;iy++) {
      for(int ix=Lbx-Nb; ix <= Ubx+Nb ;ix++) {
        ebc(iz,iy,ix,2) = 0.5*(ueb(iz,iy,ix,2) + ueb(iz-1,iy,ix,2));
        ebc(iz,iy,ix,5) = 0.5*(ueb(iz,iy,ix,5) + ueb(iz-1,iy,ix,5));
      }
    }
  }
}

/// convert conservative variables to primitive variables
void Global::primitive(float64 uc[10], float64 eb[6], float64 uf[10])
{
  float64 var, up[5], ue[5], S[4];

  // electromagnetic energy and momentum
  S[0] = rpi4*(eb[0]*eb[0] + eb[1]*eb[1] + eb[2]*eb[2] +
               eb[3]*eb[3] + eb[4]*eb[4] + eb[5]*eb[5]) * 0.5;
  S[1] = rpi4*(eb[1]*eb[5] - eb[2]*eb[4]) * c;
  S[2] = rpi4*(eb[2]*eb[3] - eb[0]*eb[5]) * c;
  S[3] = rpi4*(eb[0]*eb[4] - eb[1]*eb[3]) * c;

  // proton
  var   = mp/(qe*mp - qp*me);
  up[0] = (qe*uc[0] - me*uc[5]          ) * var;
  up[1] = (qe*uc[1] - me*uc[6] - qe*S[1]) * var;
  up[2] = (qe*uc[2] - me*uc[7] - qe*S[2]) * var;
  up[3] = (qe*uc[3] - me*uc[8] - qe*S[3]) * var;
  up[4] = (qe*uc[4] - me*uc[9] - qe*S[0]) * var;
  rhd_primitive(up, &uf[0], mp, c, G);

  // electron
  var   = me/(qp*me - qe*mp);
  ue[0] = (qp*uc[0] - mp*uc[5]          ) * var;
  ue[1] = (qp*uc[1] - mp*uc[6] - qp*S[1]) * var;
  ue[2] = (qp*uc[2] - mp*uc[7] - qp*S[2]) * var;
  ue[3] = (qp*uc[3] - mp*uc[8] - qp*S[3]) * var;
  ue[4] = (qp*uc[4] - mp*uc[9] - qp*S[0]) * var;
  rhd_primitive(ue, &uf[5], me, c, G);
}

/// convert primitive variables to conservative variables
void Global::conservative(float64 uf[10], float64 eb[6], float64 uc[10])
{
  float64 S[4], gp, ge, wp, we, p;
  float64 *up = &uf[0];
  float64 *ue = &uf[5];

  gp = lorentz(up[1], up[2], up[3]);
  ge = lorentz(ue[1], ue[2], ue[3]);

  // electromagnetic energy and momentum
  S[0] = rpi4*(eb[0]*eb[0] + eb[1]*eb[1] + eb[2]*eb[2] +
               eb[3]*eb[3] + eb[4]*eb[4] + eb[5]*eb[5]) * 0.5;
  S[1] = rpi4*(eb[1]*eb[5] - eb[2]*eb[4]) * c;
  S[2] = rpi4*(eb[2]*eb[3] - eb[0]*eb[5]) * c;
  S[3] = rpi4*(eb[0]*eb[4] - eb[1]*eb[3]) * c;

  // sum of two fluids
  wp = up[0]*mp*cc + G*up[4];
  we = ue[0]*me*cc + G*ue[4];
  p  = up[4] + ue[4];
  uc[0] = mp*gp*up[0] + me*ge*ue[0];
  uc[1] =(wp*gp*up[1] + we*ge*ue[1] + S[1])*rcc;
  uc[2] =(wp*gp*up[2] + we*ge*ue[2] + S[2])*rcc;
  uc[3] =(wp*gp*up[3] + we*ge*ue[3] + S[3])*rcc;
  uc[4] = wp*gp*gp  + we*ge*ge - p + S[0];

  // difference of two fluids
  wp *= qmp;
  we *= qme;
  p   = qmp*up[4] + qme*ue[4];
  uc[5] = qp*gp*up[0] + qe*ge*ue[0];
  uc[6] =(wp*gp*up[1] + we*ge*ue[1])*rcc;
  uc[7] =(wp*gp*up[2] + we*ge*ue[2])*rcc;
  uc[8] =(wp*gp*up[3] + we*ge*ue[3])*rcc;
  uc[9] = wp*gp*gp + we*ge*ge - p;
}

/// calculate source term (right-hand side) for conservative variables
void Global::rhs(float64 dt, float64 uf[10], float64 eb[6], float64 rhs[10])
{
  float64 U[4], J[4];
  float64 gp, ge, wp0, wpp, wpe, ro0, rop, roe;
  float64 *up = &uf[0];
  float64 *ue = &uf[5];

  gp = lorentz(up[1], up[2], up[3]);
  ge = lorentz(ue[1], ue[2], ue[3]);

  wpp  = qmp*qp*up[0];
  wpe  = qme*qe*ue[0];
  U[0] = wpp*gp    + wpe*ge;
  U[1] = wpp*up[1] + wpe*ue[1];
  U[2] = wpp*up[2] + wpe*ue[2];
  U[3] = wpp*up[3] + wpe*ue[3];

  rop  = qp*up[0];
  roe  = qe*ue[0];
  J[0] = rop*gp    + roe*ge;
  J[1] = rop*up[1] + roe*ue[1];
  J[2] = rop*up[2] + roe*ue[2];
  J[3] = rop*up[3] + roe*ue[3];

  wp0  = wpp + wpe;
  ro0  = (U[0]*J[0] - (U[1]*J[1] + U[2]*J[2] + U[3]*J[3])*rc*rc) / wp0;

  rhs[0] = 0.0;
  rhs[1] = 0.0;
  rhs[2] = 0.0;
  rhs[3] = 0.0;
  rhs[4] = 0.0;

  rhs[5] = 0.0;
  rhs[6] = dt*(U[0]*eb[0] + (U[2]*eb[5] - U[3]*eb[4])*rc);
  rhs[7] = dt*(U[0]*eb[1] + (U[3]*eb[3] - U[1]*eb[5])*rc);
  rhs[8] = dt*(U[0]*eb[2] + (U[1]*eb[4] - U[2]*eb[3])*rc);
  rhs[9] = dt*(U[1]*eb[0] + U[2]*eb[1] + U[3]*eb[2]);

  // resistivity
  rhs[6] -= dt*eta*(wp0*J[1] - ro0*U[1]);
  rhs[7] -= dt*eta*(wp0*J[2] - ro0*U[2]);
  rhs[8] -= dt*eta*(wp0*J[3] - ro0*U[3]);
  rhs[9] -= dt*eta*(wp0*J[0] - ro0*U[0]) * c;
}

/// calculate flux in x direction
void Global::xflux(float64 uf[10], float64 eb[6], float64 fx[10])
{
  float64 S[4], gp, ge, wp, we, p;
  float64 *up = &uf[0];
  float64 *ue = &uf[5];

  gp = lorentz(up[1], up[2], up[3]);
  ge = lorentz(ue[1], ue[2], ue[3]);

  // electromagnetic part
  S[0] = rpi4*(eb[1]*eb[5] - eb[2]*eb[4]) * c;
  S[1] = rpi4*(-eb[0]*eb[0] + eb[1]*eb[1] + eb[2]*eb[2]
               -eb[3]*eb[3] + eb[4]*eb[4] + eb[5]*eb[5]) * 0.5;
  S[2] =-rpi4*(eb[0]*eb[1] + eb[3]*eb[4]);
  S[3] =-rpi4*(eb[0]*eb[2] + eb[3]*eb[5]);

  // sum of two fluids
  wp = up[0]*mp*cc + G*up[4];
  we = ue[0]*me*cc + G*ue[4];
  p  = up[4] + ue[4];
  fx[0] = mp*up[1]*up[0] + me*ue[1]*ue[0];
  fx[1] =(wp*up[1]*up[1] + we*ue[1]*ue[1])*rcc + p + S[1];
  fx[2] =(wp*up[1]*up[2] + we*ue[1]*ue[2])*rcc     + S[2];
  fx[3] =(wp*up[1]*up[3] + we*ue[1]*ue[3])*rcc     + S[3];
  fx[4] = wp*gp*up[1] + we*ge*ue[1]                + S[0];

  // difference of two fluids
  wp *= qmp;
  we *= qme;
  p   = qmp*up[4] + qme*ue[4];
  fx[5] = qp*up[1]*up[0] + qe*ue[1]*ue[0];
  fx[6] =(wp*up[1]*up[1] + we*ue[1]*ue[1])*rcc + p;
  fx[7] =(wp*up[1]*up[2] + we*ue[1]*ue[2])*rcc;
  fx[8] =(wp*up[1]*up[3] + we*ue[1]*ue[3])*rcc;
  fx[9] = wp*gp*up[1] + we*ge*ue[1];
}

/// calculate flux in y direction
void Global::yflux(float64 uf[10], float64 eb[6], float64 fy[10])
{
  float64 S[4], gp, ge, wp, we, p;
  float64 *up = &uf[0];
  float64 *ue = &uf[5];

  gp = lorentz(up[1], up[2], up[3]);
  ge = lorentz(ue[1], ue[2], ue[3]);

  // electromagnetic part
  S[0] = rpi4*(eb[2]*eb[3] - eb[0]*eb[5]) * c;
  S[1] =-rpi4*(eb[1]*eb[0] + eb[4]*eb[3]);
  S[2] = rpi4*(+eb[0]*eb[0] - eb[1]*eb[1] + eb[2]*eb[2]
               +eb[3]*eb[3] - eb[4]*eb[4] + eb[5]*eb[5]) * 0.5;
  S[3] =-rpi4*(eb[1]*eb[2] + eb[4]*eb[5]);

  // sum of two fluids
  wp = up[0]*mp*cc + G*up[4];
  we = ue[0]*me*cc + G*ue[4];
  p  = up[4] + ue[4];
  fy[0] = mp*up[2]*up[0] + me*ue[2]*ue[0];
  fy[1] =(wp*up[2]*up[1] + we*ue[2]*ue[1])*rcc     + S[1];
  fy[2] =(wp*up[2]*up[2] + we*ue[2]*ue[2])*rcc + p + S[2];
  fy[3] =(wp*up[2]*up[3] + we*ue[2]*ue[3])*rcc     + S[3];
  fy[4] = wp*gp*up[2] + we*ge*ue[2]                + S[0];

  // difference of two fluids
  wp *= qmp;
  we *= qme;
  p   = qmp*up[4] + qme*ue[4];
  fy[5] = qp*up[2]*up[0] + qe*ue[2]*ue[0];
  fy[6] =(wp*up[2]*up[1] + we*ue[2]*ue[1])*rcc;
  fy[7] =(wp*up[2]*up[2] + we*ue[2]*ue[2])*rcc + p;
  fy[8] =(wp*up[2]*up[3] + we*ue[2]*ue[3])*rcc;
  fy[9] = wp*gp*up[2] + we*ge*ue[2];
}

/// calculate flux in z direction
void Global::zflux(float64 uf[10], float64 eb[6], float64 fz[10])
{
  float64 S[4], gp, ge, wp, we, p;
  float64 *up = &uf[0];
  float64 *ue = &uf[5];

  gp = lorentz(up[1], up[2], up[3]);
  ge = lorentz(ue[1], ue[2], ue[3]);

  // electromagnetic part
  S[0] = rpi4*(eb[0]*eb[4] - eb[1]*eb[3]) * c;
  S[1] =-rpi4*(eb[2]*eb[0] + eb[5]*eb[3]);
  S[2] =-rpi4*(eb[2]*eb[1] + eb[5]*eb[4]);
  S[3] = rpi4*(+eb[0]*eb[0] + eb[1]*eb[1] - eb[2]*eb[2]
               +eb[3]*eb[3] + eb[4]*eb[4] - eb[5]*eb[5]) * 0.5;

  // sum of two fluids
  wp = up[0]*mp*cc + G*up[4];
  we = ue[0]*me*cc + G*ue[4];
  p  = up[4] + ue[4];
  fz[0] = mp*up[3]*up[0] + me*ue[3]*ue[0];
  fz[1] =(wp*up[3]*up[1] + we*ue[3]*ue[1])*rcc     + S[1];
  fz[2] =(wp*up[3]*up[2] + we*ue[3]*ue[2])*rcc     + S[2];
  fz[3] =(wp*up[3]*up[3] + we*ue[3]*ue[3])*rcc + p + S[3];
  fz[4] = wp*gp*up[3] + we*ge*ue[3]                + S[0];

  // difference of two fluids
  wp *= qmp;
  we *= qme;
  p   = qmp*up[4] + qme*ue[4];
  fz[5] = qp*up[3]*up[0] + qe*ue[3]*ue[0];
  fz[6] =(wp*up[3]*up[1] + we*ue[3]*ue[1])*rcc;
  fz[7] =(wp*up[3]*up[2] + we*ue[3]*ue[2])*rcc;
  fz[8] =(wp*up[3]*up[3] + we*ue[3]*ue[3])*rcc + p;
  fz[9] = wp*gp*up[3] + we*ge*ue[3];
}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
