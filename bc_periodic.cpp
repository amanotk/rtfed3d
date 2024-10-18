// -*- C++ -*-

///
/// @brief Implementation of periodic boundary condition via MPI
///
/// $Id: bc_periodic.cpp,v 99bd1b55086c 2015/09/02 19:54:43 amano $
///
#include "bc_periodic.hpp"

template <class T_array, class T_shape>
void PeriodicBoundary::pack(float64 *buf, int &pos, T_array &array,
                            T_shape &Lb, T_shape &Ub, bool is_boundary)
{
  if(is_boundary) return;

  for(int iz=Lb[0]; iz <= Ub[0] ;iz++) {
    for(int iy=Lb[1]; iy <= Ub[1] ;iy++) {
      for(int ix=Lb[2]; ix <= Ub[2] ;ix++) {
        for(int ic=Lb[3]; ic <= Ub[3] ;ic++) {
          buf[pos] = array(iz,iy,ix,ic);
          pos++;
        }
      }
    }
  }
}

template <class T_array, class T_shape>
void PeriodicBoundary::unpack(float64 *buf, int &pos, T_array &array,
                              T_shape &Lb, T_shape &Ub, bool is_boundary)
{
  if(is_boundary) return;

  for(int iz=Lb[0]; iz <= Ub[0] ;iz++) {
    for(int iy=Lb[1]; iy <= Ub[1] ;iy++) {
      for(int ix=Lb[2]; ix <= Ub[2] ;ix++) {
        for(int ic=Lb[3]; ic <= Ub[3] ;ic++) {
          array(iz,iy,ix,ic) = buf[pos];
          pos++;
        }
      }
    }
  }
}

PeriodicBoundary::PeriodicBoundary(const int Nz, const int Ny, const int Nx,
                                   const int Nb)
{
  int Mm, Mx, My, Mz, Lbx, Lby, Lbz, Ubx, Uby, Ubz;

  calc_grid_bounds(Nz, Nb, Mz, Lbz, Ubz);
  calc_grid_bounds(Ny, Nb, My, Lby, Uby);
  calc_grid_bounds(Nx, Nb, Mx, Lbx, Ubx);
  Mm = limiter::max(Mz*My, Mz*Mx, My*Mx);

  m_lb = Lbz, Lby, Lbx;
  m_ub = Ubz, Uby, Ubx;

  // allocate buffer
  m_bufnum = Mm * 16*Nb;
  for(int dir=0; dir < 3 ;dir++) {
    m_bufsnd[dir][0] = new float64 [m_bufnum];
    m_bufsnd[dir][1] = new float64 [m_bufnum];
    m_bufrcv[dir][0] = new float64 [m_bufnum];
    m_bufrcv[dir][1] = new float64 [m_bufnum];
  }
}

PeriodicBoundary::~PeriodicBoundary()
{
  for(int dir=0; dir < 3 ;dir++) {
    delete [] m_bufsnd[dir][0];
    delete [] m_bufsnd[dir][1];
    delete [] m_bufrcv[dir][0];
    delete [] m_bufrcv[dir][1];
  }
}

void PeriodicBoundary::set_field(Global &g, T_vector &eb, int nb)
{
  const int Nb = (nb < 0) ? g.Nb : nb;
  int neighbor[3][2];
  int tag[3][2];

  mpiutil::getNeighbor(neighbor);
  mpiutil::getTag(tag);

  // directional send/recv
  for(int dir=0; dir < 3 ;dir++) {
    bool is_lower = mpiutil::isLowerBoundary(dir);
    bool is_upper = mpiutil::isUpperBoundary(dir);
    MPI_Request request[4];
    int pos;
    float64 *buf;

    blitz::TinyVector<int,4> eb_lb = eb.lbound();
    blitz::TinyVector<int,4> eb_ub = eb.ubound();

    //
    // to lower
    //
    pos = 0;
    buf = m_bufsnd[dir][0];
    // pack eb
    eb_lb[dir] = m_lb[dir];
    eb_ub[dir] = m_lb[dir]+Nb-1;
    pack(buf, pos, eb, eb_lb, eb_ub, is_lower);
    // send
    MPI_Isend(m_bufsnd[dir][0], pos, MPI_REAL8,
              neighbor[dir][0], tag[dir][0], MPI_COMM_WORLD, &request[0]);
    // receive
    MPI_Irecv(m_bufrcv[dir][1], m_bufnum, MPI_REAL8,
              neighbor[dir][1], tag[dir][0], MPI_COMM_WORLD, &request[3]);

    //
    // to upper
    //
    pos = 0;
    buf = m_bufsnd[dir][1];
    // pack eb
    eb_lb[dir] = m_ub[dir]-Nb+1;
    eb_ub[dir] = m_ub[dir];
    pack(buf, pos, eb, eb_lb, eb_ub, is_upper);
    // send
    MPI_Isend(m_bufsnd[dir][1], pos, MPI_REAL8,
              neighbor[dir][1], tag[dir][1], MPI_COMM_WORLD, &request[1]);
    // receive
    MPI_Irecv(m_bufrcv[dir][0], m_bufnum, MPI_REAL8,
              neighbor[dir][0], tag[dir][1], MPI_COMM_WORLD, &request[2]);

    // wait
    MPI_Waitall(4, request, MPI_STATUSES_IGNORE);

    //
    // from lower
    //
    pos = 0;
    buf = m_bufrcv[dir][0];
    // unpack eb
    eb_lb[dir] = m_lb[dir]-Nb;
    eb_ub[dir] = m_lb[dir]-1;
    unpack(buf, pos, eb, eb_lb, eb_ub, is_lower);

    //
    // from upper
    //
    pos = 0;
    buf = m_bufrcv[dir][1];
    // unpack eb
    eb_lb[dir] = m_ub[dir]+1;
    eb_ub[dir] = m_ub[dir]+Nb;
    unpack(buf, pos, eb, eb_lb, eb_ub, is_upper);
  }
}

void PeriodicBoundary::set_fluid(Global &g, T_vector &uf, T_vector &eb,
                                 int nb)
{
  const int Nb = (nb < 0) ? g.Nb : nb;
  int neighbor[3][2];
  int tag[3][2];

  mpiutil::getNeighbor(neighbor);
  mpiutil::getTag(tag);

  // directional send/recv
  for(int dir=0; dir < 3 ;dir++) {
    bool is_lower = mpiutil::isLowerBoundary(dir);
    bool is_upper = mpiutil::isUpperBoundary(dir);
    MPI_Request request[4];
    int pos;
    float64 *buf;

    blitz::TinyVector<int,4> uf_lb = uf.lbound();
    blitz::TinyVector<int,4> uf_ub = uf.ubound();
    blitz::TinyVector<int,4> eb_lb = eb.lbound();
    blitz::TinyVector<int,4> eb_ub = eb.ubound();

    //
    // to lower
    //
    pos = 0;
    buf = m_bufsnd[dir][0];
    // pack uf
    uf_lb[dir] = m_lb[dir];
    uf_ub[dir] = m_lb[dir]+Nb-1;
    pack(buf, pos, uf, uf_lb, uf_ub, is_lower);
    // pack eb
    eb_lb[dir] = m_lb[dir];
    eb_ub[dir] = m_lb[dir]+Nb-1;
    pack(buf, pos, eb, eb_lb, eb_ub, is_lower);
    // send
    MPI_Isend(m_bufsnd[dir][0], pos, MPI_REAL8,
              neighbor[dir][0], tag[dir][0], MPI_COMM_WORLD, &request[0]);
    // receive
    MPI_Irecv(m_bufrcv[dir][1], m_bufnum, MPI_REAL8,
              neighbor[dir][1], tag[dir][0], MPI_COMM_WORLD, &request[3]);

    //
    // to upper
    //
    pos = 0;
    buf = m_bufsnd[dir][1];
    // pack uf
    uf_lb[dir] = m_ub[dir]-Nb+1;
    uf_ub[dir] = m_ub[dir];
    pack(buf, pos, uf, uf_lb, uf_ub, is_upper);
    // pack eb
    eb_lb[dir] = m_ub[dir]-Nb+1;
    eb_ub[dir] = m_ub[dir];
    pack(buf, pos, eb, eb_lb, eb_ub, is_upper);
    // send
    MPI_Isend(m_bufsnd[dir][1], pos, MPI_REAL8,
              neighbor[dir][1], tag[dir][1], MPI_COMM_WORLD, &request[1]);
    // receive
    MPI_Irecv(m_bufrcv[dir][0], m_bufnum, MPI_REAL8,
              neighbor[dir][0], tag[dir][1], MPI_COMM_WORLD, &request[2]);

    // wait
    MPI_Waitall(4, request, MPI_STATUSES_IGNORE);

    //
    // from lower
    //
    pos = 0;
    buf = m_bufrcv[dir][0];
    // unpack uf
    uf_lb[dir] = m_lb[dir]-Nb;
    uf_ub[dir] = m_lb[dir]-1;
    unpack(buf, pos, uf, uf_lb, uf_ub, is_lower);
    // unpack eb
    eb_lb[dir] = m_lb[dir]-Nb;
    eb_ub[dir] = m_lb[dir]-1;
    unpack(buf, pos, eb, eb_lb, eb_ub, is_lower);

    //
    // from upper
    //
    pos = 0;
    buf = m_bufrcv[dir][1];
    // unpack uf
    uf_lb[dir] = m_ub[dir]+1;
    uf_ub[dir] = m_ub[dir]+Nb;
    unpack(buf, pos, uf, uf_lb, uf_ub, is_upper);
    // unpack eb
    eb_lb[dir] = m_ub[dir]+1;
    eb_ub[dir] = m_ub[dir]+Nb;
    unpack(buf, pos, eb, eb_lb, eb_ub, is_upper);
  }
}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
