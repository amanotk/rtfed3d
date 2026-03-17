// -*- C++ -*-

///
/// MPI utlility module for three dimensional domain decomposition
///
/// Author: Takanobu AMANO <amanot@stelab.nagoya-u.ac.jp>
/// $Id$
///
#include "mpiutils.hpp"

// definition of static variables
mpiutils *mpiutils::instance;
int mpiutils::tag[3][2] = { {0, 1}, {2, 3}, {4, 5} };

//
// pack data; return position after pack
//
int mpiutils::pack(void *data, int count, MPI_Datatype dtype,
                   void *buf, int bufsize, int pos)
{
  MPI_Pack(data, count, dtype, buf, bufsize, &pos, MPI_COMM_WORLD);
  return pos;
}

//
// unpack data; return position after unpack
//
int mpiutils::unpack(void *data, int count, MPI_Datatype dtype,
                     void *buf, int bufsize, int pos)
{
  MPI_Unpack(buf, bufsize, &pos, data, count, dtype, MPI_COMM_WORLD);
  return pos;
}

//
// begin boundary exchange
//
void mpiutils::bc_exchange_begin(void *buf0, void *buf1, void *buf2,
                                 int dsize, int count[3],
                                 MPI_Request req[12])
{
  bc_exchange_dir_begin(0, buf0, dsize, count[0], &req[4*0]);
  bc_exchange_dir_begin(1, buf1, dsize, count[1], &req[4*1]);
  bc_exchange_dir_begin(2, buf2, dsize, count[2], &req[4*2]);
}

//
// begin non-blocking directional boundary exchange
//
void mpiutils::bc_exchange_dir_begin(int dir, void *buf,
                                     int dsize, int count,
                                     MPI_Request req[4])
{
  MPI_Comm comm = MPI_COMM_WORLD;

  int size  = dsize*count;
  int lower = instance->m_nb_dim[dir][0];
  int upper = instance->m_nb_dim[dir][1];
  char *sndbuf_l = static_cast<char*>(buf) + size * 0;
  char *sndbuf_u = static_cast<char*>(buf) + size * 1;
  char *rcvbuf_l = static_cast<char*>(buf) + size * 2;
  char *rcvbuf_u = static_cast<char*>(buf) + size * 3;

  MPI_Isend(sndbuf_l, size, MPI_BYTE, lower, tag[dir][0], comm, &req[0]);
  MPI_Isend(sndbuf_u, size, MPI_BYTE, upper, tag[dir][1], comm, &req[1]);
  MPI_Irecv(rcvbuf_l, size, MPI_BYTE, lower, tag[dir][1], comm, &req[2]);
  MPI_Irecv(rcvbuf_u, size, MPI_BYTE, upper, tag[dir][0], comm, &req[3]);
}

//
// begin boundary exchange
//
void mpiutils::bc_exchange_begin(mpiutils::Buffer *buf0,
                                 mpiutils::Buffer *buf1,
                                 mpiutils::Buffer *buf2)
{
  bc_exchange_dir_begin(0, buf0);
  bc_exchange_dir_begin(1, buf1);
  bc_exchange_dir_begin(2, buf2);
}

//
// begin non-blocking directional boundary exchange
//
void mpiutils::bc_exchange_dir_begin(int dir, mpiutils::Buffer *buf)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int lower = instance->m_nb_dim[dir][0];
  int upper = instance->m_nb_dim[dir][1];

  MPI_Isend(buf->data[0], buf->count[0], MPI_BYTE, lower, tag[dir][0],
            comm, &buf->request[0]);
  MPI_Isend(buf->data[1], buf->count[1], MPI_BYTE, upper, tag[dir][1],
            comm, &buf->request[1]);
  MPI_Irecv(buf->data[2], buf->count[2], MPI_BYTE, lower, tag[dir][1],
            comm, &buf->request[2]);
  MPI_Irecv(buf->data[3], buf->count[3], MPI_BYTE, upper, tag[dir][0],
            comm, &buf->request[3]);
}

//
// wait for all MPI requests to complete
//
void mpiutils::wait(MPI_Request req[], int n)
{
  MPI_Waitall(n, req, MPI_STATUSES_IGNORE);
}

#ifdef __MAIN__

//
// main program for checking mpiutils usage
//
int main(int argc, char **argv)
{
  int period[3] = {1, 0, 0};

  // initialize with domain decomposition
  mpiutils::initialize(&argc, &argv, period, true);

  // show info
  mpiutils::info(std::cerr);

  // how to use utility functions
  std::cerr << "time in sec = "
            << mpiutils::getTime() << std::endl
            << "filename    = "
            << mpiutils::getFilename("test", "dat") << std::endl;;

  // test directional send/recv
  {
    const int dir = 0;
    const int N = 3;
    int buf[4*N];
    MPI_Request req[4];

    // clear buffer
    memset(buf, -1, 4*N*sizeof(int));

    // fill send buffer by coordinate
    mpiutils::getCoord(&buf[0*N]);
    mpiutils::getCoord(&buf[1*N]);

    mpiutils::bc_exchange_dir_begin(dir, buf, sizeof(int), N, req);
    mpiutils::wait(req, 4);

    // show results
    std::cerr << "--- results of directional send/recv ---" << std::endl;
    std::cerr << tfm::format("send lower => [%2d,%2d,%2d]\n"
                             "send upper => [%2d,%2d,%2d]\n",
                             buf[0*N+0], buf[0*N+1], buf[0*N+2],
                             buf[1*N+0], buf[1*N+1], buf[1*N+2]);
    std::cerr << tfm::format("recv lower => [%2d,%2d,%2d]\n"
                             "recv upper => [%2d,%2d,%2d]\n",
                             buf[2*N+0], buf[2*N+1], buf[2*N+2],
                             buf[3*N+0], buf[3*N+1], buf[3*N+2]);
  }

  // test package send/recv
  {
    const int N0 = 3;
    const int N1 = 3;
    const int N2 = 3;
    int count[3] = {N0, N1, N2};
    int buf0[4*N0];
    int buf1[4*N1];
    int buf2[4*N2];
    MPI_Request req[12];

    // clear buffer
    memset(buf0, -1, 4*N0*sizeof(int));
    memset(buf1, -1, 4*N1*sizeof(int));
    memset(buf2, -1, 4*N2*sizeof(int));

    // fill send buffer by coordinate
    mpiutils::getCoord(&buf0[0*N0]);
    mpiutils::getCoord(&buf0[1*N0]);
    mpiutils::getCoord(&buf1[0*N1]);
    mpiutils::getCoord(&buf1[1*N1]);
    mpiutils::getCoord(&buf2[0*N2]);
    mpiutils::getCoord(&buf2[1*N2]);

    mpiutils::bc_exchange_begin(buf0, buf1, buf2, sizeof(int), count, req);
    mpiutils::wait(req, 12);

    // show results
    std::cerr << "--- results of send/recv in dir. 0 ---" << std::endl;
    std::cerr << tfm::format("send lower => [%2d,%2d,%2d]\n"
                             "send upper => [%2d,%2d,%2d]\n",
                             buf0[0*N0+0], buf0[0*N0+1], buf0[0*N0+2],
                             buf0[1*N0+0], buf0[1*N0+1], buf0[1*N0+2]);
    std::cerr << tfm::format("recv lower => [%2d,%2d,%2d]\n"
                             "recv upper => [%2d,%2d,%2d]\n",
                             buf0[2*N0+0], buf0[2*N0+1], buf0[2*N0+2],
                             buf0[3*N0+0], buf0[3*N0+1], buf0[3*N0+2]);
    std::cerr << "--- results of send/recv in dir. 1 ---" << std::endl;
    std::cerr << tfm::format("send lower => [%2d,%2d,%2d]\n"
                             "send upper => [%2d,%2d,%2d]\n",
                             buf1[0*N1+0], buf1[0*N1+1], buf1[0*N1+2],
                             buf1[1*N1+0], buf1[1*N1+1], buf1[1*N1+2]);
    std::cerr << tfm::format("recv lower => [%2d,%2d,%2d]\n"
                             "recv upper => [%2d,%2d,%2d]\n",
                             buf1[2*N1+0], buf1[2*N1+1], buf1[2*N1+2],
                             buf1[3*N1+0], buf1[3*N1+1], buf1[3*N1+2]);
    std::cerr << "--- results of send/recv in dir. 2 ---" << std::endl;
    std::cerr << tfm::format("send lower => [%2d,%2d,%2d]\n"
                             "send upper => [%2d,%2d,%2d]\n",
                             buf2[0*N2+0], buf2[0*N2+1], buf2[0*N2+2],
                             buf2[1*N2+0], buf2[1*N2+1], buf2[1*N2+2]);
    std::cerr << tfm::format("recv lower => [%2d,%2d,%2d]\n"
                             "recv upper => [%2d,%2d,%2d]\n",
                             buf2[2*N2+0], buf2[2*N2+1], buf2[2*N2+2],
                             buf2[3*N2+0], buf2[3*N2+1], buf2[3*N2+2]);
  }

  // finalize
  mpiutils::finalize();

  return 0;
}

#endif

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
