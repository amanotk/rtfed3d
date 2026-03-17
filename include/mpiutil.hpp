// -*- C++ -*-
#ifndef _MPIUTIL_HPP_
#define _MPIUTIL_HPP_

///
/// MPI utlility module for three dimensional domain decomposition
///
/// Author: Takanobu AMANO <amano@eps.s.u-tokyo.ac.jp>
/// $Id: mpiutil.hpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
#define MPICH_IGNORE_CXX_SEEK
#include "common.hpp"
#include "cmdline.hpp"
#include <mpi.h>
using namespace common;

/// template for Singleton class
template <class T>
class Singleton
{
private:
  Singleton(const Singleton &);
  Singleton& operator=(const Singleton &);

protected:
  Singleton() {};
  virtual ~Singleton() {};

public:
  static T* getInstance()
  {
    static T instance;
    return &instance;
  }
};

/// stream buffer mimicking "tee" command
class teebuf : public std::streambuf
{
private:
  std::streambuf *m_sb1;
  std::streambuf *m_sb2;

  virtual int overflow(int c)
  {
    if (c == EOF) {
      return !EOF;
    } else {
      int const r1 = m_sb1->sputc(c);
      int const r2 = m_sb2->sputc(c);
      return r1 == EOF || r2 == EOF ? EOF : c;
    }
  }

  virtual int sync()
  {
    int const r1 = m_sb1->pubsync();
    int const r2 = m_sb2->pubsync();
    return r1 == 0 && r2 == 0 ? 0 : -1;
  }

public:
  teebuf(std::streambuf *sb1, std::streambuf *sb2)
    : m_sb1(sb1), m_sb2(sb2)
  {
  }
};

///
/// @class mpiutil mpiutil.hpp
/// @brief utility class for domain decomposition using MPI
///
class mpiutil : public Singleton<mpiutil>
{
  friend class Singleton<mpiutil>;

private:
  // cartesian topology communicator
  MPI_Comm m_cart;

  // number of grids
  int m_shape[3];

  // process, rank, neighbors, tag
  int m_nprocess;           ///< number of process
  int m_thisrank;           ///< rank of current process
  int m_proc_dim[3];        ///< number of process in each dir.
  int m_nb_dim[3][2];       ///< neighbor process
  int m_coord[3];           ///< carteisan topology coordinate
  int m_tag[3][2];          ///< tag
  int *m_rank[3];           ///< rank array for each dir.

  // for stdout/stderr
  std::string     m_outf;   ///< dummy standard output file
  std::string     m_errf;   ///< dummy standard error file
  std::ofstream  *m_out;    ///< dummy standard output
  std::ofstream  *m_err;    ///< dummy standard error
  std::streambuf *m_errbuf; ///< buffer of original cerr
  std::streambuf *m_outbuf; ///< buffer of original cout
  teebuf         *m_outtee; ///< buffer for replicating cout and file
  teebuf         *m_errtee; ///< buffer for replicating cerr and file

  mpiutil() {};
  ~mpiutil() {};

  // remain undefined
  mpiutil(const mpiutil &);
  mpiutil& operator=(const mpiutil &);

public:
  static void decompose(int shape[3], int domain[3])
  {
    mpiutil *instance = getInstance();

    for(int r=0; r < 3; r++) {
      domain[r] = shape[r] == 1 ? 1 : domain[r];
    }

    MPI_Dims_create(instance->m_nprocess, 3, domain);

    for(int r=0; r < 3 ;r++) {
      // check consistency
      if( shape[r] % domain[r] != 0 ) {
        std::cerr << "Error: invalid domain decomposition" << std::endl;
        std::exit(-1);
      }

      shape[r]  = shape[r] / domain[r];
    }
  }

  /// initialize MPI call
  static void initialize(int *argc, char*** argv,
                         int shape[3], int domain[3])
  {
    int period[3] = {1, 1, 1};
    initialize(argc, argv, shape, domain, period);
  }

  /// initialize MPI call
  static void initialize(int *argc, char*** argv,
                         int shape[3], int domain[3], int period[3])
  {
    mpiutil *instance = getInstance();

    // initialize
    MPI_Init(argc, argv);

    MPI_Comm_size(MPI_COMM_WORLD, &instance->m_nprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &instance->m_thisrank);

    // open dummy standard output stream
    instance->m_outf   = tfm::format("%s_PE%04d.stdout", argv[0][0],
                                     instance->m_thisrank);
    instance->m_out    = new std::ofstream(instance->m_outf.c_str());
    instance->m_outtee = new teebuf(std::cout.rdbuf(),
                                    instance->m_out->rdbuf());

    // open dummy standard error stream
    instance->m_errf   = tfm::format("%s_PE%04d.stderr", argv[0][0],
                                     instance->m_thisrank);
    instance->m_err    = new std::ofstream(instance->m_errf.c_str());
    instance->m_errtee = new teebuf(std::cerr.rdbuf(),
                                    instance->m_err->rdbuf());

    if( instance->m_thisrank == 0 ) {
      // stdout/stderr are replicated for rank==0
      instance->m_outbuf = std::cout.rdbuf(instance->m_outtee);
      instance->m_errbuf = std::cerr.rdbuf(instance->m_errtee);
    } else {
      instance->m_outbuf = std::cout.rdbuf(instance->m_out->rdbuf());
      instance->m_errbuf = std::cerr.rdbuf(instance->m_err->rdbuf());
    }

    // domain decomposition
    decompose(shape, domain);
    instance->m_shape[0] = shape[0];
    instance->m_shape[1] = shape[1];
    instance->m_shape[2] = shape[2];

    // three-dimensional domain decomposition
    instance->m_proc_dim[0] = domain[0];
    instance->m_proc_dim[1] = domain[1];
    instance->m_proc_dim[2] = domain[2];

    // cartesian topology
    MPI_Cart_create(MPI_COMM_WORLD, 3, domain, period, 0, &instance->m_cart);

    // set neighbors
    MPI_Cart_shift(instance->m_cart, 0, +1,
                   &instance->m_nb_dim[0][0], &instance->m_nb_dim[0][1]);
    MPI_Cart_shift(instance->m_cart, 1, +1,
                   &instance->m_nb_dim[1][0], &instance->m_nb_dim[1][1]);
    MPI_Cart_shift(instance->m_cart, 2, +1,
                   &instance->m_nb_dim[2][0], &instance->m_nb_dim[2][1]);

    // coordinate
    MPI_Cart_coords(instance->m_cart, instance->m_thisrank, 3,
                    instance->m_coord);

    // rank array
    for(int dir=0; dir < 3 ;dir++) {
      int coord[3] =
        { instance->m_coord[0],
          instance->m_coord[1],
          instance->m_coord[2]
        };
      // allocate
      instance->m_rank[dir] = new int [instance->m_proc_dim[dir]];
      // set rank
      for(int i=0; i < instance->m_proc_dim[dir] ;i++) {
        coord[dir] = i;
        MPI_Cart_rank(instance->m_cart, coord, &instance->m_rank[dir][i]);
      }
    }

    // tag
    instance->m_tag[0][0] = 0;
    instance->m_tag[0][1] = 1;
    instance->m_tag[1][0] = 2;
    instance->m_tag[1][1] = 3;
    instance->m_tag[2][0] = 4;
    instance->m_tag[2][1] = 5;
  }

  /// finalize MPI call
  static void finalize()
  {
    mpiutil *instance = getInstance();

    // close dummy stndard output
    instance->m_out->flush();
    instance->m_out->close();
    std::cout.rdbuf(instance->m_outbuf);
    delete instance->m_outtee;
    delete instance->m_out;

    // close dummy standard error
    instance->m_err->flush();
    instance->m_err->close();
    std::cerr.rdbuf(instance->m_errbuf);
    delete instance->m_errtee;
    delete instance->m_err;

    // delete rank array
    delete [] instance->m_rank[0];
    delete [] instance->m_rank[1];
    delete [] instance->m_rank[2];

    // finalize
    MPI_Finalize();
  }

  /// return if the lower boundary or not
  static bool isLowerBoundary(int dir)
  {
    mpiutil *instance = getInstance();
    return (instance->m_nb_dim[dir][0] == MPI_PROC_NULL);
  }

  /// return if the lower boundary or not
  static bool isUpperBoundary(int dir)
  {
    mpiutil *instance = getInstance();
    return (instance->m_nb_dim[dir][1] == MPI_PROC_NULL);
  }

  /// get global shape
  static void getGlobalShape(int shape[3])
  {
    mpiutil *instance = getInstance();
    shape[0] = instance->m_shape[0] * instance->m_proc_dim[0];
    shape[1] = instance->m_shape[1] * instance->m_proc_dim[1];
    shape[2] = instance->m_shape[2] * instance->m_proc_dim[2];
  }

  /// get global offset
  static void getGlobalOffset(int offset[3])
  {
    mpiutil *instance = getInstance();
    offset[0] = instance->m_shape[0] * instance->m_coord[0];
    offset[1] = instance->m_shape[1] * instance->m_coord[1];
    offset[2] = instance->m_shape[2] * instance->m_coord[2];
  }

  /// get local spatial range
  static void getLocalRange(int dir, float64 dh, float64 rng[3])
  {
    mpiutil *instance = getInstance();
    rng[0] = dh * instance->m_shape[dir] * instance->m_coord[dir];
    rng[1] = dh * instance->m_proc_dim[dir] + rng[0];
    rng[2] = rng[1] - rng[0];
  }

  /// get global spatial range
  static void getGlobalRange(int dir, float64 dh, float64 rng[3])
  {
    mpiutil *instance = getInstance();
    rng[0] = 0.0;
    rng[1] = dh * instance->m_shape[dir] * instance->m_proc_dim[dir];
    rng[2] = rng[1] - rng[0];
  }

  /// get number of process
  static int getNProcess()
  {
    mpiutil *instance = getInstance();
    return instance->m_nprocess;
  }

  /// get rank
  static int getThisRank()
  {
    mpiutil *instance = getInstance();
    return instance->m_thisrank;
  }

  /// get rank array
  static void getRank(int *rank[3])
  {
    mpiutil *instance = getInstance();
    rank[0] = instance->m_rank[0];
    rank[1] = instance->m_rank[1];
    rank[2] = instance->m_rank[2];
  }

  /// get decomposition
  static void getDomain(int domain[3])
  {
    mpiutil* instance = getInstance();
    domain[0] = instance->m_proc_dim[0];
    domain[1] = instance->m_proc_dim[1];
    domain[2] = instance->m_proc_dim[2];
  }

  /// get coordinate
  static void getCoord(int coord[3])
  {
    mpiutil* instance = getInstance();
    coord[0] = instance->m_coord[0];
    coord[1] = instance->m_coord[1];
    coord[2] = instance->m_coord[2];
  }

  /// get neighbor
  static void getNeighbor(int neighbor[3][2])
  {
    mpiutil* instance = getInstance();
    for(int dir=0; dir < 3 ;dir++)
      for(int i=0; i < 2 ;i++)
        neighbor[dir][i] = instance->m_nb_dim[dir][i];
  }

  /// get tags
  static void getTag(int tag[3][2])
  {
    mpiutil* instance = getInstance();
    for(int dir=0; dir < 3 ;dir++)
      for(int i=0; i < 2 ;i++)
        tag[dir][i] = instance->m_tag[dir][i];
  }

  /// get filename with PE identifier
  static std::string getFilename(std::string prefix, std::string ext)
  {
    mpiutil* instance = getInstance();
    std::string filename =
      tfm::format("%s-%03d-%03d-%03d.%s",
                  prefix,
                  instance->m_coord[0],
                  instance->m_coord[1],
                  instance->m_coord[2],
                  ext);
    return filename;
  }

  /// flush
  static void flush()
  {
    mpiutil* instance = getInstance();

    instance->m_out->flush();
    instance->m_err->flush();
  }

  /// show debugging information
  static void info(std::ostream &out)
  {
    mpiutil* instance = getInstance();

    out << tfm::format("\n"
                       " <<< INFO: mpiutil >>>"
                       "\n"
                       "Number of Process   : %4d\n"
                       "This Rank           : %4d\n"
                       "Temporary std::cout : %s\n"
                       "Temporary std::cerr : %s\n",
                       instance->m_nprocess,
                       instance->m_thisrank,
                       instance->m_outf.c_str(),
                       instance->m_errf.c_str());

    out << "Info for domain decomposition:\n";
    for(int dir=0; dir < 3 ; dir++) {
      int cl[3], cc[3], cu[3];
      // self
      MPI_Cart_coords(instance->m_cart, instance->m_thisrank, 3, cc);
      // lower
      if ( instance->m_nb_dim[dir][0] != MPI_PROC_NULL ) {
        MPI_Cart_coords(instance->m_cart, instance->m_nb_dim[dir][0], 3, cl);
      } else {
        cl[0] = cl[1] = cl[2] = -1;
      }
      // upper
      if ( instance->m_nb_dim[dir][1] != MPI_PROC_NULL ) {
        MPI_Cart_coords(instance->m_cart, instance->m_nb_dim[dir][1], 3, cu);
      } else {
        cu[0] = cu[1] = cu[2] = -1;
      }
      out << tfm::format(" neighbor in dir %1d: "
                         "[%2d,%2d,%2d] <= [%2d,%2d,%2d] => [%2d,%2d,%2d]\n",
                         dir,
                         cl[0], cl[1], cl[2],
                         cc[0], cc[1], cc[2],
                         cu[0], cu[1], cu[2]);
    }
    out << std::endl;
  }
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
