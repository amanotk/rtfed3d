// -*- C++ -*-
#ifndef _MPIUTILS_HPP_
#define _MPIUTILS_HPP_

///
/// MPI utlility module for three dimensional domain decomposition
///
/// Author: Takanobu AMANO <amano@eps.s.u-tokyo.ac.jp>
/// $Id$
///
#define MPICH_IGNORE_CXX_SEEK
#include "common.hpp"
#include "cmdline.hpp"
#include <mpi.h>
using namespace common;

//
// configuration of MPI library
//
#if   defined (OPEN_MPI)
//
// Open MPI
//
// Version 1.3.2 or below does not define MPI::INTEGER8 constant,
// which is MPI 2 standard.
//
#if OMPI_HAVE_FORTRAN_INTEGER8
namespace MPI
{
const Datatype INTEGER8 = MPI_INTEGER8;
}
#endif

#elif defined (MPICH2)
//
// MPICH 2
//
// MPICH_IGNORE_CXX_SEEK should be defined before including mpi.h
// for use of MPI 2 standard.
//

#endif

/// MPI Domain Decomposition
class MpiDomain
{
public:
  int domain[3];

  friend std::ostream& operator<<(std::ostream &os,
                                  const MpiDomain &dd)
  {
    tfm::format(os, "MpiDomain = [%4d,%4d,%4d]",
                dd.domain[0], dd.domain[1], dd.domain[2]);
    return os;
  }
};


/// MpiDomain object reader
class MpiDomainReader
{
public:
  MpiDomain operator()(const std::string &str){
    MpiDomain dom;
    std::stringstream ss(str);

    bool status = true;

    for(int i=0; i < 3 ;i++) {
      std::string nstr;

      if( !std::getline(ss, nstr, ',') ) {
        status = false;
        break;
      }

      if( nstr.empty() ) {
        status = false;
        break;
      }

      int n = std::atoi(nstr.c_str());
      if( n == 0 ) {
        status = false;
        break;
      }

      dom.domain[i] = n;
    }

    // report error and exit
    if( !status ) {
      std::cerr << "failed to parse domain" << std::endl;
      exit(-1);
    }

    return dom;
  }
};

///
/// @class mpiutils mpiutils.hpp
/// @brief utility class for domain decomposition using MPI
///
/// This class is an implementation of singleton, having only one instance
/// during execution of a program. You should first initialize() before any MPI
/// calls and finalize() at the end. These corresponds to MPI_Init() and
/// MPI_Finalize() respectively, with some additional procedures. Methods
/// implemented here are all static, meaning that you can use them like global
/// functions.
///
/// Stream buffers std::cout and std::cerr are redirected to local files for
/// each PEs and by default, these are concatenated to original buffer when
/// finalize() is called.
///
class mpiutils
{
private:
  static mpiutils *instance; ///< unique instance
  static int tag[3][2];      ///< tag for MPI communications

  // cartesian topology communicator
  MPI_Comm m_cart;

  // process, rank and neighbors
  int m_nprocess;         ///< number of process
  int m_thisrank;         ///< rank of current process
  int m_proc_dim[3];      ///< number of process in each dir.
  int m_rank_dim[3];      ///< rank of current process in eash dir.
  int m_nb_dim[3][2];     ///< neighbor process
  int m_coord[3];         ///< carteisan topology coordinate

  // for stdout/stderr
  bool            m_concat; ///< flag for concatenate cerr/cout
  std::string     m_outf;   ///< dummy standard output file
  std::string     m_errf;   ///< dummy standard error file
  std::ofstream  *m_out;    ///< dummy standard output
  std::ofstream  *m_err;    ///< dummy standard error
  std::streambuf *m_errbuf; ///< buffer of original cerr
  std::streambuf *m_outbuf; ///< buffer of original cout

  // remain undefined
  mpiutils();
  mpiutils(const mpiutils &);
  mpiutils& operator=(const mpiutils &);

  /// constructor
  mpiutils(int *argc, char*** argv, int period[3], bool concat)
  {
    int domain[3];
    get_domain(*argc, *argv, domain);

    // initialize
    MPI_Init(argc, argv);

    MPI_Comm_size(MPI_COMM_WORLD, &m_nprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_thisrank);
    m_concat   = concat;

    // check consistency between argument and nprocess
    if( domain[0]*domain[1]*domain[2] != m_nprocess ) {
      if( m_thisrank == 0 ) {
        std::cerr <<
          tfm::format("Error in mpiutils: invalid number of PEs (=%4d)\n"
                      "current domain decomposition = "
                      "[%3d, %3d, %3d] ===> expected PEs = %4d\n",
                      m_nprocess, domain[0], domain[1], domain[2],
                      (domain[0]*domain[1]*domain[2]));
      }
      // exit with error
      MPI_Finalize();
      exit(-1);
    }

    // three-dimensional domain decomposition
    {
      m_proc_dim[0] = domain[0];
      m_proc_dim[1] = domain[1];
      m_proc_dim[2] = domain[2];
      // cartesian topology
      MPI_Cart_create(MPI_COMM_WORLD, 3, domain, period, 0, &m_cart);
      MPI_Cart_get(m_cart, 3, domain, period, m_rank_dim);
      // set neighbors
      MPI_Cart_shift(m_cart, 0, +1, &m_nb_dim[0][0], &m_nb_dim[0][1]);
      MPI_Cart_shift(m_cart, 1, +1, &m_nb_dim[1][0], &m_nb_dim[1][1]);
      MPI_Cart_shift(m_cart, 2, +1, &m_nb_dim[2][0], &m_nb_dim[2][1]);
      // coordinate
      MPI_Cart_coords(m_cart, m_thisrank, 3, m_coord);
    }

    // open dummy standard error stream
    m_errf   = tfm::format("%s_PE%04d.stderr", argv[0], m_thisrank);
    m_err    = new std::ofstream(m_errf.c_str());
    m_errbuf = std::cerr.rdbuf(m_err->rdbuf());

    // open dummy standard output stream
    m_outf   = tfm::format("%s_PE%04d.stdout", argv[0], m_thisrank);
    m_out    = new std::ofstream(m_outf.c_str());
    m_outbuf = std::cout.rdbuf(m_out->rdbuf());
  }

  /// destructor
  ~mpiutils()
  {
    // close dummy stndard output
    m_out->flush();
    m_out->close();
    std::cout.rdbuf(m_outbuf);
    delete m_out;

    // close dummy standard error
    m_err->flush();
    m_err->close();
    std::cerr.rdbuf(m_errbuf);
    delete m_err;

    if( m_concat ) {
      // concatenate output stream
      concatenate_stream(m_outf, std::cout, 0, "stdout");
      concatenate_stream(m_errf, std::cerr, 1, "stderr");
    }

    // finalize
    MPI_Finalize();

    // set null
    instance = 0;
  }

  // parse command line arguments to get domain decomposition
  void get_domain(int argc, char**argv, int domain[3])
  {
    cmdline::parser parser;

    parser.add<MpiDomain>("domain", 'd', "MPI domain decomposition",
                       true, MpiDomain(), MpiDomainReader());
    parser.parse_check(argc, argv);

    MpiDomain d = parser.get<MpiDomain>("domain");
    domain[0] = d.domain[0];
    domain[1] = d.domain[1];
    domain[2] = d.domain[2];
  }

  /// concatenate sorted stream
  void concatenate_stream(std::string &filename, std::ostream &dst,
                          const int tag, const char *label)
  {
    MPI_Request r;
    int dummy = 0;
    int nprocess = instance->m_nprocess;
    int thisrank = instance->m_thisrank;

    // recieve from the previous
    if( thisrank > 0 ) {
      MPI_Irecv(&dummy, 1, MPI_INT, thisrank-1, tag, MPI_COMM_WORLD, &r);
      MPI_Wait(&r, MPI_STATUS_IGNORE);
    }

    // output
    std::ofstream f(filename.c_str(), std::ios::in);
    dst << tfm::format("--- begin %s from PE =%4d ---\n", label, thisrank);
    if( f.rdbuf()->in_avail() != 0 ) dst << f.rdbuf();
    dst << tfm::format("--- end   %s from PE =%4d ---\n", label, thisrank);
    dst.flush();

    // remove file
    std::remove(filename.c_str());

    // send to the next
    if( thisrank == nprocess-1 ) return;
    MPI_Isend(&dummy, 1, MPI_INT, thisrank+1, tag, MPI_COMM_WORLD, &r);
  }

public:
  ///
  /// MPI send/recv buffer container
  ///
  class Buffer
  {
  public:
    int  size;
    int  count[4];
    char *data[4];
    MPI_Request request[4];

    Buffer(int bufsize) : size(bufsize)
    {
      for(int i=0; i < 4 ;i++) {
        data[i]  = new char [size];
        count[i] = size;
      }
    }

    ~Buffer()
    {
      for(int i=0; i < 4 ;i++) {
        delete [] data[i];
      }
    }
  };

  /// initialize MPI call
  static mpiutils* initialize(int *argc, char*** argv,
                              int period[3], bool concat=true)
  {
    if( instance == 0 ) {
      // create instance
      instance = new mpiutils(argc, argv, period, concat);
    }
    return instance;
  }

  /// finalize MPI call
  static void finalize()
  {
    if( instance != 0 ) delete instance;
  }

  /// get unique instance
  static mpiutils* getInstance()
  {
    if( instance ) return instance;
    // report error
    std::cerr << "Error in getInstance() !"
              << std::endl
              << "===> initialization method should be called in advance."
              << std::endl;
    return 0;
  }

  /// get number of process
  static int getNProcess()
  {
    return instance->m_nprocess;
  }

  /// get rank
  static int getThisRank()
  {
    return instance->m_thisrank;
  }

  /// get coordinate
  static void getCoord(int coord[3])
  {
    coord[0] = instance->m_coord[0];
    coord[1] = instance->m_coord[1];
    coord[2] = instance->m_coord[2];
  }

  /// get neighbors
  static void getNeighbors(int neighbors[3][2])
  {
    for(int dir=0; dir < 3 ;dir++)
      for(int i=0; i < 2 ;i++)
        neighbors[dir][i] = instance->m_nb_dim[dir][i];
  }

  /// get wall clock time
  static double getTime()
  {
    return MPI_Wtime();
  }

  /// get filename with PE identifier
  static std::string getFilename(std::string prefix, std::string ext)
  {
    std::string filename =
      tfm::format("%s-%03d-%03d-%03d.%s",
                  prefix,
                  instance->m_rank_dim[0],
                  instance->m_rank_dim[1],
                  instance->m_rank_dim[2],
                  ext);
    return filename;
  }

  /// show debugging information
  static void info(std::ostream &out)
  {
    out << tfm::format("\n"
                       " <<< INFO: mpiutils >>>"
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

  /// pack data; return position after pack
  static int pack(void *data, int count, MPI_Datatype dtype,
                  void *buf, int bufsize, int pos);

  /// unpack data; return position after unpack
  static int unpack(void *data, int count, MPI_Datatype dtype,
                    void *buf, int bufsize, int pos);

  /// begin boundary exchange with non-blocking send/recv
  static void bc_exchange_begin(void *buf0, void *buf1, void *buf2,
                                int dsize, int count[3],
                                MPI_Request req[12]);

  /// begin directional boundary exchange with non-blocking send/recv
  static void bc_exchange_dir_begin(int dir, void *buf,
                                    int dsize, int count,
                                    MPI_Request req[4]);

  /// begin boundary exchange with non-blocking send/recv
  static void bc_exchange_begin(mpiutils::Buffer *buf0,
                                mpiutils::Buffer *buf1,
                                mpiutils::Buffer *buf2);

  /// begin directional boundary exchange with non-blocking send/recv
  static void bc_exchange_dir_begin(int dir, mpiutils::Buffer *buf);

  /// wait requests
  static void wait(MPI_Request req[], int n);
};


// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
