// -*- C++ -*-
#ifndef _MERSENNE_TWISTER_HPP_
#define _MERSENNE_TWISTER_HPP_

///
/// Random Number Generator
///
/// Author: Takanobu AMANO <amano@eps.s.u-tokyo.ac.jp>
/// $Id$
///
#include <iostream>
#include <cmath>

///
/// @class MersenneTwister MersenneTwister.hpp
/// @brief Uniform Random Number Generator with MersenneTwister Algorithm
///
/// This code is just a port of original source code written in C
/// by Takuji Nishimura and Makoto Matsumoto.
/// Static variables of original code are replaced by member variables.
/// @see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
///
class MersenneTwister
{
private:
  // local constants
  enum {
    N = 624,
    M = 397,
    MATRIX_A   = 0x9908b0dfUL,
    UPPER_MASK = 0x80000000UL,
    LOWER_MASK = 0x7fffffffUL
  };

  MersenneTwister& operator=(const MersenneTwister&); // disable assignment

protected:
  unsigned long *mt;   ///< size N array
  unsigned int   mti;  ///< should be initialized to N+1

  /// @name original routines taken from C source code
  //@{
  /// initializes mt[N] with a seed
  void init_genrand(unsigned long s)
  {
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
      mt[mti] =
        (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
      // See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
      // In the previous versions, MSBs of the seed affect
      // only MSBs of the array mt[].
      // 2002/01/09 modified by Makoto Matsumoto
      mt[mti] &= 0xffffffffUL;
      // for >32 bit machines
    }
  }

  ///
  /// initialize by an array with array-length
  /// init_key is the array for initializing keys
  /// key_length is its length
  /// slight change for C++, 2004/2/26
  ///
  void init_by_array(const unsigned long init_key[],
                     const unsigned int key_length)
  {
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
        + init_key[j] + j; // non linear
      mt[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines
      i++; j++;
      if (i>=N) { mt[0] = mt[N-1]; i=1; }
      if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
        - i; // non linear
      mt[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines
      i++;
      if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; // MSB is 1; assuring non-zero initial array
  }

  /// generates a random number on [0,0xffffffff]-interval
  unsigned long genrand_int32(void)
  {
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    // mag01[x] = x * MATRIX_A  for x=0,1

    if (mti >= N) { // generate N words at one time
      int kk;

      if (mti == N+1)   // if init_genrand() has not been called,
        init_genrand(5489UL); // a default initial seed is used

      for (kk=0;kk<N-M;kk++) {
        y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
        mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      for (;kk<N-1;kk++) {
        y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
        mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
      mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

      mti = 0;
    }

    y = mt[mti++];

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
  }

  /// generates a random number on [0,0x7fffffff]-interval
  long genrand_int31(void)
  {
    return (long)(genrand_int32()>>1);
  }

  /// generates a random number on [0,1]-real-interval
  double genrand_real1(void)
  {
    return genrand_int32()*(1.0/4294967295.0);
    // divided by 2^32-1
  }

  /// generates a random number on [0,1)-real-interval
  double genrand_real2(void)
  {
    return genrand_int32()*(1.0/4294967296.0);
    // divided by 2^32
  }

  /// generates a random number on (0,1)-real-interval
  double genrand_real3(void)
  {
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
    // divided by 2^32
  }

  /// generates a random number on [0,1) with 53-bit resolution
  double genrand_res53(void)
  {
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
  }
  //@}

  /// memory allocation
  void alloc()
  {
    mt = new unsigned long [N];
  }

public:
  /// default constructor
  MersenneTwister() : mti(N+1)
  {
    // default seed (taken from original code)
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
    unsigned long length=4;

    alloc();
    init_by_array(init, length);
  }

  /// constructor with seed
  MersenneTwister(const unsigned long seed)
  {
    alloc();
    init_genrand(seed);
  }

  /// constructor with seed
  MersenneTwister(const unsigned long init[], const int length) : mti(N+1)
  {
    // initialize random number generator
    alloc();
    init_by_array(init, length);
  }

  /// copy constructor
  MersenneTwister(const MersenneTwister& obj)
  {
    alloc();
    // copy data
    for(unsigned int i=0; i < N ;i++) {
      mt[i] = obj.mt[i];
    }
    mti = obj.mti;
  }

  /// destructor
  ~MersenneTwister()
  {
    delete [] mt;
  }

  /// set seed by integer
  void setSeed(const unsigned long seed)
  {
    init_genrand(seed);
  }

  /// set seed by array
  void setSeed(const unsigned long init[], const int length)
  {
    init_by_array(init, length);
  }

  /// @name uniform random number generation routines
  //@{
  unsigned long rand32() { return genrand_int32(); }
  long rand31() { return genrand_int31(); }
  double rand()  { return genrand_real3(); }
  double rand1() { return genrand_real1(); }
  double rand2() { return genrand_real2(); }
  double rand3() { return genrand_real3(); }
  //@}

  /// normal random number (Box-Muller method)
  double normal(const double sigma=1.0)
  {
    return M_SQRT2*sigma*std::sqrt(-std::log(rand()))*cos(2*M_PI*rand());
  }
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
