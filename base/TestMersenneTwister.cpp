// -*- C++ -*-

///
/// @file TestMersenneTwister.cpp
/// @brief Test code for MersenneTwister class
///
/// This code demonstrates how to use MersenneTwister object.
///
/// @todo add test routine for normal random numbers
///
/// Author: Takanobu AMANO <amanot@stelab.nagoya-u.ac.jp>
/// $Id$
///
#include <iostream>
#include <string>
#include <iomanip>
#include <ctime>
#include "boost/format.hpp"
#include "config.hpp"
#include "MersenneTwister.hpp"

using namespace std;

int main()
{
  {
    // create MersenneTwister by a default constructor
    // this should always produce the same random numbers
    MersenneTwister mt;
    cout << setw(60) << setfill('-') << "" << endl;
    for(int i=0; i < 10 ;i++) {
      cout << boost::format("uniform random number : %4d ===> %10.5f\n")
        % i % mt.rand();
    }
  }

  {
    // create MersenneTwister by a constructor with seed array
    // this should also produce the same random numbers
    unsigned long init[4] = { 1, 2, 3, 4 };
    MersenneTwister mt(init, 4);
    cout << setw(60) << setfill('-') << "" << endl;
    for(int i=0; i < 10 ;i++) {
      cout << boost::format("uniform random number : %4d ===> %10.5f\n")
        % i % mt.rand();
    }
  }

  {
    // create MersenneTwister by a constructor with seed
    // the seed is obtained by the current time, this means that this produces
    // different random numbers for each execution
    MersenneTwister mt(static_cast<unsigned long>(time(0)));
    cout << setw(60) << setfill('-') << "" << endl;
    for(int i=0; i < 10 ;i++) {
      cout << boost::format("uniform random number : %4d ===> %10.5f\n")
        % i % mt.rand();
    }
  }

  return 0;
}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
