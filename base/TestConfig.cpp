// -*- C++ -*-

///
/// @file TestConfig.cpp
/// @brief Test configuration
///
/// Show configuration given by config.hpp
///
/// Author: Takanobu AMANO <amanot@stelab.nagoya-u.ac.jp>
/// $Id$
///
#include <iostream>
#include "config.hpp"

using namespace std;

int main()
{
  cout << "--- show type info ---" << endl;
  cout << "size of   int32 (in byte) : " << sizeof(int32) << endl;
  cout << "size of  uint32 (in byte) : " << sizeof(uint32) << endl;
  cout << "size of   int64 (in byte) : " << sizeof(int64) << endl;
  cout << "size of  uint64 (in byte) : " << sizeof(uint64) << endl;
  cout << "size of float32 (in byte) : " << sizeof(float32) << endl;
  cout << "size of float64 (in byte) : " << sizeof(float64) << endl;

  cout << "--- show format string ---" << endl;
  cout << "32bit integer : " << PRId32 << endl;
  cout << "64bit integer : " << PRId64 << endl;

  return 0;
}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
