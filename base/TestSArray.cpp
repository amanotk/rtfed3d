// -*- C++ -*-

///
/// @file TestSArray.cpp
/// @brief Test code for SArray?D classes
///
/// This code demonstrates how to use SArray?D object.
///
/// Author: Takanobu AMANO <amanot@stelab.nagoya-u.ac.jp>
/// $Id$
///
#define __MAIN__
#include "boost/format.hpp"
#include "SArray.hpp"
#include "MersenneTwister.hpp"

using namespace std;
static MersenneTwister mt;

// return uniform integer random number in the range [min, max]
int rand(int min, int max)
{
  double r = mt.rand();
  return static_cast<int>(r * (max - min)) + min;
}

int main()
{
  { // 1D array
    const int N1 = 10;
    SArray1D<int,N1> a1;

    cout << boost::format("\n"
                          "*** Testing SArray1D ...\n"
                          "    - ndim = %5d\n"
                          "    - size = %5d\n"
                          "    - dims = {%3d}\n"
                          "    ===> ")
      % a1.ndim % a1.size % a1.dims[0];

    for(int i1=0; i1 < N1 ;i1++) {
      a1.data[i1] = rand(0, 100);
    }

    int *ptr = &a1.data[0];
    bool status = true;
    for(int i1=0; i1 < N1 ;i1++) {
      if( a1.data[i1] != ptr[i1] )
        status = false;
    }
    if( status ) {
      cout << "works fine !" << endl;
    } else {
      cout << "does not work !" << endl;
    }
  }

  { // 2D array
    const int N1 = 4;
    const int N2 = 5;
    SArray2D<int,N1,N2> a2;

    cout << boost::format("\n"
                          "*** Testing SArray2D ...\n"
                          "    - ndim = %5d\n"
                          "    - size = %5d\n"
                          "    - dims = {%3d,%3d}\n"
                          "    ===> ")
      % a2.ndim % a2.size % a2.dims[0] % a2.dims[1];

    for(int i1=0; i1 < N1 ;i1++) {
      for(int i2=0; i2 < N2 ;i2++) {
        a2.data[i1][i2] = rand(0, 100);
      }
    }

    int *ptr = &a2.data[0][0];
    bool status = true;
    for(int i1=0; i1 < N1; i1++) {
      for(int i2=0; i2 < N2 ;i2++) {
        if( a2.data[i1][i2] != ptr[i1*N2+i2] )
          status = false;
      }
    }
    if( status ) {
      cout << "works fine !" << endl;
    } else {
      cout << "does not work !" << endl;
    }
  }

  { // 3D array
    const int N1 = 3;
    const int N2 = 2;
    const int N3 = 4;
    SArray3D<int,N1,N2,N3> a3;

    cout << boost::format("\n"
                          "*** Testing SArray3D ...\n"
                          "    - ndim = %5d\n"
                          "    - size = %5d\n"
                          "    - dims = {%3d,%3d,%3d}\n"
                          "    ===> ")
      % a3.ndim % a3.size % a3.dims[0] % a3.dims[1] % a3.dims[2];

    for(int i1=0; i1 < N1 ;i1++) {
      for(int i2=0; i2 < N2 ;i2++) {
        for(int i3=0; i3 < N3 ;i3++) {
          a3.data[i1][i2][i3] = rand(0, 100);
        }
      }
    }

    int *ptr = &a3.data[0][0][0];
    bool status = true;
    for(int i1=0; i1 < N1; i1++) {
      for(int i2=0; i2 < N2 ;i2++) {
        for(int i3=0; i3 < N3 ;i3++) {
          if( a3.data[i1][i2][i3] != ptr[(i1*N2+i2)*N3+i3] )
            status = false;
        }
      }
    }
    if( status ) {
      cout << "works fine !" << endl;
    } else {
      cout << "does not work !" << endl;
    }
  }

  { // 4D array
    const int N1 = 2;
    const int N2 = 3;
    const int N3 = 4;
    const int N4 = 2;
    SArray4D<int,N1,N2,N3,N4> a4;

    cout << boost::format("\n"
                          "*** Testing SArray4D ...\n"
                          "    - ndim = %5d\n"
                          "    - size = %5d\n"
                          "    - dims = {%3d,%3d,%3d,%3d}\n"
                          "    ===> ")
      % a4.ndim % a4.size % a4.dims[0] % a4.dims[1] % a4.dims[2] % a4.dims[3];

    for(int i1=0; i1 < N1 ;i1++) {
      for(int i2=0; i2 < N2 ;i2++) {
        for(int i3=0; i3 < N3 ;i3++) {
          for(int i4=0; i4 < N4 ;i4++) {
            a4.data[i1][i2][i3][i4] = rand(0, 100);
          }
        }
      }
    }

    int *ptr = &a4.data[0][0][0][0];
    bool status = true;
    for(int i1=0; i1 < N1; i1++) {
      for(int i2=0; i2 < N2 ;i2++) {
        for(int i3=0; i3 < N3 ;i3++) {
          for(int i4=0; i4 < N4 ;i4++) {
            if( a4.data[i1][i2][i3][i4] != ptr[((i1*N2+i2)*N3+i3)*N4+i4] )
              status = false;
          }
        }
      }
    }
    if( status ) {
      cout << "works fine !" << endl;
    } else {
      cout << "does not work !" << endl;
    }
  }

  { // 5D array
    const int N1 = 2;
    const int N2 = 3;
    const int N3 = 2;
    const int N4 = 3;
    const int N5 = 4;
    SArray5D<int,N1,N2,N3,N4,N5> a5;

    cout << boost::format("\n"
                          "*** Testing SArray5D ...\n"
                          "    - ndim = %5d\n"
                          "    - size = %5d\n"
                          "    - dims = {%3d,%3d,%3d,%3d,%3d}\n"
                          "    ===> ")
      % a5.ndim % a5.size % a5.dims[0] % a5.dims[1] % a5.dims[2] % a5.dims[3]
      % a5.dims[4];

    for(int i1=0; i1 < N1 ;i1++) {
      for(int i2=0; i2 < N2 ;i2++) {
        for(int i3=0; i3 < N3 ;i3++) {
          for(int i4=0; i4 < N4 ;i4++) {
            for(int i5=0; i5 < N5 ;i5++) {
              a5.data[i1][i2][i3][i4][i5] = rand(0, 100);
            }
          }
        }
      }
    }

    int *ptr = &a5.data[0][0][0][0][0];
    bool status = true;
    for(int i1=0; i1 < N1; i1++) {
      for(int i2=0; i2 < N2 ;i2++) {
        for(int i3=0; i3 < N3 ;i3++) {
          for(int i4=0; i4 < N4 ;i4++) {
            for(int i5=0; i5 < N5 ;i5++) {
              if( a5.data[i1][i2][i3][i4][i5] !=
                  ptr[(((i1*N2+i2)*N3+i3)*N4+i4)*N5+i5] )
                status = false;
            }
          }
        }
      }
    }
    if( status ) {
      cout << "works fine !" << endl;
    } else {
      cout << "does not work !" << endl;
    }
  }

  { // 6D array
    const int N1 = 4;
    const int N2 = 2;
    const int N3 = 3;
    const int N4 = 3;
    const int N5 = 2;
    const int N6 = 3;
    SArray6D<int,N1,N2,N3,N4,N5,N6> a6;

    cout << boost::format("\n"
                          "*** Testing SArray6D ...\n"
                          "    - ndim = %5d\n"
                          "    - size = %5d\n"
                          "    - dims = {%3d,%3d,%3d,%3d,%3d,%3d}\n"
                          "    ===> ")
      % a6.ndim % a6.size % a6.dims[0] % a6.dims[1] % a6.dims[2] % a6.dims[3]
      % a6.dims[4] % a6.dims[5];

    for(int i1=0; i1 < N1 ;i1++) {
      for(int i2=0; i2 < N2 ;i2++) {
        for(int i3=0; i3 < N3 ;i3++) {
          for(int i4=0; i4 < N4 ;i4++) {
            for(int i5=0; i5 < N5 ;i5++) {
              for(int i6=0; i6 < N6 ;i6++) {
                a6.data[i1][i2][i3][i4][i5][i6] = rand(0, 100);
              }
            }
          }
        }
      }
    }

    int *ptr = &a6.data[0][0][0][0][0][0];
    bool status = true;
    for(int i1=0; i1 < N1; i1++) {
      for(int i2=0; i2 < N2 ;i2++) {
        for(int i3=0; i3 < N3 ;i3++) {
          for(int i4=0; i4 < N4 ;i4++) {
            for(int i5=0; i5 < N5 ;i5++) {
              for(int i6=0; i6 < N6 ;i6++) {
                if( a6.data[i1][i2][i3][i4][i5][i6] !=
                    ptr[((((i1*N2+i2)*N3+i3)*N4+i4)*N5+i5)*N6+i6] )
                  status = false;
              }
            }
          }
        }
      }
    }
    if( status ) {
      cout << "works fine !" << endl;
    } else {
      cout << "does not work !" << endl;
    }
  }

  return 0;
}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
