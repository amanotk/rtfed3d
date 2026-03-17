// -*- C++ -*-

///
/// @file TestNArray.cpp
/// @brief Test code for NArray<T,Rank> class
///
/// This code demonstrates how to use NArray<T,Rank> object.
///
/// Author: Takanobu AMANO <amanot@stelab.nagoya-u.ac.jp>
/// $Id$
///
#include "boost/format.hpp"
#include "NArray.hpp"
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
    NArray<int,1> a1(N1);

    for(int i=0; i < N1 ;i++) a1.data[i] = i+1;

    cout << "----- 1D Array -----" << endl;
    int ptr = 0;
    bool status = true;
    for(int i=0; i < N1 ;i++) {
      int x = a1.array[i];
      int y = a1(i);
      int z = a1.data[ptr];
      ptr++;
      if( x != y || x != z ) status = false;
#ifdef DEBUG_NARRAY
      // show debugging output
      cout << boost::format("array[%2d] = %8d <=> %8d\n") % i % x % y;
#endif
    }
    if( status ) {
      cout << "===> works fine !" << endl;
    } else {
      cout << "===> does not work !" << endl;
    }
  }


  { // 2D array
    const int N1 = 4;
    const int N2 = 5;
    NArray<int,2> a2(N1, N2);

    for(int i=0; i < N1*N2 ;i++) a2.data[i] = i+1;

    cout << "----- 2D Array -----" << endl;
    int ptr = 0;
    bool status = true;
    for(int i=0; i < N1; i++) {
      for(int j=0; j < N2 ;j++) {
        int x = a2.array[i][j];
        int y = a2(i,j);
        int z = a2.data[ptr];
        ptr++;
        if( x != y || x != z ) status = false;
#ifdef DEBUG_NARRAY
        // show debugging output
        cout << boost::format("array[%2d][%2d] = %8d <=> %8d\n")
          % i % j % x % y;
#endif
      }
    }
    if( status ) {
      cout << "===> works fine !" << endl;
    } else {
      cout << "===> does not work !" << endl;
    }
  }

  { // 3D array
    const int N1 = 3;
    const int N2 = 2;
    const int N3 = 4;
    NArray<int,3> a3(N1, N2, N3);

    for(int i=0; i < N1*N2*N3 ;i++) a3.data[i] = i+1;

    cout << "----- 3D Array -----" << endl;
    int ptr = 0;
    bool status = true;
    for(int i=0; i < N1; i++) {
      for(int j=0; j < N2 ;j++) {
        for(int k=0; k < N3 ;k++) {
          int x = a3.array[i][j][k];
          int y = a3(i,j,k);
          int z = a3.data[ptr];
          ptr++;
          if( x != y || x != z ) status = false;
#ifdef DEBUG_NARRAY
          // show debugging output
          cout << boost::format("array[%2d][%2d][%2d] = %8d <=> %8d\n")
            % i % j % k % x % y;
#endif
        }
      }
    }
    if( status ) {
      cout << "===> works fine !" << endl;
    } else {
      cout << "===> does not work !" << endl;
    }
  }

  { // 4D array
    const int N1 = 2;
    const int N2 = 3;
    const int N3 = 4;
    const int N4 = 2;
    NArray<int,4> a4(N1, N2, N3, N4);

    for(int i=0; i < N1*N2*N3*N4 ;i++) a4.data[i] = rand(0, 100);

    cout << "----- 4D Array -----" << endl;
    int ptr = 0;
    bool status = true;
    for(int i=0; i < N1; i++) {
      for(int j=0; j < N2 ;j++) {
        for(int k=0; k < N3 ;k++) {
          for(int l=0; l < N4 ;l++) {
            int x = a4.array[i][j][k][l];
            int y = a4(i,j,k,l);
            int z = a4.data[ptr];
            ptr++;
            if( x != y || x != z ) status = false;
#ifdef DEBUG_NARRAY
            // show debugging output
            cout << boost::format("array[%2d][%2d][%2d][%2d] = %8d <=> %8d\n")
              % i % j % k % l % x % y;
#endif
          }
        }
      }
    }
    if( status ) {
      cout << "===> works fine !" << endl;
    } else {
      cout << "===> does not work !" << endl;
    }
  }

  { // 5D array
    const int N1 = 2;
    const int N2 = 3;
    const int N3 = 2;
    const int N4 = 3;
    const int N5 = 4;
    NArray<int,5> a5(N1, N2, N3, N4, N5);

    for(int i=0; i < N1*N2*N3*N4*N5 ;i++) a5.data[i] = rand(0, 100);

    cout << "----- 5D Array -----" << endl;
    int ptr = 0;
    bool status = true;
    for(int i=0; i < N1; i++) {
      for(int j=0; j < N2 ;j++) {
        for(int k=0; k < N3 ;k++) {
          for(int l=0; l < N4 ;l++) {
            for(int m=0; m < N5 ;m++) {
              int x = a5.array[i][j][k][l][m];
              int y = a5(i,j,k,l,m);
              int z = a5.data[ptr];
              ptr++;
              if( x != y || x != z ) status = false;
#ifdef DEBUG_NARRAY
              // show debugging output
              cout << boost::format("array[%2d][%2d][%2d][%2d][%2d]"
                                    " = %8d <=> %8d\n")
                % i % j % k % l % m % x % y;
#endif
            }
          }
        }
      }
    }
    if( status ) {
      cout << "===> works fine !" << endl;
    } else {
      cout << "===> does not work !" << endl;
    }
  }

  { // 6D array
    const int N1 = 4;
    const int N2 = 2;
    const int N3 = 3;
    const int N4 = 3;
    const int N5 = 2;
    const int N6 = 3;
    NArray<int,6> a6(N1, N2, N3, N4, N5, N6);

    for(int i=0; i < N1*N2*N3*N4*N5*N6 ;i++) a6.data[i] = rand(0, 100);

    cout << "----- 6D Array -----" << endl;
    int ptr = 0;
    bool status = true;
    for(int i=0; i < N1; i++) {
      for(int j=0; j < N2 ;j++) {
        for(int k=0; k < N3 ;k++) {
          for(int l=0; l < N4 ;l++) {
            for(int m=0; m < N5 ;m++) {
              for(int n=0; n < N6 ;n++) {
                int x = a6.array[i][j][k][l][m][n];
                int y = a6(i,j,k,l,m,n);
                int z = a6.data[ptr];
                ptr++;
                if( x != y || x != z ) status = false;
#ifdef DEBUG_NARRAY
                // show debugging output
                cout << boost::format("array[%2d][%2d][%2d][%2d][%2d][%2d]"
                                      " = %8d <=> %8d\n")
                  % i % j % k % l % m % n % x % y;
#endif
              }
            }
          }
        }
      }
    }
    if( status ) {
      cout << "===> works fine !" << endl;
    } else {
      cout << "===> does not work !" << endl;
    }
  }

  return 0;
}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
