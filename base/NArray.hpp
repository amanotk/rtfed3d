// -*- C++ -*-
#ifndef _NARRAY_HPP_
#define _NARRAY_HPP_

///
/// Multidimensional Array Container
///
/// Author: Takanobu AMANO <amanot@stelab.nagoya-u.ac.jp>
/// $Id$
///
#include "config.hpp"

///
/// @class NArray NArray.hpp
/// @brief A Multidimensional Array Container Object
///
/// A multidimensional array object with automatic memory management.
/// This object always provides a contiguous memory block for multidimensional
/// array. This object allows direct access to its internal data for
/// high performance computing. It is the users responsibility to treat
/// internal pointers in appropriate ways.
///
/// @todo
/// - add support for the data pointer referring to non-internal pointer
///
template <class T, int Rank> class NArray;

/// @brief 1-dimensional array (partial specialization)
template <class T>
class NArray<T,1>
{
private:
  typedef NArray<T,1> T_array;

  // convert pointer to 1D array
  void reshape()
  {
    array = data;
  }

  // remain undefined
  //@{
  NArray& operator=(const T_array &array);
  NArray(const T_array &array);
  //@}

public:
  T* RESTRICT data;
  T* RESTRICT array;
  uint64 shape[1];
  uint64 stride[1];

  // constructor
  NArray(const uint64 n1)
  {
    shape[0] = n1;
    stride[0] = 1;

    if( shape[0] == 0 ) {
      data  = 0;
      array = 0;
    } else {
      data = new T [ shape[0] ];
      reshape();
    }
  }

  // destructor
  ~NArray()
  {
    if( data != 0 ) {
      delete [] data;
      data  = 0;
      array = 0;
    }
  }

  /// get total number of element
  uint64 getSize()
  {
    return shape[0];
  }

  /// access operator
  //{
  const T& RESTRICT operator()(int i1) const
  {
    return data[i1*stride[0]];
  }
  T& RESTRICT operator()(int i1)
  {
    return data[i1*stride[0]];
  }
  //}
};

/// @brief 2-dimensional array (partial specialization)
template <class T>
class NArray<T,2>
{
private:
  typedef NArray<T,2> T_array;

  // convert pointer to 2D array
  void reshape()
  {
    array = new T* [shape[0]];

    for(uint64 i=0; i < shape[0] ;i++) {
      array[i] = data + i*shape[1];
    }
  }

  // remain undefined
  //@{
  NArray& operator=(const T_array &array);
  NArray(const T_array &array);
  //@}

public:
  T*  RESTRICT data;
  T** RESTRICT array;
  uint64 shape[2];
  uint64 stride[2];

  // constructor
  NArray(const uint64 n1, const uint64 n2)
  {
    uint64 size = n1*n2;
    shape[0] = n1;
    shape[1] = n2;
    stride[0] = n2;
    stride[1] = 1;

    if( size == 0 ) {
      data  = 0;
      array = 0;
    } else {
      data = new T [ size ];
      reshape();
    }
  }

  // destructor
  ~NArray()
  {
    if( data != 0 ) {
      delete [] data;
      delete [] array;
      data  = 0;
      array = 0;
    }
  }

  /// get total number of element
  uint64 getSize()
  {
    return shape[0]*shape[1];
  }

  /// access operator
  //{
  const T& RESTRICT operator()(int i1, int i2) const
  {
    return data[i1*stride[0] + i2*stride[1]];
  }
  T& RESTRICT operator()(int i1, int i2)
  {
    return data[i1*stride[0] + i2*stride[1]];
  }
  //}
};

/// @brief 3-dimensional array (partial specialization)
template <class T>
class NArray<T,3>
{
private:
  typedef NArray<T,3> T_array;

  // convert pointer to 3D array
  void reshape()
  {
    array = new T** [shape[0]];

    for(uint64 i=0; i < shape[0] ;i++) {
      array[i] = new T* [shape[1]];
      for(uint64 j=0; j < shape[1] ;j++) {
        array[i][j] = data + i*shape[1]*shape[2] + j*shape[2];
      }
    }
  }

  // remain undefined
  //@{
  NArray& operator=(const T_array &array);
  NArray(const T_array &array);
  //@}

public:
  T*   RESTRICT data;
  T*** RESTRICT array;
  uint64 shape[3];
  uint64 stride[3];

  // constructor
  NArray(const uint64 n1, const uint64 n2, const uint64 n3)
  {
    uint64 size = n1*n2*n3;
    shape[0] = n1;
    shape[1] = n2;
    shape[2] = n3;
    stride[0] = n3*n2;
    stride[1] = n3;
    stride[2] = 1;

    if( size == 0 ) {
      data  = 0;
      array = 0;
    } else {
      data = new T [ size ];
      reshape();
    }
  }

  // destructor
  ~NArray()
  {
    if( data != 0 ) {
      // delete data
      delete [] data;
      // delete array
      for(uint64 i=0; i < shape[0] ;i++) {
        delete [] array[i];
      }
      delete [] array;
      // clear
      data  = 0;
      array = 0;
    }
  }

  /// get total number of element
  uint64 getSize()
  {
    return shape[0]*shape[1]*shape[2];
  }

  /// access operator
  //{
  const T& RESTRICT operator()(int i1, int i2, int i3) const
  {
    return data[i1*stride[0] + i2*stride[1] + i3*stride[2]];
  }
  T& RESTRICT operator()(int i1, int i2, int i3)
  {
    return data[i1*stride[0] + i2*stride[1] + i3*stride[2]];
  }
  //}
};

/// @brief 4-dimensional array (partial specialization)
template <class T>
class NArray<T,4>
{
private:
  typedef NArray<T,4> T_array;

  // convert pointer to 4D array
  void reshape()
  {
    array = new T*** [shape[0]];

    for(uint64 i=0; i < shape[0] ;i++) {
      array[i] = new T** [shape[1]];
      for(uint64 j=0; j < shape[1] ;j++) {
        array[i][j] = new T* [shape[2]];
        for(uint64 k=0; k < shape[2] ; k++) {
          array[i][j][k] = data +
            i*shape[1]*shape[2]*shape[3] +
            j*shape[2]*shape[3] +
            k*shape[3];
        }
      }
    }
  }

  // remain undefined
  //@{
  NArray& operator=(const T_array &array);
  NArray(const T_array &array);
  //@}

public:
  T*    RESTRICT data;
  T**** RESTRICT array;
  uint64 shape[4];
  uint64 stride[4];

  // constructor
  NArray(const uint64 n1, const uint64 n2, const uint64 n3,
         const uint64 n4)
  {
    uint64 size = n1*n2*n3*n4;
    shape[0] = n1;
    shape[1] = n2;
    shape[2] = n3;
    shape[3] = n4;
    stride[0] = n4*n3*n2;
    stride[1] = n4*n3;
    stride[2] = n4;
    stride[3] = 1;

    if( size == 0 ) {
      data  = 0;
      array = 0;
    } else {
      data = new T [ size ];
      reshape();
    }
  }

  // destructor
  ~NArray()
  {
    if( data != 0 ) {
      // delte data
      delete [] data;
      // delte array
      for(uint64 i=0; i < shape[0] ;i++) {
        for(uint64 j=0; j < shape[1] ;j++) {
          delete [] array[i][j];
        }
        delete [] array[i];
      }
      delete [] array;
      // clear
      data = 0;
      array = 0;
    }
  }

  /// get total number of element
  uint64 getSize()
  {
    return shape[0]*shape[1]*shape[2]*shape[3];
  }

  /// access operator
  //{
  const T& RESTRICT operator()(int i1, int i2, int i3,
                               int i4) const
  {
    return data[i1*stride[0] + i2*stride[1] + i3*stride[2] +
                i4*stride[3]];
  }
  T& RESTRICT operator()(int i1, int i2, int i3,
                         int i4)
  {
    return data[i1*stride[0] + i2*stride[1] + i3*stride[2] +
                i4*stride[3]];
  }
  //}
};

/// @brief 5-dimensional array (partial specialization)
template <class T>
class NArray<T,5>
{
private:
  typedef NArray<T,5> T_array;

  // convert pointer to 5D array
  void reshape()
  {
    array = new T**** [shape[0]];

    for(uint64 i=0; i < shape[0] ;i++) {
      array[i] = new T*** [shape[1]];
      for(uint64 j=0; j < shape[1] ;j++) {
        array[i][j] = new T** [shape[2]];
        for(uint64 k=0; k < shape[2] ; k++) {
          array[i][j][k] = new T* [shape[3]];
          for(uint64 l=0; l < shape[3] ;l++) {
            array[i][j][k][l] = data +
              i*shape[1]*shape[2]*shape[3]*shape[4] +
              j*shape[2]*shape[3]*shape[4] +
              k*shape[3]*shape[4] +
              l*shape[4];
          }
        }
      }
    }
  }

  // remain undefined
  //@{
  NArray& operator=(const T_array &array);
  NArray(const T_array &array);
  //@}

public:
  T*     RESTRICT data;
  T***** RESTRICT array;
  uint64 shape[5];
  uint64 stride[5];

  // constructor
  NArray(const uint64 n1, const uint64 n2, const uint64 n3,
         const uint64 n4, const uint64 n5)
  {
    uint64 size = n1*n2*n3*n4*n5;
    shape[0] = n1;
    shape[1] = n2;
    shape[2] = n3;
    shape[3] = n4;
    shape[4] = n5;
    stride[0] = n5*n4*n3*n2;
    stride[1] = n5*n4*n3;
    stride[2] = n5*n4;
    stride[3] = n5;
    stride[4] = 1;

    if( size == 0 ) {
      data  = 0;
      array = 0;
    } else {
      data = new T [ size ];
      reshape();
    }
  }

  // destructor
  ~NArray()
  {
    if( data != 0 ) {
      // delete data
      delete [] data;
      // delete array
      for(uint64 i=0; i < shape[0] ;i++) {
        for(uint64 j=0; j < shape[1] ;j++) {
          for(uint64 k=0; k < shape[2] ;k++) {
            delete [] array[i][j][k];
          }
          delete [] array[i][j];
        }
        delete [] array[i];
      }
      delete [] array;
      // clear
      data = 0;
      array = 0;
    }
  }

  /// get total number of element
  uint64 getSize()
  {
    return shape[0]*shape[1]*shape[2]*shape[3]*shape[4];
  }

  /// access operator
  //{
  const T& RESTRICT operator()(int i1, int i2, int i3,
                               int i4, int i5) const
  {
    return data[i1*stride[0] + i2*stride[1] + i3*stride[2] +
                i4*stride[3] + i5*stride[4]];
  }
  T& RESTRICT operator()(int i1, int i2, int i3,
                         int i4, int i5)
  {
    return data[i1*stride[0] + i2*stride[1] + i3*stride[2] +
                i4*stride[3] + i5*stride[4]];
  }
  //}
};

/// @brief 6-dimensional array (partial specialization)
template <class T>
class NArray<T,6>
{
private:
  typedef NArray<T,6> T_array;

  // convert pointer to 6D array
  void reshape()
  {
    array = new T***** [shape[0]];

    for(uint64 i=0; i < shape[0] ;i++) {
      array[i] = new T**** [shape[1]];
      for(uint64 j=0; j < shape[1] ;j++) {
        array[i][j] = new T*** [shape[2]];
        for(uint64 k=0; k < shape[2] ; k++) {
          array[i][j][k] = new T** [shape[3]];
          for(uint64 l=0; l < shape[3] ;l++) {
            array[i][j][k][l] = new T* [shape[4]];
            for(uint64 m=0; m < shape[4] ;m++) {
              array[i][j][k][l][m] = data +
                i*shape[1]*shape[2]*shape[3]*shape[4]*shape[5] +
                j*shape[2]*shape[3]*shape[4]*shape[5] +
                k*shape[3]*shape[4]*shape[5] +
                l*shape[4]*shape[5] +
                m*shape[5];
            }
          }
        }
      }
    }
  }

  // remain undefined
  //@{
  NArray& operator=(const T_array &array);
  NArray(const T_array &array);
  //@}

public:
  T*      RESTRICT data;
  T****** RESTRICT array;
  uint64 shape[6];
  uint64 stride[6];

  // constructor
  NArray(const uint64 n1, const uint64 n2, const uint64 n3,
         const uint64 n4, const uint64 n5, const uint64 n6)
  {
    uint64 size = n1*n2*n3*n4*n5*n6;
    shape[0] = n1;
    shape[1] = n2;
    shape[2] = n3;
    shape[3] = n4;
    shape[4] = n5;
    shape[5] = n6;
    stride[0] = n6*n5*n4*n3*n2;
    stride[1] = n6*n5*n4*n3;
    stride[2] = n6*n5*n4;
    stride[3] = n6*n5;
    stride[4] = n6;
    stride[5] = 1;

    if( size == 0 ) {
      data  = 0;
      array = 0;
    } else {
      data = new T [ size ];
      reshape();
    }
  }

  // destructor
  ~NArray()
  {
    if( data != 0 ) {
      // delete data
      delete [] data;
      // delete array
      for(uint64 i=0; i < shape[0] ;i++) {
        for(uint64 j=0; j < shape[1] ;j++) {
          for(uint64 k=0; k < shape[2] ;k++) {
            for(uint64 l=0; l < shape[3] ;l++) {
              delete [] array[i][j][k][l];
            }
            delete [] array[i][j][k];
          }
          delete [] array[i][j];
        }
        delete [] array[i];
      }
      delete [] array;
      // clear
      data = 0;
      array = 0;
    }
  }

  /// get total number of element
  uint64 getSize()
  {
    return shape[0]*shape[1]*shape[2]*shape[3]*shape[4]*shape[5];
  }

  /// access operator
  //{
  const T& RESTRICT operator()(int i1, int i2, int i3,
                               int i4, int i5, int i6) const
  {
    return data[i1*stride[0] + i2*stride[1] + i3*stride[2] +
                i4*stride[3] + i5*stride[4] + i6*stride[5]];
  }
  T& RESTRICT operator()(int i1, int i2, int i3,
                         int i4, int i5, int i6)
  {
    return data[i1*stride[0] + i2*stride[1] + i3*stride[2] +
                i4*stride[3] + i5*stride[4] + i6*stride[5]];
  }
  //}
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
