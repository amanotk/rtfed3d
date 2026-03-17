// -*- C++ -*-
#ifndef _SARRAY_HPP_
#define _SARRAY_HPP_

///
/// Multidimensional Static Array Container
///
/// Author: Takanobu AMANO <amano@eps.s.u-tokyo.ac.jp>
/// $Id$
///

///
/// @class SArray1D SArray.hpp
/// @brief A 1D Static Array Container Object
///
template <class T, size_t N0>
class SArray1D
{
public:
  static const size_t ndim;
  static const size_t size;
  static const size_t dims[1];
  T data[N0];
  T* ptr;

  SArray1D() { ptr = &data[0]; }
};

#ifdef __MAIN__
template <class T, size_t N0>
const size_t SArray1D<T,N0>::ndim = 1;
template <class T, size_t N0>
const size_t SArray1D<T,N0>::size = N0;
template <class T, size_t N0>
const size_t SArray1D<T,N0>::dims[] = {N0};
#endif

///
/// @class SArray2D SArray.hpp
/// @brief A 2D Static Array Container Object
///
template <class T, size_t N0, size_t N1>
class SArray2D
{
public:
  static const size_t ndim;
  static const size_t size;
  static const size_t dims[2];
  T data[N0][N1];
  T* ptr;

  SArray2D() { ptr = &data[0][0]; }
};

#ifdef __MAIN__
template <class T, size_t N0, size_t N1>
const size_t SArray2D<T,N0,N1>::ndim = 2;
template <class T, size_t N0, size_t N1>
const size_t SArray2D<T,N0,N1>::size = N0*N1;
template <class T, size_t N0, size_t N1>
const size_t SArray2D<T,N0,N1>::dims[] = {N0, N1};
#endif

///
/// @class SArray3D SArray.hpp
/// @brief A 3D Static Array Container Object
///
template <class T, size_t N0, size_t N1, size_t N2>
class SArray3D
{
public:
  static const size_t ndim;
  static const size_t size;
  static const size_t dims[3];
  T data[N0][N1][N2];
  T* ptr;

  SArray3D() { ptr = &data[0][0][0]; }
};

#ifdef __MAIN__
template <class T, size_t N0, size_t N1, size_t N2>
const size_t SArray3D<T,N0,N1,N2>::ndim = 3;
template <class T, size_t N0, size_t N1, size_t N2>
const size_t SArray3D<T,N0,N1,N2>::size = N0*N1*N2;
template <class T, size_t N0, size_t N1, size_t N2>
const size_t SArray3D<T,N0,N1,N2>::dims[] = {N0, N1, N2};
#endif

///
/// @class SArray4D SArray.hpp
/// @brief A 4D Static Array Container Object
///
template <class T, size_t N0, size_t N1, size_t N2, size_t N3>
class SArray4D
{
public:
  static const size_t ndim;
  static const size_t size;
  static const size_t dims[4];
  T data[N0][N1][N2][N3];
  T* ptr;

  SArray4D() { ptr = &data[0][0][0][0]; }
};

#ifdef __MAIN__
template <class T, size_t N0, size_t N1, size_t N2, size_t N3>
const size_t SArray4D<T,N0,N1,N2,N3>::ndim = 4;
template <class T, size_t N0, size_t N1, size_t N2, size_t N3>
const size_t SArray4D<T,N0,N1,N2,N3>::size = N0*N1*N2*N3;
template <class T, size_t N0, size_t N1, size_t N2, size_t N3>
const size_t SArray4D<T,N0,N1,N2,N3>::dims[] = {N0, N1, N2, N3};
#endif

///
/// @class SArray5D SArray.hpp
/// @brief A 5D Static Array Container Object
///
template <class T, size_t N0, size_t N1, size_t N2, size_t N3, size_t N4>
class SArray5D
{
public:
  static const size_t ndim;
  static const size_t size;
  static const size_t dims[5];
  T data[N0][N1][N2][N3][N4];
  T* ptr;

  SArray5D() { ptr = &data[0][0][0][0][0]; }
};

#ifdef __MAIN__
template <class T, size_t N0, size_t N1, size_t N2, size_t N3, size_t N4>
const size_t SArray5D<T,N0,N1,N2,N3,N4>::ndim = 5;
template <class T, size_t N0, size_t N1, size_t N2, size_t N3, size_t N4>
const size_t SArray5D<T,N0,N1,N2,N3,N4>::size = N0*N1*N2*N3*N4;
template <class T, size_t N0, size_t N1, size_t N2, size_t N3, size_t N4>
const size_t SArray5D<T,N0,N1,N2,N3,N4>::dims[] = {N0, N1, N2, N3, N4};
#endif

///
/// @class SArray6D SArray.hpp
/// @brief A 6D Static Array Container Object
///
template <class T, size_t N0, size_t N1, size_t N2, size_t N3, size_t N4,
          size_t N5>
class SArray6D
{
public:
  static const size_t ndim;
  static const size_t size;
  static const size_t dims[6];
  T data[N0][N1][N2][N3][N4][N5];
  T* ptr;

  SArray6D() { ptr = &data[0][0][0][0][0][0]; }
};

#ifdef __MAIN__
template <class T, size_t N0, size_t N1, size_t N2, size_t N3, size_t N4,
          size_t N5>
const size_t SArray6D<T,N0,N1,N2,N3,N4,N5>::ndim = 6;
template <class T, size_t N0, size_t N1, size_t N2, size_t N3, size_t N4,
          size_t N5>
const size_t SArray6D<T,N0,N1,N2,N3,N4,N5>::size = N0*N1*N2*N3*N4*N5;
template <class T, size_t N0, size_t N1, size_t N2, size_t N3, size_t N4,
          size_t N5>
const size_t SArray6D<T,N0,N1,N2,N3,N4,N5>::dims[] = {N0, N1, N2, N3, N4, N5};
#endif

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
