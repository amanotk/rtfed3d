// -*- C++ -*-
#ifndef _HDF5UTIL_HPP_
#define _HDF5UTIL_HPP_

///
/// @brief HDF5 Utility
///
/// $Id: hdf5util.hpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
#include <mpi.h>
#include <hdf5.h>

namespace hdf5util
{

void create_file(const char *filename);

hid_t open_file(const char *filename);

void close_file(hid_t file);

void put_attribute_type(hid_t dest, const char *name, hid_t type,
                        void *data);

void put_attribute_type(hid_t dest, const char *name, hid_t type,
                        const int n, void *data);

template <class T>
void put_attribute(hid_t dest, const char *name, T *data);

template <class T>
void put_attribute(hid_t dest, const char *name, const int n, T *data);

void extend_dimension(hid_t dest, const char *name);

void create_dataset(hid_t dest, const char *name, hid_t type,
                    const int rank, hsize_t dims[], hsize_t maxsh[]);

void create_dataset(hid_t dest, const char *name, hid_t type,
                    const int rank, hsize_t dims[], hsize_t mdim[],
                    hsize_t chunk[]);

void write_dataset(hid_t dest, const char *name,
                   const int rank, hsize_t dims[], hsize_t count[],
                   hsize_t loffset[], hsize_t goffset[],
                   void *data);
}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
