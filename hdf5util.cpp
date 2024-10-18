// -*- C++ -*-

///
/// @brief Implementation of HDF5 Utility
///
/// $Id: hdf5util.cpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
#include "hdf5util.hpp"

namespace hdf5util
{

const int MAXDIMS = 32;

void create_file(const char *filename)
{
  hid_t file, prop;

  // set property
  prop = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(prop, MPI_COMM_WORLD, MPI_INFO_NULL);

  // create a new file
  file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, prop);

  H5Pclose(prop);
  H5Fclose(file);
}

hid_t open_file(const char *filename)
{
  hid_t file, prop;

  prop = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(prop, MPI_COMM_WORLD, MPI_INFO_NULL);

  file = H5Fopen(filename, H5F_ACC_RDWR, prop);

  H5Pclose(prop);

  return file;
}

void close_file(hid_t file)
{
  H5Fclose(file);
}

void put_attribute_type(hid_t dest, const char *name, hid_t type, void *data)
{
  hid_t scalar = H5Screate(H5S_SCALAR);
  hid_t attr   = H5Acreate(dest, name, type, scalar, H5P_DEFAULT, H5P_DEFAULT);

  H5Awrite(attr, type, data);

  H5Sclose(scalar);
  H5Aclose(attr);
}

void put_attribute_type(hid_t dest, const char *name, hid_t type,
                        const int n, void *data)
{
  hsize_t dims[1] = {n};
  hid_t array = H5Screate_simple(1, dims, NULL);
  hid_t attr  = H5Acreate(dest, name, type, array, H5P_DEFAULT, H5P_DEFAULT);

  H5Awrite(attr, type, data);

  H5Sclose(array);
  H5Aclose(attr);
}

template <>
void put_attribute(hid_t dest, const char *name, int *data)
{
  put_attribute_type(dest, name, H5T_NATIVE_INT, data);
}

template <>
void put_attribute(hid_t dest, const char *name, float *data)
{
  put_attribute_type(dest, name, H5T_NATIVE_FLOAT, data);
}

template <>
void put_attribute(hid_t dest, const char *name, double *data)
{
  put_attribute_type(dest, name, H5T_NATIVE_DOUBLE, data);
}

template <>
void put_attribute(hid_t dest, const char *name, const int n, int *data)
{
  put_attribute_type(dest, name, H5T_NATIVE_INT, n, data);
}

template <>
void put_attribute(hid_t dest, const char *name, const int n, float *data)
{
  put_attribute_type(dest, name, H5T_NATIVE_FLOAT, n, data);
}

template <>
void put_attribute(hid_t dest, const char *name, const int n, double *data)
{
  put_attribute_type(dest, name, H5T_NATIVE_DOUBLE, n, data);
}

void extend_dimension(hid_t dest, const char *name)
{
  hsize_t dims[MAXDIMS], mdim[MAXDIMS];

  hid_t dataset = H5Dopen(dest, name, H5P_DEFAULT);
  hid_t space   = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(space, dims, mdim);

  dims[0]++;
  H5Dset_extent(dataset, dims);

  H5Sclose(space);
  H5Dclose(dataset);
}

void create_dataset(hid_t dest, const char *name, hid_t type,
                    const int rank, hsize_t dims[], hsize_t mdim[])
{
  // chunk size == dims
  create_dataset(dest, name, type, rank, dims, mdim, dims);
}

void create_dataset(hid_t dest, const char *name, hid_t type,
                    const int rank, hsize_t dims[], hsize_t mdim[],
                    hsize_t chunk[])
{
  hid_t space   = H5Screate_simple(rank, dims, mdim);
  hid_t plist   = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dtype   = H5Tcopy(type);

  H5Pset_chunk(plist, rank, chunk);
  H5Tset_order(dtype, H5T_ORDER_LE);

  hid_t dataset = H5Dcreate(dest, name, dtype, space,
                            H5P_DEFAULT, plist, H5P_DEFAULT);

  H5Sclose(space);
  H5Pclose(plist);
  H5Tclose(dtype);
  H5Dclose(dataset);
}

void write_dataset(hid_t dest, const char *name,
                   const int rank, hsize_t dims[], hsize_t count[],
                   hsize_t loffset[], hsize_t goffset[],
                   void *data)
{
  hid_t dataset = H5Dopen(dest, name, H5P_DEFAULT);
  hid_t dspace  = H5Dget_space(dataset);
  hid_t mspace  = H5Screate_simple(rank, dims, NULL);
  hid_t dtype   = H5Dget_type(dataset);
  hid_t plist   = H5Pcreate(H5P_DATASET_XFER);

  // destination dataspace
  H5Sselect_hyperslab(dspace, H5S_SELECT_SET, goffset, NULL, count, NULL);

  // source dataspace
  H5Sselect_hyperslab(mspace, H5S_SELECT_SET, loffset, NULL, count, NULL);

  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dataset, dtype, mspace, dspace, plist, data);

  H5Pclose(plist);
  H5Tclose(dtype);
  H5Sclose(dspace);
  H5Sclose(mspace);
  H5Dclose(dataset);
}

}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
