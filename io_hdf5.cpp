// -*- C++ -*-

///
/// @brief Implementation of HDF5 I/O
///
/// $Id: io_hdf5.cpp,v 641a96a86755 2015/09/14 19:35:09 amano $
///
#include "io_hdf5.hpp"

void HDF5IO::write_parameter(Global &g, std::string filename)
{
  int offset[3], shape[3];
  hsize_t dims[3], count[3], loffset[3], goffset[3], gshape[3];

  hid_t file = hdf5util::open_file(filename.c_str());

  mpiutil::getGlobalOffset(offset);
  mpiutil::getGlobalShape(shape);

  // number of grids etc.
  hdf5util::put_attribute(file, "Nb", &g.Nb);
  hdf5util::put_attribute(file, "Nz", &shape[0]);
  hdf5util::put_attribute(file, "Ny", &shape[1]);
  hdf5util::put_attribute(file, "Nx", &shape[2]);

  // grid size
  hdf5util::put_attribute(file, "delx", &g.delx);
  hdf5util::put_attribute(file, "dely", &g.dely);
  hdf5util::put_attribute(file, "delz", &g.delz);

  // physical parameters
  hdf5util::put_attribute(file, "c", &g.c);
  hdf5util::put_attribute(file, "gam", &g.gam);
  hdf5util::put_attribute(file, "eta", &g.eta);

  // electron
  hdf5util::put_attribute(file, "qe", &g.qe);
  hdf5util::put_attribute(file, "me", &g.me);

  // proton
  hdf5util::put_attribute(file, "qp", &g.qp);
  hdf5util::put_attribute(file, "mp", &g.mp);

  // grid
  dims[0]    = g.Mz;
  dims[1]    = g.My;
  dims[2]    = g.Mx;
  count[0]   = g.Nz;
  count[1]   = g.Ny;
  count[2]   = g.Nx;
  loffset[0] = g.Lbz;
  loffset[1] = g.Lby;
  loffset[2] = g.Lbx;
  goffset[0] = offset[0];
  goffset[1] = offset[1];
  goffset[2] = offset[2];
  gshape[0]  = shape[0];
  gshape[1]  = shape[1];
  gshape[2]  = shape[2];

  // z
  hdf5util::create_dataset(file, "zig", H5T_NATIVE_DOUBLE,
                           1, &gshape[0], &gshape[0]);
  hdf5util::write_dataset(file, "zig", 1, &dims[0], &count[0],
                          &loffset[0], &goffset[0], g.zig.data());

  // y
  hdf5util::create_dataset(file, "yig", H5T_NATIVE_DOUBLE,
                           1, &gshape[1], &gshape[1]);
  hdf5util::write_dataset(file, "yig", 1, &dims[1], &count[1],
                          &loffset[1], &goffset[1], g.yig.data());

  // x
  hdf5util::create_dataset(file, "xig", H5T_NATIVE_DOUBLE,
                           1, &gshape[2], &gshape[2]);
  hdf5util::write_dataset(file, "xig", 1, &dims[2], &count[2],
                          &loffset[2], &goffset[2], g.xig.data());

  // close
  hdf5util::close_file(file);
}

void HDF5IO::write_field(Global &g, std::string filename, bool init)
{
  static int Nwc;
  hid_t file, root;
  hsize_t dims[5], mdim[5], count[5], loffset[5], goffset[5], gshape[5];

  int offset[3], shape[3];
  mpiutil::getGlobalOffset(offset);
  mpiutil::getGlobalShape(shape);

  if( init ) {
    //
    // create file and write parameters
    //
    hdf5util::create_file(filename.c_str());
    write_parameter(g, filename);
  }

  // open file and root group
  file = hdf5util::open_file(filename.c_str());
  root = H5Gopen(file, "/", 0);

  if( init ) {
    Nwc = 1;

    //
    // create dataset
    //
    dims[0] = 1;
    mdim[0] = H5S_UNLIMITED;
    hdf5util::create_dataset(root, "physt", H5T_NATIVE_DOUBLE, 1, dims, mdim);
    hdf5util::create_dataset(root, "delt", H5T_NATIVE_DOUBLE, 1, dims, mdim);

    dims[1] = shape[0];
    dims[2] = shape[1];
    dims[3] = shape[2];
    dims[4] = 6;
    for(int k=1; k < 5 ;k++) mdim[k] = dims[k];
    hdf5util::create_dataset(root, "eb", H5T_NATIVE_DOUBLE, 5, dims, mdim);

    dims[1] = shape[0];
    dims[2] = shape[1];
    dims[3] = shape[2];
    dims[4] = 10;
    for(int k=1; k < 5 ;k++) mdim[k] = dims[k];
    hdf5util::create_dataset(root, "uf", H5T_NATIVE_DOUBLE, 5, dims, mdim);
  } else {
    // extend time dimension
    hdf5util::extend_dimension(root, "physt");
    hdf5util::extend_dimension(root, "delt");
    hdf5util::extend_dimension(root, "eb");
    hdf5util::extend_dimension(root, "uf");
  }

  // time and time step
  dims[0]    = 1;
  count[0]   = 1;
  loffset[0] = 0;
  goffset[0] = Nwc-1;
  hdf5util::write_dataset(root, "physt", 1, dims, count,
                          loffset, goffset, &g.physt);
  hdf5util::write_dataset(root, "delt", 1, dims, count,
                          loffset, goffset, &g.delt);

  // write data
  dims[0]    = 1;
  dims[1]    = g.Mz;
  dims[2]    = g.My;
  dims[3]    = g.Mx;
  dims[4]    = 0; // to be overwritten
  count[0]   = 1;
  count[1]   = g.Nz;
  count[2]   = g.Ny;
  count[3]   = g.Nx;
  count[4]   = 0; // to be overwritten
  loffset[0] = 0;
  loffset[1] = g.Lbz;
  loffset[2] = g.Lby;
  loffset[3] = g.Lbx;
  loffset[4] = 0;
  goffset[0] = Nwc-1;
  goffset[1] = offset[0];
  goffset[2] = offset[1];
  goffset[3] = offset[2];
  goffset[4] = 0;

  dims[4]  = 6;
  count[4] = 6;
  hdf5util::write_dataset(root, "eb", 5, dims, count,
                          loffset, goffset, g.ebc.data());

  dims[4]  = 10;
  count[4] = 10;
  hdf5util::write_dataset(root, "uf", 5, dims, count,
                          loffset, goffset, g.uf.data());

  Nwc++;

  H5Gclose(root);
  hdf5util::close_file(file);
}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
