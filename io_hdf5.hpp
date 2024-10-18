// -*- C++ -*-
#ifndef _IO_HDF5_HPP_
#define _IO_HDF5_HPP_

///
/// @brief HDF5 I/O
///
/// $Id: io_hdf5.hpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
#include "io.hpp"
#include "hdf5.h"
#include "hdf5util.hpp"

///
/// HDF5 I/O
///
class HDF5IO : public BaseIO
{

public:
  HDF5IO() {}
  virtual ~HDF5IO() {}

  virtual void write_parameter(Global &g, std::string filename);
  virtual void write_field(Global &g, std::string filename,
                           bool init=false);
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
