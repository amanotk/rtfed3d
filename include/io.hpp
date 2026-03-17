// -*- C++ -*-
#ifndef _IO_HPP_
#define _IO_HPP_

///
/// @breif Disk I/O
///
/// $Id: io.hpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
#include "global.hpp"

///
/// Base I/O class
///
class BaseIO
{
public:
  BaseIO() {}
  virtual ~BaseIO() {}

  virtual void write_field(Global &g, std::string filename,
                           bool init=false) = 0;
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
