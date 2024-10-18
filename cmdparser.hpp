// -*- C++ -*-
#ifndef _CMDPARSER_HPP_
#define _CMDPARSER_HPP_

///
/// @brief Common Command-line Parser
///
/// $Id: cmdparser.hpp,v 8d17220e7571 2015/09/02 13:32:37 amano $
///
#include "cmdline.hpp"

///
/// Common command-line Parser class
///
class CmdParser : public cmdline::parser
{
public:
  CmdParser(std::string config="")
  {
    add_default(config);
  }

  void add_default(std::string config)
  {
    this->add<std::string>("config", 'c', "configuration file",
                            false, config);
    this->add<int>("xdomain", 'x', "#PE in x direction", false, 0);
    this->add<int>("ydomain", 'y', "#PE in y direction", false, 0);
    this->add<int>("zdomain", 'z', "#PE in z direction", false, 0);
  }
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
