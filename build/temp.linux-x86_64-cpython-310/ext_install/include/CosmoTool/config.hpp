/*+
This is CosmoTool (./src/config.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

guilhem.lavaux@gmail.com

This software is a computer program whose purpose is to provide a toolbox for cosmological
data analysis (e.g. filters, generalized Fourier transforms, power spectra, ...)

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
+*/

#ifndef __COSMOTOOL_CONFIG_HPP
#define __COSMOTOOL_CONFIG_HPP

#include <string>
#include <stdint.h>
#include <exception>
#include <stdexcept>
#include <cstring>

namespace CosmoTool
{

#define NUMDIMS 3
#define NUMCUBES 8

  /**
   * Base type to specity at what precision we
   * must achieve computations.
   */
  typedef double ComputePrecision;
  /**
   * Coordinate type (should be a 3-array).
   */
  typedef double Coordinates[NUMDIMS];

  /* 
   * Single precision coordinates.
   */
  typedef float FCoordinates[NUMDIMS];
  
  /**
   * This function is used whenever one needs a general
   * conversion between mass and luminosity (or the opposite).
   * It should take a "mass" (or luminosity) in input, a unit is
   * given to convert this mass into solar units. The output should
   * be the "luminosity" (or mass), in solar units.
   */
  typedef double (*BiasFunction)(double mass, double unit);
  
  /**
   * Function to copy the coordinates "a" into "b".
   */
  inline void copyCoordinates(const Coordinates& a, Coordinates& b)
  {
    memcpy(b, a, sizeof(a));
  }

  /**
   * Base exception class for all exceptions handled by
   * this library.
   */
  class Exception : public std::runtime_error
  {
  public:
    Exception(const std::string& mess)
      : std::runtime_error(mess), msg(mess), msgok(true) {}
    Exception()
      : std::runtime_error("No message"), msgok(false) {}

    virtual ~Exception() throw () {}

    const char *getMessage() const { return msgok ? msg.c_str() : "No message"; };
    virtual const char *what() const throw () { return msgok ? msg.c_str() : "What 'what' ?"; };

  private:
    std::string msg;
    bool msgok;
  };

  /**
   * Exception raised when an invalid argument has been
   * passed to a function of the library.
   */
  class InvalidArgumentException : public Exception
  {
  public:
    InvalidArgumentException(const std::string& mess)
      : Exception(mess) {}
    InvalidArgumentException()
      : Exception() {}
  };

  /**
   */
  class InvalidRangeException : public Exception
  {
  public:
    InvalidRangeException(const std::string& mess)
      : Exception(mess) {}
    InvalidRangeException()
      : Exception() {}
  };
  
  /**
   */
  class NoSuchFileException : public Exception
  {
  public:
    NoSuchFileException(const std::string& mess)
      : Exception(mess) {}
    NoSuchFileException()
      : Exception() {}
  };

  /**
   */
  class InvalidFileFormatException : public Exception
  {
  public:
    InvalidFileFormatException(const std::string& mess)
      : Exception(mess) {}
    InvalidFileFormatException()
      : Exception() {}
  };

  class EndOfFileException: public Exception
  {
  public:
    EndOfFileException(const std::string& mess)
      : Exception(mess) {}
    EndOfFileException()
      : Exception() {}
  };

  class FilesystemFullException: public Exception
  {
  public:
    FilesystemFullException(const std::string& mess)
      : Exception(mess) {}
    FilesystemFullException()
      : Exception() {}
  };
}

#endif
