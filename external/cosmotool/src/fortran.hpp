/*+
This is CosmoTool (./src/fortran.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __COSMO_FORTRAN_HPP
#define __COSMO_FORTRAN_HPP

#include <string>
#include <stdint.h>
#include <inttypes.h>
#include <iostream>
#include "config.hpp"

namespace CosmoTool
{
  class InvalidUnformattedAccess : public Exception
  {
  public:
    InvalidUnformattedAccess()
       : Exception("Invalid unformatted fortran file format") {}
  };

  class FortranTypes
  {
  public:
    enum Ordering {
      LittleEndian, BigEndian
    };

    enum CheckpointSize {
      Check_32bits, Check_64bits
    };
  };

  class UnformattedRead: public FortranTypes
  {
  public:

    UnformattedRead(const std::string& fname);
    UnformattedRead(const char *fname);
    ~UnformattedRead();

    // Todo implement primitive description
    void setOrdering(Ordering o);
    void setCheckpointSize(CheckpointSize cs);

    uint64_t getBlockSize() const { return checkPointRef; }
    
    void beginCheckpoint(bool bufferRecord = false);
    void endCheckpoint(bool autodrop = false);
    
    double readReal64();
    float readReal32();
    uint32_t readUint32();
    int32_t readInt32();
    int64_t readInt64();

    void skip(int64_t off);

    int64_t position() const;
    void seek(int64_t pos);

    void readOrderedBuffer(void *buffer, int size);
  protected:
    bool swapOrdering;
    CheckpointSize cSize;
    uint64_t checkPointRef;
    uint64_t checkPointAccum;
    std::ifstream *f;
    uint8_t *recordBuffer;

  };

 class UnformattedWrite: public FortranTypes
  {
  public:

    UnformattedWrite(const std::string& fname);
    UnformattedWrite(const char *fname);
    ~UnformattedWrite();

    // Todo implement primitive description
    void setOrdering(Ordering o);
    void setCheckpointSize(CheckpointSize cs);
    
    void beginCheckpoint();
    void endCheckpoint();
    
    void writeReal64(double d);
    void writeReal32(float f);
    void writeInt32(int32_t i);
    void writeInt64(int64_t i);
    void writeInt8(int8_t c);

    void writeOrderedBuffer(void *buffer, int size);
  protected:
    bool swapOrdering;
    CheckpointSize cSize;
    std::streamoff checkPointRef;
    uint64_t checkPointAccum;
    std::ofstream *f;

  };

};

#endif
