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

    UnformattedRead(const std::string& fname)
      throw (NoSuchFileException);
    UnformattedRead(const char *fname)
      throw (NoSuchFileException);
    ~UnformattedRead();

    // Todo implement primitive description
    void setOrdering(Ordering o);
    void setCheckpointSize(CheckpointSize cs);
    
    void beginCheckpoint()
      throw (InvalidUnformattedAccess,EndOfFileException);
    void endCheckpoint(bool autodrop = false)
      throw (InvalidUnformattedAccess);
    
    double readReal64()
      throw (InvalidUnformattedAccess);
    float readReal32()
      throw (InvalidUnformattedAccess);
    int32_t readInt32()
      throw (InvalidUnformattedAccess);
    int64_t readInt64()
      throw (InvalidUnformattedAccess);

    void skip(int64_t off)
      throw (InvalidUnformattedAccess);

  protected:
    bool swapOrdering;
    CheckpointSize cSize;
    uint64_t checkPointRef;
    uint64_t checkPointAccum;
    std::ifstream *f;

    void readOrderedBuffer(void *buffer, int size)
      throw (InvalidUnformattedAccess);
  };

 class UnformattedWrite: public FortranTypes
  {
  public:

    UnformattedWrite(const std::string& fname)
      throw (NoSuchFileException);
    UnformattedWrite(const char *fname)
      throw (NoSuchFileException);
    ~UnformattedWrite();

    // Todo implement primitive description
    void setOrdering(Ordering o);
    void setCheckpointSize(CheckpointSize cs);
    
    void beginCheckpoint()
      throw (FilesystemFullException,InvalidUnformattedAccess);
    void endCheckpoint()
      throw (FilesystemFullException,InvalidUnformattedAccess);
    
    void writeReal64(double d)
      throw (FilesystemFullException);
    void writeReal32(float f)
      throw (FilesystemFullException);
    void writeInt32(int32_t i)
      throw (FilesystemFullException);
    void writeInt64(int64_t i)
      throw (FilesystemFullException);
    void writeInt8(int8_t c)
      throw (FilesystemFullException);

    void writeOrderedBuffer(void *buffer, int size)
      throw(FilesystemFullException);
  protected:
    bool swapOrdering;
    CheckpointSize cSize;
    std::streamoff checkPointRef;
    uint64_t checkPointAccum;
    std::ofstream *f;

  };

};

#endif
