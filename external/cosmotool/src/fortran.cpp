/*+
This is CosmoTool (./src/fortran.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#include "config.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include "fortran.hpp"

using namespace std;
using namespace CosmoTool;

UnformattedRead::UnformattedRead(const string& fname)
{
  f = new ifstream(fname.c_str());
  if (!*f)
    throw NoSuchFileException();

  swapOrdering = false;
  cSize = Check_32bits;
  checkPointRef = checkPointAccum = 0;
  recordBuffer = 0;
}

UnformattedRead::UnformattedRead(const char *fname)
{
  f = new ifstream(fname);
  if (!*f)
    throw NoSuchFileException();
  swapOrdering = false;
  cSize = Check_32bits;
  checkPointRef = checkPointAccum = 0;
  recordBuffer = 0;
}


UnformattedRead::~UnformattedRead()
{
  if (recordBuffer != 0)
    delete[] recordBuffer;
  delete f;
}

int64_t UnformattedRead::position() const
{
 return f->tellg();
}

void UnformattedRead::seek(int64_t pos)
{
  f->clear();
  f->seekg(pos, istream::beg);
  checkPointRef = checkPointAccum = 0;
}

// Todo implement primitive description
void UnformattedRead::setOrdering(Ordering o)
{
  if (o == LittleEndian)
    {
    }
}

void UnformattedRead::setCheckpointSize(CheckpointSize cs)
{
  cSize = cs;
}
    
void UnformattedRead::skip(int64_t off)
{
  if (checkPointAccum == 0 && checkPointRef == 0)
    {
      // We are not in a checked block
      f->seekg(off, ios::cur);
      return;
    }
  if (off < 0)
    throw InvalidUnformattedAccess();
  
  if ((checkPointAccum+off) > checkPointRef)
       throw InvalidUnformattedAccess();

  f->seekg(off, ios::cur);
  checkPointAccum += off;
}

void UnformattedRead::beginCheckpoint(bool bufferRecord)
{
  if (checkPointAccum != 0)
    throw InvalidUnformattedAccess();
  
  checkPointRef = (cSize == Check_32bits) ? 4 : 8;
  checkPointAccum = 0;  

  checkPointRef = (cSize == Check_32bits) ? readUint32() : readInt64();
  checkPointAccum = 0;

  if (f->eof())
    throw EndOfFileException();

  if (bufferRecord) {
    std::cout << "Use fast I/O mode" << std::endl;
    recordBuffer = new uint8_t[checkPointRef];
    f->read(reinterpret_cast<char *>(recordBuffer), checkPointRef);
  }
}

void UnformattedRead::endCheckpoint(bool autodrop)
{
  if (recordBuffer != 0) {
    delete[] recordBuffer;
    recordBuffer = 0;
  }

  if (checkPointRef != checkPointAccum)
    {
      if (!autodrop || checkPointAccum > checkPointRef) {
	throw InvalidUnformattedAccess();
      }
      f->seekg(checkPointRef-checkPointAccum, ios::cur);
    }

  int64_t oldCheckPoint = checkPointRef;

  checkPointRef = (cSize == Check_32bits) ? 4 : 8;
  checkPointAccum = 0;
  checkPointRef = (cSize == Check_32bits) ? readUint32() : readInt64();
  
  if (oldCheckPoint != checkPointRef) {
    throw InvalidUnformattedAccess();
  }

  checkPointAccum = checkPointRef = 0;
}

void UnformattedRead::readOrderedBuffer(void *buffer, int size)
{
  if ((checkPointAccum+(uint64_t)size) > checkPointRef)
       throw InvalidUnformattedAccess();

  if (recordBuffer != 0) {
    memcpy(buffer, &recordBuffer[checkPointAccum], size);
  } else {
    f->read((char *)buffer, size);
  }
  
  if (swapOrdering)
    {
      char *cb = (char *)buffer;
      for (int i = 0; i < size/2; i++)
	swap(cb[i], cb[size-i-1]);
    }
  checkPointAccum += size;
}
  
double UnformattedRead::readReal64()
{
  union
  {
    char b[8];
    double d;
  } a;

  readOrderedBuffer(&a, 8);

  return a.d;
}

float UnformattedRead::readReal32()
{
  union
  {
    char b[4];
    float f;
  } a;

  readOrderedBuffer(&a, 4);

  return a.f;
}

uint32_t UnformattedRead::readUint32()
{
  union
  {
    char b[4];
    uint32_t i;
  } a;

  readOrderedBuffer(&a, 4);

  return a.i;
}

int32_t UnformattedRead::readInt32()
{
  union
  {
    char b[4];
    int32_t i;
  } a;

  readOrderedBuffer(&a, 4);

  return a.i;
}

int64_t UnformattedRead::readInt64()
{
  union
  {
    char b[8];
    int64_t i;
  } a;

  readOrderedBuffer(&a, 8);

  return a.i;
}

//// UnformattedWrite

UnformattedWrite::UnformattedWrite(const string& fname)
{
  f = new ofstream(fname.c_str());
  if (!*f)
    throw NoSuchFileException();

  swapOrdering = false;
  cSize = Check_32bits;
  checkPointRef = checkPointAccum = 0;
}

UnformattedWrite::UnformattedWrite(const char *fname)
{
  f = new ofstream(fname);
  if (!*f)
    throw NoSuchFileException();
  swapOrdering = false;
  cSize = Check_32bits;
  checkPointRef = checkPointAccum = 0;
}


UnformattedWrite::~UnformattedWrite()
{
  delete f;
}

// Todo implement primitive description
void UnformattedWrite::setOrdering(Ordering o)
{
  if (o == LittleEndian)
    {
    }
}

void UnformattedWrite::setCheckpointSize(CheckpointSize cs)
{
  cSize = cs;
}
    
void UnformattedWrite::beginCheckpoint()
{
  if (checkPointAccum != 0)
    throw InvalidUnformattedAccess();
  
  checkPointRef = f->tellp();
  if (cSize == Check_32bits)
    writeInt32(0);
  else
    writeInt64(0);

  checkPointAccum = 0;

  if (!*f)
    throw FilesystemFullException();
}

void UnformattedWrite::endCheckpoint()
{
  if (checkPointAccum == 0)
    throw InvalidUnformattedAccess();

  streampos curPos = f->tellp();

  int64_t deltaPos = curPos-checkPointRef;

  deltaPos -= (cSize == Check_32bits) ? 4 : 8;
  // This is a sanity check.
  if (checkPointAccum != deltaPos)
    throw InvalidUnformattedAccess();

  uint64_t saveAccum = checkPointAccum;

  f->seekp(checkPointRef);
  if (cSize == Check_32bits)
    writeInt32(saveAccum);
  else
    writeInt64(saveAccum);

  f->seekp(curPos);
  if (cSize == Check_32bits)
    writeInt32(saveAccum);
  else
    writeInt64(saveAccum);

  checkPointAccum = checkPointRef = 0;
}

void UnformattedWrite::writeOrderedBuffer(void *buffer, int size)
{
  f->write((char *)buffer, size);
  
  if (swapOrdering)
    {
      char *cb = (char *)buffer;
      for (int i = 0; i < size/2; i++)
	swap(cb[i], cb[size-i-1]);
    }
  checkPointAccum += size;

  if (!*f)
    throw FilesystemFullException();
}
  
void UnformattedWrite::writeReal64(double d)
{
  union
  {
    char b[8];
    double d;
  } a;

  a.d = d;

  writeOrderedBuffer(&a, 8);
}

void UnformattedWrite::writeReal32(float f)
{
  union
  {
    char b[4];
    float f;
  } a;
  
  a.f = f;

  writeOrderedBuffer(&a, 4);
}

void UnformattedWrite::writeInt32(int32_t i)
{
  union
  {
    char b[4];
    int32_t i;
  } a;

  a.i = i;
  writeOrderedBuffer(&a, 4);
}

void UnformattedWrite::writeInt64(int64_t i)
{
  union
  {
    char b[8];
    int64_t i;
  } a;

  a.i = i;
  writeOrderedBuffer(&a, 8);
}

void UnformattedWrite::writeInt8(int8_t i)
{
  union { char b; int8_t i; } a;

  a.i = i;
  writeOrderedBuffer(&a, 1);
}
