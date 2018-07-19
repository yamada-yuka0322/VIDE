/*+
This is CosmoTool (./src/loadGadget.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#include <cmath>
#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include "load_data.hpp"
#include "loadGadget.hpp"
#include "fortran.hpp"

using namespace CosmoTool;
using namespace std;


void loadGadgetHeader(UnformattedRead *f, GadgetHeader& h, SimuData *data, int id)
{
  f->beginCheckpoint();
  for (int i = 0; i < 6; i++)
    h.npart[i] = f->readInt32();
  for (int i = 0; i < 6; i++)    
    h.mass[i] = f->readReal64();
  data->time = h.time = f->readReal64();
  h.redshift = f->readReal64();
  h.flag_sfr = f->readInt32();
  h.flag_feedback = f->readInt32();
  for (int i = 0; i < 6; i++)
    h.npartTotal[i] = f->readInt32();
  h.flag_cooling = f->readInt32();
  h.num_files = f->readInt32();
  data->BoxSize = h.BoxSize = f->readReal64();
  data->Omega_M = h.Omega0 = f->readReal64();
  data->Omega_Lambda = h.OmegaLambda = f->readReal64();
  data->Hubble = h.HubbleParam = f->readReal64();
  (int)f->readInt32(); // stellarage
  (int)f->readInt32(); // metals
  for (int i = 0; i < 6; i++)
    h.npartTotal[i] |= ((unsigned long)f->readInt32()) << 32;
  (int)f->readInt32(); // entropy instead of u
  h.flag_doubleprecision = f->readInt32();  

  f->endCheckpoint(true);
  
  ssize_t NumPart = 0, NumPartTotal = 0;
  for(int k=0; k<6; k++)
    {
      NumPart += h.npart[k];
      NumPartTotal += (id < 0) ? h.npart[k] : h.npartTotal[k];
    }
  data->NumPart = NumPart;
  data->TotalNumPart = NumPartTotal;
}

template<typename T>
T myRead64(UnformattedRead *f) { return f->readReal64(); }

template<typename T>
T myRead32(UnformattedRead *f) { return f->readReal32(); }

template<typename T>
T myReadI64(UnformattedRead *f) { return f->readInt64(); }

template<typename T>
T myReadI32(UnformattedRead *f) { return f->readInt32(); }

struct BlockInfo {
   int64_t position, size;
};

SimuData *CosmoTool::loadGadgetMulti(const char *fname, int id,
				     int loadflags, int GadgetFormat, 
				     SimuFilter filter)
{
  SimuData *data;
  int p, n;
  UnformattedRead *f;
  GadgetHeader h;
  float velmul;
  boost::function0<double> readToDouble;
  boost::function0<float> readToSingle;
  boost::function0<int64_t> readID;
  long float_size;

  if (GadgetFormat > 2) {
     cerr << "Unknown gadget format" << endl;
     return 0;
  }

  if (id >= 0) {
    int k = snprintf(0, 0, "%s.%d", fname, id)+1;
    char *out_fname = new char[k];
    snprintf(out_fname, k, "%s.%d", fname, id);

    f = new UnformattedRead(out_fname);
    if (f == 0)
      return 0;

    delete[] out_fname;

  } else {

    f = new UnformattedRead(fname);
    if (f == 0)
      return 0;

  }

  typedef std::map<std::string, BlockInfo> BlockMap;
  BlockMap blockTable;
  
  if (GadgetFormat == 2) {
    int64_t startBlock = 0;
    char block[5];
    uint32_t blocksize;

    try {
      while (true) {
        f->beginCheckpoint();
        f->readOrderedBuffer(block, 4);
        block[4] = 0;
        blocksize = f->readUint32();
        f->endCheckpoint();
        blockTable[block].position = f->position();
        blockTable[block].size = blocksize;
        f->skip(blocksize);
      }
    } catch (EndOfFileException&) {}

    f->seek(startBlock);
  }

  data = new SimuData;
  if (data == 0) {
    delete f;
    return 0;
  }

  ssize_t NumPart = 0, NumPartTotal = 0;
#define ENSURE2(name,out_sz) { \
    int64_t sz; \
    if (GadgetFormat == 2) { \
     BlockMap::iterator iter = blockTable.find(name); \
     if (iter == blockTable.end()) { \
        std::cerr << "GADGET2: Cannot find block named '" << name << "'" << endl; \
        if (data) delete data; \
        delete f; \
        return 0; \
     } \
     f->seek(iter->second.position); \
     sz = iter->second.size; \
     out_sz = sz;\
    } else if (GadgetFormat==1) { \
      int64_t oldpos = f->position(); \
      f->beginCheckpoint(); \
      out_sz = f->getBlockSize(); \
      f->endCheckpoint(true); \
      f->seek(oldpos); \
    } \
\
  } 
#define ENSURE(name) ENSURE2(name,sz);
  
  try
    {
      ENSURE("HEAD");
      loadGadgetHeader(f, h, data, id);

      velmul = sqrt(h.time);
      
      if (h.flag_doubleprecision) {
        //cout << "Gadget format with double precision" << endl;
        readToDouble = boost::bind(myRead64<double>, f);
        readToSingle = boost::bind(myRead64<float>, f);
        float_size = sizeof(double);
      } else {
        //cout << "Gadget format with single precision" << endl;
        readToDouble = boost::bind(myRead32<double>, f);
        readToSingle = boost::bind(myRead32<float>, f);
        float_size = sizeof(float);
      }
      NumPart = data->NumPart;
      NumPartTotal = data->TotalNumPart;
    }
  catch (const InvalidUnformattedAccess& e)
    {
      cerr << "Invalid format while reading header" << endl;
      delete data;
      delete f;
      return 0;
    }


  if (loadflags & NEED_TYPE)
    {
      int p = 0;
 
      data->type = new int[data->NumPart];
      for (int k = 0; k < 6; k++)
        for (int n = 0; n < h.npart[k]; n++,p++)
            data->type[p] = k;
    }

  if (loadflags & NEED_POSITION) {
      
    for (int i = 0; i < 3; i++) {
        data->Pos[i] = new float[data->NumPart];
        if (data->Pos[i] == 0) {
        	delete data;
    	return 0;
    }
  }
    
    try
      {
        ENSURE("POS ");
	f->beginCheckpoint(true); // Use more memory but faster I/O
        if (f->getBlockSize() != NumPart*float_size*3) {
          // Check that single would work
          if (f->getBlockSize() == NumPart*sizeof(float)*3) {
            // Change to single
            cout << "Change back to single. Buggy header." << endl;
            readToDouble = boost::bind(myRead32<double>, f);
            readToSingle = boost::bind(myRead32<float>, f);
            float_size = sizeof(float);
          }
        }
	for(int k = 0, p = 0; k < 6; k++) {
	  for(int n = 0; n < h.npart[k]; n++) {
//            if ((n%1000000)==0) cout << n << endl;
	    data->Pos[0][p] = readToSingle();
	    data->Pos[1][p] = readToSingle();
	    data->Pos[2][p] = readToSingle();
	    p++;
	  }
	}
	f->endCheckpoint();
      }
    catch (const InvalidUnformattedAccess& e)
      {
        cerr << "Invalid format while reading positions" << endl;
        delete f;
        delete data;
        return 0;
      }
    
  } else {
    long float_size = (h.flag_doubleprecision) ? sizeof(double) : sizeof(float);
    // Skip positions
    f->skip(NumPart * 3 * float_size + 2*4);
  }
  
  if (loadflags & NEED_VELOCITY) {
    for (int i = 0; i < 3; i++)
      {
	data->Vel[i] = new float[data->NumPart];
	if (data->Vel[i] == 0) 
	  {
	    delete f;
	    delete data;
	    return 0;
	  }
      }
    
    try
      {
        ENSURE("VEL ");
	f->beginCheckpoint(true);
	for(int k = 0, p = 0; k < 6; k++) {
	  for(int n = 0; n < h.npart[k]; n++) {
	    // THIS IS GADGET 1
//            if ((n%1000000)==0) cout << n << endl;
	    data->Vel[0][p] = readToSingle()*velmul;
	    data->Vel[1][p] = readToSingle()*velmul;
	    data->Vel[2][p] = readToSingle()*velmul;
	    p++;
	  }
	}
	f->endCheckpoint();
      }
    catch (const InvalidUnformattedAccess& e)
      {
	cerr << "Invalid format while reading velocities" << endl;
	delete f;
	delete data;
	return 0;
      }
    
    // THE VELOCITIES ARE IN PHYSICAL COORDINATES
///    // TODO: FIX THE UNITS OF THESE FUNKY VELOCITIES !!!
  } else {
    // Skip velocities
    f->skip(NumPart*3*float_size+2*4);
  }

  // Skip ids
  if (loadflags & NEED_GADGET_ID) {
    try
      {
        int64_t idSize;
        ENSURE2("ID  ", idSize);

        if (idSize / data->NumPart == 8) {
          readID = boost::bind(myReadI64<int64_t>, f);
        } else
        if (idSize / data->NumPart == 4) {
          readID = boost::bind(myReadI32<int64_t>, f);
        } else {
          throw InvalidUnformattedAccess();
        }

	f->beginCheckpoint(true);
	data->Id = new int64_t[data->NumPart];
	if (data->Id == 0)
	  {
	    delete f;
	    delete data;
	    return 0;
	  }
	
	for(int k = 0, p = 0; k < 6; k++)
	  {
	    for(int n = 0; n < h.npart[k]; n++)
	      {
		data->Id[p] = readID();
		p++;
	      }
	  }
	f->endCheckpoint();
      }
    catch (const InvalidUnformattedAccess& e)
      {
	cerr << "Invalid unformatted access while reading ID" << endl;
	delete f;
	delete data;
	return 0;
      }
  } else {
    f->skip(2*4);
    for (int k = 0; k < 6; k++)
      f->skip(h.npart[k]*4);
  }

  if (loadflags & NEED_MASS) {
    bool do_load = false;

    for (int k = 0; k < 6; k++)
      {
         do_load = do_load || ((h.mass[k] == 0)&&(h.npart[k]>0));
      }

    try
      {
        long l = 0;
	if (do_load) {
          ENSURE("MASS");
          f->beginCheckpoint();
        }
        data->Mass = new float[NumPart];
        for (int k = 0; k < 6; k++)
	  {
            if (h.mass[k] == 0) {
              for(int n = 0; n < h.npart[k]; n++)
	        {
//            if ((n%1000000)==0) cout << n << endl;
                  data->Mass[l++] = readToSingle();
                }
            } else {
              for(int n = 0; n < h.npart[k]; n++)
	        {
//            if ((n%1000000)==0) cout << n << endl;
                  data->Mass[l++] = h.mass[k];
                }
            }
          }
        if (do_load)
          f->endCheckpoint();
      }
    catch (const InvalidUnformattedAccess& e)
      {
	cerr << "Invalid unformatted access while reading ID" << endl;
	delete f;
	delete data;
	return 0;
      }
    catch (const EndOfFileException& e)
      {
        for (int k = 0; k < 6; k++)
          cerr << "mass[" << k << "] = " << h.mass[k] << endl;
      }
  } else {
    f->skip(2*4);
    for (int k = 0; k < 6; k++)
      if (h.mass[k] == 0)
        f->skip(h.npart[k]*4);
  }

  delete f;

  return data;
}

#undef ENSURE


void CosmoTool::writeGadget(const std::string& fname, SimuData *data, int GadgetFormat)
{
  UnformattedWrite *f;
  int npart[6], npartTotal[6];
  float mass[6];

  if (data->Pos[0] == 0 || data->Vel[0] == 0 || data->Id == 0) {
    cerr << "Invalid simulation data array" << endl;
    return;
  }

  f = new UnformattedWrite(fname);
  if (f == 0) {
    cerr << "Cannot create file" << endl;
    return;
  }

  for (int i = 0; i < 6; i++)
    {
      npart[i] = npartTotal[i] = 0;
      mass[i] = 0;
    }
  mass[1] = 1.0;
  npart[1] = data->NumPart;
  npartTotal[1] = data->TotalNumPart;
  
  f->beginCheckpoint();
  for (int i = 0; i < 6; i++)
    f->writeInt32(npart[i]);
  for (int i = 0; i < 6; i++)
    f->writeReal64(mass[i]);

  f->writeReal64(data->time);
  f->writeReal64(1/data->time-1);
  f->writeInt32(0);
  f->writeInt32(0);

  for (int i = 0; i < 6; i++)
    f->writeInt32(npartTotal[i]);
  f->writeInt32(0);
  f->writeInt32(1);
  f->writeReal64(data->BoxSize);
  f->writeReal64(data->Omega_M);
  f->writeReal64(data->Omega_Lambda);
  f->writeReal64(data->Hubble);
  char buf[100] = { 0, };
  f->writeOrderedBuffer(buf, 96);
  f->endCheckpoint();
 
  f->beginCheckpoint();
  for(int n = 0; n < data->NumPart; n++) {
    for (int k = 0; k < 3; k++)
      f->writeReal32(data->Pos[k][n]);
  }
  f->endCheckpoint();

  float velmul = 1.0;
  velmul = sqrt(data->time);

  f->beginCheckpoint();
  for(int n = 0; n < data->NumPart; n++) {
    for (int k = 0; k < 3; k++)
      f->writeReal32(data->Vel[k][n]/velmul);
  }
  f->endCheckpoint();

  f->beginCheckpoint();
  for(int n = 0; n < data->NumPart; n++)
    {
      f->writeInt32(data->Id[n]);
    }
  f->endCheckpoint();
  delete f;
}

  
