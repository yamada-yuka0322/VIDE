#include <cmath>
#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
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
  f->endCheckpoint(true);
  
  long NumPart = 0, NumPartTotal = 0;
  for(int k=0; k<6; k++)
    {
      NumPart += h.npart[k];
      NumPartTotal += (id < 0) ? h.npart[k] : h.npartTotal[k];
    }
  data->NumPart = NumPart;
  data->TotalNumPart = NumPartTotal;
}

SimuData *CosmoTool::loadGadgetMulti(const char *fname, int id,
				     int loadflags, int GadgetFormat, 
				     SimuFilter filter)
{
  SimuData *data;
  int p, n;
  UnformattedRead *f;
  GadgetHeader h;
  float velmul;

  if (id >= 0) {
    int k = snprintf(0, 0, "%s.%d", fname, id)+1;
    char *out_fname = new char[k];
    snprintf(out_fname, k, "%s.%d", fname, id);

    f = new UnformattedRead(out_fname);
    if (f == 0)
      return 0;

    delete out_fname;

  } else {

    f = new UnformattedRead(fname);
    if (f == 0)
      return 0;

  }

  data = new SimuData;
  if (data == 0) {
    delete f;
    return 0;
  }

  long NumPart = 0, NumPartTotal = 0;
  
  try
    {
      loadGadgetHeader(f, h, data, id);

      if (GadgetFormat == 1)
        velmul = sqrt(h.time);
      else if (GadgetFormat == 2)
        velmul = 1/(h.time);
      else {
        cerr << "unknown gadget format" << endl;
	abort();
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
	f->beginCheckpoint();
	for(int k = 0, p = 0; k < 6; k++) {
	  for(int n = 0; n < h.npart[k]; n++) {
	    data->Pos[0][p] = f->readReal32();
	    data->Pos[1][p] = f->readReal32();
	    data->Pos[2][p] = f->readReal32();
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
    // Skip positions
    f->skip(NumPart * 3 * sizeof(float) + 2*4);
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
	f->beginCheckpoint();
	for(int k = 0, p = 0; k < 6; k++) {
	  for(int n = 0; n < h.npart[k]; n++) {
	    // THIS IS GADGET 1
	    data->Vel[0][p] = f->readReal32()*velmul;
	    data->Vel[1][p] = f->readReal32()*velmul;
	    data->Vel[2][p] = f->readReal32()*velmul;
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
    f->skip(NumPart*3*sizeof(float)+2*4);
  }

  // Skip ids
  if (loadflags & NEED_GADGET_ID) {
    try
      {
	f->beginCheckpoint();
	data->Id = new long[data->NumPart];
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
		data->Id[p] = f->readInt32();
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

  delete f;

  return data;
}



void CosmoTool::writeGadget(const char *fname, SimuData *data, int GadgetFormat)
{
  UnformattedWrite *f;
  int npart[6];
  float mass[6];

  if (data->Pos[0] == 0 || data->Vel[0] == 0 || data->Id == 0)
    return;

  f = new UnformattedWrite(fname);
  if (f == 0)
    return;

  for (int i = 0; i < 6; i++)
    {
      npart[i] = 0;
      mass[i] = 0;
    }
  mass[1] = 1.0;

  npart[1] = data->NumPart;
  
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
    f->writeInt32(npart[i]);
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
  if (GadgetFormat == 1)
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

  
