/* Reads in FLASH v3 files in HDF5 format */

#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "load_data.hpp"
#include "loadFlash.hpp"
#include "h5_readFlash.hpp"
#include "H5Cpp.h"

using namespace CosmoTool;
using namespace std;
using namespace H5;

SimuData *CosmoTool::loadFlashMulti(const char *fname, int id, int loadflags)
{
  SimuData *data;
  int p, n;
  H5File *fileID;
  H5std_string filename;
  //char filename[81];
  double lbox, time, hubble, omegam, omegalambda, redshift;
  int npart;

  const double kpc2cm = 3.08568025e21;
  const double km2cm = 1.e5;
  const double hubble2cm = 3.240779270005e-18;

  if (id != 0)
    throw NoSuchFileException();

  data = new SimuData;
  if (data == 0) {
    return 0;
  }

  filename = fname;
  try {
    H5File file (filename, H5F_ACC_RDONLY);

    // simulation info
    h5_read_flash3_header_info(&file, &time, &redshift);
    data->time    = 1/(1+redshift);    

    h5_read_runtime_parameters(&file, &lbox, &npart, &hubble, &omegam, &omegalambda);
    data->TotalNumPart = data->NumPart = npart;
    data->Hubble = hubble/hubble2cm;
    data->BoxSize = lbox/kpc2cm*data->Hubble;
    data->Omega_M = omegam;
    data->Omega_Lambda = omegalambda;

    // particle data
    if (loadflags& NEED_POSITION) {
      for (int i = 0; i < 3; i++) {
        data->Pos[i] = new float[data->NumPart];
        if (data->Pos[i] == 0) {
	  delete data;
	return 0;
      }
    } }


    if (loadflags &NEED_VELOCITY) {
    for (int i = 0; i < 3; i++) {
      data->Vel[i] = new float[data->NumPart];
      if (data->Vel[i] == 0) {
	delete data;
	return 0;
      }
    } }
 
    if (loadflags & NEED_GADGET_ID) {
    data->Id = new int[data->NumPart];
    if (data->Id == 0) {
      delete data;
      return 0;
    }
    }

    int offset = 0;

    if (loadflags &(NEED_GADGET_ID|NEED_POSITION|NEED_VELOCITY))
    h5_read_flash3_particles(&file, &npart, &npart, &offset, 
			     data->Pos[0], data->Pos[1], data->Pos[2], 
			     data->Vel[0], data->Vel[1], data->Vel[2],
			     data->Id); 

    for (int i = 0; i < 3; i++) {
      for (int n = 0; n < data->NumPart; n++) {
	if (loadflags& NEED_POSITION) data->Pos[i][n] = data->Pos[i][n] * data->Hubble / kpc2cm;
	if (loadflags&NEED_VELOCITY) data->Vel[i][n] = data->Vel[i][n] * data->Hubble / km2cm;
      }
    }

    file.close();
  } catch (const FileIException& e) {
    throw NoSuchFileException();
  }

  return data;
}


