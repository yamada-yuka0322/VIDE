/*+
This is CosmoTool (./src/loadFlash.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
    data->Id = new int64_t[data->NumPart];
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


