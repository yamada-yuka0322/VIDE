/* This file contains the functions that read the data from the HDF5 file
 * The functions accept the PARAMESH data through arguments.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "hdf5_flash.h"
#include "H5Cpp.h"

using namespace H5;

void h5_read_runtime_parameters
     (H5File* file,    /* file handle */
      double* LBox,
      int* numPart,
      double* hubble,
      double* omegam,
      double* omegalambda);

void h5_read_flash3_particles (H5File* file,
                        int* totalparticles,
                        int* localnp,
                        int* particle_offset,
                        float *pos1,
                        float *pos2,
                        float *pos3,
                        float *vel1,
                        float *vel2,
                        float *vel3,
                        int    *id);

void h5_read_flash3_header_info(H5File* file,
				double* time,                   /* simulation time */
				double *redshift);              /* simulation redshift */
