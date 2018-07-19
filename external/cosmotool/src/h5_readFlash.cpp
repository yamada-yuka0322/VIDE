/*+
This is CosmoTool (./src/h5_readFlash.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
/*+
!!

This file has been developped by P. M. Sutter

!!

*/
/* This file contains the functions that read the data from the HDF5 file
 * The functions accept the PARAMESH data through arguments.
 */

#include "h5_readFlash.hpp"

using namespace H5;

/* indices of attributes in fof_particle_type */
int iptag_out = 0;
int ipx_out   = 0;
int ipy_out   = 1;
int ipz_out   = 2;
int ipvx_out   = 3;
int ipvy_out   = 4;
int ipvz_out   = 5;

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/*  n*_runtime_parameters should be set by the caller to 
    the maximum number of runtime parameters to read. 
*/

void h5_read_runtime_parameters
     (H5File* file,    /* file handle */
      double* LBox,    // box size
      int* numPart,
      double *hubble,
      double *omegam,
      double *omegalambda)       // number of particles
{

 int MAX_PARM = 200;

 int nreal_runtime_parameters, nint_runtime_parameters, nstr_runtime_parameters;

 char real_runtime_parameter_names[MAX_PARM][RUNTIME_PARAMETER_STRING_SIZE];
 char int_runtime_parameter_names[MAX_PARM][RUNTIME_PARAMETER_STRING_SIZE];
 char str_runtime_parameter_names[MAX_PARM][RUNTIME_PARAMETER_STRING_SIZE];

 double real_runtime_parameter_values[MAX_PARM];
 int    int_runtime_parameter_values[MAX_PARM];
 char   str_runtime_parameter_values[MAX_PARM][RUNTIME_PARAMETER_STRING_SIZE];

  int rank;
  hsize_t dimens_1d, maxdimens_1d;

  real_runtime_params_t *real_rt_parms;
  int_runtime_params_t *int_rt_parms;
  str_runtime_params_t *str_rt_parms;
  log_runtime_params_t *log_rt_parms;

  double omegarad;

  int i;

  nint_runtime_parameters = MAX_PARM;
  nreal_runtime_parameters = MAX_PARM;

  /* integer runtime parameters */
  int_rt_parms = (int_runtime_params_t *) malloc(nint_runtime_parameters * sizeof(int_runtime_params_t)); 

  rank = 1;
  DataSet dataset = file->openDataSet("integer runtime parameters");
  
  IntType int_rt_type = dataset.getIntType();  
  //int_rt_type = H5Dget_type(dataset);

  DataSpace dataspace = dataset.getSpace();
  //dataspace = H5Dget_space(dataset);

  int ndims = dataspace.getSimpleExtentDims(&dimens_1d, NULL);
  //H5Sget_simple_extent_dims(dataspace, &dimens_1d, &maxdimens_1d);

  /* don't read in more than we can handle */
  if (nint_runtime_parameters < dimens_1d) {
    dimens_1d = nint_runtime_parameters;
  } else {
    nint_runtime_parameters = dimens_1d;
  }
  DataSpace memspace(rank, &dimens_1d);
  //memspace = H5Screate_simple(rank, &dimens_1d, NULL);

  dataset.read(int_rt_parms, int_rt_type, memspace, dataspace,
		    H5P_DEFAULT);
  //status = H5Dread(dataset, int_rt_type, memspace, dataspace,
	//	    H5P_DEFAULT, int_rt_parms);


  for (i = 0; i < nint_runtime_parameters; i++) {
    strncpy(int_runtime_parameter_names[i], int_rt_parms[i].name, RUNTIME_PARAMETER_STRING_SIZE);
    int_runtime_parameter_values[i] = int_rt_parms[i].value;
  }

  free(int_rt_parms);
  memspace.close();
  dataspace.close();
  dataset.close();
  //H5Sclose(dataspace);
  //H5Dclose(dataset);

  /* done with int runtime parameters */

  /* reals */

  real_rt_parms = (real_runtime_params_t *) malloc(nreal_runtime_parameters * sizeof(real_runtime_params_t)); 
  
  rank = 1;
  dataset = file->openDataSet("real runtime parameters");
  //dataset = H5Dopen(*file_identifier, "real runtime parameters"); 
  
  dataspace = dataset.getSpace();
  FloatType real_rt_type = dataset.getFloatType();  
  ndims = dataspace.getSimpleExtentDims(&dimens_1d, NULL);
  //dataspace = H5Dget_space(dataset);
  //real_rt_type = H5Dget_type(dataset);
  //H5Sget_simple_extent_dims(dataspace, &dimens_1d, &maxdimens_1d);

  /* don't read in more than we can handle */
  if (nreal_runtime_parameters < dimens_1d) {
    dimens_1d = nreal_runtime_parameters;
  } else {
    nreal_runtime_parameters = dimens_1d;
  }
   memspace = DataSpace(rank, &dimens_1d);
  //memspace = H5Screate_simple(rank, &dimens_1d, NULL);

  dataset.read(real_rt_parms, real_rt_type, memspace, dataspace,
		    H5P_DEFAULT);
  //status = H5Dread(dataset, real_rt_type, memspace, dataspace,
//		    H5P_DEFAULT, real_rt_parms);


  for (i = 0; i < nreal_runtime_parameters; i++) {
    strncpy(real_runtime_parameter_names[i], real_rt_parms[i].name, RUNTIME_PARAMETER_STRING_SIZE);
    real_runtime_parameter_values[i] = real_rt_parms[i].value;
  }

  free(real_rt_parms);
  memspace.close();
  dataspace.close();
  dataset.close();
  //H5Sclose(dataspace);
  //H5Dclose(dataset);

  /* done with reals */

  // grab the data we want
  for (i = 0; i < nreal_runtime_parameters; i++) {
    if (strncmp(real_runtime_parameter_names[i],"xmax",4) == 0 ) {
      *LBox = real_runtime_parameter_values[i];
    }
    if (strncmp(real_runtime_parameter_names[i],"hubbleconstant", 14) == 0 ) {
      *hubble = real_runtime_parameter_values[i];      
    }
    if (strncmp(real_runtime_parameter_names[i],"omegamatter", 11) == 0 ) {
      *omegam = real_runtime_parameter_values[i];
    }
    if (strncmp(real_runtime_parameter_names[i],"omegaradiation", 11) == 0 ) {
      omegarad = real_runtime_parameter_values[i];
    }
    if (strncmp(real_runtime_parameter_names[i],"cosmologicalconstant", 20) == 0 ) {
      *omegalambda = real_runtime_parameter_values[i];
    }
  }
  
  for (i = 0; i < nint_runtime_parameters; i++) {
    if (strncmp(int_runtime_parameter_names[i],"pt_numx",7) == 0 ) {
      *numPart = int_runtime_parameter_values[i];
      *numPart = *numPart * *numPart * *numPart;
    }
  }
}

     
/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
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
                        int64_t    *id)
{

  herr_t   status;
  hsize_t  dimens_1d, maxdims_1d;
  hsize_t  start_1d, stride_1d, count_1d;
  int rank;

  int      numProps, i, p;
  int      numPartBuffer = 5000, sizePartBuffer, sizePart;
  int      pstack, poffset, pcount;
  int      iptag, ipx, ipy, ipz, ipvx, ipvy, ipvz;
  char     *propName;
  double   *partBuffer;
  char    part_names[50][OUTPUT_PROP_LENGTH];
  int     string_size;  
         
//  char part_names[NPART_PROPS][OUTPUT_PROP_LENGTH];
  hsize_t dimens_2d[2], maxdimens_2d[2];
  hsize_t start_2d[2], count_2d[2], stride_2d[2];

  /* skip this routine if no particles to read */
  if ((*localnp) == 0) {
    return;
  }

 /* first determine how many particle properties are
     present in the input data file, and determine which of these
     are the properties we are interested in */
  DataSet dataset = file->openDataSet("particle names");
  DataSpace dataspace = dataset.getSpace();
  //dataset = H5Dopen(*file_identifier, "particle names");
  //dataspace = H5Dget_space(dataset);

  int ndims = dataspace.getSimpleExtentDims(dimens_2d, NULL);
  //H5Sget_simple_extent_dims(dataspace, dimens_2d, maxdimens_2d);

  //total number of particle properties
  numProps = dimens_2d[0];
  
  string_size = OUTPUT_PROP_LENGTH;
  StrType string_type = H5Tcopy(H5T_C_S1);
  string_type.setSize(string_size);
  //status = H5Tset_size(string_type, string_size);
  
  rank = 2;

  start_2d[0] = 0;
  start_2d[1] = 0;

  stride_2d[0] = 1;
  stride_2d[1] = 1;

  count_2d[0] = dimens_2d[0];
  count_2d[1] = dimens_2d[1];

  dataspace.selectHyperslab(H5S_SELECT_SET, count_2d, start_2d); 
  //status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d,
  //                      stride_2d, count_2d, NULL);

  DataSpace memspace(rank, dimens_2d);
  //memspace = H5Screate_simple(rank, dimens_2d, NULL);

  dataset.read(part_names, string_type, H5S_ALL, H5S_ALL, H5P_DEFAULT);
  //status = H5Dread(dataset, string_type, H5S_ALL, H5S_ALL,
  //              H5P_DEFAULT, part_names);

 
  string_type.close();
  memspace.close();
  dataspace.close();
  dataset.close(); 
  //H5Tclose(string_type);
  //H5Sclose(memspace);
  //H5Sclose(dataspace);
  //H5Dclose(dataset);

  for (i=0;i<numProps;i++) {
    if (strncmp(part_names[i], "tag", 3) == 0)  { iptag = i+1; }
    if (strncmp(part_names[i], "posx", 4) == 0) { ipx = i+1; }
    if (strncmp(part_names[i], "posy", 4) == 0) { ipy = i+1; }
    if (strncmp(part_names[i], "posz", 4) == 0) { ipz = i+1; }
    if (strncmp(part_names[i], "velx", 4) == 0) { ipvx = i+1; }
    if (strncmp(part_names[i], "vely", 4) == 0) { ipvy = i+1; }
    if (strncmp(part_names[i], "velz", 4) == 0) { ipvz = i+1; }
  }

  if ((iptag < 0) || (ipx < 0) || (ipy < 0) || (ipz < 0) || (ipvx < 0) || 
       (ipvy < 0) || (ipvz < 0) ) {
    printf("One or more required particle attributes not found in file!\n");
    return;
  }

  //printf("iptag = %d, ipx = %d, ipy = %d, ipz = %d\n", iptag, ipx, ipy, ipz);
  //printf("ipvx = %d, ipvy = %d, ipvz = %d\n", ipvx, ipvy, ipvz);

  //read particles
  dataset = file->openDataSet("tracer particles");
  //dataset = H5Dopen(*file_identifier, "tracer particles");
  
  FloatType datatype = dataset.getFloatType();
  //datatype = H5Dget_type(dataset);

  /* establish read-in particle buffer */
  sizePart = numProps*(sizeof(double));
  sizePartBuffer = numPartBuffer * sizePart;
  partBuffer = (double *)malloc(sizePartBuffer);

  dataspace = dataset.getSpace();
  ndims = dataspace.getSimpleExtentDims(dimens_2d, NULL);
  //dataspace = H5Dget_space(dataset);
  //H5Sget_simple_extent_dims(dataspace, dimens_2d, maxdimens_2d);
 
  /*insert particle properties (numPartBuffer) particles at a time*/
  pstack = (*localnp);
  poffset = 0;
  if (pstack > numPartBuffer) {
   pcount = numPartBuffer;
  }
  else {
   pcount = pstack;
  }

  while ( pstack > 0) {
    rank       = 2;
    maxdimens_2d[0] = (hsize_t) (*totalparticles);
    maxdimens_2d[1] = (hsize_t) (numProps);

    start_2d[0]  = (hsize_t) (*particle_offset + poffset);
    start_2d[1]  = (hsize_t) 0;

    stride_2d[0] = 1;
    stride_2d[1] = 1;

    count_2d[0]  = (hsize_t) (pcount);
    count_2d[1]  = (hsize_t) (numProps);

    dimens_2d[0] = (pcount);
    dimens_2d[1] = (numProps);

    dataspace.selectHyperslab(H5S_SELECT_SET, count_2d, start_2d); 
    //status     = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d,
    //                                 stride_2d, count_2d, NULL);

    memspace = DataSpace(rank, dimens_2d);
    //memspace   = H5Screate_simple(rank, dimens_2d, maxdimens_2d);

    /* read data from the dataset */
   dataset.read(partBuffer, datatype, memspace, dataspace, H5P_DEFAULT);
   //status = H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, partBuffer);

    /* convert buffer into particle struct */

    if (id) {
    for(p=0; p < (pcount); p++) {
      id[p+poffset] = (int64_t)  *(partBuffer+iptag-1+p*numProps);
    } }

    if (pos1 && pos2 && pos3) {
    for(p=0; p < (pcount); p++) {
      pos1[p+poffset]  = (float) *(partBuffer+ipx-1+p*numProps);
      pos2[p+poffset]  = (float) *(partBuffer+ipy-1+p*numProps);
      pos3[p+poffset]  = (float) *(partBuffer+ipz-1+p*numProps);
    }
    }
   
    
    if (vel1 && vel2 && vel3) {
    for(p=0; p < (pcount); p++) {
      vel1[p+poffset]  = (float) *(partBuffer+ipvx-1+p*numProps);
      vel2[p+poffset]  = (float) *(partBuffer+ipvy-1+p*numProps);
      vel3[p+poffset]  = (float) *(partBuffer+ipvz-1+p*numProps);
    }
    }

    memspace.close();
    //status = H5Sclose(memspace);
   /* advance buffer */
    pstack = pstack - pcount;
    poffset = poffset + pcount;
    if (pstack > numPartBuffer) {
      pcount = numPartBuffer;
    }
    else {
      pcount = pstack;
    }

  } /* end while */


  datatype.close();
  dataspace.close();
  dataset.close();
  //status = H5Tclose(datatype);
  //status = H5Sclose(dataspace);
  //status = H5Dclose(dataset);
}
 

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

void h5_read_flash3_header_info(H5File* file,
				double* time,                  /* simulation time */
				double* redshift)              /* redshift of checkpoint */
{

  herr_t status;

  int file_version;

  hid_t sp_type, si_type;

  hsize_t dimens_1d, maxdimens_1d;
  hid_t string_type;
  real_list_t *real_list;
  int* num_real, num_int;
  int MAX_SCALARS = 100;
  char real_names[MAX_SCALARS][MAX_STRING_LENGTH];
  double real_values[MAX_SCALARS];
  char int_names[MAX_SCALARS][MAX_STRING_LENGTH];
  int int_values[MAX_SCALARS];
  int i;

  H5std_string DATASET_NAME;

  string_type = H5Tcopy(H5T_C_S1);  
  H5Tset_size(string_type, MAX_STRING_LENGTH);

  DataSet dataset = file->openDataSet("real scalars");
  DataSpace dataspace = dataset.getSpace();

  /* read extent of 'dataspace' (i.e. # of name/value pairs) into 'dimens_1d' */
  int ndims = dataspace.getSimpleExtentDims(&dimens_1d, NULL);

  if (dimens_1d > MAX_SCALARS) {
    printf("Error: reading more than MAX_SCALARS runtime parameters in checkpoint file!\n");
  }

  /* malloc a pointer to a list of real_list_t's */
  real_list = (real_list_t *) malloc(dimens_1d * sizeof(real_list_t));

  // create a new simple dataspace of 1 dimension and size of 'dimens_1d' 
  DataSpace memspace(1, &dimens_1d);

  // create an empty vessel sized to hold one real_list_t's worth of data 
  CompType real_list_type( sizeof(real_list_t) );

  // subdivide the empty vessel into its component sections (name and value) 
  real_list_type.insertMember(
          "name",
          HOFFSET(real_list_t, name),
          string_type);

  real_list_type.insertMember(
          "value",
          HOFFSET(real_list_t, value),
          PredType::NATIVE_DOUBLE);

  // read the data into 'real_list' 
  dataset.read( real_list, real_list_type, memspace, dataspace,
                H5P_DEFAULT);


  if (status < 0) {
    printf("Error readingruntime parameterss from data file\n");
  }

  for (i = 0; i < dimens_1d; i++) {
    strncpy(real_names[i], real_list[i].name, MAX_STRING_LENGTH);
    real_values[i] = real_list[i].value;

    if (strncmp(real_names[i],"time",4) == 0 ) {
      *time = real_values[i];
    }
    if (strncmp(real_names[i],"redshift",8) == 0 ) {
      *redshift = real_values[i];
    }
  }

  free(real_list);
  real_list_type.close();
  memspace.close();
  dataspace.close();
  dataset.close();
}
