/* general header file for the HDF 5 IO in FLASH */


#ifndef _HDF5_FLASH_H
#define _HDF5_FLASH_H

/* pull in some basic FLASH information */

//#include "flash_defines.fh"

/* define an integer file format version number, that is stored
   in the output files.  This way, people can check this number
   before reading, and compare against a published format to know
   what is stored in the file.  In theory, this number should be
   incremented anytime a change is made to the file format */

/* File format history:

   1 -- original version

   2 -- added build records:
         "FLASH build date"
         "FLASH build directory"
         "FLASH build machine"
         "FLASH setup call"

        added the "file format version" record

   3 -- added the "run comment record" (why was this
        not done long ago?)

   4 -- added extrema attributes to the variable records

   5 -- redshift included

   6 -- added the Module data to the attributes of "/"

   7 -- make build info attributes on "/"
*/

#define FILE_FORMAT_VERSION 7

#define RUNTIME_PARAMETER_STRING_SIZE 80

#define TIMER_NAME_STRING_SIZE 30

#define MAX_STRING_LENGTH 80

#define LIST_STRING_SIZE 80

#define OUTPUT_PROP_LENGTH 24

typedef struct real_list_t {
  char name[LIST_STRING_SIZE];
  double value;
} real_list_t;

typedef struct int_runtime_params_t {
  char name[RUNTIME_PARAMETER_STRING_SIZE];
  int value;
} int_runtime_params_t;

typedef struct real_runtime_params_t {
  char name[RUNTIME_PARAMETER_STRING_SIZE];
  double value;
} real_runtime_params_t;

typedef struct str_runtime_params_t {
  char value[RUNTIME_PARAMETER_STRING_SIZE];
  char name[RUNTIME_PARAMETER_STRING_SIZE];
} str_runtime_params_t;

typedef struct log_runtime_params_t {
  int value;
  char name[RUNTIME_PARAMETER_STRING_SIZE];
} log_runtime_params_t;


#define MAX_TIMER_PARENTS 20
#define MAX_TIMER_CALL_STACK_DEPTH 20
typedef struct timer_data_t {
  char name[TIMER_NAME_STRING_SIZE];
  double t_value[MAX_TIMER_PARENTS];
  int    t_counts[MAX_TIMER_PARENTS];
  int    t_on[MAX_TIMER_PARENTS];
  int    t_stacks[MAX_TIMER_PARENTS][MAX_TIMER_CALL_STACK_DEPTH];
  int    t_num_parents;
  int    t_stack_sizes[MAX_TIMER_PARENTS];
} full_timer_data_t;


typedef struct sim_params_t {
  int total_blocks;
  int nsteps;
  int nxb;
  int nyb;
  int nzb;
  double time; 
  double timestep;
  double redshift;
    
} sim_params_t;

typedef struct sim_params_sp_t {
  int total_blocks;
  int nsteps;
  int nxb;
  int nyb;
  int nzb;
  float time; 
  float timestep;
  float redshift;
    
} sim_params_sp_t;

typedef struct sim_info_t {
  int file_format_version;
  char setup_call[400];
  char file_creation_time[MAX_STRING_LENGTH];
  char flash_version[MAX_STRING_LENGTH];
  char build_date[MAX_STRING_LENGTH];
  char build_dir[MAX_STRING_LENGTH];
  char build_machine[MAX_STRING_LENGTH];
  char cflags[400];
  char fflags[400];
  char setup_time_stamp[MAX_STRING_LENGTH];
  char build_time_stamp[MAX_STRING_LENGTH];
} sim_info_t;

/* define some particle property constants */

#if FLASH_NUMBER_OF_INT_PARTICLE_PROPS > 0
#define NUMINTPROPS  2*((FLASH_NUMBER_OF_INT_PARTICLE_PROPS+1)/2)
#else
#define NUMINTPROPS  2
#endif

#if FLASH_NUMBER_OF_REAL_PARTICLE_PROPS > 0
#define NUMREALPROPS FLASH_NUMBER_OF_REAL_PARTICLE_PROPS
#else
#define NUMREALPROPS 1
#endif

#define NUMACTUALINTPROPS FLASH_NUMBER_OF_INT_PARTICLE_PROPS
#define NUMACTUALREALPROPS FLASH_NUMBER_OF_REAL_PARTICLE_PROPS


/* set the dimension and grid variables -- the variable N_DIM is set 
   in the compile line */


/* mdim is the maximum dimension -- this is set in tree.fh */
#define MDIM 3
#define MGID 15

/* 3-d problem */
#if N_DIM == 3 

#define NDIM  3

#define NGID 15

#define k2d 1
#define k3d 1


/* 2-d problem */
#elif N_DIM == 2

#define NDIM  2

#define NGID 9

#define k2d 1
#define k3d 0


/* 1-d problem */
#else

#define NDIM 1

#define NGID 5

#define k2d 0
#define k3d 0

#endif


#endif

