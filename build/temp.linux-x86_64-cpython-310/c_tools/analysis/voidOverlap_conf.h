/** @file voidOverlap_conf.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.5
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef VOIDOVERLAP_CONF_H
#define VOIDOVERLAP_CONF_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef VOIDOVERLAP_CONF_PACKAGE
/** @brief the program name (used for printing errors) */
#define VOIDOVERLAP_CONF_PACKAGE "voidOverlap"
#endif

#ifndef VOIDOVERLAP_CONF_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define VOIDOVERLAP_CONF_PACKAGE_NAME "voidOverlap"
#endif

#ifndef VOIDOVERLAP_CONF_VERSION
/** @brief the program version */
#define VOIDOVERLAP_CONF_VERSION "0"
#endif

/** @brief Where the command line options are stored */
struct voidOverlap_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * configFile_arg;	/**< @brief Configuration file path.  */
  char * configFile_orig;	/**< @brief Configuration file path original value given at command line.  */
  const char *configFile_help; /**< @brief Configuration file path help description.  */
  char * partFile1_arg;	/**< @brief Particle file for catalog 1.  */
  char * partFile1_orig;	/**< @brief Particle file for catalog 1 original value given at command line.  */
  const char *partFile1_help; /**< @brief Particle file for catalog 1 help description.  */
  char * volFile1_arg;	/**< @brief Volume file for catalog 1.  */
  char * volFile1_orig;	/**< @brief Volume file for catalog 1 original value given at command line.  */
  const char *volFile1_help; /**< @brief Volume file for catalog 1 help description.  */
  char * voidFile1_arg;	/**< @brief Void info file for catalog 1.  */
  char * voidFile1_orig;	/**< @brief Void info file for catalog 1 original value given at command line.  */
  const char *voidFile1_help; /**< @brief Void info file for catalog 1 help description.  */
  char * infoFile1_arg;	/**< @brief Extra info file for catalog 1.  */
  char * infoFile1_orig;	/**< @brief Extra info file for catalog 1 original value given at command line.  */
  const char *infoFile1_help; /**< @brief Extra info file for catalog 1 help description.  */
  char * zoneFile1_arg;	/**< @brief Zone file for catalog 1.  */
  char * zoneFile1_orig;	/**< @brief Zone file for catalog 1 original value given at command line.  */
  const char *zoneFile1_help; /**< @brief Zone file for catalog 1 help description.  */
  char * zonePartFile1_arg;	/**< @brief Zone-particle file for catalog 1.  */
  char * zonePartFile1_orig;	/**< @brief Zone-particle file for catalog 1 original value given at command line.  */
  const char *zonePartFile1_help; /**< @brief Zone-particle file for catalog 1 help description.  */
  char * centerFile1_arg;	/**< @brief Barycenter file for catalog 1.  */
  char * centerFile1_orig;	/**< @brief Barycenter file for catalog 1 original value given at command line.  */
  const char *centerFile1_help; /**< @brief Barycenter file for catalog 1 help description.  */
  char * shapeFile1_arg;	/**< @brief Shape file for catalog 1.  */
  char * shapeFile1_orig;	/**< @brief Shape file for catalog 1 original value given at command line.  */
  const char *shapeFile1_help; /**< @brief Shape file for catalog 1 help description.  */
  char * partFile2_arg;	/**< @brief Particle file for catalog 2.  */
  char * partFile2_orig;	/**< @brief Particle file for catalog 2 original value given at command line.  */
  const char *partFile2_help; /**< @brief Particle file for catalog 2 help description.  */
  char * volFile2_arg;	/**< @brief Volume file for catalog 2.  */
  char * volFile2_orig;	/**< @brief Volume file for catalog 2 original value given at command line.  */
  const char *volFile2_help; /**< @brief Volume file for catalog 2 help description.  */
  char * voidFile2_arg;	/**< @brief Void info file for catalog 2.  */
  char * voidFile2_orig;	/**< @brief Void info file for catalog 2 original value given at command line.  */
  const char *voidFile2_help; /**< @brief Void info file for catalog 2 help description.  */
  char * infoFile2_arg;	/**< @brief Extra info file for catalog 2.  */
  char * infoFile2_orig;	/**< @brief Extra info file for catalog 2 original value given at command line.  */
  const char *infoFile2_help; /**< @brief Extra info file for catalog 2 help description.  */
  char * zoneFile2_arg;	/**< @brief Zone file for catalog 2.  */
  char * zoneFile2_orig;	/**< @brief Zone file for catalog 2 original value given at command line.  */
  const char *zoneFile2_help; /**< @brief Zone file for catalog 2 help description.  */
  char * zonePartFile2_arg;	/**< @brief Zone-particle file for catalog 2.  */
  char * zonePartFile2_orig;	/**< @brief Zone-particle file for catalog 2 original value given at command line.  */
  const char *zonePartFile2_help; /**< @brief Zone-particle file for catalog 2 help description.  */
  char * centerFile2_arg;	/**< @brief Barycenter file for catalog 2.  */
  char * centerFile2_orig;	/**< @brief Barycenter file for catalog 2 original value given at command line.  */
  const char *centerFile2_help; /**< @brief Barycenter file for catalog 2 help description.  */
  char * shapeFile2_arg;	/**< @brief Shape file for catalog 2.  */
  char * shapeFile2_orig;	/**< @brief Shape file for catalog 2 original value given at command line.  */
  const char *shapeFile2_help; /**< @brief Shape file for catalog 2 help description.  */
  int isObservation_flag;	/**< @brief We are working with observational data (default=off).  */
  const char *isObservation_help; /**< @brief We are working with observational data help description.  */
  char * outfile_arg;	/**< @brief Output file.  */
  char * outfile_orig;	/**< @brief Output file original value given at command line.  */
  const char *outfile_help; /**< @brief Output file help description.  */
  int useID_flag;	/**< @brief Use unique catalog ID to match voids; otherwise use volumes (default=off).  */
  const char *useID_help; /**< @brief Use unique catalog ID to match voids; otherwise use volumes help description.  */
  double overlapFrac_arg;	/**< @brief threshold fraction of voronoi radius to count as matched (default='0.25').  */
  char * overlapFrac_orig;	/**< @brief threshold fraction of voronoi radius to count as matched original value given at command line.  */
  const char *overlapFrac_help; /**< @brief threshold fraction of voronoi radius to count as matched help description.  */
  char * periodic_arg;	/**< @brief Set of edges which are periodic (default='xy').  */
  char * periodic_orig;	/**< @brief Set of edges which are periodic original value given at command line.  */
  const char *periodic_help; /**< @brief Set of edges which are periodic help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int configFile_given ;	/**< @brief Whether configFile was given.  */
  unsigned int partFile1_given ;	/**< @brief Whether partFile1 was given.  */
  unsigned int volFile1_given ;	/**< @brief Whether volFile1 was given.  */
  unsigned int voidFile1_given ;	/**< @brief Whether voidFile1 was given.  */
  unsigned int infoFile1_given ;	/**< @brief Whether infoFile1 was given.  */
  unsigned int zoneFile1_given ;	/**< @brief Whether zoneFile1 was given.  */
  unsigned int zonePartFile1_given ;	/**< @brief Whether zonePartFile1 was given.  */
  unsigned int centerFile1_given ;	/**< @brief Whether centerFile1 was given.  */
  unsigned int shapeFile1_given ;	/**< @brief Whether shapeFile1 was given.  */
  unsigned int partFile2_given ;	/**< @brief Whether partFile2 was given.  */
  unsigned int volFile2_given ;	/**< @brief Whether volFile2 was given.  */
  unsigned int voidFile2_given ;	/**< @brief Whether voidFile2 was given.  */
  unsigned int infoFile2_given ;	/**< @brief Whether infoFile2 was given.  */
  unsigned int zoneFile2_given ;	/**< @brief Whether zoneFile2 was given.  */
  unsigned int zonePartFile2_given ;	/**< @brief Whether zonePartFile2 was given.  */
  unsigned int centerFile2_given ;	/**< @brief Whether centerFile2 was given.  */
  unsigned int shapeFile2_given ;	/**< @brief Whether shapeFile2 was given.  */
  unsigned int isObservation_given ;	/**< @brief Whether isObservation was given.  */
  unsigned int outfile_given ;	/**< @brief Whether outfile was given.  */
  unsigned int useID_given ;	/**< @brief Whether useID was given.  */
  unsigned int overlapFrac_given ;	/**< @brief Whether overlapFrac was given.  */
  unsigned int periodic_given ;	/**< @brief Whether periodic was given.  */

} ;

/** @brief The additional parameters to pass to parser functions */
struct voidOverlap_conf_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure voidOverlap_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure voidOverlap_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *voidOverlap_info_purpose;
/** @brief the usage string of the program */
extern const char *voidOverlap_info_usage;
/** @brief all the lines making the help output */
extern const char *voidOverlap_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int voidOverlap_conf (int argc, char **argv,
  struct voidOverlap_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use voidOverlap_conf_ext() instead
 */
int voidOverlap_conf2 (int argc, char **argv,
  struct voidOverlap_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int voidOverlap_conf_ext (int argc, char **argv,
  struct voidOverlap_info *args_info,
  struct voidOverlap_conf_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int voidOverlap_conf_dump(FILE *outfile,
  struct voidOverlap_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int voidOverlap_conf_file_save(const char *filename,
  struct voidOverlap_info *args_info);

/**
 * Print the help
 */
void voidOverlap_conf_print_help(void);
/**
 * Print the version
 */
void voidOverlap_conf_print_version(void);

/**
 * Initializes all the fields a voidOverlap_conf_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void voidOverlap_conf_params_init(struct voidOverlap_conf_params *params);

/**
 * Allocates dynamically a voidOverlap_conf_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized voidOverlap_conf_params structure
 */
struct voidOverlap_conf_params *voidOverlap_conf_params_create(void);

/**
 * Initializes the passed voidOverlap_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void voidOverlap_conf_init (struct voidOverlap_info *args_info);
/**
 * Deallocates the string fields of the voidOverlap_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void voidOverlap_conf_free (struct voidOverlap_info *args_info);

/**
 * The config file parser (deprecated version)
 * @param filename the name of the config file
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use voidOverlap_conf_config_file() instead
 */
int voidOverlap_conf_configfile (const char *filename,
  struct voidOverlap_info *args_info,
  int override, int initialize, int check_required);

/**
 * The config file parser
 * @param filename the name of the config file
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int voidOverlap_conf_config_file (const char *filename,
  struct voidOverlap_info *args_info,
  struct voidOverlap_conf_params *params);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int voidOverlap_conf_required (struct voidOverlap_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* VOIDOVERLAP_CONF_H */
