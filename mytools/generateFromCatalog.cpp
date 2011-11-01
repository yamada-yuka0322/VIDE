#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <boost/format.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "generateFromCatalog_conf.h"
#include "contour_pixels.hpp"

using namespace std;
using boost::format;

struct NYU_Data
{
  int index;
  int sector;
  int region;
  double ra, dec;
  double cz;
  double fgotten;
  double phi_z;
};

void loadData(const string& fname, vector<NYU_Data> & data)
{
  ifstream f(fname.c_str());

  while (!f.eof())
    {
      NYU_Data d;
      f >> d.index >> d.sector >> d.region >> d.ra >> d.dec >> d.cz >> d.fgotten >> d.phi_z;
      data.push_back(d);
    }
}


int main(int argc, char **argv)
{
  generateFromCatalog_info args_info;
  generateFromCatalog_conf_params args_params;
 
  generateFromCatalog_conf_init(&args_info);
  generateFromCatalog_conf_params_init(&args_params);
  
  args_params.check_required = 0;
  if (generateFromCatalog_conf_ext (argc, argv, &args_info, &args_params))
    return 1;
  
  if (!args_info.configFile_given)
    {
      if (generateFromCatalog_conf_required (&args_info, GENERATEFROMCATALOG_CONF_PACKAGE))
        return 1;
    }
  else
    {
      args_params.check_required = 1;
      args_params.initialize = 0;
      if (generateFromCatalog_conf_config_file (args_info.configFile_arg,
					 &args_info,
					 &args_params))
	return 1;
    }
  
  generateFromCatalog_conf_print_version();

  cout << "Loading NYU data..." << endl;
  vector<NYU_Data> data;
  Healpix_Map<float> mask;
  vector<int> pixel_list;
  loadData(args_info.catalog_arg, data);

  cout << "Loading mask..." << endl;
  read_Healpix_map_from_fits(args_info.mask_arg, mask);

  computeContourPixels(mask,pixel_list);

  return 0;
}
