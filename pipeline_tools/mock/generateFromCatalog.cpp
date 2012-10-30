#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <boost/format.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "generateFromCatalog_conf.h"
#include "contour_pixels.hpp"
#include <netcdfcpp.h>
#include <CosmoTool/fortran.hpp>

using namespace std;
using boost::format;
using namespace CosmoTool;

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

struct Position
{
  double xyz[3];
};

struct ParticleData
{
  vector<int> id_gal;
  int id_mask;
  vector<Position> pos;
  double box[3][2];
  double Lmax;
};

typedef vector<NYU_Data> NYU_VData;

void loadData(const string& fname, NYU_VData & data)
{
  ifstream f(fname.c_str());

  while (!f.eof())
    {
      NYU_Data d;
      f >> d.index >> d.sector >> d.region >> d.ra >> d.dec >> d.cz >> d.fgotten >> d.phi_z;
      data.push_back(d);
    }
}

void generateGalaxiesInCube(NYU_VData& data, ParticleData& output_data)
{
  double d2r = M_PI/180;

  output_data.pos.resize(data.size());
  output_data.id_gal.resize(data.size());

  for (int j = 0; j < 3; j++)
    {
      output_data.box[j][0] = -INFINITY;
      output_data.box[j][1] = INFINITY;
    }

  for (int i = 0; i < data.size(); i++)
    {  
      double ra = data[i].ra*d2r, dec = data[i].dec*d2r;
      Position& p = output_data.pos[i];

      p.xyz[0] = data[i].cz*cos(ra)*cos(dec);
      p.xyz[1] = data[i].cz*sin(ra)*cos(dec);
      p.xyz[2] = data[i].cz*sin(dec);
      output_data.id_gal[i] = data[i].index;

      for (int j = 0; j < 3; j++)
	{
	  if (p.xyz[j] > output_data.box[j][0])
	    output_data.box[j][0] = p.xyz[j];
	  if (p.xyz[j] < output_data.box[j][1])
	    output_data.box[j][1] = p.xyz[j];	  
	}
    }
  cout << format("Galaxy position generated: %d galaxies") % output_data.pos.size() << endl;
  cout << format("box is %g < x < %g;  %g < y < %g;   %g < z < %g") 
    % (1e-2*output_data.box[0][1]) % (1e-2*output_data.box[0][0])
    % (1e-2*output_data.box[1][1]) % (1e-2*output_data.box[1][0])
    % (1e-2*output_data.box[2][1]) % (1e-2*output_data.box[2][0]) << endl;
  
}

static double cube(double x)
{
  return x*x*x;
}

void generateBoxMask(generateFromCatalog_info& args ,
                         Healpix_Map<float>& mask,
                         vector<int>& pixel_list,
                         NYU_VData& data,
                         ParticleData& output_data)
{
  int idx = -1;
  int insertion = 0;
  double volume = pixel_list.size()*1.0/mask.Npix()*4*M_PI;
  int numToInsert;

  idx = output_data.id_mask;

  cout << "Generate box mask..." << endl;

  double Rmax = output_data.Lmax, t = args.box_thickness_arg;
  output_data.Lmax += args.box_thickness_arg;


  volume *= Rmax*Rmax/1e6 * args.box_thickness_arg;
  numToInsert = (int)floor(volume*args.box_density_arg);

  for (int i = 0; i < numToInsert; i++)
    {
      Position p;
      bool stop_here;

      do
        {

          int p0 = (int)floor(drand48()*pixel_list.size());
          vec3 v = mask.pix2vec(pixel_list[p0]);
          double r = Rmax*pow(drand48()*(cube(1+t/Rmax)-1) + 1,1./3);

          p.xyz[0] = v.x * r;
          p.xyz[1] = v.y * r;
          p.xyz[2] = v.z * r;

          stop_here = true;
          for (int j = 0; j < 3; j++)
            {
              if (p.xyz[j] > output_data.box[j][0] ||
                  p.xyz[j] < output_data.box[j][1])
                stop_here = false;
            }
        }
      while (!stop_here);

      output_data.pos.push_back(p);
      output_data.id_gal.push_back(idx);
      insertion++;
    }
}

void generateSurfaceMask(generateFromCatalog_info& args ,
			 Healpix_Map<float>& mask, 
			 vector<int>& pixel_list, 
			 NYU_VData& data, 
			 ParticleData& output_data)
{
  // Find the first free index
  int idx = -1;
  int insertion = 0;
  double volume = pixel_list.size()*1.0/mask.Npix()*4*M_PI;
  int numToInsert;

  for (int i = 0; i < output_data.id_gal.size(); i++)
    {
      if (idx < output_data.id_gal[i])
	idx = output_data.id_gal[i]+1;
    }

  output_data.id_mask = idx;
  
  cout << "Generate surface mask..." << endl;

  double Rmax = -1;
  for (int j = 0; j < 3; j++)
    {
      Rmax = max(Rmax, max(output_data.box[j][0], -output_data.box[j][1]));
    }

  output_data.Lmax = Rmax;


  cout << format("Rmax is %g, surface volume is %g") % (Rmax/100) % (volume/(4*M_PI)) << endl;
  volume *= Rmax*Rmax*Rmax/3/1e6;
  numToInsert = (int)floor(volume*args.density_fake_arg);
  cout << format("3d volume to fill: %g (Mpc/h)^3") % volume << endl;

  cout << format("Will insert %d particles") % numToInsert << endl;

  double pct = 0;
  for (int i = 0; i < numToInsert; i++)
    {
      double new_pct = i*100./numToInsert;      

      if (new_pct-pct > 5.)
	{
	  pct = new_pct;
	  cout << format(" .. %3.0f %%") % pct << endl;
	}

      Position p;
      bool stop_here;

      do
	{

	  int p0 = (int)floor(drand48()*pixel_list.size());
	  vec3 v = mask.pix2vec(pixel_list[p0]);
	  double r = Rmax*pow(drand48(),1./3);
	  
	  p.xyz[0] = v.x * r;
	  p.xyz[1] = v.y * r;
	  p.xyz[2] = v.z * r;

	  stop_here = true;
	  for (int j = 0; j < 3; j++)
	    {
	      if (p.xyz[j] > output_data.box[j][0] ||
		  p.xyz[j] < output_data.box[j][1])
		stop_here = false;
	    }
	}
      while (!stop_here);
      
      output_data.pos.push_back(p);
      output_data.id_gal.push_back(idx);
      insertion++;
    }
  cout << format("Done. Inserted %d particles.") % insertion << endl;
}

void saveData(ParticleData& pdata)
{
  NcFile f("particles.nc", NcFile::Replace);
  
  assert(f.is_valid());

  NcDim *d = f.add_dim("space", 3);
  NcDim *p = f.add_dim("Np", pdata.pos.size());
  NcVar *v = f.add_var("particles", ncDouble, d, p);
  double *x = new double[pdata.pos.size()];

  for (int j = 0; j < 3; j++)
    {

      for (int i = 0; i < pdata.pos.size(); i++)
	x[i] = pdata.pos[i].xyz[j];

      v->put_rec(d, x, j);
    }

  v = f.add_var("id_gal", ncInt, p);
  v->put(&pdata.id_gal[0], pdata.id_gal.size());

  delete[] x;
  
}

void saveForZobov(ParticleData& pdata, const string& fname, const string& paramname)
{
  UnformattedWrite f(fname);
  static const char axis[] = { 'X', 'Y', 'Z' };
  double Lmax = pdata.Lmax;

  f.beginCheckpoint();
  f.writeInt32(pdata.pos.size());
  f.endCheckpoint();

  for (int j = 0; j < 3; j++)
    {
      cout << format("Writing %c components...") % axis[j] << endl;
      f.beginCheckpoint();
      for (uint32_t i = 0; i < pdata.pos.size(); i++)
	{
	  f.writeReal32((pdata.pos[i].xyz[j]+Lmax)/(2*Lmax));
	}
      f.endCheckpoint();
    }

  NcFile fp(paramname.c_str(), NcFile::Replace);

  fp.add_att("range_x_min", -Lmax);
  fp.add_att("range_x_max", Lmax);
  fp.add_att("range_y_min", -Lmax);
  fp.add_att("range_y_max", Lmax);
  fp.add_att("range_z_min", -Lmax);
  fp.add_att("range_z_max", Lmax);

  NcDim *NumPart_dim = fp.add_dim("numpart_dim", pdata.pos.size());
  NcVar *v = fp.add_var("particle_ids", ncInt, NumPart_dim);
  NcVar *v2 = fp.add_var("expansion", ncDouble, NumPart_dim);

  double *expansion_fac = new double[pdata.pos.size()];

  for (int i = 0; i <  pdata.pos.size(); i++)
    expansion_fac[i] = 1.0;

  v->put(&pdata.id_gal[0], pdata.id_gal.size());
  v2->put(expansion_fac, pdata.pos.size());

  delete[] expansion_fac;
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
  Healpix_Map<float> o_mask;
  vector<int> pixel_list;
  ParticleData output_data;

  loadData(args_info.catalog_arg, data);
  

  cout << "Loading mask..." << endl;
  read_Healpix_map_from_fits(args_info.mask_arg, o_mask);

  Healpix_Map<float> mask;

  mask.SetNside(128, RING);
  mask.Import(o_mask);

  computeContourPixels(mask,pixel_list);

  // We compute a cube holding all the galaxies + the survey surface mask

  generateGalaxiesInCube(data, output_data);
  generateSurfaceMask(args_info, mask, pixel_list, data, output_data);
  computeFilledPixels(mask,pixel_list);
  generateBoxMask(args_info, mask, pixel_list, data, output_data);
  
  saveForZobov(output_data, args_info.output_arg, args_info.params_arg);
  //  saveData(output_data);

  return 0;
}
