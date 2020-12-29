/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/mock/generateFromCatalog.cpp
    Copyright (C) 2010-2014 Guilhem Lavaux
    Copyright (C) 2011-2014 P. M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
+*/



#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <boost/format.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "generateFromCatalog_conf.h"
#include "contour_pixels.hpp"
#include <netcdf>
#include <CosmoTool/fortran.hpp>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>

#define LIGHT_SPEED 299792.458

using namespace std;
using boost::format;
using namespace CosmoTool;
using namespace netCDF;

struct NYU_Data
{
  int index;
  int sector;
  int region;
  double ra, dec;
  double cz;
  double fgotten;
  double phi_z;
  long uniqueID;
};

struct Position
{
  double xyz[3];
};

struct ParticleData
{
  vector<int> id_gal;
  vector<double> ra;
  vector<double> dec;
  vector<double> redshift;
  vector<double> catalogID;
  vector<long> uniqueID;
  int id_mask;
  // PMS
  int mask_index;
  // END PMS
  vector<Position> pos;
  double box[3][2];
  double Lmax;
};

typedef vector<NYU_Data> NYU_VData;

// this defines the expansion function that we will integrate
// Laveaux & Wandelt (2012) Eq. 24
struct my_expan_params { double Om; double w0; double wa; };
double expanFun (double z, void * p) {
  struct my_expan_params * params = (struct my_expan_params *)p;
  double Om = (params->Om);
  double w0 = (params->w0);
  double wa = (params->wa);

  //const double h0 = 1.0;
  const double h0 = 0.71;
  double ez;

  double wz = w0 + wa*z/(1+z);

  ez = Om*pow(1+z,3) + (1.-Om);
  //ez = Om*pow(1+z,3) + pow(h0,2) * (1.-Om)*pow(1+z,3+3*wz);

  ez = sqrt(ez);
  //ez = sqrt(ez)/h0;

  ez = 1./ez;

  return  ez;
}

void loadData(const string& fname, NYU_VData & data)
{
  ifstream f(fname.c_str());
   
  while (!f.eof())
    {
      NYU_Data d;
      f >> d.index >> d.sector >> d.region >> d.ra >> d.dec >> d.cz >> d.fgotten >> d.phi_z;
      // Maubert - avoid double counting of last particle in array if EOF is after a blank line
      if (!f) 
      { 
      continue; 
      } 
      // End Maubert
      d.uniqueID = d.index;
      data.push_back(d);
    }
}

void generateGalaxiesInCube(NYU_VData& data, ParticleData& output_data, 
                            bool useComoving, double omegaM)
{
  double d2r = M_PI/180;

  gsl_function expanF;
  expanF.function = &expanFun;
  struct my_expan_params expanParams;
  double maxZ = 2.0, z, result, error, *dL, *redshifts;
  int numZ = 1000, iZ;
  size_t nEval;

  expanParams.Om = omegaM;
  expanParams.w0 = -1.0;
  expanParams.wa = 0.0;
  expanF.params = &expanParams;

  dL = (double *) malloc(numZ * sizeof(double));
  redshifts = (double *) malloc(numZ * sizeof(double));

  for (iZ = 0; iZ < numZ; iZ++) {
    z = iZ * maxZ/numZ;
    //gsl_integration_qng(&expanF, 0.0, z, 1.e-6, 1.e-6, &result, &error, &nEval);
    dL[iZ] = result;
    redshifts[iZ] = z;
  }

  gsl_interp *interp = gsl_interp_alloc(gsl_interp_linear, numZ); 
  gsl_interp_init(interp, redshifts, dL, numZ);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  output_data.pos.resize(data.size());
  output_data.id_gal.resize(data.size());
  output_data.ra.resize(data.size());
  output_data.dec.resize(data.size());
  output_data.redshift.resize(data.size());
  output_data.uniqueID.resize(data.size());

  for (int j = 0; j < 3; j++)
    {
      output_data.box[j][0] = -INFINITY;
      output_data.box[j][1] = INFINITY;
    }

  for (int i = 0; i < data.size(); i++)
    {  
      double ra = data[i].ra*d2r, dec = data[i].dec*d2r;
      Position& p = output_data.pos[i];

      if (useComoving) {
        //double pos = gsl_interp_eval(interp, redshifts, dL, data[i].cz, acc);
        // Maubert - Lower boundary in redshift set to 0 to be consistent with pruneVoids (was 1.e-6 before).
        gsl_integration_qng(&expanF, 0.0, data[i].cz/LIGHT_SPEED,
                                        1.e-6, 
                                        1.e-6, &result, &error, &nEval);
        double Dc = result*LIGHT_SPEED;
        p.xyz[0] = Dc*cos(ra)*cos(dec);
        p.xyz[1] = Dc*sin(ra)*cos(dec);
        p.xyz[2] = Dc*sin(dec);
      } else {
        p.xyz[0] = data[i].cz*cos(ra)*cos(dec);
        p.xyz[1] = data[i].cz*sin(ra)*cos(dec);
        p.xyz[2] = data[i].cz*sin(dec);
      }
//printf("CREATE %e %e\n", data[i].cz, sqrt(p.xyz[0]*p.xyz[0] + p.xyz[1]*p.xyz[1] + p.xyz[2]*p.xyz[2]));
      output_data.id_gal[i] = data[i].index;
      output_data.ra[i] = ra;
      output_data.dec[i] = dec;
      output_data.redshift[i] = data[i].cz;
      output_data.uniqueID[i] = data[i].uniqueID;

      for (int j = 0; j < 3; j++)
	{
	  if (p.xyz[j] > output_data.box[j][0])
	    output_data.box[j][0] = p.xyz[j];
	  if (p.xyz[j] < output_data.box[j][1])
	    output_data.box[j][1] = p.xyz[j];	  
	}
//printf("INSERT GAL %d %e %e %e\n", output_data.id_gal[i], p.xyz[0], p.xyz[1], p.xyz[2]);
    }

  // normalize box
  float left = 1.e99;
  float right = -1.e99;
  for (int j = 0; j < 3; j++) {
    if (output_data.box[j][1] < left) left = output_data.box[j][1];
    if (output_data.box[j][0] > right) right = output_data.box[j][0];
  }
  for (int j = 0; j < 3; j++) {
    output_data.box[j][1] = left;
    output_data.box[j][0] = right;
  }

  cout << format("Galaxy position generated: %d galaxies") % output_data.pos.size() << endl;
  cout << format("box is %g < x < %g;  %g < y < %g;   %g < z < %g") 
    % (1e-2*output_data.box[0][1]) % (1e-2*output_data.box[0][0])
    % (1e-2*output_data.box[1][1]) % (1e-2*output_data.box[1][0])
    % (1e-2*output_data.box[2][1]) % (1e-2*output_data.box[2][0]) << endl;
 
  gsl_interp_free(interp); 
}

void generateSurfaceMask(generateFromCatalog_info& args ,
			 Healpix_Map<float>& mask, 
			 vector<int>& pixel_list, 
			 vector<int>& full_mask_list, 
			 NYU_VData& data, 
			 ParticleData& output_data,
			 bool useComoving, 
			 double omegaM)
{

  //Maubert - Needed for comobile distance in mock_sphere
  gsl_function expanF;
  expanF.function = &expanFun;
  struct my_expan_params expanParams;
  expanParams.Om = omegaM;
  expanParams.w0 = -1.0;
  expanParams.wa = 0.0;
  expanF.params = &expanParams;
  
  double result, error ;
  size_t nEval;
  //End Maubert - Needed for comobile distance in mock_sphere

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
 
// PMS
  output_data.mask_index = output_data.id_gal.size();
// END PMS 
  cout << "Generate surface mask..." << endl;

  double Rmax = -1;
  for (int j = 0; j < 3; j++)
    {
      Rmax = max(Rmax, max(output_data.box[j][0], -output_data.box[j][1]));
    }

  output_data.Lmax = Rmax;

// PMS - write a small text file with galaxy position
  FILE *fp;
  fp = fopen("galaxies.txt", "w");
  for (int i = 0; i < data.size(); i++) {  
    Position& p = output_data.pos[i];
    fprintf(fp, "%e %e %e\n", 
            (p.xyz[0]), 
            (p.xyz[1]), 
            (p.xyz[2]));
  }
  fclose(fp);
// END PMS

  cout << format("Rmax is %g, surface volume is %g") % (Rmax/100) % (volume/(4*M_PI)) << endl;
  volume *= Rmax*Rmax*Rmax/3/1e6;
  numToInsert = (int)floor(volume*args.density_fake_arg);
  cout << format("3d volume to fill: %g (Mpc/h)^3") % volume << endl;

  cout << format("Will insert %d particles") % numToInsert << endl;

  fp = fopen("mock_galaxies.txt", "w");

  double pct = 0;
  for (int i = 0; i < numToInsert; i++) {
    double new_pct = i*100./numToInsert;      

    if (new_pct-pct > 5.) {
	    pct = new_pct;
	    cout << format(" .. %3.0f %%") % pct << endl;
 	  }

    Position p;
    bool stop_here;

    do {
	    int p0 = (int)floor(drand48()*pixel_list.size());
	    vec3 v = mask.pix2vec(pixel_list[p0]);
	    double r = Rmax*pow(drand48(),1./3);
	  
	    p.xyz[0] = v.x * r;
	    p.xyz[1] = v.y * r;
	    p.xyz[2] = v.z * r;

	    stop_here = true;
	    for (int j = 0; j < 3; j++) {
	      if (p.xyz[j] > output_data.box[j][0] ||
		        p.xyz[j] < output_data.box[j][1])
		      stop_here = false;
	    }
	  }
    while (!stop_here);
      
// PMS : write mock galaxies to a small file for plotting
    fprintf(fp, "%e %e %e\n", 
            (p.xyz[0]), 
            (p.xyz[1]), 
            (p.xyz[2]));
// END PMS
    output_data.pos.push_back(p);
    output_data.id_gal.push_back(idx);
    output_data.ra.push_back(-1);
    output_data.dec.push_back(-1);
    output_data.redshift.push_back(-1);
    output_data.uniqueID.push_back(-1);
//printf("INSERT MOCK %d %e %e %e\n", idx, p.xyz[0], p.xyz[1], p.xyz[2]);
    insertion++;
  }

  fclose(fp);

  // PMS
  // TEST - insert mock galaxies along box edge
  fp = fopen("mock_boundary.txt", "w");
  double dx[3];
  dx[0] = output_data.box[0][1] - output_data.box[0][0];
  dx[1] = output_data.box[1][1] - output_data.box[1][0];
  dx[2] = output_data.box[2][1] - output_data.box[2][0];

  int nPart = 100;
// TEST
  for (int iDir = 0; iDir < 0; iDir++) {
  for (int iFace = 0; iFace < 0; iFace++) {
  //for (int iDir = 0; iDir < 3; iDir++) {
  //for (int iFace = 0; iFace < 2; iFace++) {

    int iy = (iDir + 1) % 3;
    int iz = (iDir + 2) % 3;

    for (int i = 0; i < nPart; i++) {
    for (int j = 0; j < nPart; j++) {
      Position p;

      p.xyz[iDir] = output_data.box[iDir][iFace];
      p.xyz[iy] = i * dx[iy]/nPart + output_data.box[iy][0];
      p.xyz[iz] = j * dx[iz]/nPart + output_data.box[iz][0];

      output_data.pos.push_back(p);
      output_data.id_gal.push_back(idx);
      output_data.ra.push_back(-1);
      output_data.dec.push_back(-1);
      output_data.redshift.push_back(-1);
      output_data.uniqueID.push_back(-1);
      insertion++;

      fprintf(fp, "%e %e %e\n", 
              (p.xyz[0]),
              (p.xyz[1]), 
              (p.xyz[2]));
    } 
    }
  } 
  } 
  fclose(fp);
  // END PMS

   // PMS
  // TEST - insert mock galaxies along spheres of survey redshift boundaries
  fp = fopen("mock_sphere.txt", "w");
  //Maubert - insert mock galaxies according to useComoving specification
  // Added & compute rmin & rmax out of loop
  double rmin ;
  double rmax ; 
  if (useComoving) {
        gsl_integration_qng(&expanF, 0.0, args.zMin_arg,
                                        1.e-6, 
                                        1.e-6, &result, &error, &nEval);
        rmin = result* LIGHT_SPEED;
        
        gsl_integration_qng(&expanF, 0.0, args.zMax_arg,
                                        1.e-6, 
                                        1.e-6, &result, &error, &nEval);
        rmax = result* LIGHT_SPEED;
  } else {
         rmin = args.zMin_arg * LIGHT_SPEED;
         rmax = args.zMax_arg * LIGHT_SPEED;
  }
    
  //for (int q = 0; q < 0; q++) {
  for (int q = 0; q < full_mask_list.size(); q++) {
    vec3 v = mask.pix2vec(full_mask_list[q]);
	  
    Position p;


    if (rmin > 0.) {
      p.xyz[0] = v.x * rmin;
      p.xyz[1] = v.y * rmin;
      p.xyz[2] = v.z * rmin;
      output_data.pos.push_back(p);
      output_data.id_gal.push_back(idx);
      output_data.ra.push_back(-1);
      output_data.dec.push_back(-1);
      output_data.redshift.push_back(-1);
      output_data.uniqueID.push_back(-1);
      insertion++;
      fprintf(fp, "%e %e %e\n", 
              (p.xyz[0]),
              (p.xyz[1]), 
              (p.xyz[2]));
    }


    p.xyz[0] = v.x * rmax;
    p.xyz[1] = v.y * rmax;
    p.xyz[2] = v.z * rmax;
    output_data.pos.push_back(p);
    output_data.id_gal.push_back(idx);
    output_data.ra.push_back(-1);
    output_data.dec.push_back(-1);
    output_data.redshift.push_back(-1);
    output_data.uniqueID.push_back(-1);
    insertion++;
    fprintf(fp, "%e %e %e\n", 
            (p.xyz[0]),
            (p.xyz[1]), 
            (p.xyz[2]));
  }
  fclose(fp);
  // END PMS

  cout << format("Done. Inserted %d particles.") % insertion << endl;
}

void saveData(ParticleData& pdata)
{
  NcFile f("particles.nc", NcFile::replace);

  NcDim d = f.addDim("space", 3);
  NcDim p = f.addDim("Np", pdata.pos.size());
  NcVar v = f.addVar("particles", ncDouble, {d, p});
  double *x = new double[pdata.pos.size()];

  for (int j = 0; j < 3; j++)
    {

      for (int i = 0; i < pdata.pos.size(); i++)
	x[i] = pdata.pos[i].xyz[j];

      v.putVar({size_t(j), 0}, {1, pdata.pos.size()}, x);
    }

  v = f.addVar("id_gal", ncInt, std::vector<NcDim>({p}));
  v.putVar(&pdata.id_gal[0]);

  delete[] x;
  
}

void saveForZobov(ParticleData& pdata, const string& fname, const string& paramname)
{
  UnformattedWrite f(fname);
  static const char axis[] = { 'X', 'Y', 'Z' };
  double Lmax = pdata.Lmax;
  double r2d = 180./M_PI;

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

  cout << format("Writing RA...")  << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < pdata.pos.size(); i++) {
	  f.writeReal32(pdata.ra[i]*r2d);
	}
  f.endCheckpoint();
   
  cout << format("Writing Dec...")  << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < pdata.pos.size(); i++) {
	  f.writeReal32(pdata.dec[i]*r2d);
	}
  f.endCheckpoint();
    
  cout << format("Writing Redshift...")  << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < pdata.pos.size(); i++) {
	  f.writeReal32(pdata.redshift[i]);
	}
  f.endCheckpoint();
   
  cout << format("Writing Unique ID...")  << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < pdata.pos.size(); i++) {
	  f.writeInt64(pdata.uniqueID[i]);
	}
  f.endCheckpoint();
   
  NcFile fp(paramname.c_str(), NcFile::replace);

  fp.putAtt("range_x_min", ncDouble, -Lmax/100.);
  fp.putAtt("range_x_max", ncDouble, Lmax/100.);
  fp.putAtt("range_y_min", ncDouble, -Lmax/100.);
  fp.putAtt("range_y_max", ncDouble, Lmax/100.);
  fp.putAtt("range_z_min", ncDouble, -Lmax/100.);
  fp.putAtt("range_z_max", ncDouble, Lmax/100.);
  fp.putAtt("mask_index", ncInt, pdata.mask_index); // PMS
  fp.putAtt("is_observation", ncInt, 1); // PMS

  int nOutputPart = pdata.mask_index;
  //int nOutputPart = pdata.pos.size();

  NcDim NumPart_dim = fp.addDim("numpart_dim", nOutputPart);
  NcVar v = fp.addVar("particle_ids", ncInt, NumPart_dim);
  //NcVar v2 = fp.addVar("expansion", ncDouble, NumPart_dim);

  //double *expansion_fac = new double[pdata.pos.size()];

  //for (int i = 0; i <  pdata.pos.size(); i++)
  //  expansion_fac[i] = 1.0;

  v.putVar({0}, {size_t(nOutputPart)}, &pdata.id_gal[0]);
  //v2->put(expansion_fac, pdata.pos.size());

  //delete[] expansion_fac;

/*
  FILE *infoFile = fopen("sample_info.txt", "w");
  fprintf(infoFile, "x_min = %f\n", -Lmax/100.);  
  fprintf(infoFile, "x_max = %f\n", Lmax/100.);
  fprintf(infoFile, "y_min = %f\n", -Lmax/100.);  
  fprintf(infoFile, "y_max = %f\n", Lmax/100.);  
  fprintf(infoFile, "z_min = %f\n", -Lmax/100.);
  fprintf(infoFile, "z_max = %f\n", Lmax/100.);
  fprintf(infoFile, "mask_index = %d\n", pdata.mask_index);
  fprintf(infoFile, "total_particles = %d\n", pdata.pos.size());
  fclose(infoFile);
*/

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

  cout << "Loading data " << args_info.catalog_arg << "..." << endl;
  vector<NYU_Data> data;
  Healpix_Map<float> o_mask;
  vector<int> pixel_list;
  vector<int> full_mask_list;
  ParticleData output_data;

  loadData(args_info.catalog_arg, data);
  

  cout << "Loading mask..." << endl;
  read_Healpix_map_from_fits(args_info.mask_arg, o_mask);

  Healpix_Map<float> mask;

  mask.SetNside(128, RING);
  mask.Import(o_mask);

  computeContourPixels(mask,pixel_list);
  computeMaskPixels(mask,full_mask_list);

  // We compute a cube holding all the galaxies + the survey surface mask

  cout << "Placing galaxies..." << endl;
  generateGalaxiesInCube(data, output_data, args_info.useComoving_flag, 
                         args_info.omegaM_arg);
  generateSurfaceMask(args_info, mask, pixel_list, full_mask_list,
                      data, output_data,args_info.useComoving_flag, 
                         args_info.omegaM_arg);
  
  saveForZobov(output_data, args_info.output_arg, args_info.params_arg);
  //  saveData(output_data);
  
  // PMS
  FILE *fp = fopen("mask_index.txt", "w");
  fprintf(fp, "%d", output_data.mask_index);
  fclose(fp);

  fp = fopen("total_particles.txt", "w");
  fprintf(fp, "%d", output_data.pos.size());
  fclose(fp);
  printf("Done!\n");
  // END PMS
  return 0;
}
