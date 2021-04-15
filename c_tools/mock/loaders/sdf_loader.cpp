/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/mock/loaders/sdf_loader.cpp
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



#include <iostream>
#include <boost/format.hpp>
#include <vector>
#include <cassert>
#include <string>
#include "sdfloader_internal.hpp"
#include "simulation_loader.hpp"
#include <libSDF/SDF.h>

#undef SDFgetfloatOrDie
#undef SDFgetdoubleOrDie
#undef SDFgetintOrDie

using boost::format;

using namespace std;
using namespace CosmoTool;

static const double one_kpc = 3.08567802e16; /* km */
static const double one_Gyr = 3.1558149984e16; /* sec */

static void SDFgetfloatOrDie(SDF *sdfp, const char *name, float *v)
{
   if( SDFgetfloat(sdfp, name, v) ) {
     cerr << format("SDFgetfloat(%s) failed") % name << endl;
     abort();
   }
}

static void SDFgetdoubleOrDie(SDF *sdfp, const char *name, double *v)
{
   if( SDFgetdouble(sdfp, name, v) ) {
     cerr << format("SDFgetdouble(%s) failed") % name << endl;
     abort();
   }
}

static void SDFgetintOrDie(SDF *sdfp, const char *name, int *v)
{
   if( SDFgetint(sdfp, name, v) ) {
     cerr << format("SDFgetint(%s) failed") % name << endl;
     abort();
   }
}

class SDFLoader: public SimulationLoader
{
private:
  int load_flags;
  bool onefile;
  int _num_files;
  double unitMpc;
  SimuData *sdf_header;
  SDF *sdfp;
  SimulationPreprocessor *preproc;
  int num_splitting;
public:
  SDFLoader(SDF *sdf_file, SimuData *header, int flags, int num_splitting,
            SimulationPreprocessor *p)
    : sdfp(sdf_file), load_flags(flags),
      sdf_header(header), preproc(p)
  {
    this->num_splitting = num_splitting;
  }
  
  ~SDFLoader()
  {
    delete sdf_header;
  }
  
  SimuData *getHeader() {
    return sdf_header;
  }
  
  int num_files() {
    return num_splitting;
  }

  int64_t getStart(int id)
  {
    return sdf_header->TotalNumPart * int64_t(id) / num_splitting;
  }

  int64_t getNumberInSplit(int id)
  {
    return getStart(id+1)-getStart(id);
  }

  void rescaleParticles(SimuData *d,
                        double rescale_position,
                        double rescale_velocity)
  {
#define INFINITY std::numeric_limits<float>::max()
    float shift = 0.5*d->BoxSize;
    rescale_position /= d->time;
    float min_pos[3] = {INFINITY,INFINITY, INFINITY}, max_pos[3] = {-INFINITY,-INFINITY,-INFINITY};

    if (d->Pos[0] != 0)
      {
        for (int k = 0; k < 3; k++)
          {
            for (int64_t i = 0; i < d->NumPart; i++)
              {
                d->Pos[k][i] = (d->Pos[k][i]*rescale_position + shift);
                min_pos[k] = std::min(min_pos[k], d->Pos[k][i]);
                max_pos[k] = std::max(max_pos[k], d->Pos[k][i]);
              }
            cout << boost::format("min[%d] = %g, max[%d] = %g") % k % min_pos[k] % k %max_pos[k] << endl;
          }
      }
    if (d->Vel[0] != 0)
      {
        for (int k = 0; k < 3; k++)
          {
            for (int64_t i = 0; i < d->NumPart; i++)
              {
                d->Vel[k][i] *= rescale_velocity;
              }
          }
      }
  }


  void enforceBoxSize(SimuData *d)
  {
    for (int k = 0; k < 3; k++)
      {
        for (int64_t i = 0; i < d->NumPart; i++)
         {
           while (d->Pos[k][i]<0)
             d->Pos[k][i] += d->BoxSize;
           while (d->Pos[k][i]>=d->BoxSize)
             d->Pos[k][i] -= d->BoxSize;
         }
       }
  }

  SimuData *loadFile(int id) {
    SimuData *d;
    int64_t starts[7];
    int nobjs[7];
    int strides[7];
    void *addrs[7];
    const char *names[7];
    int ret;

    if (id >= num_splitting)
      return 0;

    d = new SimuData;
    d->BoxSize = sdf_header->BoxSize;
    d->time = sdf_header->time;
    d->Hubble = sdf_header->Hubble;
    d->Omega_M = sdf_header->Omega_M;
    d->Omega_Lambda = sdf_header->Omega_Lambda;
    
    d->NumPart = getNumberInSplit(id);

    int64_t numPartToLoad = getNumberInSplit(id);
    int64_t start = getStart(id);
    int k = 0;
#define SETUP_READ(name, vec) \
    names[k] = name, starts[k] = start, nobjs[k] = numPartToLoad, strides[k] = sizeof(vec[0]), addrs[k] = vec, k++;

    if (load_flags & NEED_POSITION)
      {
        const char *name_template[3] = { "x", "y", "z" };
        for (int q = 0; q < 3; q++)
          {
            d->Pos[q] = new float[numPartToLoad];
            SETUP_READ(name_template[q], d->Pos[q]);
          }
      }
    if (load_flags & NEED_VELOCITY)
      {
        const char *name_template[3] = { "vx", "vy", "vz" };
        for (int q = 0; q < 3; q++)
          {
            d->Vel[q] = new float[numPartToLoad];
            SETUP_READ(name_template[q], d->Vel[q]);
          }
      }
    if (load_flags & NEED_GADGET_ID)
      {
        d->Id = new int64_t[numPartToLoad];
        SETUP_READ("ident", d->Id);
      }
#undef SETUP_READ

    ret = SDFseekrdvecsarr(sdfp, k, (char**)names, starts, nobjs, addrs, strides);
    if (ret != 0)
      {
        cerr << format("Splitting %d/%d, SDF read failure: %s") % id % num_splitting % SDFerrstring << endl;
        abort();
      }

    if (load_flags & (NEED_POSITION | NEED_VELOCITY))
      rescaleParticles(d, d->Hubble*1e-5, one_kpc/one_Gyr);

    enforceBoxSize(d);
    applyTransformations(d);
    basicPreprocessing(d, preproc);

    return d;
  }
};

SimulationLoader *sdfLoader(const std::string& snapshot, int flags, 
                            int num_splitting,
                            SimulationPreprocessor *p)
{
  SDF *sdfp;
  int fileversion;
  SimuData *hdr;
  int64_t gnobj;

  sdfp = SDFopen(0, snapshot.c_str());
  if (sdfp == 0)
    {
      return 0;
    }
   cout << "Loading SDF with artificial splitting " << num_splitting << endl;
  hdr = new SimuData;
  
  SDFgetintOrDefault(sdfp, "version", &fileversion, 1);
  double h0;
  if (fileversion == 1)
    {
      SDFgetfloatOrDie(sdfp, "Omega0", &hdr->Omega_M);
      SDFgetfloatOrDie(sdfp, "Lambda_prime", &hdr->Omega_Lambda);
      SDFgetfloatOrDie(sdfp, "H0", &hdr->Hubble);
      hdr->Hubble *= 1000.*(one_kpc/one_Gyr);
      h0 = hdr->Hubble / 100.;
    }
  else
    {
      float Or, Of;
      SDFgetfloatOrDie(sdfp, "Omega0_m", &hdr->Omega_M);
      SDFgetfloatOrDie(sdfp, "Omega0_lambda", &hdr->Omega_Lambda);
      SDFgetfloatOrDie(sdfp, "Omega0_r", &Or);
      SDFgetfloatOrDie(sdfp, "Omega0_fld", &Of);

      hdr->Omega_M += Or;
      hdr->Omega_Lambda += Of;
      SDFgetfloatOrDie(sdfp, "h_100", &hdr->Hubble);
      hdr->Hubble *= 100.;
      h0 = hdr->Hubble / 100.; 
    }

  if (SDFhasname("R0", sdfp))
    {
      double R0;

      SDFgetdoubleOrDie(sdfp, "R0", &R0);
      hdr->BoxSize = 2.0*0.001*R0*h0;
    }

  if (SDFhasname("Rx", sdfp))
    {
      double R0;

      SDFgetdoubleOrDie(sdfp, "Rx", &R0);
      hdr->BoxSize = 2.0*0.001*R0*h0;
    }

  if (SDFhasname("redshift", sdfp))
    {
      double redshift;

      SDFgetdoubleOrDie(sdfp, "redshift", &redshift);
      hdr->time = 1/(1+redshift);
    }
   
  if (SDFgetint64(sdfp, "npart", &gnobj))
    {
      gnobj = SDFnrecs("x", sdfp);
      cerr << format("[Warning] No 'npart' found in SDF file '%s', guessing npart=%ld from SDFnrecs") % snapshot % gnobj << endl;
    }
  hdr->NumPart = hdr->TotalNumPart = gnobj;

  return new SDFLoader(sdfp, hdr, flags, num_splitting, p);
}


void sdfLoaderInit(int& argc, char **& argv)
{
}
