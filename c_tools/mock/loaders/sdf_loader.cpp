/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/mock/loaders/sdf_loader.cpp
    Copyright (C) 2010-2013 Guilhem Lavaux
    Copyright (C) 2011-2013 P. M. Sutter

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
#include <libsdf/mpmy.h>
#include <libsdf/SDF.h>
#include <libsdf/error.h>
#include <libsdf/cosmo.h>

using boost::format;

using namespace std;
using namespace CosmoTool;

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
    float shift = 0.5*d->BoxSize;
    rescale_position /= d->time;

    if (d->Pos[0] != 0)
      {
        for (int k = 0; k < 3; k++)
          {
            for (int64_t i = 0; i < d->NumPart; i++)
              {
                d->Pos[k][i] = (d->Pos[k][i]*rescale_position + shift);
              }
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
    *d = *sdf_header;
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
        const char *name_template[3] = { "x", "y", "z" };
        for (int q = 0; q < 3; q++)
          {
            d->Vel[q] = new float[numPartToLoad];
            SETUP_READ(name_template[q], d->Vel[q]);
          }
      }
    if (load_flags & NEED_GADGET_ID)
      {
        d->Id = new long[numPartToLoad];
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
      rescaleParticles(d, 0.001*d->Hubble, one_kpc/one_Gyr);

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
  
  SDFgetintOrDefault(sdfp, "version", &fileversion, 1);
  if (fileversion == 1)
    {
      SDFgetfloatOrDie(sdfp, "Omega0", &hdr->Omega_M);
      SDFgetfloatOrDie(sdfp, "Lambda_prime", &hdr->Omega_Lambda);
      SDFgetfloatOrDie(sdfp, "H0", &hdr->Hubble);
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
      SDFgetfloatOrDie(sdfp, "H0", &hdr->Hubble);
      
    }
  double h0;
  h0 = hdr->Hubble*10.0*(one_kpc/one_Gyr);

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
  MPMY_Init(&argc, &argv);
}
