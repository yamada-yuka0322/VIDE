/*+
This is CosmoTool (./src/load_data.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "load_data.hpp"

using namespace CosmoTool;

//#define LARGE_CONTROL
//#define LITTLE_ENDIAN

#define NEW(t,n)  ((t *)malloc(sizeof(t)*n))
#define SKIP(f) fread(&dummy,sizeof(dummy),1,f);
#define WRITE_DUM(f) fwrite(&dummy, sizeof(dummy),1,f);

static int dummy;

void CosmoTool::writeGadget(GadgetData *data, const char *fname)
{
  FILE *f;
  int k, n, p;

  f = fopen(fname, "w");
  if (f == NULL) {
    fprintf(stderr, "Cannot write gadget to file %s\n", fname);
    return;
  }

  dummy = 256;
  WRITE_DUM(f);
  fwrite(&data->header, sizeof(data->header), 1, f);
  WRITE_DUM(f);

  dummy = sizeof(float)*3*data->NumPart;
  WRITE_DUM(f);
  for(k=0,p=0;k<5;k++) {
    for(n=0;n<data->header.npart[k];n++) {
      fwrite(&data->particles[p].Pos[0], sizeof(float), 3, f);
      p++;
    }
  }
  WRITE_DUM(f);

  dummy = sizeof(float)*3*data->NumPart;  
  WRITE_DUM(f);
  for(k=0,p=0;k<6;k++) {
    for(n=0;n<data->header.npart[k];n++) {
      fwrite(&data->particles[p].Vel[0], sizeof(float), 3, f);
      p++;
    }
  }
  WRITE_DUM(f);

  dummy = sizeof(int)*data->NumPart;
  WRITE_DUM(f);
  for(k=0,p=0;k<6;k++)
    {
      for(n=0;n<data->header.npart[k];n++)
	{
	  fwrite(&data->particles[p].Id, sizeof(int), 1, f);
	  p++;
	}
    }
  WRITE_DUM(f);

  if(data->ntot_withmasses>0) {
    dummy = sizeof(float)*data->NumPart;
    WRITE_DUM(f);
  }
  for(k=0, p=0; k<6; k++)
    {
      for(n=0;n<data->header.npart[k];n++)
	{
	  if(data->header.mass[k]==0)
	    fwrite(&data->particles[p].Mass, sizeof(float), 1, f);
	  p++;
	}
    }
  if(data->ntot_withmasses>0)
    WRITE_DUM(f);      

  if(data->header.npart[0]>0) {
    dummy = data->header.npart[0]*sizeof(float);
    WRITE_DUM(f);
    for(n=0, p=0; n<data->header.npart[0];p++,n++) {
      fwrite(&data->particles[p].U, sizeof(float), 1, f);
    }
    WRITE_DUM(f);
    
    WRITE_DUM(f);
    for(n=0, p=0; n<data->header.npart[0];p++,n++) {
      fwrite(&data->particles[p].Rho, sizeof(float), 1, f);
    }
    WRITE_DUM(f);
    
    if(data->header.flag_cooling) {
      WRITE_DUM(f);
      for(n=0, p=0; n<data->header.npart[0];p++,n++) {
	fwrite(&data->particles[p].Ne, sizeof(float), 1, f);
      }
      WRITE_DUM(f);
    }
  }
  
  fclose(f);
}

GadgetData *CosmoTool::loadGadget(const char *fname)
{
  FILE *f;
  GadgetData *data;
  int p, k, n;

  f = fopen(fname, "r");
  if (f == NULL)
    return NULL;
  
  data = NEW(GadgetData, 1);
  SKIP(f);
  fread(&data->header, sizeof(data->header), 1, f);
  SKIP(f);
  
  for(k=0, data->ntot_withmasses=0; k<5; k++) {
    if(data->header.mass[k]==0)
      data->ntot_withmasses+= data->header.npart[k];
  }

  for(k=0, data->NumPart=0; k<5; k++)
    data->NumPart+= data->header.npart[k];
  
  data->particles = NEW(ParticleState, data->NumPart);

  SKIP(f);
  for(k=0,p=0;k<5;k++) {
    for(n=0;n<data->header.npart[k];n++) {
      fread(&data->particles[p].Pos[0], sizeof(float), 3, f);
      p++;
    }
  }
  SKIP(f);
  
  SKIP(f);
  for(k=0,p=0;k<6;k++) {
    for(n=0;n<data->header.npart[k];n++) {
      fread(&data->particles[p].Vel[0], sizeof(float), 3, f);
      p++;
    }
  }
  SKIP(f);
  

  SKIP(f);
  for(k=0,p=0;k<6;k++)
    {
      for(n=0;n<data->header.npart[k];n++)
	{
	  fread(&data->particles[p].Id, sizeof(int), 1, f);
	  p++;
	}
    }
  SKIP(f);

  if(data->ntot_withmasses>0)
    SKIP(f);
  for(k=0, p=0; k<6; k++)
    {
      for(n=0;n<data->header.npart[k];n++)
	{
	  data->particles[p].Type=k;

	  if(data->header.mass[k]==0)
	    fread(&data->particles[p].Mass, sizeof(float), 1, f);
	  else
	    data->particles[p].Mass= data->header.mass[k];
	  p++;
	}
    }
  if(data->ntot_withmasses>0)
    SKIP(f);      

  if(data->header.npart[0]>0)
    {
      SKIP(f);
      for(n=0, p=0; n<data->header.npart[0];p++,n++) {
	  fread(&data->particles[p].U, sizeof(float), 1, f);
	}
      SKIP(f);

      SKIP(f);
      for(n=0, p=0; n<data->header.npart[0];p++,n++) {
	fread(&data->particles[p].Rho, sizeof(float), 1, f);
      }
      SKIP(f);

      if(data->header.flag_cooling)
	{
	  SKIP(f);
	  for(n=0, p=0; n<data->header.npart[0];p++,n++)
	    {
	      fread(&data->particles[p].Ne, sizeof(float), 1, f);
	    }
	  SKIP(f);
	}
      else
	for(n=0, p=0; n<data->header.npart[0];p++,n++)
	  {
	    data->particles[p].Ne= 1.0;
	  }
    }


  fclose(f);

  return data;
}

void CosmoTool::freeGadget(GadgetData *data)
{
  free(data->particles);
  free(data);
}

void CosmoTool::writePersoSet(ParticleSet *set, const char *fname)
{
  FILE *f;
  int i;

  f = fopen(fname, "w");
  if (f == NULL) {
    perror("writePersoSet");
    return;
  }

  fwrite(&set->header, sizeof(set->header), 1, f);
  fwrite(set->Npart, sizeof(set->Npart[0]), set->header.Ntypes, f);
  
  for (i=0;i<set->header.Ntypes;i++)
    fwrite(set->particles[i], sizeof(ParticleState), set->Npart[i], f);

  fclose(f);
}

ParticleSet *CosmoTool::loadPersoSet(const char *fname)
{
  ParticleSet *set;
  FILE *f;
  int i;

  f = fopen(fname, "r");
  if (f == NULL) {
    perror("loadPersoSet");
    return NULL;
  }

  set = NEW(ParticleSet, 1);
  fread(&set->header, sizeof(set->header), 1, f);

  set->Npart = NEW(int, set->header.Ntypes);
  fread(set->Npart, sizeof(set->Npart[0]), set->header.Ntypes, f);;
  
  set->particles = NEW(ParticleState *, set->header.Ntypes);
  for (i=0;i<set->header.Ntypes;i++) {
    set->particles[i] = NEW(ParticleState, set->Npart[i]);
    fread(set->particles[i], sizeof(ParticleState), set->Npart[i], f);
  }

  fclose(f);

  return set;
}

void CosmoTool::freePersoSet(ParticleSet *set)
{
  int i;

  for (i=0;i<set->header.Ntypes;i++) {
    free(set->particles[i]);
  }
  if (set->Npart != NULL) {
    free(set->particles);
    free(set->Npart);
  }
}

#ifdef WANT_MAIN
int main(int argc, char **argv) {
  GadgetData *data;
  FILE *plot;
  int i;
  double bl;
  int N;
  double rms;

  if (argc < 3) {
      fprintf(stderr, "Usage: %s [GADGET DATA FILE] [BOXSIZE] [N PARTIC]\n", argv[0]);
      return -1;
  }

  plot = fopen("plot", "w");
  
  bl = atof(argv[2]);
  data = loadGadget(argv[1]);

  printf("Redshift: %lg\n", data->header.redshift);
  rms = 0;
  N = atoi(argv[3]);
  for (i=0;i<data->NumPart;i++) {
    if (i == data->header.npart[0])
	    fprintf(plot,"\n\n");

    fprintf(plot, "%f %f %f\n", data->particles[i].Pos[0], data->particles[i].Pos[1], data->particles[i].Pos[2]);


    /* Compute the RMS */
    {
      /* First find the nearest grid node. */
      int k;
      int x;
      double dx;

      for (k=0;k<3;k++) {
	x = data->particles[i].Pos[k] / bl * N;
	dx = data->particles[i].Pos[k]-x*bl/N;
	rms += dx*dx;
      }
    }
  }

  printf("delta rms = %e\n", sqrt(rms/data->NumPart));
  freeGadget(data);

  fclose(plot);

  return 0;
}
#endif

#define LEN0 200.0

GadgetData *CosmoTool::loadSimulationData(const char *fname)
{
  GadgetData *gd = NEW(GadgetData, 1);
  FILE *f;
  int lineNo;
  char line[1024];
  int i;
  int j;
  
  gd->header.BoxSize = LEN0;

  f = fopen(fname, "r");
  lineNo = 0;
  while (!feof(f))
    {
      fgets(line, sizeof(line), f);
      lineNo++;
    }
  lineNo--;
  rewind(f);

  gd->NumPart = lineNo;
  gd->particles = NEW(ParticleState, lineNo);

  i = 0;
  while (!feof(f))
    {
      fgets(line, sizeof(line), f);
      int r = sscanf(line, "%*d %*d %f %f %f %f %f %f %f %f %f %*f %*f %*f %f %f %f",
		     &gd->particles[i].Pos[0], &gd->particles[i].Pos[1], &gd->particles[i].Pos[2],
		     &gd->particles[i].Init[0], &gd->particles[i].Init[1], &gd->particles[i].Init[2],
		     &gd->particles[i].Vel[0], &gd->particles[i].Vel[1], &gd->particles[i].Vel[2],
		     &gd->particles[i].VelInit[0], &gd->particles[i].VelInit[1], &gd->particles[i].VelInit[2]
		     );
      if (r != 12)
	{
	  printf("line %d: '%s'\n", i, line);
	  printf("returned r=%d\n", r);
	  abort();
	}
      assert(r == 12);
      for (j = 0; j < 3; j++)
	{
	  gd->particles[i].Vel[j] *= 100.0 * LEN0 / (0.9641010);
	  gd->particles[i].VelInit[j] *= 100.0 * 1/71. * LEN0 / (0.9641010);
	  gd->particles[i].Pos[j] *= LEN0;
	  gd->particles[i].Init[j] *= LEN0;
	}

      gd->particles[i].Type = 0;
      gd->particles[i].Mass = 1.0;
      gd->particles[i].Id = i;

      i++;
    }
  fclose(f);

  return gd;
}

#ifndef LITTLE_ENDIAN 
#define read_buf(b, n) \
{ \
  int k; \
  control_size -= n; \
  for (k = (n-1); k >= 0; k--) \
    fread(&b[k], 1, 1, infile); \
}
#else
#define read_buf(b, n) \
{ \
  int k; \
  control_size -= n; \
  for (k = 0; k < n; k++) \
    fread(&b[k], 1, 1, infile); \
}
#endif

#define read_int(i) \
{ \
  char *o = (char*)&(i); \
  read_buf(o, 4); \
}

#define read_real(f) \
{ \
  char *o = (char*)&(f); \
  read_buf(o, 4); \
}

#define read_characters(c, n) { \
  int k; \
  control_size -= n; \
  fread(c, 1, n, outfile); \
}

#define push_dummy_control(id) \
{ int control_size = 0;

#define pop_dummy_control() }

#if defined(LARGE_CONTROL) && defined(LITTLE_ENDIAN)
#define push_control(id) \
{ \
  int control_size = 0; \
  int control_size2 = 0; \
  char *intbuf = (char*)&control_size; \
  fread(&control_size, 8, 1, infile);

#define pop_control(id) \
  fread(&control_size2, 8, 1, infile); \
  assert(control_size == 0); \
}
#elif !defined(LARGE_CONTROL) && defined(LITTLE_ENDIAN)
#define push_control(id) \
{ \
  int control_size = 0; \
  int control_size2 = 0; \
  char *intbuf = (char*)&control_size; \
  fread(&control_size, 4, 1, infile);

#define pop_control(id) \
  fread(&control_size2, 4, 1, infile); \
  assert(control_size == 0); \
}

#elif defined(LARGE_CONTROL) && !defined(LITTLE_ENDIAN)
#define push_control(id) \
{ \
  int control_size = 0; \
  int control_size2 = 0; \
  char *intbuf = (char*)&control_size; \
  fread(&control_size, 8, 1, infile);

#define pop_control(id) \
  fread(&control_size2, 8, 1, infile); \
  assert(control_size == 0); \
}

#elif !defined(LARGE_CONTROL) && !defined(LITTLE_ENDIAN)
#define push_control(id) \
{ \
  int control_size = 0; \
  int control_size2 = 0; \
  char *intbuf = (char*)&control_size; \
  fread(&control_size, 4, 1, infile);

#define pop_control(id) \
  fread(&control_size2, 4, 1, infile); \
  assert(control_size == 0); \
}

#endif

GadgetData *CosmoTool::loadHydra(const char *fname)
{
  GadgetData *gd = NEW(GadgetData, 1);
  FILE *f;
  int version0, version1, version2;
  int irun, nobj, ngas, ndark, intl, nlmx, perr;
  float dtnorm, sft0, sftmin, sftmax;
  int pad3;
  float h100, box100, zmet0;
  int lcool;
  float rmnorm0;
  int pad4, pad5;
  float tstart, omega0, xlambda0, h0t0, rcen, rmax2;
  float rmbary;
  int j;
  float atime;

  f = fopen(fname, "r");
#define infile f
  push_control(0);
  read_int(version0);
  read_int(version1);
  read_int(version2);
  pop_control(0);

  if (version0 != 4)
    {
      fclose(f);
      return NULL;
    }
  push_control(1);
  
  for (j = 0; j < 200; j++)
    {
      int mydummy;
      read_int(mydummy);
    }
  for (j = 0; j < 5; j++)
  {
	  float mydummy;
	  read_real(mydummy);
  }
  read_real(atime);
  gd->header.time = atime;
  gd->header.redshift = 1/atime - 1;

  for (j = 6; j < 100; j++)
    {
      int mydummy;
      read_int(mydummy);
    }
  read_int(irun);
  read_int(nobj);
  read_int(ngas);
  read_int(ndark);
  read_int(intl);
  read_int(nlmx);
  read_int(perr);
  read_real(dtnorm);
  read_real(sft0);
  read_real(sftmin);
  read_real(sftmax);
  read_int(pad3);
  read_real(h100);
  read_real(box100);
  read_real(zmet0);
  read_int(lcool);
  read_real(rmbary);
  read_real(rmnorm0);
  read_int(pad4);
  read_int(pad5);
  read_real(tstart);
  read_real(omega0);
  read_real(xlambda0);
  read_real(h0t0);
  read_real(rcen);
  read_real(rmax2);
  for (j = 0; j < 74; j++)
    {
      int mydummy;
      read_int(mydummy);
    }
  pop_control(1);

  gd->header.npart[1] = ndark;
  gd->header.npart[0] = ngas;
  gd->header.num_files = 1;
  gd->header.flag_cooling = lcool;
  gd->header.BoxSize = box100 * 1000;
  gd->header.HubbleParam = h100;
  gd->header.Omega0 = omega0;
  gd->header.OmegaLambda = xlambda0;

  push_control(2);
  for (j = 0; j < nobj; j++)
    {
      int mydummy;
      read_int(mydummy);
    }
  pop_control(2);

  gd->NumPart = nobj;
  gd->ntot_withmasses = nobj;
  gd->particles = NEW(ParticleState, nobj);

  push_control(3);
  for (j = 0; j < nobj; j++)
    {
      float rm;
      gd->particles[j].Id = j;
      read_real(gd->particles[j].Mass);
    }
  pop_control(3);

  push_control(4);
  for (j = 0; j < nobj; j++)
    {
      int k;
      for (k = 0; k < 3; k++)
	{
	  read_real(gd->particles[j].Pos[k]);
	  gd->particles[j].Pos[k] *= gd->header.BoxSize;
	}
    }
  pop_control(4);

  push_control(5);
  for (j = 0; j < nobj; j++)
    {
      int k;
      for (k = 0; k < 3; k++)
	{
	  read_real(gd->particles[j].Vel[k]);
	  gd->particles[j].Vel[k] *= 100.0 * box100  / h0t0 * atime;
	}
    }
  pop_control(5);

  fclose(f);
#undef infile

  return gd;
}
