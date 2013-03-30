#include "libqhull/qhull_a.h"
#include "voz.h"

#define NUMCPU 2
#define DL for (d=0;d<3;d++)
#define BF 1e30

#define QHULL_MAX_PARTICLES (1L<<24-1)

int posread(char *posfile, float ***p, float fact);

int main(int argc, char *argv[]) {
  int i, np;
  float **rfloat, rtemp[3];
  FILE *pos, *scr;
  char *posfile, scrfile[200], systemstr[90], *suffix, *outDir, *vobozPath;
  float xmin,xmax,ymin,ymax,zmin,zmax;
  
  int isitinbuf;
  char isitinmain, d;
  int numdiv;
  int p; 
  int nvp, nvpall, nvpbuf, nvpmin, nvpmax, nvpbufmin, nvpbufmax; /* yes, the insurance */
  float width, width2, totwidth, totwidth2, bf, s, g;
  float border, boxsize;
  float c[3];
  int numGuards;
  int b[3];
  int numThreads;
  int mockIndex;

  if (argc != 10) {
    printf("Wrong number of arguments.\n");
    printf("arg1: position file\n");
    printf("arg2: buffer size (default 0.1)\n");
    printf("arg3: box size\n");
    printf("arg4: number of divisions (default 2)\n");
    printf("arg5: suffix describing this run\n");
    printf("arg6: number of parallel threads\n");
    printf("arg7: location of voboz executables\n");
    printf("arg8: output directory\n");
    printf("arg9: index of mock galaxies\n\n");
    exit(0);
  }
  posfile = argv[1];
  suffix = argv[2];
  if (sscanf(suffix,"%f",&border) != 1) {
    printf("That's no border size; try again.\n");
    exit(0);
  }
  suffix = argv[3];
  if (sscanf(suffix,"%f",&boxsize) != 1) {
    printf("That's no boxsize; try again.\n");
    exit(0);
  }
  suffix = argv[4];
  if (sscanf(suffix,"%d",&numdiv) != 1) {
    printf("That's no number of divisions; try again.\n");
    exit(0);
  }
  if (numdiv < 2) {
    printf("Cannot have a number of divisions less than 2.  Resetting to 2:\n");
    numdiv = 2;
  }

  suffix = argv[5];
  vobozPath = argv[6];
  if (sscanf(vobozPath,"%d",&numThreads) != 1) {
    printf("That's no number of threads; try again.\n");
    exit(0);
  }
  vobozPath = argv[7];
  outDir = argv[8];
  if (sscanf(argv[9],"%d",&mockIndex) != 1) {
    printf("That's no mock galaxy index; try again.\n");
    exit(0);
  }


  /* Read the position file */
  np = posread(posfile,&rfloat,1./boxsize);
  /* Boxsize should be the range in r, yielding a range 0-1 */

  width = 1./(float)numdiv;
  width2 = 0.5*width;
  if (border > 0.) bf = border;
  else bf = 0.1;

  /* In units of 0-1, the thickness of each subregion's buffer*/
  totwidth = width+2.*bf;
  totwidth2 = width2 + bf;
  
  s = width/(float)NGUARD;
  if ((bf*bf - 2.*s*s) < 0.) {
    printf("Not enough guard points for given border.\nIncrease guards to >= %f\n.",
	   totwidth/sqrt(0.5*bf*bf));
    printf("bf = %f\n",bf);
    exit(0);
  }
  g = (bf / 2.)*(1. + sqrt(1 - 2.*s*s/(bf*bf)));
  printf("s = %f, bf = %f, g = %f.\n",s,bf,g);
  
  nvpmax = 0; nvpbufmax = 0; nvpmin = np; nvpbufmin = np;
  
  for (b[0] = 0; b[0] < numdiv; b[0]++) {
   c[0] = ((float)b[0]+0.5)*width;
   for (b[1] = 0; b[1] < numdiv; b[1]++) {
    c[1] = ((float)b[1]+0.5)*width;
    for (b[2] = 0; b[2] < numdiv; b[2]++) {
      c[2] = ((float)b[2]+0.5)*width;

      nvp = 0; /* Number of particles excluding buffer */
      nvpbuf = 0; /* Number of particles to tesselate, including
		     buffer */
      xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
      for (i=0; i<np; i++) {
	isitinbuf = 1; isitinmain = 1;
	for (d=0; d<3; d++) {
	  rtemp[d] = rfloat[i][d] - c[d];
	  if (rtemp[d] > 0.5) rtemp[d] --;
	  if (rtemp[d] < -0.5) rtemp[d] ++;
	  isitinbuf = isitinbuf && (fabs(rtemp[d]) < totwidth2);
	  isitinmain = isitinmain && (fabs(rtemp[d]) <= width2);
	}
	if (isitinbuf) {
	  nvpbuf++;
	}
	if (isitinmain) {
	  nvp++;
	  if (rtemp[0] < xmin) xmin = rtemp[0];
	  if (rtemp[0] > xmax) xmax = rtemp[0];
	  if (rtemp[1] < ymin) ymin = rtemp[1];
	  if (rtemp[1] > ymax) ymax = rtemp[1];
	  if (rtemp[2] < zmin) zmin = rtemp[2];
	  if (rtemp[2] > zmax) zmax = rtemp[2];
	}
      }
      if (nvp > nvpmax) nvpmax = nvp;
      if (nvpbuf > nvpbufmax) nvpbufmax = nvpbuf;
      if (nvp < nvpmin) nvpmin = nvp;
      if (nvpbuf < nvpbufmin) nvpbufmin = nvpbuf;

      printf("b=(%d,%d,%d), c=(%f,%f,%f), nvp=%d, nvpbuf=%d\n",
	     b[0],b[1],b[2],c[0],c[1],c[2],nvp,nvpbuf);
    }
   }
  }
  printf("Nvp range: %d,%d\n",nvpmin,nvpmax);
  printf("Nvpbuf range: %d,%d\n",nvpbufmin,nvpbufmax);

  numGuards = 6*(NGUARD+1)*(NGUARD+1);
  if (nvpbufmax+numGuards >= QHULL_MAX_PARTICLES)
    {
      printf("Too many particles to tesselate per division (Qhull will overflow). Increase divisions or reduce buffer size.");
      fflush(stdout);
      exit(1);
    }

  /* Output script file */
  sprintf(scrfile,"scr%s",suffix);
  printf("Writing script file to %s.\n",scrfile);fflush(stdout);
  scr = fopen(scrfile,"w");
  if (scr == NULL) {
    printf("Problem opening script file.\n");
    fflush(stdout);
    exit(1);
  }
  fprintf(scr,"#!/bin/bash -f\n");
  p = 0;
  for (b[0]=0;b[0]<numdiv; b[0]++) {
   for (b[1] = 0; b[1] < numdiv; b[1]++) {
    for (b[2] = 0; b[2] < numdiv; b[2]++) {
// TEST
      fprintf(scr,"%s/voz1b1 %s %f %f,%f,%f %s %d %d %d %d %s&\n",
        vobozPath,
	      posfile,border,boxsize,boxsize,boxsize,suffix,numdiv,b[0],b[1],b[2], 
        outDir);
// END TEST
      p++;
      if ((p == numThreads)) { fprintf(scr, "wait\n"); p = 0; }
    }
   }
  }
  fprintf(scr,"wait\n");
  fprintf(scr,"%s/voztie %d %s %s %d\n", vobozPath, numdiv,suffix, outDir, mockIndex);
  fclose(scr);

  sprintf(systemstr,"chmod u+x %s",scrfile);
  system(systemstr);

  return(0);
}
