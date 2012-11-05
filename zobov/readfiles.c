#include <stdio.h>
#include <stdlib.h>

#define DL for (d=0;d<3;d++) /* Dimension loop */
#define BF 1e30

/* Positions */
/* Returns number of particles read */
int posread(char *posfile, float ***p, float fact) {

  FILE *pos;
  int np,dum,d,i;
  float xmin,xmax,ymin,ymax,zmin,zmax;
  float *ptemp;

  pos = fopen(posfile, "r");
  if (pos == NULL) {
    printf("Unable to open position file %s\n\n",posfile);
    exit(0);
  }
  /* Fortran77 4-byte headers and footers */
  /* Delete "dum" statements if you don't need them */

  /* Read number of particles */
  fread(&dum,1,4,pos); fread(&np,1, sizeof(int),pos); fread(&dum,1,4,pos);

  /* Allocate the arrays */
  (*p) = (float **)malloc(np*sizeof(float *));
  ptemp = (float *)malloc(np*sizeof(float));

  printf("np = %d\n",np);

  /* Fill the arrays */
  fread(&dum,1,4,pos); 
  fread(ptemp,np,4,pos);
  for (i=0; i<np; i++) {
    (*p)[i] = (float *)malloc(3*sizeof(float));
    if ((*p)[i] == NULL) {
      printf("Unable to allocate particle array in readfiles!\n");
      fflush(stdout);
      exit(0);
    }
    (*p)[i][0] = ptemp[i];
  }
  fread(&dum,1,4,pos); 
  fread(&dum,1,4,pos); 
  fread(ptemp,np,4,pos);
  for (i=0; i<np; i++) (*p)[i][1] = ptemp[i];
  fread(&dum,1,4,pos); 
  fread(&dum,1,4,pos); 
  fread(ptemp,np,4,pos);
  for (i=0; i<np; i++) (*p)[i][2] = ptemp[i];
  fread(&dum,1,4,pos); 

  fclose(pos);
  free(ptemp);

  /* Get into physical units (Mpc/h) */
  

  printf("%f\n",fact);fflush(stdout);
  for (i=0; i<np; i++) DL (*p)[i][d] *= fact;

  /* Test range -- can comment out */
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np;i++) {
    if ((*p)[i][0]<xmin) xmin = (*p)[i][0]; if ((*p)[i][0]>xmax) xmax = (*p)[i][0];
    if ((*p)[i][1]<ymin) ymin = (*p)[i][1]; if ((*p)[i][1]>ymax) ymax = (*p)[i][1];
    if ((*p)[i][2]<zmin) zmin = (*p)[i][2]; if ((*p)[i][2]>zmax) zmax = (*p)[i][2];
  }
  printf("np: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);

  return(np);
}

/* Velocities */
/* Returns number of particles read */
/* Used in voboz, but not zobov, which doesn't use velocities */
int velread(char *velfile, float ***v, float fact) {

  FILE *vel;
  int np,dum,d,i;
  float xmin,xmax,ymin,ymax,zmin,zmax;

  vel = fopen(velfile, "r");
  if (vel == NULL) {
    printf("Unable to open velocity file %s\n\n",velfile);
    exit(0);
  }
  /* Fortran77 4-byte headers and footers */
  /* Delete "dum" statements if you don't need them */

  /* Read number of particles */
  fread(&dum,1,4,vel); fread(&np,1, sizeof(int),vel); fread(&dum,1,4,vel);

  /* Allocate the arrays */
  (*v) = (float **)malloc(3*sizeof(float*));
  for (i=0;i<3;i++) (*v)[i] = (float *)malloc(np*sizeof(float));

  /* Fill the arrays */
  fread(&dum,1,4,vel); fread((*v)[0],np,4,vel); fread(&dum,1,4,vel);
  fread(&dum,1,4,vel); fread((*v)[1],np,4,vel); fread(&dum,1,4,vel);
  fread(&dum,1,4,vel); fread((*v)[2],np,4,vel); fread(&dum,1,4,vel);

  fclose(vel);

  /* Convert from code units into physical units (km/sec) */
  
  for (i=0; i<np; i++) DL (*v)[d][i] *= fact;

  /* Test range -- can comment out */
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np;i++) {
    if ((*v)[0][i] < xmin) xmin = (*v)[0][i]; if ((*v)[0][i] > xmax) xmax = (*v)[0][i];
    if ((*v)[1][i] < ymin) ymin = (*v)[1][i]; if ((*v)[1][i] > ymax) ymax = (*v)[1][i];
    if ((*v)[2][i] < zmin) zmin = (*v)[2][i]; if ((*v)[2][i] > zmax) zmax = (*v)[2][i];
  }
  printf("vx: %f,%f; vy: %f,%f; vz: %f,%f\n",xmin,xmax, ymin,ymax, zmin,zmax);
  
  return(np);
}
