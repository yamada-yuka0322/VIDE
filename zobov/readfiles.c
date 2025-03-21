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


int readPosAndIntensity(const char* posfile, float*** pos, float** intensity, int* numElements) {
    FILE* file = fopen(posfile, "rb");
    if (file == NULL) {
        fprintf(stderr, "Unable to open file %s\n", posfile);
        return -1;
    }

    // データの個数を読み取る
    int np;
    fread(&np, sizeof(int), 1, file);
    *numElements = np; // データの個数をセット

    // pos データのメモリを動的に確保
    *pos = (float**)malloc(np * sizeof(float*));
    if (*pos == NULL) {
        fprintf(stderr, "Memory allocation failed for pos\n");
        fclose(file);
        return -1;
    }

    float* temp = (float*)malloc(np * 3 * sizeof(float)); // 位置データを一時的に格納
    float* _temp = (float*)malloc(np * sizeof(float)); // 位置データを一時的に格納
    if (temp == NULL) {
        fprintf(stderr, "Memory allocation failed for temp\n");
        fclose(file);
        free(*pos);
        return -1;
    }

    // pos データの読み込み
    fread(temp, sizeof(float), np * 3, file);
    for (int i = 0; i < np; i++) {
        (*pos)[i] = (float*)malloc(3 * sizeof(float));
        if ((*pos)[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for pos[%d]\n", i);
            fclose(file);
            free(temp);
            for (int j = 0; j < i; j++) {
                free((*pos)[j]);
            }
            free(*pos);
            return -1;
        }
        (*pos)[i][0] = temp[i];
        (*pos)[i][1] = temp[np + i];
        (*pos)[i][2] = temp[2 * np + i];
    }

    free(temp); // temp メモリ解放

    // intensity データのメモリを動的に確保
    *intensity = (float*)malloc(np * sizeof(float));
    fread(_temp, sizeof(float), np , file);
    if (*intensity == NULL) {
        fprintf(stderr, "Memory allocation failed for intensity\n");
        fclose(file);
        for (int i = 0; i < np; i++) {
          intensity[i] = _temp[i];
          free((*pos)[i]);
        }
        free(*pos);
        return -1;
    }

    free(_temp);

    fclose(file); // ファイルをクローズ

    return 0; // 成功
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
