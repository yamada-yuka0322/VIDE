#include "libqhull/qhull_a.h"
#include "voz.h"

#define DL for (d=0;d<2;d++)
#define BF 1e30

#define MAX(a,b) ( ((a) < (b)) ? (a) : (b) )

int delaunadj_2D (coordT *points, int nvp, int nvpbuf, int nvpall, PARTADJ **adjs);
int vorvol_2D (coordT *deladjs, coordT *points, pointT *intpoints, int numpoints, float *area);
int posread(char *posfile, float ***p, float fact);
int readPosAndIntensity(char *posfile, float ***p, float **intensity, float fact);

int main(int argc, char *argv[]) {
  int exitcode;
  pid_t i, j, np;
  float **r, *weight;
  coordT rtemp[2], *parts;
  coordT deladjs[2*MAXVERVER], points[2*MAXVERVER];
  pointT intpoints[2*MAXVERVER];
  FILE *pos, *out;
  char *posfile, outfile[200], *suffix, *outDir;
  PARTADJ *adjs;
  float *vols;
  realT predict, xmin,xmax,ymin,ymax,zmin,zmax, weightmin, weightmax;
  pid_t *orig;
  
  int isitinbuf;
  char isitinmain, d;
  int numdiv;
  pid_t nvp, nvpall, nvpbuf;
  realT width, width2, totwidth, totwidth2, bf, s, g;
  float border, boxsize;
  realT c[2];
  int b[2];
  double totalvol;

  if (argc != 9) {
    printf("Wrong number of arguments.\n");
    printf("arg1: position file\n");
    printf("arg2: border size\n");
    printf("arg3: boxsize\n");
    printf("arg4: suffix\n");
    printf("arg5: number of divisions\n");
    printf("arg6-7: b[0-1]\n\n");
    printf("arg8: output directory\n");
    exit(0);
  }
  posfile = argv[1];
  if (sscanf(argv[2],"%f",&border) != 1) {
    printf("That's no border size; try again.\n");
    exit(0);
  }
  if (sscanf(argv[3],"%f",&boxsize) != 1) {
    printf("That's no boxsize; try again.\n");
    exit(0);
  }
  suffix = argv[4];
  if (sscanf(argv[5],"%d",&numdiv) != 1) {
    printf("%s is no number of divisions; try again.\n",argv[5]);
    exit(0);
  }
  if (numdiv == 1) {
    printf("Only using one division; should only use for an isolated segment.\n");
  }
  if (numdiv < 1) {
    printf("Cannot have a number of divisions less than 1.  Resetting to 1.\n");
    numdiv = 1;
  }
  if (sscanf(argv[6],"%d",&b[0]) != 1) {
    printf("That's no b index; try again.\n");
    exit(0);
  }
  if (sscanf(argv[7],"%d",&b[1]) != 1) {
    printf("That's no b index; try again.\n");
    exit(0);
  }
  outDir = argv[8];
  
  /* Boxsize should be the range in r, yielding a range 0-1 */
  np = readPosAndIntensity(posfile,&r,&weight,1./boxsize);
  np = posread(posfile,&r,1./boxsize);
  printf("%d particles\n",np);fflush(stdout);

  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF; weightmin = BF; weightmax = -BF;
  for (i=0; i<np;i++) {
    if (r[i][0]<xmin) xmin = r[i][0]; if (r[i][0]>xmax) xmax = r[i][0];
    if (r[i][1]<ymin) ymin = r[i][1]; if (r[i][1]>ymax) ymax = r[i][1];
    if (r[i][2]<zmin) zmin = r[i][2]; if (r[i][2]>zmax) zmax = r[i][2];
    if (weight[i]<weightmin) weightmin = weight[i]; if (weight[i]>weightmax) weightmax = weight[i];
  }
  printf("np: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);
  printf("np: %d, weight: %f,%f;\n",np,weightmin,weightmax); fflush(stdout);

  width = 1./(float)numdiv;
  width2 = 0.5*width;
  if (border > 0.) bf = border;
  else bf = 0.1;
      /* In units of 0-1, the thickness of each subregion's buffer*/
  totwidth = width+2.*bf;
  totwidth2 = width2 + bf;
  
  s = width/(float)NGUARD;
  if ((bf*bf - 2.*s*s) < 0.) {
    printf("bf = %f, s = %f.\n",bf,s);
    printf("Not enough guard points for given border.\nIncrease guards to >= %f\n.",
	   sqrt(2.)*width/bf);
    exit(0);
  }
  g = (bf / 2.)*(1. + sqrt(1 - 2.*s*s/(bf*bf)));
  printf("s = %f, bf = %f, g = %f.\n",s,bf,g);

  fflush(stdout);

  adjs = (PARTADJ *)malloc(np*sizeof(PARTADJ));
  if (adjs == NULL) {
    printf("Unable to allocate adjs\n");
    exit(1);
  }
  
  DL c[d] = ((float)b[d])*width;
  printf("c: %f,%f\n",c[0],c[1]);
  /* Assign temporary array*/
  nvpbuf = 0; /* Number of particles to tesselate, including
		 buffer */
  nvp = 0; /* Without the buffer */
  for (i=0; i<np; i++) {
    isitinbuf = 1;
    isitinmain = 1;
    DL {
      rtemp[d] = (double)r[i][d] - (double)c[d];
      if (rtemp[d] > 0.5) rtemp[d] --;
      if (rtemp[d] < -0.5) rtemp[d] ++;
      isitinbuf = isitinbuf && (fabs(rtemp[d]) < totwidth2);
      isitinmain = isitinmain && (fabs(rtemp[d]) <= width2);
    }
  
    if (isitinbuf) nvpbuf++;
    if (isitinmain) nvp++;
  }
  
  nvpbuf += 4*(NGUARD+1); /* number of guard
					points */

  parts = (coordT *)malloc(2*nvpbuf*sizeof(coordT));
  orig = (pid_t *)malloc(nvpbuf*sizeof(pid_t));

  if (parts == NULL) {
    printf("Unable to allocate parts\n");
    fflush(stdout);
    exit(1);
  }
  if (orig == NULL) {
    printf("Unable to allocate orig\n");
    fflush(stdout);
    exit(1);
  }

  nvp = 0; nvpall = 0; /* nvp = number of particles without buffer */
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np; i++) {
    isitinmain = 1;
    DL {
      rtemp[d] = (realT)r[i][d] - (realT)c[d];
      if (rtemp[d] > 0.5) rtemp[d] --;
      if (rtemp[d] < -0.5) rtemp[d] ++;
      isitinmain = isitinmain && (fabs(rtemp[d]) <= width2);
    }
    if (isitinmain) {
      parts[2*nvp] = rtemp[0];
      parts[2*nvp+1] = rtemp[1];
      orig[nvp] = i;
      nvp++;
      if (rtemp[0] < xmin) xmin = rtemp[0];
      if (rtemp[0] > xmax) xmax = rtemp[0];
      if (rtemp[1] < ymin) ymin = rtemp[1];
      if (rtemp[1] > ymax) ymax = rtemp[1];
    }
  }
  printf("nvp = %d\n",nvp);
  printf("x: %f,%f; y: %f,%f\n",xmin,xmax,ymin,ymax);
  nvpbuf = nvp;
  for (i=0; i<np; i++) {
    isitinbuf = 1;
    isitinmain = 1;

    DL {
      rtemp[d] = (realT)r[i][d] - (realT)c[d];
      if (rtemp[d] > 0.5) rtemp[d] --;
      if (rtemp[d] < -0.5) rtemp[d] ++;
      isitinbuf = isitinbuf && (fabs(rtemp[d])<totwidth2);
      isitinmain = isitinmain && (fabs(rtemp[d]) <= width2);
    }
    if (isitinbuf && !isitinmain) {
      /*printf("%3.3f ",sqrt(rtemp[0]*rtemp[0] + rtemp[1]*rtemp[1] +
	rtemp[2]*rtemp[2]));
	printf("|%2.2f,%2.2f,%2.2f,%f,%f",r[i][0],r[i][1],r[i][2],width2,totwidth2);*/
      parts[2*nvpbuf] = rtemp[0];
      parts[2*nvpbuf+1] = rtemp[1];
      orig[nvpbuf] = i;

      nvpbuf++;
      if (rtemp[0] < xmin) xmin = rtemp[0];
      if (rtemp[0] > xmax) xmax = rtemp[0];
      if (rtemp[1] < ymin) ymin = rtemp[1];
      if (rtemp[1] > ymax) ymax = rtemp[1];
    }
  }
  printf("nvpbuf = %d\n",nvpbuf);
  printf("x: %f,%f; y: %f,%f\n",xmin,xmax,ymin,ymax);
  nvpall = nvpbuf;
  predict = pow(width+2.*bf,3)*(float)np;
  printf("There should be ~ %f points; there are %d\n",predict,nvpbuf);

  for (i=0;i<np;i++) free(r[i]);
  free(r);
  
  /* Add guard points */
  for (i = 0; i < NGUARD + 1; i++) {
    parts[2 * nvpall]     = -width2 + (realT)i * s;
    parts[2 * nvpall + 1] = -width2 - g;
    nvpall++;
    parts[2 * nvpall]     = -width2 + (realT)i * s;
    parts[2 * nvpall + 1] = width2 + g;
    nvpall++;
  }

for (j = 0; j < NGUARD + 1; j++) {
    parts[2 * nvpall]     = -width2 - g;
    parts[2 * nvpall + 1] = -width2 + (realT)j * s;
    nvpall++;

    parts[2 * nvpall]     = width2 + g;
    parts[2 * nvpall + 1] = -width2 + (realT)j * s;
    nvpall++;
}
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=nvpbuf;i<nvpall;i++) {
    if (parts[2*i] < xmin) xmin = parts[2*i];
    if (parts[2*i] > xmax) xmax = parts[2*i];
    if (parts[2*i+1] < ymin) ymin = parts[2*i+1];
    if (parts[2*i+1] > ymax) ymax = parts[2*i+1];
  }
  
  printf("Added guard points to total %d points (should be %d)\n",nvpall,
	 nvpbuf + 4*(NGUARD+1));
  printf("x: %f,%f; y: %f,%f\n",xmin,xmax,ymin,ymax);
  
  /* Do tesselation*/
  printf("File read.  Tessellating ...\n"); fflush(stdout);
  exitcode = delaunadj_2D(parts, nvp, nvpbuf, nvpall, &adjs);
  if (exitcode != 0)
   {
     printf("Error while tesselating. Stopping here."); fflush(stdout);
     exit(1);
   }
  
  /* Now calculate volumes*/
  printf("Now finding volumes ...\n"); fflush(stdout);
  vols = (float *)malloc(nvp*sizeof(float));
  if (vols == NULL) {
    printf("Unable to allocate vols\n");
    exit(1);
  }
  
  for (i=0; i<nvp; i++) { /* Just the original particles
			     Assign adjacency coordinate array*/
    /* Volumes */
    for (j = 0; j < adjs[i].nadj; j++)
      DL {
	deladjs[2*j + d] = parts[2*adjs[i].adj[j]+d] - parts[2*i+d];
	if (deladjs[2*j+d] < -0.5) deladjs[2*j+d]++;
	if (deladjs[2*j+d] > 0.5) deladjs[2*j+d]--;
      }
    
    exitcode = vorvol_2D(deladjs, points, intpoints, adjs[i].nadj, &(vols[i]));
    vols[i] *= (float)np;
    if (i % 1000 == 0)
      printf("%d: %d, volume: %f, weight: %f\n",i,adjs[i].nadj,vols[i],weight[i]);
  }

  /* Get the adjacencies back to their original values */

  for (i=0; i<nvp; i++)
    for (j = 0; j < adjs[i].nadj; j++)
      adjs[i].adj[j] = orig[adjs[i].adj[j]];
  
  totalvol = 0.;
  for (i=0;i<nvp; i++) {
    totalvol += (double)vols[i];
  }
  printf("Average volume = %g\n",totalvol/(float)nvp);
 
  /* Now the output!
     First number of particles */
  sprintf(outfile,"%s/part.%s.%02d.%02d",outDir,suffix,b[0],b[1]);

  printf("Output to %s\n\n",outfile);
  out = fopen(outfile,"w");
  if (out == NULL) {
    printf("Unable to open %s\n",outfile);
    exit(1);
  }
  fwrite(&np,1, sizeof(int),out);
  fwrite(&nvp,1, sizeof(int),out);
  printf("nvp = %d\n",nvp);

  /* Tell us where the original particles were */
  fwrite(orig,sizeof(pid_t),nvp,out);
  /* Volumes*/
  fwrite(vols,sizeof(float),nvp,out);
  /* Weights*/
  fwrite(weight,sizeof(float),nvp,out);
  /* Adjacencies */
  for (i=0;i<nvp;i++) {
    fwrite(&(adjs[i].nadj),1,sizeof(pid_t),out);
    if (adjs[i].nadj > 0)
      fwrite(adjs[i].adj,adjs[i].nadj,sizeof(pid_t),out);
    else printf("0");
  }
  fclose(out);
  
  return(0);
}
