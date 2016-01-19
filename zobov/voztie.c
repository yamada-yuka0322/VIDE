#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "voz.h"

int compare_int(const void *a, const void *b)
{
  const int *ia = (const int *)a;
  const int *ib = (const int *)b;
  if ((*ia) < (*ib))
    return -1;
  else if ((*ia) > (*ib))
    return 1;
  else
    return 0;
}

int main(int argc, char *argv[]) {

  FILE *part, *adj, *vol;
  char partfile[200], *suffix, adjfile[200], volfile[200], *outDir;
  float *vols, volstemp;
  int *cnt_adj;
  
  PARTADJ *adjs;

  int numdiv,np,np2,na;

  pid_t i,j,k,p,nout;
  pid_t nvp,npnotdone,nvpmax, nvpsum, *orig;
  double avgnadj, avgvol;

  int numRemoved = 0;
  
  // PMS
  int mockIndex;
  // END PMS

  if (argc != 5) {
    printf("Wrong number of arguments.\n");
    printf("arg1: number of divisions (default 2)\n");
    printf("arg2: suffix describing this run\n");
    printf("arg3: file directory\n");
    printf("arg4: Beginning index of mock particles\n\n");
    exit(0);
  }
  if (sscanf(argv[1],"%d",&numdiv) != 1) {
    printf("That's no number of divisions; try again.\n");
    exit(0);
  }
  if (numdiv < 1) {
    printf("Cannot have a number of divisions less than 1.  Resetting to 1:\n");
    numdiv = 1;
  }

  suffix = argv[2];
  outDir = argv[3];
  if (sscanf(argv[4],"%d",&mockIndex) != 1) {
    printf("That's no mock galaxy index; try again.\n");
    exit(0);
  }
  
  np = -1; nvpmax = -1; nvpsum = 0;

  for (i = 0; i < numdiv; i++) {
   for (j = 0; j < numdiv; j++) {
    for (k = 0; k < numdiv; k++) {
      sprintf(partfile,"%s/part.%s.%02d.%02d.%02d",outDir,suffix,i,j,k);
      part = fopen(partfile,"r");
      if (part == NULL) {
	printf("Unable to open file %s.\n\n",partfile);
	exit(0);
      }
      fread(&np2,1,sizeof(int),part);
      fread(&nvp,1,sizeof(int),part);
      nvpsum += nvp;
      if (np == -1)
	np = np2;
      else 
	if (np2 != np) {
	  printf("Incompatible total particle numbers: %d,%d\n\n",np,np2);
	  exit(0);
	}
      if (nvp > nvpmax) nvpmax = nvp;
      fclose(part);
    }
   }
  }

  printf("We have %d particles to tie together.\n",np); fflush(stdout);
  printf("The maximum number of particles in a file is %d.\n",nvpmax);

  // PMS
  if (mockIndex == -1) mockIndex = np;
  // END PMS

  cnt_adj = (int *)malloc(np*sizeof(int));
  if (cnt_adj == NULL)
    printf("Could not allocate cnt_adj.\n");

  adjs = (PARTADJ *)malloc(np*sizeof(PARTADJ));
  if (adjs == NULL) printf("Couldn't allocate adjs.\n");
  vols = (float *)malloc(np*sizeof(float));
  if (vols == NULL) printf("Couldn't allocate vols.\n");
  orig = (pid_t *)malloc(nvpmax*sizeof(pid_t));
  if (orig == NULL) printf("Couldn't allocate orig.\n");
  if ((cnt_adj==NULL) || (vols == NULL) || (orig == NULL) || (adjs == NULL)) {
    printf("Not enough memory to allocate. Exiting.\n");
    exit(0);
  }
  for (p=0;p<np;p++)
    {
      vols[p] = -1.;
      adjs[p].nadj = 0;
      adjs[p].adj = 0;
    }

  nvpsum = 0;
  for (i = 0; i < numdiv; i++) {
   for (j = 0; j < numdiv; j++) {
    for (k = 0; k < numdiv; k++) {
      sprintf(partfile,"%s/part.%s.%02d.%02d.%02d",outDir,suffix,i,j,k);
      part = fopen(partfile,"r");
      if (part == NULL) {
	printf("Unable to open file %s.\n\n",partfile);
	exit(0);
      }
      fread(&np2,1,sizeof(pid_t),part);
      fread(&nvp,1,sizeof(pid_t),part);
      /*printf("nvp = %d\n",nvp);fflush(stdout);*/

      nvpsum += nvp;

      fread(orig,nvp,sizeof(pid_t),part);
      for (p=0;p<nvp;p++) {
	fread(&volstemp,1,sizeof(float),part);
	if (vols[orig[p]] > -1.)
	  if (fabs(vols[orig[p]]-volstemp)/volstemp > 1.5e-3 && orig[p] < mockIndex) {
	    printf("Inconsistent volumes for p. %d: (%10g,%10g)!\n",
		   orig[p],vols[orig[p]],volstemp);
	    //exit(0);
	  }
	vols[orig[p]] = volstemp;
      }
      
      for (p=0;p<nvp;p++) {
	fread(&na,1,sizeof(int),part);
        pid_t pid = orig[p];
	if (na > 0) {
          pid_t *previous_adj = adjs[pid].adj;
          assert(adjs[pid].nadj == 0 || adjs[pid].nadj == na);
	  adjs[pid].nadj = na;
          adjs[pid].adj = (pid_t *)malloc(na*sizeof(pid_t));
	  if (adjs[pid].adj == 0) {
	    printf("Couldn't allocate adjs[orig[%d]].adj.\n",p);
	    exit(0);
	  }
	  fread(adjs[pid].adj,na,sizeof(pid_t),part);
          if (previous_adj != 0)
           {
             qsort(previous_adj, na, sizeof(pid_t), compare_int);
             qsort(adjs[pid].adj, na, sizeof(pid_t), compare_int);
             if (memcmp(previous_adj, adjs[pid].adj, sizeof(pid_t)*na) != 0)
               {
                 printf("Inconsistent adjacencies between two divisions.\n");
                 abort();
               }
           }
	} else {
	  printf("0"); fflush(stdout);
	}
      }
      fclose(part);
      printf("%d ",k);
    }
   }
  }
  printf("\n");

  // PMS : remove mock galaxies and anything adjacent to a mock galaxy
  printf("\nRemoving mock galaxies...\n"); 
 
  // completely unlink mock particles
  for (i = mockIndex; i < np; i++) {
    vols[i] = 1.e-29;
    adjs[i].nadj = 0;
  }
  
  // unlink particles adjacent to mock galaxies
  for (i = 0; i < mockIndex; i++) {
    for (j = 0; j < adjs[i].nadj; j++) {
      if (adjs[i].adj[j] > mockIndex) {
//printf("KILLING %d\n", i);
        vols[i] = 1.e-29;
        adjs[i].nadj = 0;
        numRemoved++;
        break;
      }
    }
  }

  // update all other adjacencies
  for (i = 0; i < mockIndex; i++) {
    
    int numAdjSaved = 0;
    for (j = 0; j < adjs[i].nadj; j++) {

      if ( adjs[adjs[i].adj[j]].nadj != 0) {
        adjs[i].adj[numAdjSaved] = adjs[i].adj[j];
        numAdjSaved++;
      }

    }
    adjs[i].nadj = numAdjSaved; 
  }

  printf("Removed %d mock galaxies and %d adjacent galaxies.\n", np-mockIndex, 
                                                                 numRemoved); 
  printf("There are %d galaxies remaining.\n", mockIndex-numRemoved);

  // END PMS

  npnotdone = 0; avgnadj = 0.; avgvol = 0.;
  for (p=0;p<np;p++) {
    // PMS
    if (vols[p] == 1.e-29) continue;
    // END PMS
    if (vols[p] == -1.) npnotdone++;
    avgnadj += (double)(adjs[p].nadj);
    avgvol += (double)(vols[p]);
  }
  if (npnotdone > 0)
    printf("%d particles not done!\n", npnotdone);
  printf("%d particles done more than once.\n",nvpsum-np);
  avgnadj /= (double)np;
  avgvol /= (double)np;
  printf("Average # adjacencies = %lf (%f for Poisson)\n",avgnadj,
	 48.*3.141593*3.141593/35.+2.);
  printf("Average volume = %lf\n",avgvol);
    
  /* Now the output! */

  sprintf(adjfile,"%s/adj%s.dat",outDir,suffix);
  sprintf(volfile,"%s/vol%s.dat",outDir,suffix);

  printf("Outputting to%s, %s\n\n", adjfile, volfile);

  adj = fopen(adjfile,"w");
  if (adj == NULL) {
    printf("Unable to open %s\n",adjfile);
    exit(0);
  }
  fwrite(&mockIndex,1, sizeof(int),adj);
  /* Adjacencies: first the numbers of adjacencies, 
     and the number we're actually going to write per particle */

  // Recount the number of adjacencies after merge
  for(i=0;i<mockIndex;i++)
    cnt_adj[i] = adjs[i].nadj;

  for (i=0;i<mockIndex;i++)
    fwrite(&cnt_adj[i],1,sizeof(int),adj);
    
  /* Now the lists of adjacencies (without double counting) */
  for (i=0;i<mockIndex;i++) {
      nout = 0;
      for (j=0;j<adjs[i].nadj; j++) 
         if (adjs[i].adj[j] > i)
           nout++;
      fwrite(&nout,1,sizeof(int),adj);      
      for (j=0;j<adjs[i].nadj; j++) {
        pid_t id = adjs[i].adj[j];
	if (id > i) 
	  fwrite(&id,1,sizeof(pid_t),adj);
      }
    }

  fclose(adj);
  
  /* Volumes */
  vol = fopen(volfile,"w");
  if (vol == NULL) {
    printf("Unable to open %s\n",volfile);
    exit(0);
  }
// PMS
  fwrite(&mockIndex,1, sizeof(int),vol);
  fwrite(vols,sizeof(float),mockIndex,vol);
  //fwrite(&np,1, sizeof(int),vol);
  //fwrite(vols,sizeof(float),np,vol);
// END PMS

  fclose(vol);

  free(vols);
  free(cnt_adj);
  free(adjs);
  free(orig);

  return(0);
}
