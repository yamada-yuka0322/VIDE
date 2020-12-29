/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/stacking/pruneVoids.cpp
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



// Reads in the void catalog and removes any void that could potentially
//   be affected by a mock particle. It does this by computing the longest
//   particle distance within each void and comparing it to the distance
//   of the nearest mock particle. If the void could potentially by rotated
//   to include this particle, we throw out the void.

// This is placed here instead of using the edgeAvoidance option in 
//   stackVoidsZero so that we can optionally filter the entire
//   catalog at once before the stacking phase. This is useful
//   for producing a "clean" void catalog for other purposes.

#include "gsl/gsl_math.h"
#include "gsl/gsl_eigen.h"
#include "string.h"
#include "ctype.h"
#include "stdlib.h"
#include <math.h>
#include <stdio.h>
#include <netcdf>
#include "pruneVoids_conf.h"
#include <vector>
#include "assert.h"
#include "voidTree.hpp"
#include "loadZobov.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>

#define LIGHT_SPEED 299792.458
#define MPC2Z 100./LIGHT_SPEED
#define Z2MPC LIGHT_SPEED/100.

#define CENTRAL_VOID 1
#define EDGE_VOID 2

using namespace std;

typedef struct partStruct {
  float x, y, z, vol;
  int nadj, ncnt;
  int *adj;
} PART;

typedef struct zoneStruct {
  int numPart;
  int *partIDs;
} ZONE2PART;

typedef struct voidZoneStruct {
  int numZones;
  int *zoneIDs;
} VOID2ZONE;

typedef struct voidStruct {
  float vol, coreDens, zoneVol, densCon, voidProb, radius;
  float rescaledCoreDens;
  int voidID, numPart, numZones, coreParticle, zoneNumPart;
  float maxRadius, nearestMock, centralDen, redshift, redshiftInMpc;
  float nearestMockFromCore, nearestGalFromCore;
  float nearestEdge;
  float center[3], macrocenter[3];
  int accepted;
  int voidType;
  int parentID, numChildren, level;
  bool isLeaf, hasHighCentralDen;
  gsl_vector *eval;
  gsl_matrix *evec;
  float ellip;
} VOID;

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

void openFiles(string outputDir, string sampleName, 
               string prefix, string dataPortion, 
               int mockIndex, int numKept,
               FILE** fpZobov, FILE** fpCenters, 
               FILE** fpCentersNoCut, 
               FILE** fpBarycenter, FILE** fpDistances, FILE** fpShapes, 
               FILE** fpSkyPositions);

void closeFiles(FILE* fpZobov, FILE* fpCenters, 
                FILE* fpCentersNoCut, 
                FILE* fpBarycenter, FILE* fpDistances, FILE* fpShapes, 
                FILE* fpSkyPositions);

void outputVoids(string outputDir, string sampleName, string prefix, 
                string dataPortion, int mockIndex, 
                vector<VOID> voids, 
                bool isObservation, double *boxLen,
                bool doTrim, bool doCentralDenCut);

int main(int argc, char **argv) {

  // initialize arguments
  pruneVoids_info args;
  pruneVoids_conf_params args_params;

  pruneVoids_conf_init(&args);
  pruneVoids_conf_params_init(&args_params);

  args_params.check_required = 0;
  if (pruneVoids_conf_ext (argc, argv, &args, &args_params))
    return 1;

  if (!args.configFile_given) {
    if (pruneVoids_conf_required (&args,
                                           PRUNEVOIDS_CONF_PACKAGE))
      return 1;
  } else {
    args_params.check_required = 1;
    args_params.initialize = 0;
    if (pruneVoids_conf_config_file (args.configFile_arg,
                                              &args,
                                              &args_params))
    return 1;
  }

  // initialize cosmology integrator and interpolator
  gsl_function expanF;
  expanF.function = &expanFun;
  struct my_expan_params expanParams;
  expanParams.Om = args.omegaM_arg;
  expanParams.w0 = -1.0;
  expanParams.wa = 0.0;
  expanF.params = &expanParams;
  double result, error;
  size_t nEval;

  int iZ, numZ = 4000;
  double maxZ = 5.0, z, *dL, *redshifts;
  dL = (double *) malloc(numZ * sizeof(double));  
  redshifts = (double *) malloc(numZ * sizeof(double));
  for (iZ = 0; iZ < numZ; iZ++) {
    z = iZ * maxZ/numZ;
    gsl_integration_qng(&expanF, 0.0, z, 1.e-6, 1.e-6, &result, &error, &nEval);
    dL[iZ] = result*LIGHT_SPEED/100.;
    //printf("HERE %e %e\n", z, dL[iZ]);
    redshifts[iZ] = z;
  }
  gsl_interp *interp = gsl_interp_alloc(gsl_interp_linear, 4000);
  gsl_interp_init(interp, dL, redshifts, numZ);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();

 
  int i, p, p2, numPartTot, numZonesTot, dummy, iVoid;
  int numVoids, mockIndex;
  double tolerance;
  FILE *fp;
  PART *part, *voidPart;
  ZONE2PART *zones2Parts;
  VOID2ZONE *void2Zones;
  vector<VOID> voids;
  float *temp, junk, voidVol;
  int junkInt, voidID, numPart, numZones, zoneID, partID, maxNumPart;
  int coreParticle, zoneNumPart;
  float coreDens, zoneVol, densCon, voidProb, dist[3], dist2, minDist, maxDist;
  float centralRad, centralDen;
  double nearestEdge, redshift;
  char line[500], junkStr[10];
  string outputDir, sampleName, dataPortion, prefix;
  int mask_index;
  double ranges[3][2], boxLen[3], mul; 
  double volNorm, radius;
  int clock1, clock2, clock3, clock4;
  double interval;
  int periodicX=0, periodicY=0, periodicZ=0;
  string dataPortions[2];

  gsl_eigen_symmv_workspace *eigw = gsl_eigen_symmv_alloc(3);

  numVoids = args.numVoids_arg;
  mockIndex = args.mockIndex_arg;
  tolerance = args.tolerance_arg;

  clock1 = clock();
  printf("Pruning parameters: %f %f %f %s\n", args.zMin_arg, 
                                             args.zMax_arg,
                                             args.rMin_arg,
                                             args.periodic_arg);

  // check for periodic box
  periodicX = 0;
  periodicY = 0;
  periodicZ = 0;
  if (!args.isObservation_flag) {
    if ( strchr(args.periodic_arg, 'x') != NULL) {
      periodicX = 1;
      printf("Will assume x-direction is periodic.\n");
    }
    if ( strchr(args.periodic_arg, 'y') != NULL) {
      periodicY = 1;
      printf("Will assume y-direction is periodic.\n");
    }
    if ( strchr(args.periodic_arg, 'z') != NULL) {
      periodicZ = 1;
      printf("Will assume z-direction is periodic.\n");
    }
  }

  // load box size
  printf("\n Getting info...\n");
  netCDF::NcFile f_info(args.extraInfo_arg, netCDF::NcFile::read);
  f_info.getAtt("range_x_min").getValues(&ranges[0][0]);
  f_info.getAtt("range_x_max").getValues(&ranges[0][1]);
  f_info.getAtt("range_y_min").getValues(&ranges[1][0]);
  f_info.getAtt("range_y_max").getValues(&ranges[1][1]);
  f_info.getAtt("range_z_min").getValues(&ranges[2][0]);
  f_info.getAtt("range_z_max").getValues(&ranges[2][1]);

  printf(" Range xmin %e\n", ranges[0][0]);
  printf(" Range xmax %e\n", ranges[0][1]);
  printf(" Range ymin %e\n", ranges[1][0]);
  printf(" Range ymax %e\n", ranges[1][1]);
  printf(" Range zmin %e\n", ranges[2][0]);
  printf(" Range zmax %e\n", ranges[2][1]);

  boxLen[0] = ranges[0][1] - ranges[0][0];
  boxLen[1] = ranges[1][1] - ranges[1][0];
  boxLen[2] = ranges[2][1] - ranges[2][0];
 
  // read in all particle positions
  clock3 = clock();
  printf("\n Loading particles...\n");
  fp = fopen(args.partFile_arg, "r");
  fread(&dummy, 1, 4, fp); 
  fread(&numPartTot, 1, 4, fp);
  fread(&dummy, 1, 4, fp);

  part = (PART *) malloc(numPartTot * sizeof(PART));
  temp = (float *) malloc(numPartTot * sizeof(float));
 
  volNorm = numPartTot/(boxLen[0]*boxLen[1]*boxLen[2]);
  printf(" VOL NORM = %f\n", volNorm);

  printf(" CENTRAL DEN = %f\n", args.maxCentralDen_arg);

  fread(&dummy, 1, 4, fp);
  fread(temp, numPartTot, 4, fp);
  mul = ranges[0][1] - ranges[0][0];
  for (p = 0; p < numPartTot; p++)
    part[p].x = mul*temp[p];
  fread(&dummy, 1, 4, fp);
  fread(&dummy, 1, 4, fp);
  fread(temp, numPartTot, 4, fp);
  mul = ranges[1][1] - ranges[1][0];
  for (p = 0; p < numPartTot; p++) 
    part[p].y = mul*temp[p];
  fread(&dummy, 1, 4, fp);
  fread(&dummy, 1, 4, fp);
  fread(temp, numPartTot, 4, fp);
  mul = ranges[2][1] - ranges[2][0];
  for (p = 0; p < numPartTot; p++) 
    part[p].z = mul*temp[p];

  if (!args.isObservation_flag) {
    for (p = 0; p < numPartTot; p++)  {
      part[p].x += ranges[0][0];
      part[p].y += ranges[1][0];
      part[p].z += ranges[2][0];
    }
  }
  fclose(fp); 

  clock4 = clock();
  interval = 1.*(clock4 - clock3)/CLOCKS_PER_SEC;
  printf(" Read %d particles (%.2f sec)...\n", numPartTot, interval);
  

  if (mockIndex == -1) mockIndex = numPartTot;
  
  // read in desired voids
  clock3 = clock();
  printf(" Loading voids...\n");
  fp = fopen(args.voidDesc_arg ,"r");
  fgets(line, sizeof(line), fp);
  sscanf(line, "%d %s %d %s", &junkInt, junkStr, &junkInt, junkStr);
  fgets(line, sizeof(line), fp);

  voids.resize(numVoids);
  i = 0;
  while (fgets(line, sizeof(line), fp) != NULL) {
    sscanf(line, "%d %d %d %f %f %d %d %f %d %f %f\n", &iVoid, &voidID, 
           &coreParticle, &coreDens, &zoneVol, &zoneNumPart, &numZones, 
           &voidVol, &numPart, &densCon, &voidProb);

    i++;
    voids[i-1].coreParticle = coreParticle;
    voids[i-1].zoneNumPart = zoneNumPart;
    voids[i-1].coreDens = coreDens;
    voids[i-1].zoneVol = zoneVol;
    voids[i-1].voidID = voidID;
    voids[i-1].vol = voidVol;
    voids[i-1].numPart = numPart;
    voids[i-1].numZones = numZones;
    voids[i-1].densCon = densCon;
    voids[i-1].voidProb = voidProb;

    voids[i-1].radius = pow(voidVol/volNorm*3./4./M_PI, 1./3.);
    voids[i-1].accepted = 1;

    voids[i-1].isLeaf = true;
    voids[i-1].hasHighCentralDen = false;
    voids[i-1].numChildren = 0;
    voids[i-1].parentID = -1;  
    voids[i-1].level = 0;  

    voids[i-1].eval = gsl_vector_alloc(3);
    voids[i-1].evec = gsl_matrix_alloc(3,3);
    voids[i-1].ellip = 0;
  }
  fclose(fp);

  // load up the zone membership for each void
  printf(" Loading void-zone membership info...\n");
  fp = fopen(args.void2Zone_arg, "r");
  fread(&numZonesTot, 1, 4, fp);
  void2Zones = (VOID2ZONE *) malloc(numZonesTot * sizeof(VOID2ZONE));

  for (iZ = 0; iZ < numZonesTot; iZ++) {
    fread(&numZones, 1, 4, fp);

    void2Zones[iZ].numZones = numZones;
    void2Zones[iZ].zoneIDs = (int *) malloc(numZones * sizeof(int));

    for (p = 0; p < numZones; p++) {
      fread(&void2Zones[iZ].zoneIDs[p], 1, 4, fp);
    }
  }
  fclose(fp);

  // now the particles-zone
  printf(" Loading particle-zone membership info...\n");
  fp = fopen(args.zone2Part_arg, "r");
  fread(&dummy, 1, 4, fp);
  fread(&numZonesTot, 1, 4, fp);
 
  zones2Parts = (ZONE2PART *) malloc(numZonesTot * sizeof(ZONE2PART));
  for (iZ = 0; iZ < numZonesTot; iZ++) {
    fread(&numPart, 1, 4, fp);
    
    zones2Parts[iZ].numPart = numPart;
    zones2Parts[iZ].partIDs = (int *) malloc(numPart * sizeof(int));

    for (p = 0; p < numPart; p++) {
      fread(&zones2Parts[iZ].partIDs[p], 1, 4, fp);
    }
  }

  // and finally volumes
  printf(" Loading particle volumes...\n");
  fp = fopen(args.partVol_arg, "r");
  fread(&mask_index, 1, 4, fp);
  if (mask_index != mockIndex) {
    printf("NON-MATCHING MOCK INDICES!? %d %d\n", mask_index, mockIndex);
    exit(-1);
  } 
  for (p = 0; p < mask_index; p++) {
    fread(&temp[0], 1, 4, fp);
    part[p].vol = temp[0];
  }
  fclose(fp);

/*
  // and finally finally adjacencies
  printf(" Loading particle adjacencies...\n");
  fp = fopen(args.partAdj_arg, "r");
  fread(&mask_index, 1, 4, fp);
  if (mask_index != mockIndex) {
    printf("NON-MATCHING MOCK INDICES!? %d %d\n", mask_index, mockIndex);
    exit(-1);
  } 
  int tempInt;
  for (p = 0; p < mask_index; p++) {
    fread(&tempInt, 1, 4, fp);
    part[p].nadj = tempInt;
    part[p].ncnt = 0;
    if (part[p].nadj > 0)
      part[p].adj = (int *) malloc(part[p].nadj * sizeof(int));
  }
  for (p = 0; p < mask_index; p++) {
    fread(&tempInt, 1, 4, fp);
    int nin = tempInt;
    if (nin > 0) {
      for (int nAdj = 0; nAdj < nin; nAdj++) {
        int tempAdj;
        fread(&tempAdj, 1, 4, fp);

        // this bit has been readjusted just in case we are 
        //   accidentally still linking to mock particles

        //if (tempAdj < mask_index) {
          assert(p < tempAdj);
          //if (part[p].ncnt == part[p].nadj) {
          //  printf("OVERFLOW %d\n", p);
          //} else if (part[tempAdj].ncnt == part[tempAdj].nadj) {
          //  printf("OVERFLOW %d\n", tempAdj);
          //} else {
            part[p].adj[part[p].ncnt] = tempAdj;
            part[p].ncnt++;
          if (tempAdj < mask_index) {
            part[tempAdj].adj[part[tempAdj].ncnt] = p;
            part[tempAdj].ncnt++;
          }
          //}
        //} 
      }
//printf("ADJ %d %d %d %d %d\n", p, nin, part[p].nadj, nAdj, tempInt);
    }
  }
  fclose(fp);
*/

  clock4 = clock();
  interval = 1.*(clock4 - clock3)/CLOCKS_PER_SEC;
  printf(" Read voids (%.2f sec)...\n", interval);

  // load voids *again* using Guilhem's code so we can get tree
  clock3 = clock();
  //if (!args.isObservation_flag) {
    printf(" Re-loading voids and building tree..\n");
    ZobovRep zobovCat;
    if (!loadZobov(args.voidDesc_arg, args.zone2Part_arg,
                   args.void2Zone_arg, 
                   0, zobovCat)) {
       printf("Error loading catalog!\n");
       return -1;
    }
    VoidTree *tree;
    tree = new VoidTree(zobovCat);
    zobovCat.allZones.erase(zobovCat.allZones.begin(), zobovCat.allZones.end());
 
    // copy tree information to our own data structures
    for (iVoid = 0; iVoid < numVoids; iVoid++) {
      voidID = voids[iVoid].voidID;
      voids[iVoid].parentID = tree->getParent(voidID);
      voids[iVoid].numChildren = tree->getChildren(voidID).size();

      // compute level in tree
      int level = 0;
      int parentID = tree->getParent(voidID);
      while (parentID != -1) {
        level++;
        parentID = tree->getParent(parentID);
      }
      voids[iVoid].level = level;
    }
  //} // end re-load
  clock4 = clock();
  interval = 1.*(clock4 - clock3)/CLOCKS_PER_SEC;
  printf(" Re-read voids (%.2f sec)...\n", interval);

  // check boundaries
  printf(" Computing void properties...\n");

  maxNumPart = 0;
  for (iVoid = 0; iVoid < numVoids; iVoid++) {
    if (voids[iVoid].numPart > maxNumPart) maxNumPart = voids[iVoid].numPart;
  }
  
  voidPart = (PART *) malloc(maxNumPart * sizeof(PART));
  
  for (iVoid = 0; iVoid < numVoids; iVoid++) {

    voidID = voids[iVoid].voidID;
    printf("  DOING %d (of %d) %d %d %f\n", iVoid, numVoids, voidID, 
                                           voids[iVoid].numPart, 
                                           voids[iVoid].radius);

    voids[iVoid].center[0] = part[voids[iVoid].coreParticle].x;
    voids[iVoid].center[1] = part[voids[iVoid].coreParticle].y;
    voids[iVoid].center[2] = part[voids[iVoid].coreParticle].z;
 
    // first load up particles into a buffer 
    clock3 = clock();
    i = 0;
    for (iZ = 0; iZ < void2Zones[voidID].numZones; iZ++) {
      zoneID = void2Zones[voidID].zoneIDs[iZ];

      for (p = 0; p < zones2Parts[zoneID].numPart; p++) {
        partID = zones2Parts[zoneID].partIDs[p];

        if (partID > mask_index || 
            (part[partID].vol < 1.e-27 && part[partID].vol > 0.)) {
           printf("BAD PART!? %d %d %e", partID, mask_index, part[partID].vol);
           exit(-1);
        }

        voidPart[i].x = part[partID].x;
        voidPart[i].y = part[partID].y;
        voidPart[i].z = part[partID].z;
        voidPart[i].vol = part[partID].vol;

/*
        // testing for edge contamination
        if (part[partID].vol < 1.e-27)  {
           printf("CONTAMINATED!! %d %d\n", iVoid, partID);
        } else {
           //printf("NORMAL!! %d %d %e\n", iVoid, partID, part[partID].vol);
        }
        for (int iAdj = 0; iAdj < part[partID].ncnt; iAdj++) {
          if (part[partID].adj[iAdj] > mockIndex) {
            printf("CONTAMINATED!! %d %d %d\n", iVoid, partID, iAdj);
          } 
        }
*/
        i++;
      }
    }

    clock4 = clock();
    interval = 1.*(clock4 - clock3)/CLOCKS_PER_SEC;
    //printf("   %.2f for buffer\n", interval);

    // compute macrocenters
    clock3 = clock();
    double weight = 0.;
    voids[iVoid].macrocenter[0] = 0.;
    voids[iVoid].macrocenter[1] = 0.;
    voids[iVoid].macrocenter[2] = 0.;

    for (p = 0; p < voids[iVoid].numPart; p++) {
      dist[0] = voidPart[p].x - voids[iVoid].center[0];
      dist[1] = voidPart[p].y - voids[iVoid].center[1];
      dist[2] = voidPart[p].z - voids[iVoid].center[2];

      if (periodicX && fabs(dist[0]) > boxLen[0]/2.) 
          dist[0] = dist[0] - copysign(boxLen[0], dist[0]);
      if (periodicY && fabs(dist[1]) > boxLen[1]/2.) 
          dist[1] = dist[1] - copysign(boxLen[1], dist[1]);
      if (periodicZ && fabs(dist[2]) > boxLen[2]/2.) 
          dist[2] = dist[2] - copysign(boxLen[2], dist[2]);

      voids[iVoid].macrocenter[0] += voidPart[p].vol*(dist[0]);
      voids[iVoid].macrocenter[1] += voidPart[p].vol*(dist[1]);
      voids[iVoid].macrocenter[2] += voidPart[p].vol*(dist[2]);
      weight += voidPart[p].vol;
    }
    voids[iVoid].macrocenter[0] /= weight;
    voids[iVoid].macrocenter[1] /= weight;
    voids[iVoid].macrocenter[2] /= weight;
    voids[iVoid].macrocenter[0] += voids[iVoid].center[0];
    voids[iVoid].macrocenter[1] += voids[iVoid].center[1];
    voids[iVoid].macrocenter[2] += voids[iVoid].center[2];

    if (periodicX) {
      if (voids[iVoid].macrocenter[0] > ranges[0][1])
        voids[iVoid].macrocenter[0] = voids[iVoid].macrocenter[0] - boxLen[0];
      if (voids[iVoid].macrocenter[0] < ranges[0][0])
        voids[iVoid].macrocenter[0] = boxLen[0] + voids[iVoid].macrocenter[0];
    }
    if (periodicY) {
      if (voids[iVoid].macrocenter[1] > ranges[1][1])
        voids[iVoid].macrocenter[1] = voids[iVoid].macrocenter[1] - boxLen[1];
      if (voids[iVoid].macrocenter[1] < ranges[1][0])
        voids[iVoid].macrocenter[1] = boxLen[1] + voids[iVoid].macrocenter[1];
    }
    if (periodicZ) {
      if (voids[iVoid].macrocenter[2] > ranges[2][1])
        voids[iVoid].macrocenter[2] = voids[iVoid].macrocenter[2] - boxLen[2];
      if (voids[iVoid].macrocenter[2] < ranges[2][0])
        voids[iVoid].macrocenter[2] = boxLen[2] + voids[iVoid].macrocenter[2];
    }
    clock4 = clock();
    interval = 1.*(clock4 - clock3)/CLOCKS_PER_SEC;
    //printf("   %.2f for macrocenter\n", interval);

    // compute central density
    clock3 = clock();
    centralRad = voids[iVoid].radius/args.centralRadFrac_arg;
    centralDen = 0.;
    int numCentral = 0;
    for (p = 0; p < voids[iVoid].numPart; p++) {
      dist[0] = fabs(voidPart[p].x - voids[iVoid].macrocenter[0]);
      dist[1] = fabs(voidPart[p].y - voids[iVoid].macrocenter[1]);
      dist[2] = fabs(voidPart[p].z - voids[iVoid].macrocenter[2]);

      if (periodicX) dist[0] = fmin(dist[0], boxLen[0]-dist[0]);
      if (periodicY) dist[1] = fmin(dist[1], boxLen[1]-dist[1]);
      if (periodicZ) dist[2] = fmin(dist[2], boxLen[2]-dist[2]);

      dist2 = pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2);
      if (sqrt(dist2) < centralRad) numCentral += 1;
    }
    voids[iVoid].centralDen = numCentral / (volNorm*4./3. * M_PI * 
                                            pow(centralRad, 3.));

    clock4 = clock();
    interval = 1.*(clock4 - clock3)/CLOCKS_PER_SEC;
    //printf("   %.2f for central density\n", interval);

    //coreParticle = voids[iVoid].coreParticle;
    //voids[iVoid].rescaledCoreDens = voids[iVoid].coreDens*(pow(1.*mockIndex/numPartTot,3));
    //  // compute distance from core to nearest mock
    //  minDist = 1.e99;
    //  for (p = mockIndex; p < numPartTot; p++) {
    //    dist[0] = part[coreParticle].x - part[p].x;
    //    dist[1] = part[coreParticle].y - part[p].y;
    //    dist[2] = part[coreParticle].z - part[p].z;
    //
    //    dist2 = pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2);
    //    if (dist2 < minDist) minDist = dist2;
    //  }
    //  voids[iVoid].nearestMockFromCore = sqrt(minDist);
    // 
    //  // compute distance from core to nearest mock
    //  minDist = 1.e99;
    //  for (p = 0; p < mockIndex; p++) {
    //    dist[0] = part[coreParticle].x - part[p].x;
    //    dist[1] = part[coreParticle].y - part[p].y;
    //    dist[2] = part[coreParticle].z - part[p].z;
    //
    //    dist2 = pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2);
    //    if (dist2 < minDist && dist2 > 1.e-10) minDist = dist2;
    //  }
    //  voids[iVoid].nearestGalFromCore = sqrt(minDist);

    // compute maximum extent
/*
    if (args.isObservation_flag) {
      maxDist = 0.;
      for (p = 0; p < voids[iVoid].numPart; p++) {
      for (p2 = p; p2 < voids[iVoid].numPart; p2++) {
  
        dist[0] = voidPart[p].x - voidPart[p2].x;
        dist[1] = voidPart[p].y - voidPart[p2].y;
        dist[2] = voidPart[p].z - voidPart[p2].z;

        dist2 = pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2);
        if (dist2 > maxDist) maxDist = dist2;
      }
      }
      voids[iVoid].maxRadius = sqrt(maxDist)/2.;
    } else {
*/

      clock3 = clock();
      maxDist = 0.;
      for (p = 0; p < voids[iVoid].numPart; p++) {
  
        dist[0] = fabs(voidPart[p].x - voids[iVoid].macrocenter[0]);
        dist[1] = fabs(voidPart[p].y - voids[iVoid].macrocenter[1]);
        dist[2] = fabs(voidPart[p].z - voids[iVoid].macrocenter[2]);

        if (periodicX) dist[0] = fmin(dist[0], boxLen[0]-dist[0]);
        if (periodicY) dist[1] = fmin(dist[1], boxLen[1]-dist[1]);
        if (periodicZ) dist[2] = fmin(dist[2], boxLen[2]-dist[2]);

        dist2 = pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2);
        if (dist2 > maxDist) maxDist = dist2;
      }
      voids[iVoid].maxRadius = sqrt(maxDist);
      clock4 = clock();
      interval = 1.*(clock4 - clock3)/CLOCKS_PER_SEC;
      //printf("   %.2f for maximum extent\n", interval);
//    }
   
    clock3 = clock(); 
    if (args.isObservation_flag) {
      // compute distance from center to nearest mock
      minDist = 1.e99;
      for (p = mockIndex; p < numPartTot; p++) {
        dist[0] = voids[iVoid].macrocenter[0] - part[p].x;
        dist[1] = voids[iVoid].macrocenter[1] - part[p].y;
        dist[2] = voids[iVoid].macrocenter[2] - part[p].z;

        dist2 = pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2);
        if (dist2 < minDist) minDist = dist2;
      }
      voids[iVoid].nearestMock = sqrt(minDist);
    
    } else {
      voids[iVoid].nearestMock = 1.e99;
    }
 
    if (args.isObservation_flag) {
      voids[iVoid].redshiftInMpc = 
                        sqrt(pow(voids[iVoid].macrocenter[0] - boxLen[0]/2.,2) + 
                             pow(voids[iVoid].macrocenter[1] - boxLen[1]/2.,2) + 
                             pow(voids[iVoid].macrocenter[2] - boxLen[2]/2.,2));
      voids[iVoid].redshiftInMpc = voids[iVoid].redshiftInMpc;


      if (args.useComoving_flag) {
        redshift = gsl_interp_eval(interp, dL, redshifts,
                 voids[iVoid].redshiftInMpc, acc);
        //printf("HELLO %e %e\n", redshift, args.zMax_arg);
        nearestEdge = fabs((redshift-args.zMax_arg)*LIGHT_SPEED/100.); 
        voids[iVoid].redshift = redshift;
      } else {
        redshift = voids[iVoid].redshiftInMpc;
        nearestEdge = fabs(redshift-args.zMax_arg*LIGHT_SPEED/100.); 
        voids[iVoid].redshift = voids[iVoid].redshiftInMpc/LIGHT_SPEED*100.;
      }
      //nearestEdge = fmin(fabs(redshift-args.zMin_arg*LIGHT_SPEED/100.), 
      //                   fabs(redshift-args.zMax_arg*LIGHT_SPEED/100.)); 

    } else {

      voids[iVoid].redshiftInMpc = voids[iVoid].macrocenter[2];
      if (args.useComoving_flag) {
        voids[iVoid].redshift = gsl_interp_eval(interp, dL, redshifts,
                 voids[iVoid].redshiftInMpc, acc);
      } else {
        voids[iVoid].redshift = voids[iVoid].macrocenter[2]/LIGHT_SPEED*100.;
      }

      nearestEdge = 1.e99;
     
      if (!periodicX) {
        nearestEdge = fmin(nearestEdge,
                           fabs(voids[iVoid].macrocenter[0] - ranges[0][0]));
        nearestEdge = fmin(nearestEdge,
                           fabs(voids[iVoid].macrocenter[0] - ranges[0][1]));
      }
      if (!periodicY) {
        nearestEdge = fmin(nearestEdge,
                           fabs(voids[iVoid].macrocenter[1] - ranges[1][0]));
        nearestEdge = fmin(nearestEdge,
                           fabs(voids[iVoid].macrocenter[1] - ranges[1][1]));
      }
      if (!periodicZ) {
        nearestEdge = fmin(nearestEdge,
                           fabs(voids[iVoid].macrocenter[2] - ranges[2][0]));
        nearestEdge = fmin(nearestEdge,
                           fabs(voids[iVoid].macrocenter[2] - ranges[2][1]));
      }
    }

    voids[iVoid].nearestEdge = nearestEdge;

    clock4 = clock();
    interval = 1.*(clock4 - clock3)/CLOCKS_PER_SEC;
    //printf("   %.2f for nearest edge\n", interval);

    // compute eigenvalues and vectors for orientation and shape
    clock3 = clock();
    double inertia[9];
    for (int i = 0; i < 9; i++) inertia[i] = 0.;

    for (int p = 0; p < voids[iVoid].numPart; p++) {
      dist[0] = voidPart[p].x - voids[iVoid].macrocenter[0];
      dist[1] = voidPart[p].y - voids[iVoid].macrocenter[1];
      dist[2] = voidPart[p].z - voids[iVoid].macrocenter[2];

      if (periodicX && fabs(dist[0]) > boxLen[0]/2.)
          dist[0] = dist[0] - copysign(boxLen[0], dist[0]);
      if (periodicY && fabs(dist[1]) > boxLen[1]/2.)
          dist[1] = dist[1] - copysign(boxLen[1], dist[1]);
      if (periodicZ && fabs(dist[2]) > boxLen[2]/2.)
          dist[2] = dist[2] - copysign(boxLen[2], dist[2]);

      for (int i = 0; i < 3; i++) {
      for (int j = i; j < 3; j++) {
        if (i == j) inertia[i*3+j] += dist[0]*dist[0] + dist[1]*dist[1] + 
                                    dist[2]*dist[2];
        inertia[i*3+j] -= dist[i]*dist[j];      
      }
      }
    }
    inertia[1*3+0] = inertia[0*3+1];
    inertia[2*3+0] = inertia[0*3+2];
    inertia[2*3+1] = inertia[1*3+2];

    gsl_matrix_view m = gsl_matrix_view_array(inertia, 3, 3);
    gsl_eigen_symmv(&m.matrix, voids[iVoid].eval, voids[iVoid].evec, eigw);
   
    float a = sqrt(2.5*(gsl_vector_get(voids[iVoid].eval,1) +
                        gsl_vector_get(voids[iVoid].eval,2) - 
                        gsl_vector_get(voids[iVoid].eval,0)));
    float b = sqrt(2.5*(gsl_vector_get(voids[iVoid].eval,2) +
                        gsl_vector_get(voids[iVoid].eval,0) - 
                        gsl_vector_get(voids[iVoid].eval,1)));
    float c = sqrt(2.5*(gsl_vector_get(voids[iVoid].eval,0) +
                        gsl_vector_get(voids[iVoid].eval,1) - 
                        gsl_vector_get(voids[iVoid].eval,2)));
    float ca;
    float cb = c/b;

    float smallest = 1.e99;
    float largest = 0.;
    for (int i = 0; i < 3; i++) {
      if (gsl_vector_get(voids[iVoid].eval,i) < smallest) 
        smallest = gsl_vector_get(voids[iVoid].eval,i);
      if (gsl_vector_get(voids[iVoid].eval,i) > largest) 
        largest = gsl_vector_get(voids[iVoid].eval,i);
    }
    // TEST
    voids[iVoid].ellip = 1.0 - sqrt(sqrt(fabs(smallest/largest)));

    //if (a < c)  ca = a/c;
    //if (a >= c) ca = c/a;
    //voids[iVoid].ellip = fabs(1.0 - ca);

    //if (a < c)  ca = a*a/(c*c);
    //if (a >= c) ca = (c*c)/(a*a);
    //voids[iVoid].ellip = sqrt(fabs(1.0 - ca));

    clock4 = clock();
    interval = 1.*(clock4 - clock3)/CLOCKS_PER_SEC;
    //printf("   %.2f for ellipticity\n", interval);
  } // iVoid

  gsl_eigen_symmv_free(eigw);
    
  int numWrong = 0;
  int numHighDen = 0;
  int numCentral = 0;
  int numEdge = 0;
  int numNearZ = 0;
  int numAreParents = 0;
  int numTooSmall = 0;

  printf(" Picking winners and losers...\n");
  printf("  Starting with %d voids\n", (int) voids.size());

  for (iVoid = 0; iVoid < voids.size(); iVoid++) {
    voids[iVoid].accepted = 1;
  }

/*
  int j = 0;
  for (iVoid = 0; iVoid < voids.size(); iVoid++) {
    if (voids[iVoid].densCon < 1.5) {
//      voids[iVoid].accepted = -4;
    }
  }
*/

  // toss out voids that are obviously wrong
  int iGood = 0;
  for (iVoid = 0; iVoid < voids.size(); iVoid++) {
    if (voids[iVoid].densCon > 1.e4 || isnan(voids[iVoid].vol) || 
         isinf(voids[iVoid].vol)) {
      numWrong++;
    } else {
      voids[iGood++] = voids[iVoid];
    }
  }
  voids.resize(iGood);
  printf("  1st filter: rejected %d obviously bad\n", numWrong); 

  iGood = 0;
  for (iVoid = 0; iVoid < voids.size(); iVoid++) {
    if (voids[iVoid].radius < args.rMin_arg) {
      numTooSmall++;
    } else {
      voids[iGood++] = voids[iVoid];
    }
  }
  voids.resize(iGood);
  printf("  2nd filter: rejected %d too small\n", numTooSmall); 


  iGood = 0;
  for (iVoid = 0; iVoid < voids.size(); iVoid++) {
    // *always* clean out near edges since there are no mocks there
    if (tolerance*voids[iVoid].maxRadius > voids[iVoid].nearestEdge || 
        tolerance*voids[iVoid].radius > voids[iVoid].nearestEdge) {
      numNearZ++;
    } else {
      voids[iGood++] = voids[iVoid];
    }
  }
  voids.resize(iGood);
  printf("  3rd filter: rejected %d too close to high redshift boundaries\n", numNearZ); 

  numNearZ = 0;
  iGood = 0;
  for (iVoid = 0; iVoid < voids.size(); iVoid++) {
    // assume the lower z-boundary is "soft" in observations
    if (args.isObservation_flag && 
        voids[iVoid].redshift < args.zMin_arg) {
      numNearZ++;
     } else {
      voids[iGood++] = voids[iVoid];
    }
  }
  voids.resize(iGood);

  //Maubert - Uncommented this part : to be sure that voids do not cross maximum redshift asked for in zrange
  iGood = 0;
  for (iVoid = 0; iVoid < voids.size(); iVoid++) {
    // just in case
    if (args.isObservation_flag && 
        voids[iVoid].redshift > args.zMax_arg) {
      numNearZ++;
     } else {
      voids[iGood++] = voids[iVoid];
    }
  }
  voids.resize(iGood);
  // Maubert - End of Uncommented part  


  printf("  4th filter: rejected %d outside redshift boundaries\n", numNearZ); 

  // take only top-level voids
  numAreParents = 0;
  iGood = 0;
  for (iVoid = 0; iVoid < voids.size(); iVoid++) {
    if (voids[iVoid].parentID != -1) {
      numAreParents++;
      voids[iVoid].isLeaf = true;
     } else {
      voids[iVoid].isLeaf = false;
    }
  }


  for (iVoid = 0; iVoid < voids.size(); iVoid++) {
    if (voids[iVoid].centralDen > args.maxCentralDen_arg) {
      voids[iVoid].accepted = -1;
      voids[iVoid].hasHighCentralDen = true;
      numHighDen++;
    } else {
      voids[iVoid].hasHighCentralDen = false;
    }
  }

  for (iVoid = 0; iVoid < voids.size(); iVoid++) {
    if (tolerance*voids[iVoid].maxRadius < voids[iVoid].nearestMock) {
      voids[iVoid].voidType = CENTRAL_VOID;
      numCentral++;
    } else {
      voids[iVoid].voidType = EDGE_VOID;
      numEdge++;
    }
  }

  printf("  Number kept: %d (out of %d)\n", (int) voids.size(), numVoids);
  printf("   We have %d edge voids\n", numEdge); 
  printf("   We have %d central voids\n", numCentral); 
  printf("   We have %d too high central density\n", numHighDen); 
  printf("   We have %d that are not leaf nodes\n", numAreParents); 
 
  
  outputDir = string(args.outputDir_arg);
  sampleName = (string(args.sampleName_arg)+".out");

  dataPortions[0] = "central";
  dataPortions[1] = "all";

  printf(" Output fully trimmed catalog...\n");
  prefix = "";
  for (int i = 0; i < 2; i++) {
    dataPortion = dataPortions[i]; 
  
    outputVoids(outputDir, sampleName, prefix, dataPortion,
                mockIndex,
                voids,
                args.isObservation_flag, boxLen, true, true);
  }

  printf(" Output fully untrimmed catalog...\n");
  prefix = "untrimmed_";
  for (int i = 0; i < 2; i++) {
    dataPortion = dataPortions[i]; 
  
    outputVoids(outputDir, sampleName, prefix, dataPortion,
                mockIndex,
                voids,
                args.isObservation_flag, boxLen, false, false);
  }

  printf(" Output partially-trimmed catalogs...\n");
  prefix = "untrimmed_dencut_";
  for (int i = 0; i < 2; i++) {
    dataPortion = dataPortions[i]; 
  
    outputVoids(outputDir, sampleName, prefix, dataPortion,
                mockIndex,
                voids,
                args.isObservation_flag, boxLen, false, true);
  }

  prefix = "trimmed_nodencut_";
  for (int i = 0; i < 2; i++) {
    dataPortion = dataPortions[i]; 
  
    outputVoids(outputDir, sampleName, prefix, dataPortion,
                mockIndex,
                voids,
                args.isObservation_flag, boxLen, true, false);
  }





  clock2 = clock();
  printf(" Time: %f sec (for %d voids)\n", 
         (1.*clock2-clock1)/CLOCKS_PER_SEC, numVoids);
  clock1 = clock();


  printf("Done!\n");
} // end main


// ----------------------------------------------------------------------------
void openFiles(string outputDir, string sampleName, 
               string prefix, string dataPortion, 
               int mockIndex, int numKept,
               FILE** fpZobov, FILE** fpCenters, 
               FILE** fpBarycenter, FILE** fpDistances, FILE** fpShapes, 
               FILE** fpSkyPositions) {

  *fpZobov = fopen((outputDir+"/"+prefix+"voidDesc_"+dataPortion+"_"+sampleName).c_str(), "w");
  fprintf(*fpZobov, "%d particles, %d voids.\n", mockIndex, numKept);
  fprintf(*fpZobov, "Void# FileVoid# CoreParticle CoreDens ZoneVol Zone#Part Void#Zones VoidVol Void#Part VoidDensContrast VoidProb\n");

  *fpBarycenter = fopen((outputDir+"/"+prefix+"macrocenters_"+dataPortion+"_"+sampleName).c_str(), "w");

  *fpCenters = fopen((outputDir+"/"+prefix+"centers_"+dataPortion+"_"+sampleName).c_str(), "w");
  fprintf(*fpCenters, "# center x,y,z (Mpc/h), volume (normalized), radius (Mpc/h), redshift, volume (Mpc/h^3), void ID, density contrast, num part, parent ID, tree level, number of children, central density\n");

  *fpDistances = fopen((outputDir+"/"+prefix+"boundaryDistances_"+dataPortion+"_"+sampleName).c_str(), "w");

  *fpSkyPositions = fopen((outputDir+"/"+prefix+"sky_positions_"+dataPortion+"_"+sampleName).c_str(), "w");
  fprintf(*fpSkyPositions, "# RA, dec, redshift, radius (Mpc/h), void ID\n");

  *fpShapes = fopen((outputDir+"/"+prefix+"shapes_"+dataPortion+"_"+sampleName).c_str(), "w");
  fprintf(*fpShapes, "# void ID, ellip, eig(1), eig(2), eig(3), eigv(1)-x, eiv(1)-y, eigv(1)-z, eigv(2)-x, eigv(2)-y, eigv(2)-z, eigv(3)-x, eigv(3)-y, eigv(3)-z\n");
  
} // end openFiles


// ----------------------------------------------------------------------------
void closeFiles(FILE* fpZobov, FILE* fpCenters, 
                FILE* fpBarycenter, FILE* fpDistances, FILE* fpShapes, 
                FILE* fpSkyPositions) {

  fclose(fpZobov);
  fclose(fpCenters);
  fclose(fpBarycenter);
  fclose(fpDistances);
  fclose(fpShapes);
  fclose(fpSkyPositions);

} // end closeFile

// ----------------------------------------------------------------------------
void outputVoids(string outputDir, string sampleName, string prefix, 
                string dataPortion, int mockIndex,
                vector<VOID> voids, 
                bool isObservation, double *boxLen, bool doTrim,
                bool doCentralDenCut) {

    int iVoid;
    VOID outVoid;
    FILE *fp, *fpZobov, *fpCenters, *fpCentersNoCut, *fpBarycenter,
         *fpDistances, *fpShapes, *fpSkyPositions;


    openFiles(outputDir, sampleName, prefix, dataPortion,
              mockIndex, voids.size(),
              &fpZobov, &fpCenters, &fpBarycenter,
              &fpDistances, &fpShapes, &fpSkyPositions);


    for (iVoid = 0; iVoid < voids.size(); iVoid++) {
      outVoid = voids[iVoid];


      if (dataPortion == "central" && outVoid.voidType == EDGE_VOID) {
        continue;
      } 

      if (doTrim && outVoid.isLeaf) {
         continue;
       } 

      if (doCentralDenCut && outVoid.hasHighCentralDen) {
         continue;
       } 

     double outCenter[3];
     outCenter[0] = outVoid.macrocenter[0];
     outCenter[1] = outVoid.macrocenter[1];
     outCenter[2] = outVoid.macrocenter[2];

     //if (isObservation) {
     //  outCenter[0] = (outVoid.macrocenter[0]-boxLen[0]/2.)*100.;
     //  outCenter[1] = (outVoid.macrocenter[1]-boxLen[1]/2.)*100.;
     //  outCenter[2] = (outVoid.macrocenter[2]-boxLen[2]/2.)*100.;
     //}

     fprintf(fpZobov, "%d %d %d %f %f %d %d %f %d %f %f\n", 
             iVoid, 
             outVoid.voidID, 
             outVoid.coreParticle, 
             outVoid.coreDens, 
             outVoid.zoneVol, 
             outVoid.zoneNumPart, 
             outVoid.numZones, 
             outVoid.vol, 
             outVoid.numPart, 
             outVoid.densCon, 
             outVoid.voidProb);

     fprintf(fpBarycenter, "%d  %e %e %e\n", 
             outVoid.voidID,
             outVoid.macrocenter[0],
             outVoid.macrocenter[1],
             outVoid.macrocenter[2]);

     fprintf(fpDistances, "%d %e %e %e %e %e\n", 
             outVoid.voidID, 
             outVoid.nearestMock,
             outVoid.radius,
             outVoid.rescaledCoreDens,
             outVoid.nearestMockFromCore,
             outVoid.nearestGalFromCore);

       fprintf(fpCenters, "%.2f %.2f %.2f %.2f %.2f %.5f %.2f %d %f %d %d %d %d %.2f\n",
             outCenter[0],
             outCenter[1],
             outCenter[2],
             outVoid.vol,
             outVoid.radius,
             outVoid.redshift, 
             4./3.*M_PI*pow(outVoid.radius, 3),
             outVoid.voidID,
             outVoid.densCon,
             outVoid.numPart,
             outVoid.parentID,
             outVoid.level,
             outVoid.numChildren,
             outVoid.centralDen);

     double phi = atan2(outVoid.macrocenter[1]-boxLen[1]/2.,
                        outVoid.macrocenter[0]-boxLen[0]/2.);
     if (phi < 0) phi += 2.*M_PI;
     double RA = phi * 180./M_PI;
           
     double theta = acos((outVoid.macrocenter[2]-boxLen[2]/2.) / 
                          outVoid.redshiftInMpc);
     double dec = (M_PI/2. - theta) * 180./M_PI;
             
     fprintf(fpSkyPositions, "%.2f %.2f %.5f %.2f %d\n",
             RA,                  
             dec,
             outVoid.redshift,
             outVoid.radius,
             outVoid.voidID);

     fprintf(fpShapes, "%d %.6f %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e\n", 
             outVoid.voidID,
             outVoid.ellip,
             gsl_vector_get(outVoid.eval, 0),
             gsl_vector_get(outVoid.eval, 1),
             gsl_vector_get(outVoid.eval, 2),
             gsl_matrix_get(outVoid.evec, 0 ,0),
             gsl_matrix_get(outVoid.evec, 0 ,1),
             gsl_matrix_get(outVoid.evec, 0 ,2),
             gsl_matrix_get(outVoid.evec, 1 ,0),
             gsl_matrix_get(outVoid.evec, 1 ,1),
             gsl_matrix_get(outVoid.evec, 1 ,2),
             gsl_matrix_get(outVoid.evec, 2 ,0),
             gsl_matrix_get(outVoid.evec, 2 ,1),
             gsl_matrix_get(outVoid.evec, 2 ,2)
            );

  } // end iVoid

    closeFiles(fpZobov, fpCenters, fpBarycenter,
              fpDistances, fpShapes, fpSkyPositions);

} // end outputVoids
