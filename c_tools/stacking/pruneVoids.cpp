/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/stacking/pruneVoids.cpp
    Copyright (C) 2010-2013 Guilhem Lavaux
    Copyright (C) 2011-2013 Paul M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

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
#include <netcdfcpp.h>
#include "pruneVoids_conf.h"
#include <vector>

#define LIGHT_SPEED 299792.458
#define MPC2Z 100./LIGHT_SPEED
#define Z2MPC LIGHT_SPEED/100.

#define CENTRAL_VOID 1
#define EDGE_VOID 2

typedef struct partStruct {
  float x, y, z, vol;
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
  int voidID, numPart, numZones, coreParticle, zoneNumPart;
  float maxRadius, nearestMock, centralDen, redshift, redshiftInMpc;
  float nearestEdge;
  float center[3], barycenter[3];
  int accepted;
  int voidType;
  gsl_vector *eval;
  gsl_matrix *evec;
} VOID;

void outputVoid(int iVoid, VOID outVoid, FILE* fpZobov, FILE* fpCenters, 
                FILE* fpCenterNoCut, 
                FILE* fpSkyPositions, FILE* fpBarycenters, FILE* fpDistances, 
                FILE* fpShapes, bool isObservation, double *boxLen);

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

  int i, p, p2, numPartTot, numZonesTot, dummy, iVoid, iZ;
  int numVoids, mockIndex, numKept;
  double tolerance;
  FILE *fp, *fpZobovCentral, *fpZobovAll, *fpCentersCentral, *fpCentersAll,
       *fpCentersNoCutCentral, *fpCentersNoCutAll, *fpBarycenterCentral,
       *fpBarycenterAll, *fpDistancesCentral, *fpDistancesAll,
       *fpShapesCentral, *fpShapesAll, *fpSkyPositionsCentral,
       *fpSkyPositionsAll;
  PART *part, *voidPart;
  ZONE2PART *zones2Parts;
  VOID2ZONE *void2Zones;
  std::vector<VOID> voids;
  float *temp, junk, voidVol;
  int junkInt, voidID, numPart, numZones, zoneID, partID, maxNumPart;
  int coreParticle, zoneNumPart;
  float coreDens, zoneVol, densCon, voidProb, dist[3], dist2, minDist, maxDist;
  float centralRad, centralDen;
  double nearestEdge, redshift;
  char line[500], junkStr[10];
  int mask_index;
  double ranges[3][2], boxLen[3], mul; 
  double volNorm, radius;
  int clock1, clock2;
  int periodicX=0, periodicY=0, periodicZ=0;

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
  NcFile f_info(args.extraInfo_arg);
  ranges[0][0] = f_info.get_att("range_x_min")->as_double(0);
  ranges[0][1] = f_info.get_att("range_x_max")->as_double(0);
  ranges[1][0] = f_info.get_att("range_y_min")->as_double(0);
  ranges[1][1] = f_info.get_att("range_y_max")->as_double(0);
  ranges[2][0] = f_info.get_att("range_z_min")->as_double(0);
  ranges[2][1] = f_info.get_att("range_z_max")->as_double(0);

  boxLen[0] = ranges[0][1] - ranges[0][0];
  boxLen[1] = ranges[1][1] - ranges[1][0];
  boxLen[2] = ranges[2][1] - ranges[2][0];
 
  // read in all particle positions
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

  printf(" Read %d particles...\n", numPartTot);

  if (mockIndex == -1) mockIndex = numPartTot;
  
  // read in desired voids
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
  
    voids[i-1].eval = gsl_vector_alloc(3);
    voids[i-1].evec = gsl_matrix_alloc(3,3);
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
  free(temp);

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

        i++;
      }
    }

    // compute barycenters
    double weight = 0.;
    voids[iVoid].barycenter[0] = 0.;
    voids[iVoid].barycenter[1] = 0.;
    voids[iVoid].barycenter[2] = 0.;

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

      voids[iVoid].barycenter[0] += voidPart[p].vol*(dist[0]);
      voids[iVoid].barycenter[1] += voidPart[p].vol*(dist[1]);
      voids[iVoid].barycenter[2] += voidPart[p].vol*(dist[2]);
      weight += voidPart[p].vol;
    }
    voids[iVoid].barycenter[0] /= weight;
    voids[iVoid].barycenter[1] /= weight;
    voids[iVoid].barycenter[2] /= weight;
    voids[iVoid].barycenter[0] += voids[iVoid].center[0];
    voids[iVoid].barycenter[1] += voids[iVoid].center[1];
    voids[iVoid].barycenter[2] += voids[iVoid].center[2];

    if (periodicX) {
      if (voids[iVoid].barycenter[0] > boxLen[0])
        voids[iVoid].barycenter[0] = voids[iVoid].barycenter[0] - boxLen[0];
      if (voids[iVoid].barycenter[0] < 0)
        voids[iVoid].barycenter[0] = boxLen[0] + voids[iVoid].barycenter[0];
    }
    if (periodicY) {
      if (voids[iVoid].barycenter[1] > boxLen[1])
        voids[iVoid].barycenter[1] = voids[iVoid].barycenter[1] - boxLen[1];
      if (voids[iVoid].barycenter[1] < 0)
        voids[iVoid].barycenter[1] = boxLen[1] + voids[iVoid].barycenter[1];
    }
    if (periodicZ) {
      if (voids[iVoid].barycenter[2] > boxLen[2])
        voids[iVoid].barycenter[2] = voids[iVoid].barycenter[2] - boxLen[2];
      if (voids[iVoid].barycenter[2] < 0)
        voids[iVoid].barycenter[2] = boxLen[2] + voids[iVoid].barycenter[2];
    }

    // compute central density
    centralRad = voids[iVoid].radius/args.centralRadFrac_arg;
    centralDen = 0.;
    int numCentral = 0;
    for (p = 0; p < voids[iVoid].numPart; p++) {
      dist[0] = fabs(voidPart[p].x - voids[iVoid].barycenter[0]);
      dist[1] = fabs(voidPart[p].y - voids[iVoid].barycenter[1]);
      dist[2] = fabs(voidPart[p].z - voids[iVoid].barycenter[2]);

      if (periodicX) dist[0] = fmin(dist[0], boxLen[0]-dist[0]);
      if (periodicY) dist[1] = fmin(dist[1], boxLen[1]-dist[1]);
      if (periodicZ) dist[2] = fmin(dist[2], boxLen[2]-dist[2]);

      dist2 = pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2);
      if (sqrt(dist2) < centralRad) numCentral += 1;
    }
    voids[iVoid].centralDen = numCentral / (volNorm*4./3. * M_PI * 
                                            pow(centralRad, 3.));

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
      maxDist = 0.;
      for (p = 0; p < voids[iVoid].numPart; p++) {
  
        dist[0] = fabs(voidPart[p].x - voids[iVoid].barycenter[0]);
        dist[1] = fabs(voidPart[p].y - voids[iVoid].barycenter[1]);
        dist[2] = fabs(voidPart[p].z - voids[iVoid].barycenter[2]);

        if (periodicX) dist[0] = fmin(dist[0], boxLen[0]-dist[0]);
        if (periodicY) dist[1] = fmin(dist[1], boxLen[1]-dist[1]);
        if (periodicZ) dist[2] = fmin(dist[2], boxLen[2]-dist[2]);

        dist2 = pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2);
        if (dist2 > maxDist) maxDist = dist2;
      }
      voids[iVoid].maxRadius = sqrt(maxDist);
//    }
    
    if (args.isObservation_flag) {
      // compute distance from center to nearest mock
      minDist = 1.e99;
      for (p = mockIndex; p < numPartTot; p++) {
        dist[0] = voids[iVoid].barycenter[0] - part[p].x;
        dist[1] = voids[iVoid].barycenter[1] - part[p].y;
        dist[2] = voids[iVoid].barycenter[2] - part[p].z;

        dist2 = pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2);
        if (dist2 < minDist) minDist = dist2;
      }
      voids[iVoid].nearestMock = sqrt(minDist);
    
    } else {
      voids[iVoid].nearestMock = 1.e99;
    }
 
    if (args.isObservation_flag) {
      voids[iVoid].redshiftInMpc = 
                        sqrt(pow(voids[iVoid].barycenter[0] - boxLen[0]/2.,2) + 
                             pow(voids[iVoid].barycenter[1] - boxLen[1]/2.,2) + 
                             pow(voids[iVoid].barycenter[2] - boxLen[2]/2.,2));
      voids[iVoid].redshiftInMpc = voids[iVoid].redshiftInMpc;
      redshift = voids[iVoid].redshiftInMpc;
      nearestEdge = fabs(redshift-args.zMax_arg*LIGHT_SPEED/100.); 
      //nearestEdge = fmin(fabs(redshift-args.zMin_arg*LIGHT_SPEED/100.), 
      //                   fabs(redshift-args.zMax_arg*LIGHT_SPEED/100.)); 
      voids[iVoid].redshift = voids[iVoid].redshiftInMpc/LIGHT_SPEED*100.;

    } else {

      voids[iVoid].redshiftInMpc = voids[iVoid].barycenter[2];
      voids[iVoid].redshift = voids[iVoid].barycenter[2]/LIGHT_SPEED*100.;

      nearestEdge = 1.e99;
     
      if (!periodicX) {
        nearestEdge = fmin(nearestEdge,
                           fabs(voids[iVoid].barycenter[0] - ranges[0][0]));
        nearestEdge = fmin(nearestEdge,
                           fabs(voids[iVoid].barycenter[0] - ranges[0][1]));
      }
      if (!periodicY) {
        nearestEdge = fmin(nearestEdge,
                           fabs(voids[iVoid].barycenter[1] - ranges[1][0]));
        nearestEdge = fmin(nearestEdge,
                           fabs(voids[iVoid].barycenter[1] - ranges[1][1]));
      }
      if (!periodicZ) {
        nearestEdge = fmin(nearestEdge,
                           fabs(voids[iVoid].barycenter[2] - ranges[2][0]));
        nearestEdge = fmin(nearestEdge,
                           fabs(voids[iVoid].barycenter[2] - ranges[2][1]));
      }
    }

    voids[iVoid].nearestEdge = nearestEdge;

    // compute eigenvalues and vectors for orientation and shape
    double inertia[9];
    for (int p = 0; p < voids[iVoid].numPart; p++) {
      dist[0] = voidPart[p].x - voids[iVoid].barycenter[0];
      dist[1] = voidPart[p].y - voids[iVoid].barycenter[1];
      dist[2] = voidPart[p].z - voids[iVoid].barycenter[2];

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
  } // iVoid

  gsl_eigen_symmv_free(eigw);
    
  int numWrong = 0;
  int numHighDen = 0;
  int numCentral = 0;
  int numEdge = 0;
  int numNearZ = 0;
  int numTooSmall = 0;

  printf(" Picking winners and losers...\n");
  printf("  Starting with %d voids\n", voids.size());

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
    if (voids[iVoid].densCon > 1.e4) {
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
    if (tolerance*voids[iVoid].maxRadius > voids[iVoid].nearestEdge) {
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
  printf("  4th filter: rejected %d too close to low redshift boundaries\n", numNearZ); 

  for (iVoid = 0; iVoid < voids.size(); iVoid++) {
    if (voids[iVoid].centralDen > args.maxCentralDen_arg) {
      voids[iVoid].accepted = -1;
      numHighDen++;
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

  printf("  Number kept: %d (out of %d)\n", voids.size(), numVoids);
  printf("   We have %d edge voids\n", numEdge); 
  printf("   We have %d central voids\n", numCentral); 
  printf("   We have %d too high central density\n", numHighDen); 
  
  printf(" Output...\n");
  fpZobovCentral = fopen((std::string(args.outputDir_arg)+"/voidDesc_central_"+std::string(args.sampleName_arg)+".out").c_str(), "w");
  fprintf(fpZobovCentral, "%d particles, %d voids.\n", mockIndex, numKept);
  fprintf(fpZobovCentral, "Void# FileVoid# CoreParticle CoreDens ZoneVol Zone#Part Void#Zones VoidVol Void#Part VoidDensContrast VoidProb\n");

  fpZobovAll = fopen((std::string(args.outputDir_arg)+"/voidDesc_all_"+std::string(args.sampleName_arg)+".out").c_str(), "w");
  fprintf(fpZobovAll, "%d particles, %d voids.\n", mockIndex, numKept);
  fprintf(fpZobovAll, "Void# FileVoid# CoreParticle CoreDens ZoneVol Zone#Part Void#Zones VoidVol Void#Part VoidDensContrast VoidProb\n");

  fpBarycenterCentral = fopen((std::string(args.outputDir_arg)+"/barycenters_central_"+std::string(args.sampleName_arg)+".out").c_str(), "w");
  fpBarycenterAll = fopen((std::string(args.outputDir_arg)+"/barycenters_all_"+std::string(args.sampleName_arg)+".out").c_str(), "w");

  fpCentersCentral = fopen((std::string(args.outputDir_arg)+"/centers_central_"+std::string(args.sampleName_arg)+".out").c_str(), "w");
  fprintf(fpCentersCentral, "# center x,y,z (Mpc/h), volume (normalized), radius (Mpc/h), redshift, volume (Mpc/h^3), void ID, density contrast, num part\n");

  fpCentersAll = fopen((std::string(args.outputDir_arg)+"/centers_all_"+std::string(args.sampleName_arg)+".out").c_str(), "w");
  fprintf(fpCentersAll, "# center x,y,z (Mpc/h), volume (normalized), radius (Mpc/h), redshift, volume (Mpc/h^3), void ID, density contrast, num part\n");

  fpCentersNoCutCentral = fopen((std::string(args.outputDir_arg)+"/centers_nocut_central_"+std::string(args.sampleName_arg)+".out").c_str(), "w");
  fprintf(fpCentersNoCutCentral, "# center x,y,z (Mpc/h), volume (normalized), radius (Mpc/h), redshift, volume (Mpc/h^3), void ID, density contrast, num part\n");

  fpCentersNoCutAll = fopen((std::string(args.outputDir_arg)+"/centers_nocut_all_"+std::string(args.sampleName_arg)+".out").c_str(), "w");
  fprintf(fpCentersNoCutAll, "# center x,y,z (Mpc/h), volume (normalized), radius (Mpc/h), redshift, volume (Mpc/h^3), void ID, density contrast, num part\n");


  fpDistancesCentral = fopen((std::string(args.outputDir_arg)+"boundaryDistancesCentral_"+std::string(args.sampleName_arg)+".out").c_str(), "w");
  fpDistancesAll = fopen((std::string(args.outputDir_arg)+"boundaryDistancesAll_"+std::string(args.sampleName_arg)+".out").c_str(), "w");

  fpSkyPositionsCentral = fopen((std::string(args.outputDir_arg)+"/sky_positions_central_"+std::string(args.sampleName_arg)+".out").c_str(), "w");
  fprintf(fpSkyPositionsCentral, "# RA, dec, redshift, radius (Mpc/h), void ID\n");

  fpSkyPositionsAll = fopen((std::string(args.outputDir_arg)+"/sky_positions_all_"+std::string(args.sampleName_arg)+".out").c_str(), "w");
  fprintf(fpSkyPositionsAll, "# RA, dec, redshift, radius (Mpc/h), void ID\n");

  fpShapesCentral = fopen((std::string(args.outputDir_arg)+"/shapes_central_"+std::string(args.sampleName_arg)+".out").c_str(), "w");
  fprintf(fpShapesCentral, "# void ID, eig(1), eig(2), eig(3), eigv(1)-x, eiv(1)-y, eigv(1)-z, eigv(2)-x, eigv(2)-y, eigv(2)-z, eigv(3)-x, eigv(3)-y, eigv(3)-z\n");
  
  fpShapesAll = fopen((std::string(args.outputDir_arg)+"/shapes_all_"+std::string(args.sampleName_arg)+".out").c_str(), "w");
  fprintf(fpShapesAll, "# void ID, eig(1), eig(2), eig(3), eigv(1)-x, eiv(1)-y, eigv(1)-z, eigv(2)-x, eigv(2)-y, eigv(2)-z, eigv(3)-x, eigv(3)-y, eigv(3)-z\n");


  for (iVoid = 0; iVoid < voids.size(); iVoid++) {

    if (voids[iVoid].voidType == CENTRAL_VOID) {
      outputVoid(iVoid, voids[iVoid], fpZobovCentral, fpCentersCentral, 
                 fpCentersNoCutCentral, fpSkyPositionsCentral,
                 fpBarycenterCentral, fpDistancesCentral, fpShapesCentral,
                 args.isObservation_flag, boxLen);
   } 

    if (voids[iVoid].voidType == EDGE_VOID || 
        voids[iVoid].voidType == CENTRAL_VOID) {
      outputVoid(iVoid, voids[iVoid], fpZobovAll, fpCentersAll, 
                 fpCentersNoCutAll, fpSkyPositionsAll,
                 fpBarycenterAll, fpDistancesAll, fpShapesAll,
                 args.isObservation_flag, boxLen);
    }
 }

  fclose(fpZobovCentral);
  fclose(fpZobovAll);
  fclose(fpCentersCentral);
  fclose(fpCentersAll);
  fclose(fpCentersNoCutCentral);
  fclose(fpCentersNoCutAll);
  fclose(fpBarycenterCentral);
  fclose(fpBarycenterAll);
  fclose(fpDistancesCentral);
  fclose(fpDistancesAll);
  fclose(fpShapesCentral);
  fclose(fpShapesAll);
  fclose(fpSkyPositionsCentral);
  fclose(fpSkyPositionsAll);

  clock2 = clock();
  printf(" Time: %f sec (for %d voids)\n", 
         (1.*clock2-clock1)/CLOCKS_PER_SEC, numVoids);
  clock1 = clock();


  printf("Done!\n");
  return 0;
} // end main


// ----------------------------------------------------------------------------
void outputVoid(int iVoid, VOID outVoid, FILE* fpZobov, FILE* fpCenters, 
                FILE* fpCenterNoCut, FILE* fpSkyPositions, 
                FILE* fpBarycenters, FILE* fpDistances, FILE* fpShapes,
                bool isObservation, double *boxLen) {

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

     fprintf(fpBarycenters, "%d  %e %e %e\n", 
             outVoid.voidID,
             outVoid.barycenter[0],
             outVoid.barycenter[1],
             outVoid.barycenter[2]);

     fprintf(fpDistances, "%d %e\n", 
             outVoid.voidID, 
             outVoid.nearestMock);

     double outCenter[3];
     outCenter[0] = outVoid.barycenter[0];
     outCenter[1] = outVoid.barycenter[1];
     outCenter[2] = outVoid.barycenter[2];

     if (isObservation) {
       outCenter[0] = (outVoid.barycenter[0]-boxLen[0]/2.)*100.;
       outCenter[1] = (outVoid.barycenter[1]-boxLen[1]/2.)*100.;
       outCenter[2] = (outVoid.barycenter[2]-boxLen[2]/2.)*100.;
     }

     if (outVoid.accepted == 1) {
       fprintf(fpCenters, "%.2f %.2f %.2f %.2f %.2f %.5f %.2f %d %f %d\n",
             outCenter[0],
             outCenter[1],
             outCenter[2],
             outVoid.vol,
             outVoid.radius,
             outVoid.redshift, 
             4./3.*M_PI*pow(outVoid.radius, 3),
             outVoid.voidID,
             outVoid.densCon,
             outVoid.numPart);
     } 

       fprintf(fpCenterNoCut, "%.2f %.2f %.2f %.2f %.2f %.5f %.2f %d %f %d\n",
             outCenter[0],
             outCenter[1],
             outCenter[2],
             outVoid.vol,
             outVoid.radius,
             outVoid.redshift, 
             4./3.*M_PI*pow(outVoid.radius, 3),
             outVoid.voidID,
             outVoid.densCon,
             outVoid.numPart);

     fprintf(fpSkyPositions, "%.2f %.2f %.5f %.2f %d\n",
             atan((outVoid.barycenter[1]-boxLen[1]/2.) / 
                  (outVoid.barycenter[0]-boxLen[0]/2.)) * 180/M_PI + 180,  
             asin((outVoid.barycenter[2]-boxLen[2]/2.) / 
                   outVoid.redshiftInMpc) * 180/M_PI,  
             outVoid.redshift,
             outVoid.radius,
             outVoid.voidID);

     fprintf(fpShapes, "%d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n", 
             outVoid.voidID,
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

} // end outputVoid
