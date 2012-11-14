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
#include "string.h"
#include "ctype.h"
#include "stdlib.h"
#include <math.h>
#include <stdio.h>
#include <netcdfcpp.h>
#include "pruneVoids_conf.h"

#define LIGHT_SPEED 299792.458
#define MPC2Z 100./LIGHT_SPEED
#define Z2MPC LIGHT_SPEED/100.

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
  float center[3], barycenter[3];
  int accepted;
} VOID;

int main(int argc, char **argv) {

  // initialize arguments
  pruneVoids_info args_info;
  pruneVoids_conf_params args_params;

  pruneVoids_conf_init(&args_info);
  pruneVoids_conf_params_init(&args_params);

  args_params.check_required = 0;
  if (pruneVoids_conf_ext (argc, argv, &args_info, &args_params))
    return 1;

  if (!args_info.configFile_given) {
    if (pruneVoids_conf_required (&args_info,
                                           PRUNEVOIDS_CONF_PACKAGE))
      return 1;
  } else {
    args_params.check_required = 1;
    args_params.initialize = 0;
    if (pruneVoids_conf_config_file (args_info.configFile_arg,
                                              &args_info,
                                              &args_params))
    return 1;
  }

  int i, p, p2, numPartTot, numZonesTot, dummy, iVoid, iZ;
  int numVoids, mockIndex, numKept;
  double tolerance;
  FILE *fp, *fpBarycenter, *fpDistances, *fpSkyPositions, *fpInfo;
  PART *part, *voidPart;
  ZONE2PART *zones2Parts;
  VOID2ZONE *void2Zones;
  VOID *voids;
  float *temp, junk, voidVol;
  int junkInt, voidID, numPart, numZones, zoneID, partID, maxNumPart;
  int coreParticle, zoneNumPart;
  float coreDens, zoneVol, densCon, voidProb, dist[3], dist2, minDist, maxDist;
  float centralRad, centralDen;
  double nearestEdge, redshift;
  char line[500], junkStr[10];
  int mask_index;
  double ranges[2][3], boxLen[3], mul; 
  double volNorm, radius;
  int clock1, clock2;
  int periodicX=0, periodicY=0, periodicZ=0;

  numVoids = args_info.numVoids_arg;
  mockIndex = args_info.mockIndex_arg;
  tolerance = args_info.tolerance_arg;

  clock1 = clock();
  printf("Pruning parameters: %f %f %f %s\n", args_info.zMin_arg, 
                                             args_info.zMax_arg,
                                             args_info.rMin_arg,
                                             args_info.periodic_arg);

  // check for periodic box
  if (!args_info.isObservation_flag) {
    if ( strchr(args_info.periodic_arg, 'x') != NULL) {
      periodicX = 1;
      printf("Will assume x-direction is periodic.\n");
    }
    if ( strchr(args_info.periodic_arg, 'y') != NULL) {
      periodicY = 1;
      printf("Will assume y-direction is periodic.\n");
    }
    if ( strchr(args_info.periodic_arg, 'z') != NULL) {
      periodicZ = 1;
      printf("Will assume z-direction is periodic.\n");
    }
  }

  // load box size
  printf("\n Getting info...\n");
  NcFile f_info(args_info.extraInfo_arg);
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
  fp = fopen(args_info.partFile_arg, "r");
  fread(&dummy, 1, 4, fp); 
  fread(&numPartTot, 1, 4, fp);
  fread(&dummy, 1, 4, fp);

  part = (PART *) malloc(numPartTot * sizeof(PART));
  temp = (float *) malloc(numPartTot * sizeof(float));
 
  volNorm = numPartTot/(boxLen[0]*boxLen[1]*boxLen[2]);
  printf("VOL NORM = %f\n", volNorm);

  printf("CENTRAL DEN = %f\n", args_info.maxCentralDen_arg);

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

  if (!args_info.isObservation_flag) {
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
  fp = fopen(args_info.voidDesc_arg ,"r");
  fgets(line, sizeof(line), fp);
  sscanf(line, "%d %s %d %s", &junkInt, junkStr, &junkInt, junkStr);
  fgets(line, sizeof(line), fp);

  voids = (VOID *) malloc(numVoids * sizeof(VOID));
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
  }
  fclose(fp);

  // load up the zone membership for each void
  printf(" Loading void-zone membership info...\n");
  fp = fopen(args_info.void2Zone_arg, "r");
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
  fp = fopen(args_info.zone2Part_arg, "r");
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
  fp = fopen(args_info.partVol_arg, "r");
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

      if (periodicX) dist[0] = fmin(dist[0], abs(boxLen[0]-dist[0]));
      if (periodicY) dist[1] = fmin(dist[1], abs(boxLen[1]-dist[1]));
      if (periodicZ) dist[2] = fmin(dist[2], abs(boxLen[2]-dist[2]));

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

    // compute central density
    centralRad = voids[iVoid].radius/args_info.centralRadFrac_arg;
    centralRad *= centralRad;
    centralDen = 0.;
    for (p = 0; p < voids[iVoid].numPart; p++) {
      dist[0] = voidPart[p].x - voids[iVoid].barycenter[0];
      dist[1] = voidPart[p].y - voids[iVoid].barycenter[1];
      dist[2] = voidPart[p].z - voids[iVoid].barycenter[2];

      if (periodicX) dist[0] = fmin(dist[0], abs(boxLen[0]-dist[0]));
      if (periodicY) dist[1] = fmin(dist[1], abs(boxLen[1]-dist[1]));
      if (periodicZ) dist[2] = fmin(dist[2], abs(boxLen[2]-dist[2]));

      dist2 = pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2);
      if (dist2 < centralRad) centralDen += 1;
    }
    voids[iVoid].centralDen = centralDen / (4./3. * M_PI * pow(centralRad, 3./2.));

    // compute maximum extent
    //if (args_info.isObservation_flag) {
   //   maxDist = 0.;
   //   for (p = 0; p < voids[iVoid].numPart; p++) {
   //   for (p2 = p; p2 < voids[iVoid].numPart; p2++) {
  
//        dist[0] = voidPart[p].x - voidPart[p2].x;
//        dist[1] = voidPart[p].y - voidPart[p2].y;
//        dist[2] = voidPart[p].z - voidPart[p2].z;

//        dist2 = pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2);
//        if (dist2 > maxDist) maxDist = dist2;
//      }
//      }
//      voids[iVoid].maxRadius = sqrt(maxDist)/2.;
//    } else {
     maxDist = 0.;
      for (p = 0; p < voids[iVoid].numPart; p++) {
  
        dist[0] = voidPart[p].x - voids[iVoid].barycenter[0];
        dist[0] = voidPart[p].y - voids[iVoid].barycenter[1];
        dist[0] = voidPart[p].z - voids[iVoid].barycenter[2];

        if (periodicX) dist[0] = fmin(dist[0], abs(boxLen[0]-dist[0]));
        if (periodicY) dist[1] = fmin(dist[1], abs(boxLen[1]-dist[1]));
        if (periodicZ) dist[2] = fmin(dist[2], abs(boxLen[2]-dist[2]));

        dist2 = pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2);
        if (dist2 > maxDist) maxDist = dist2;
      }
      voids[iVoid].maxRadius = sqrt(maxDist);
//    }
    
    if (args_info.isObservation_flag) {
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
 
    if (args_info.isObservation_flag) {
      voids[iVoid].redshiftInMpc = 
                        sqrt(pow(voids[iVoid].barycenter[0] - boxLen[0]/2.,2) + 
                             pow(voids[iVoid].barycenter[1] - boxLen[1]/2.,2) + 
                             pow(voids[iVoid].barycenter[2] - boxLen[2]/2.,2));
      voids[iVoid].redshiftInMpc = voids[iVoid].redshiftInMpc;
      redshift = voids[iVoid].redshiftInMpc;
      nearestEdge = fmin(fabs(redshift-args_info.zMin_arg*LIGHT_SPEED/100.), 
                         fabs(redshift-args_info.zMax_arg*LIGHT_SPEED/100.)); 
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

    if (nearestEdge < voids[iVoid].nearestMock) {
      voids[iVoid].nearestMock = nearestEdge;
    }
  } // iVoid

  printf(" Picking winners and losers...\n");
  for (iVoid = 0; iVoid < numVoids; iVoid++) {
    voids[iVoid].accepted = 1;
  }

  for (iVoid = 0; iVoid < numVoids; iVoid++) {
    if (strcmp(args_info.dataPortion_arg, "edge")  == 0 &&
        tolerance*voids[iVoid].maxRadius < voids[iVoid].nearestMock) {
      voids[iVoid].accepted = 0;
    }

    if (strcmp(args_info.dataPortion_arg, "central") == 0 && 
        tolerance*voids[iVoid].maxRadius > voids[iVoid].nearestMock) {
      voids[iVoid].accepted = 0;
    }

    if (voids[iVoid].radius < args_info.rMin_arg) {
      voids[iVoid].accepted = 0;
    }

    if (voids[iVoid].centralDen > args_info.maxCentralDen_arg) {
      voids[iVoid].accepted = -1;
    }

  }

  numKept = 0;
  for (iVoid = 0; iVoid < numVoids; iVoid++) {
    if (voids[iVoid].accepted == 1) numKept++;
  }

  printf(" Number kept: %d (out of %d)\n", numKept, numVoids);
  
  printf(" Output...\n");
  fp = fopen(args_info.output_arg, "w");
  fpBarycenter = fopen(args_info.outCenters_arg, "w");
  fpInfo = fopen(args_info.outInfo_arg, "w");
  fpDistances = fopen(args_info.outDistances_arg, "w");
  fpSkyPositions = fopen(args_info.outSkyPositions_arg, "w");
  fprintf(fp, "%d particles, %d voids.\n", mockIndex, numKept);
  fprintf(fp, "see column in master void file\n");
  fprintf(fpInfo, "# center x,y,z (Mpc/h), volume (normalized), radius (Mpc/h), redshift, volume (Mpc/h^3), void ID\n");
  fprintf(fpSkyPositions, "# RA, dec, redshift, radius (Mpc/h), void ID\n");
  for (iVoid = 0; iVoid < numVoids; iVoid++) {

    if (voids[iVoid].accepted != 1) continue;

     fprintf(fp, "%d %d %d %f %f %d %d %f %d %f %f\n", 
             iVoid, 
             voids[iVoid].voidID, 
             voids[iVoid].coreParticle, 
             voids[iVoid].coreDens, 
             voids[iVoid].zoneVol, 
             voids[iVoid].zoneNumPart, 
             voids[iVoid].numZones, 
             voids[iVoid].vol, 
             voids[iVoid].numPart, 
             voids[iVoid].densCon, 
             voids[iVoid].voidProb);

     fprintf(fpBarycenter, "%d %e %e %e\n", 
             voids[iVoid].voidID,
             voids[iVoid].barycenter[0],
             voids[iVoid].barycenter[1],
             voids[iVoid].barycenter[2]);

     fprintf(fpDistances, "%d %e\n", 
             voids[iVoid].voidID, 
             voids[iVoid].nearestMock);

     double outCenter[3];
     outCenter[0] = voids[iVoid].barycenter[0];
     outCenter[1] = voids[iVoid].barycenter[1];
     outCenter[2] = voids[iVoid].barycenter[2];

     if (args_info.isObservation_flag) {
       outCenter[0] = (voids[iVoid].barycenter[0]-boxLen[0]/2.)*100.;
       outCenter[1] = (voids[iVoid].barycenter[1]-boxLen[1]/2.)*100.;
       outCenter[2] = (voids[iVoid].barycenter[2]-boxLen[2]/2.)*100.;
     }

     fprintf(fpInfo, "%.2f %.2f %.2f %.2f %.2f %.5f %.2f %d\n",
             outCenter[0],
             outCenter[1],
             outCenter[2],
             voids[iVoid].vol,
             voids[iVoid].radius,
             voids[iVoid].redshift, 
             4./3.*M_PI*pow(voids[iVoid].radius, 3),
             voids[iVoid].voidID
             );

     fprintf(fpSkyPositions, "%.2f %.2f %.5f %.2f %d\n",
             atan((voids[iVoid].barycenter[1]-boxLen[1]/2.)/(voids[iVoid].barycenter[0]-boxLen[0]/2.)) * 180/M_PI + 180,  
             asin((voids[iVoid].barycenter[2]-boxLen[2]/2.)/voids[iVoid].redshiftInMpc) * 180/M_PI,  
             voids[iVoid].redshift,
             voids[iVoid].radius,
             voids[iVoid].voidID);

  }
  fclose(fp);
  fclose(fpInfo);
  fclose(fpBarycenter);
  fclose(fpDistances);

  // print the centers catalog again but without central density cuts
  fpInfo = fopen(args_info.outNoCutInfo_arg, "w");
  fprintf(fpInfo, "# center x,y,z (km/s), volume (normalized), radius (Mpc/h), redshift, volume (Mpc/h^3), void ID\n");
  for (iVoid = 0; iVoid < numVoids; iVoid++) {

    if (voids[iVoid].accepted == 0) continue;

     double outCenter[3];
     outCenter[0] = voids[iVoid].barycenter[0];
     outCenter[1] = voids[iVoid].barycenter[1];
     outCenter[2] = voids[iVoid].barycenter[2];

     if (args_info.isObservation_flag) {
       outCenter[0] = (voids[iVoid].barycenter[0]-boxLen[0]/2.)*100.;
       outCenter[1] = (voids[iVoid].barycenter[1]-boxLen[1]/2.)*100.;
       outCenter[2] = (voids[iVoid].barycenter[2]-boxLen[2]/2.)*100.;
     }


     fprintf(fpInfo, "%.2f %.2f %.2f %.2f %.2f %.5f %.2f %d\n",
             outCenter[0],
             outCenter[1],
             outCenter[2],
             voids[iVoid].vol,
             voids[iVoid].radius,
             voids[iVoid].redshift, 
             4./3.*M_PI*pow(voids[iVoid].radius, 3),
             voids[iVoid].voidID);
  }
  fclose(fpInfo);

  clock2 = clock();
  printf(" Time: %f sec (for %d voids)\n", (1.*clock2-clock1)/CLOCKS_PER_SEC, numVoids);
  clock1 = clock();


  printf("Done!\n");
  return 0;
} // end main
