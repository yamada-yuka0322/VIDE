/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/analysis/voidOverlap.cpp
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



// =============================================================================
//
//   Program: voidOverlap
//
//   Description: Takes two void catalogs and reports the "overlap" between
//                them.
//
// ============================================================================


#include "string.h"
#include "ctype.h"
#include "stdlib.h"
#include <stdio.h>
#include <math.h>
#include "voidOverlap_conf.h"
#include <vector>
#include <netcdfcpp.h>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

typedef struct partStruct {
  float x, y, z, volume;
  float ra, dec, redshift;
  long uniqueID;
} PART;

typedef struct zoneStruct {
  int numPart;
  int *partIDs;
} ZONE2PART;

typedef struct voidZoneStruct {
  int numZones;
  int *zoneIDs;
} VOID2ZONE;

typedef struct matchProps {
  int matchID;
  float commonVolOrig, commonVolProg;
  float dist;
  float merit;
} MATCHPROPS;

typedef struct voidStruct {
  float vol, coreDens, zoneVol, densCon, voidProb, radius;
  int voidID, numPart, numZones, coreParticle, zoneNumPart;
  float maxRadius, nearestMock, centralDen, redshift, redshiftInMpc;
  float nearestEdge;
  float barycenter[3];
  float ellipticity;
  std::vector<MATCHPROPS> matches;
  int numMatches;
  int numBigMatches;
  float radiusMpc;
} VOID;

typedef struct catalog {
  int numVoids, numPartTot, numZonesTot;
  float boxLen[3];
  std::vector<PART> part;
  std::vector<ZONE2PART> zones2Parts;
  std::vector<VOID2ZONE> void2Zones;
  std::vector<VOID> voids;
} CATALOG;

void loadCatalog(const char *partFile, const char *volFile, 
                 const char *voidFile, const char *zoneFile,
                 const char *infoFile, const char *centerFile,
                 const char *shapeFile,
                 const char *zonePartFile, CATALOG& catalog);

float getDist(CATALOG& catalog1, CATALOG& catalog2, int& iVoid1, int& iVoid2,
               bool periodicX, bool periodicY, bool operiodicZ);

void sortMatches(std::vector<MATCHPROPS>& matches, 
                 CATALOG& catalog2, int mode);

// ----------------------------------------------------------------------------

int main(int argc, char **argv) {

  int p1, p2, iZ1, iZ2, iVoid1, iVoid2, iVoid, zoneID1, zoneID2, iMatch;
  int partID1, partID2;
  int voidID1, voidID2;
  bool periodicX=false, periodicY=false, periodicZ=false, match;
  float dist[3], rdist, r1, r2;
  FILE *fp;
  int closestMatchID;
  float closestMatchDist;
  float commonVolRatio;
  MATCHPROPS newMatch;
  int MAX_MATCHES = 100;
  float matchDist = 0.25;
  bool alreadyMatched;

  CATALOG catalog1, catalog2;

  // initialize arguments
  voidOverlap_info args;
  voidOverlap_conf_params params;

  voidOverlap_conf_init(&args);
  voidOverlap_conf_params_init(&params);

  params.check_required = 0;
  if (voidOverlap_conf_ext (argc, argv, &args, &params))
    return 1;

  if (!args.configFile_given) {
    if (voidOverlap_conf_required (&args, VOIDOVERLAP_CONF_PACKAGE))
      return 1;
  } else {
    params.check_required = 1;
    params.initialize = 0;
    if (voidOverlap_conf_config_file (args.configFile_arg, &args, &params))
    return 1;
  }

  loadCatalog(args.partFile1_arg, args.volFile1_arg, args.voidFile1_arg,
              args.zoneFile1_arg, args.infoFile1_arg, args.centerFile1_arg,
              args.shapeFile1_arg,
              args.zonePartFile1_arg, catalog1);

  loadCatalog(args.partFile2_arg, args.volFile2_arg, args.voidFile2_arg,
              args.zoneFile2_arg, args.infoFile2_arg, args.centerFile2_arg,
              args.shapeFile2_arg,
              args.zonePartFile2_arg, catalog2);

  // check for periodic box
    if ( strchr(args.periodic_arg, 'x') != NULL) {
      periodicX = true;
      printf("Will assume x-direction is periodic.\n");
    }
    if ( strchr(args.periodic_arg, 'y') != NULL) {
      periodicY = true;
      printf("Will assume y-direction is periodic.\n");
    }
    if ( strchr(args.periodic_arg, 'z') != NULL) {
      periodicZ = true;
      printf("Will assume z-direction is periodic.\n");
    }

  // find closest voids
  printf(" Finding nearest matches...\n");
  for (iVoid1 = 0; iVoid1 < catalog1.numVoids; iVoid1++) {
    for (iVoid2 = 0; iVoid2 < catalog2.numVoids; iVoid2++) {
      rdist = getDist(catalog1, catalog2, iVoid1, iVoid2, 
                      periodicX, periodicY, periodicZ);

      newMatch.matchID = iVoid2;
      newMatch.commonVolOrig = 0;
      newMatch.commonVolProg = 0;
      newMatch.dist = rdist; 

      //if (rdist > 1.5) continue;
      if (rdist/catalog1.voids[iVoid1].radius > 2.0) continue;

      // see if center is contained in void
      match = false;
      voidID1 = catalog1.voids[iVoid1].voidID;
      for (iZ1 = 0; iZ1 < catalog1.void2Zones[voidID1].numZones; iZ1++) {
        zoneID1 = catalog1.void2Zones[voidID1].zoneIDs[iZ1];
        for (p1 = 0; p1 < catalog1.zones2Parts[zoneID1].numPart; p1++) {
          partID1 = catalog1.zones2Parts[zoneID1].partIDs[p1];     

          dist[0] = fabs(catalog1.part[partID1].x - 
                         catalog2.voids[iVoid2].barycenter[0]);
          dist[1] = fabs(catalog1.part[partID1].y - 
                         catalog2.voids[iVoid2].barycenter[1]);
          dist[2] = fabs(catalog1.part[partID1].z - 
                         catalog2.voids[iVoid2].barycenter[2]);

          if (periodicX) dist[0] = fmin(dist[0], 1.0 - dist[0]);
          if (periodicY) dist[1] = fmin(dist[1], 1.0 - dist[1]);
          if (periodicZ) dist[2] = fmin(dist[2], 1.0 - dist[2]);

          rdist = sqrt(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2]);

          r1 = pow(3./4./M_PI*catalog1.part[partID1].volume / 
                   catalog1.numPartTot, 1./3.);
          if (rdist <= matchDist*4*r1) {
            match = true;
            break;
          }
        }
      }
       
      if (!match) continue; 
      //if (rdist/catalog1.voids[iVoid1].radius > 1.) continue;

      //if (catalog1.voids[iVoid1].matches.size() < MAX_MATCHES) {
        catalog1.voids[iVoid1].matches.push_back(newMatch);
      //} else {
      //  // find the farthest match
      //  float farthestMatchDist = 0;
      //  int   farthestMatchID = 0;
      //  for (iMatch = 0; iMatch < MAX_MATCHES; iMatch++) {
      //    if (catalog1.voids[iVoid1].matches[iMatch].dist > farthestMatchDist){
      //      farthestMatchDist = catalog1.voids[iVoid1].matches[iMatch].dist;
       //     farthestMatchID   = iMatch;
       //   }           
       // } 
       // if (rdist < farthestMatchDist)
       //   catalog1.voids[iVoid1].matches[farthestMatchID] = newMatch;
      //}
    }
  } 

  // pick the closest matches to speed up computation
  for (iVoid1 = 0; iVoid1 < catalog1.numVoids; iVoid1++) {
    int numPotential = catalog1.voids[iVoid1].matches.size();
    catalog1.voids[iVoid1].numMatches = numPotential;
    
    if (numPotential > MAX_MATCHES) {
      sortMatches(catalog1.voids[iVoid1].matches, catalog2, 1);
    }

    catalog1.voids[iVoid1].matches.resize(MAX_MATCHES); 
  }

  printf(" Determining overlap...\n");
  for (iVoid1 = 0; iVoid1 < catalog1.numVoids; iVoid1++) {
    printf("  Working on void %d of %d...\n", iVoid1+1, catalog1.numVoids); 
    voidID1 = catalog1.voids[iVoid1].voidID;

    for (iMatch = 0; iMatch < catalog1.voids[iVoid1].matches.size();iMatch++) {
      iVoid2 = catalog1.voids[iVoid1].matches[iMatch].matchID;
      voidID2 = catalog2.voids[iVoid2].voidID;

      for (iZ1 = 0; iZ1 < catalog1.void2Zones[voidID1].numZones; iZ1++) {
        zoneID1 = catalog1.void2Zones[voidID1].zoneIDs[iZ1];

        for (p1 = 0; p1 < catalog1.zones2Parts[zoneID1].numPart; p1++) {
          partID1 = catalog1.zones2Parts[zoneID1].partIDs[p1];

          alreadyMatched = false;

          for (iZ2 = 0; iZ2 < catalog2.void2Zones[voidID2].numZones; iZ2++) {
            zoneID2 = catalog2.void2Zones[voidID2].zoneIDs[iZ2];

            for (p2 = 0; p2 < catalog2.zones2Parts[zoneID2].numPart; p2++) {
              partID2 = catalog2.zones2Parts[zoneID2].partIDs[p2];     

              match = false;

              if (args.useID_flag) {
                if (catalog1.part[partID1].uniqueID == 
                    catalog2.part[partID2].uniqueID) match = true;
              } else {
                dist[0] = fabs(catalog1.part[partID1].x - 
                               catalog2.part[partID2].x);
                dist[1] = fabs(catalog1.part[partID1].y - 
                               catalog2.part[partID2].y);
                dist[2] = fabs(catalog1.part[partID1].z - 
                               catalog2.part[partID2].z);

                if (periodicX) dist[0] = fmin(dist[0], 1.0 - dist[0]);
                if (periodicY) dist[1] = fmin(dist[1], 1.0 - dist[1]);
                if (periodicZ) dist[2] = fmin(dist[2], 1.0 - dist[2]);

                rdist = sqrt(dist[0]*dist[0] + dist[1]*dist[1] + 
                             dist[2]*dist[2]);

                r1 = pow(3./4./M_PI*catalog1.part[partID1].volume / 
                         catalog1.numPartTot, 1./3.);
                r2 = pow(3./4./M_PI*catalog2.part[partID2].volume / 
                         catalog2.numPartTot, 1./3.);
                if (rdist <= matchDist*(r1+r2)) match = true;
              }
       
              if (match) {
                if (!alreadyMatched) {
                  catalog1.voids[iVoid1].matches[iMatch].commonVolOrig += 
                    catalog1.part[partID1].volume;
                } 
                catalog1.voids[iVoid1].matches[iMatch].commonVolProg += 
                  catalog2.part[partID2].volume;
                alreadyMatched = true;
              } // end if match

            } // end p2
          } // end iZ2
        } // end p1
      } // end iZ1
    } // end iVoid2

    for (iMatch = 0; iMatch < catalog1.voids[iVoid1].matches.size(); iMatch++) {
      int matchID = catalog1.voids[iVoid1].matches[iMatch].matchID;
      catalog1.voids[iVoid1].matches[iMatch].merit = 
        catalog1.voids[iVoid1].matches[iMatch].commonVolOrig * 
        catalog1.voids[iVoid1].matches[iMatch].commonVolOrig / 
        catalog1.voids[iVoid1].vol /
        catalog2.voids[iVoid1].vol;
        //catalog1.voids[iVoid1].matches[iMatch].commonVolOrig * 
        //catalog1.voids[iVoid1].matches[iMatch].commonVolProg / 
        //catalog1.voids[iVoid1].vol /
        //catalog2.voids[matchID].vol;
        //pow(catalog1.voids[iVoid1].matches[iMatch].commonVol,2) / 
        //catalog1.voids[iVoid1].vol /
        //catalog2.voids[matchID].vol;
    }
//    sortMatches(catalog1.voids[iVoid1].matches, catalog2);
//printf("BEST VOL %e\n", catalog2.voids[catalog1.voids[iVoid1].matches[0].matchID].vol/catalog1.voids[iVoid1].vol);
  } // end match finding

  printf(" Sorting matches...\n");
  for (iVoid1 = 0; iVoid1 < catalog1.numVoids; iVoid1++) {
    sortMatches(catalog1.voids[iVoid1].matches, catalog2, 0);
    if (catalog1.voids[iVoid1].matches.size() > 0 && 
        catalog1.voids[iVoid1].matches[0].merit < 1.e-6)
      sortMatches(catalog1.voids[iVoid1].matches, catalog2, 1);
  }

  // count up significant matches
  printf(" Categorizing matches...\n");
  for (iVoid1 = 0; iVoid1 < catalog1.numVoids; iVoid1++) {
    closestMatchDist = 0.;
    for (iMatch = 0; iMatch < catalog1.voids[iVoid1].matches.size(); iMatch++) {
      commonVolRatio = catalog1.voids[iVoid1].matches[iMatch].commonVolOrig / 
                     //  catalog1.voids[iVoid1].numPart;
                     catalog1.voids[iVoid1].vol;
      if (commonVolRatio > 0.2) catalog1.voids[iVoid1].numBigMatches++;
    } 
  }

  // output summary
  printf(" Output...\n");
  std::string filename;
  filename = string(args.outfile_arg);
  filename = filename.append("summary.out");
  fp = fopen(filename.c_str(), "w");
  fprintf(fp, "# void ID, radius, radius ratio, common volume ratio (to original), common volume ratio (to match), relative dist, num matches, num significant matches, match ID, merit, ellipticity ratio, density contrast\n");
  for (iVoid1 = 0; iVoid1 < catalog1.numVoids; iVoid1++) {
    int voidID = catalog1.voids[iVoid1].voidID;
    if (catalog1.voids[iVoid1].numMatches > 0) {
      iVoid2 = catalog1.voids[iVoid1].matches[0].matchID;
      float rRatio = catalog2.voids[iVoid2].radius / 
                     catalog1.voids[iVoid1].radius;
      float ellipRatio = catalog2.voids[iVoid2].ellipticity / 
                         catalog1.voids[iVoid1].ellipticity;
      commonVolRatio = catalog1.voids[iVoid1].matches[0].commonVolOrig / 
                       //catalog1.voids[iVoid1].numPart;
                       catalog1.voids[iVoid1].vol;
      float volRatio = catalog1.voids[iVoid1].matches[0].commonVolProg / 
                       catalog2.voids[iVoid2].vol;
      rdist = catalog1.voids[iVoid1].matches[0].dist;
      rdist /= catalog1.voids[iVoid1].radius;

      fprintf(fp, "%d %.4f %.4f %e %e %.4f %d %d %d %e %e %e\n", voidID, 
      //fprintf(fp, "%d %.4f %.4f %.4f %.4f %.4f %d %d %d\n", voidID, 
                                   catalog1.voids[iVoid1].radiusMpc,
                                   rRatio, 
                                   commonVolRatio,
                                   volRatio,
                                   rdist, 
                                   catalog1.voids[iVoid1].numMatches,
                                   catalog1.voids[iVoid1].numBigMatches,
                                   catalog2.voids[iVoid2].voidID,
                                   catalog1.voids[iVoid1].matches[0].merit,
                                   ellipRatio,
                                   catalog1.voids[iVoid1].densCon);

    } else {
      fprintf(fp, "%d %.4f 0.0 0.0 0.0 0.0 0.0 0 0 0 0.0 0.0\n", voidID, 
                                           catalog1.voids[iVoid1].radiusMpc);
    }
  } // end printing
  fclose(fp);

  // output detail
  printf(" Output detail...\n");
  filename = string(args.outfile_arg);
  filename = filename.append("detail.out");
  fp = fopen(filename.c_str(), "w");
  fprintf(fp, "# match ID, merit, relative dist, relative radius\n");
  for (iVoid1 = 0; iVoid1 < catalog1.numVoids; iVoid1++) {
    int voidID = catalog1.voids[iVoid1].voidID;
    fprintf(fp,"#%d (%.4f Mpc/h):\n", voidID, catalog1.voids[iVoid1].radiusMpc);
    for (iMatch = 0; iMatch < catalog1.voids[iVoid1].matches.size(); iMatch++) {
      commonVolRatio = catalog1.voids[iVoid1].matches[iMatch].commonVolOrig / 
                       //catalog1.voids[iVoid1].numPart;
                       catalog1.voids[iVoid1].vol;

      int matchID = catalog1.voids[iVoid1].matches[iMatch].matchID;
      fprintf(fp, "%d %e %.4f %.4f\n", matchID, 
                                   catalog1.voids[iVoid1].matches[iMatch].merit,
                                   catalog1.voids[iVoid1].matches[iMatch].dist/
                                     catalog1.voids[iVoid1].radius,
                                   catalog2.voids[matchID].radius/
                                     catalog1.voids[iVoid1].radius);
    }

    fprintf(fp, "\n ");
    fprintf(fp, "\n ");
  } // end printing detail
  fclose(fp);


  printf("\nDone!\n");
  return 0;
} // end main

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void loadCatalog(const char *partFile, const char *volFile, 
                 const char *voidFile, const char *zoneFile,
                 const char *infoFile, const char *centerFile,
                 const char *shapeFile,
                 const char *zonePartFile, CATALOG& catalog) {

  int i, p, numPartTot, numZonesTot, dummy, iVoid, iZ, numVolTot;
  FILE *fp;
  float *temp, junk, voidVol, coreParticle, coreDens, zoneVol, zoneNumPart;
  float densCon, voidProb, volNorm;
  long *temp2;
  int junkInt, voidID, numPart, numZones, zoneID, partID;
  char line[500], junkStr[10];
  float ranges[3][2];

  printf("Loading info...\n");
  NcFile f_info(infoFile);
  ranges[0][0] = f_info.get_att("range_x_min")->as_double(0);
  ranges[0][1] = f_info.get_att("range_x_max")->as_double(0);
  ranges[1][0] = f_info.get_att("range_y_min")->as_double(0);
  ranges[1][1] = f_info.get_att("range_y_max")->as_double(0);
  ranges[2][0] = f_info.get_att("range_z_min")->as_double(0);
  ranges[2][1] = f_info.get_att("range_z_max")->as_double(0);

  catalog.boxLen[0] = ranges[0][1] - ranges[0][0];
  catalog.boxLen[1] = ranges[1][1] - ranges[1][0];
  catalog.boxLen[2] = ranges[2][1] - ranges[2][0];
  f_info.close();
  
  // read in all particle positions
  printf("Loading particles...\n");
  fp = fopen(partFile, "r");
  fread(&dummy, 1, 4, fp); 
  fread(&numPartTot, 1, 4, fp);
  fread(&dummy, 1, 4, fp);

  catalog.part.resize(numPartTot);
  catalog.numPartTot = numPartTot;

  volNorm = numPartTot/(catalog.boxLen[0]*catalog.boxLen[1]*catalog.boxLen[2]);

  temp = (float *) malloc(numPartTot * sizeof(float));
  temp2 = (long *) malloc(numPartTot * sizeof(long));

  fread(&dummy, 1, 4, fp);
  fread(temp, numPartTot, 4, fp);
  for (p = 0; p < numPartTot; p++) 
    catalog.part[p].x = temp[p];
  fread(&dummy, 1, 4, fp);
  fread(&dummy, 1, 4, fp);
  fread(temp, numPartTot, 4, fp);
  for (p = 0; p < numPartTot; p++) 
    catalog.part[p].y = temp[p];
  fread(&dummy, 1, 4, fp);
  fread(&dummy, 1, 4, fp);
  fread(temp, numPartTot, 4, fp);
  for (p = 0; p < numPartTot; p++) 
    catalog.part[p].z = temp[p];
  fread(&dummy, 1, 4, fp);
  fread(&dummy, 1, 4, fp);
  fread(temp, numPartTot, 4, fp);
  for (p = 0; p < numPartTot; p++) 
    catalog.part[p].ra = temp[p];
  fread(&dummy, 1, 4, fp);
  fread(&dummy, 1, 4, fp);
  fread(temp, numPartTot, 4, fp);
  for (p = 0; p < numPartTot; p++) 
    catalog.part[p].dec = temp[p];
  fread(&dummy, 1, 4, fp);
  fread(&dummy, 1, 4, fp);
  fread(temp, numPartTot, 4, fp);
  for (p = 0; p < numPartTot; p++) 
    catalog.part[p].redshift = temp[p];
  fread(&dummy, 1, 4, fp);
  fread(&dummy, 1, 4, fp);
  fread(temp2, numPartTot, 8, fp);
  for (p = 0; p < numPartTot; p++) 
    catalog.part[p].uniqueID = temp2[p];
 
  free(temp2);
  fclose(fp); 

  printf(" Read %d particles...\n", catalog.numPartTot);

  // read in all particle volumes
  printf(" Loading volumes...\n");
  fp = fopen(volFile, "r");
  fread(&numVolTot, 1, 4, fp);
  fread(temp, numPartTot, 4, fp);
  for (p = 0; p < numPartTot; p++) 
    catalog.part[p].volume = temp[p];
  fclose(fp); 
  free(temp);

  // read in desired voids
  printf(" Loading voids...\n");
  fp = fopen(voidFile ,"r");
  fgets(line, sizeof(line), fp);
  sscanf(line, "%d %s %d %s", &junkInt, junkStr, &catalog.numVoids, junkStr);
  fgets(line, sizeof(line), fp);
  
  catalog.voids.resize(catalog.numVoids);
  i = 0;
  while (fgets(line, sizeof(line), fp) != NULL) {
   sscanf(line, "%d %d %d %f %f %d %d %f %d %f %f\n", &iVoid, &voidID,
           &coreParticle, &coreDens, &zoneVol, &zoneNumPart, &numZones,
           &voidVol, &numPart, &densCon, &voidProb);

    catalog.voids[i].coreParticle = coreParticle;
    catalog.voids[i].zoneNumPart = zoneNumPart;
    catalog.voids[i].coreDens = coreDens;
    catalog.voids[i].zoneVol = zoneVol;
    catalog.voids[i].voidID = voidID;
    catalog.voids[i].vol = voidVol;
    catalog.voids[i].numPart = numPart;
    catalog.voids[i].numZones = numZones;
    catalog.voids[i].densCon = densCon;
    catalog.voids[i].voidProb = voidProb;

    catalog.voids[i].radius = pow(voidVol/catalog.numPartTot*3./4./M_PI, 1./3.);
    catalog.voids[i].numMatches = 0;
    catalog.voids[i].numBigMatches = 0;
  
    catalog.voids[i].radiusMpc = pow(voidVol/volNorm*3./4./M_PI, 1./3.);

    i++;
  }
  fclose(fp);

  catalog.numVoids = i - 1;
  catalog.voids.resize(catalog.numVoids);
  printf("  Read %d voids.\n", catalog.numVoids);

  printf(" Loading barycenters\n");
  fp = fopen(centerFile, "r");
  float tempBary[3];
  float tempFloat;
  int   tempInt;
  iVoid = 0;
  //fgets(line, sizeof(line), fp);
  while (fgets(line, sizeof(line), fp) != NULL) {
    sscanf(line, "%d %f %f %f\n", 
           &tempInt, &tempBary[0], &tempBary[1], &tempBary[2]);
    //sscanf(line, "%f %f %f %f %f %f %f %d %f %d %d %d %d\n", 
    //       &tempBary[0], &tempBary[1], &tempBary[2],
    //       &tempFloat, &tempFloat, &tempFloat, &tempFloat, &tempInt,
    //       &tempFloat, &tempInt, &tempInt, &tempInt, &tempInt);

//if (iVoid < 10) printf("BARY %d %d %e %e %e\n", iVoid, catalog.voids[iVoid].voidID, tempBary[0], tempBary[1], tempBary[2]);
    tempBary[0] = (tempBary[0] - ranges[0][0])/catalog.boxLen[0];
    tempBary[1] = (tempBary[1] - ranges[1][0])/catalog.boxLen[1];
    tempBary[2] = (tempBary[2] - ranges[2][0])/catalog.boxLen[2];
    catalog.voids[iVoid].barycenter[0] = tempBary[0];
    catalog.voids[iVoid].barycenter[1] = tempBary[1];
    catalog.voids[iVoid].barycenter[2] = tempBary[2];
    iVoid++;
  }
  fclose(fp);

  printf(" Loading shapes\n");
  fp = fopen(shapeFile, "r");
  iVoid = 0;
  float tempEllip;
  fgets(line, sizeof(line), fp);
  while (fgets(line, sizeof(line), fp) != NULL) {
    sscanf(line, "%d %f %f %f %f %f %f %f %f %f %f %f %f %f\n", 
           &tempInt, &tempEllip, &tempFloat, &tempFloat, &tempFloat,
           &tempFloat, &tempFloat, &tempFloat,
           &tempFloat, &tempFloat, &tempFloat,
           &tempFloat, &tempFloat, &tempFloat);

//if (iVoid < 10) printf("SHAPE %d %d %e\n", iVoid, catalog.voids[iVoid].voidID, tempEllip);
    catalog.voids[iVoid].ellipticity = tempEllip;
    iVoid++;
  }
  fclose(fp);


  // load up the zone membership for each void
  printf(" Loading zone-void membership info...\n");
  fp = fopen(zoneFile, "r");
  fread(&catalog.numZonesTot, 1, 4, fp);

  catalog.void2Zones.resize(catalog.numZonesTot);
    
  for (iZ = 0; iZ < catalog.numZonesTot; iZ++) {
    fread(&numZones, 1, 4, fp);

    catalog.void2Zones[iZ].numZones = numZones;
    catalog.void2Zones[iZ].zoneIDs = (int *) malloc(numZones * sizeof(int));

    for (p = 0; p < numZones; p++) {
      fread(&catalog.void2Zones[iZ].zoneIDs[p], 1, 4, fp);
    }
  }
  fclose(fp);

  // now the zone membership
  printf(" Loading particle-zone membership info...\n");
  fp = fopen(zonePartFile, "r");
  fread(&dummy, 1, 4, fp);
  fread(&numZonesTot, 1, 4, fp);

  catalog.zones2Parts.resize(numZonesTot); 

  for (iZ = 0; iZ < numZonesTot; iZ++) {
    fread(&numPart, 1, 4, fp);
    
    catalog.zones2Parts[iZ].numPart = numPart;
    catalog.zones2Parts[iZ].partIDs = (int *) malloc(numPart * sizeof(int));

    for (p = 0; p < numPart; p++) {
      fread(&catalog.zones2Parts[iZ].partIDs[p], 1, 4, fp);
    }
  }

  fclose(fp);

} // end loadCatalog

// ----------------------------------------------------------------------------
float getDist(CATALOG& catalog1, CATALOG& catalog2, int& iVoid1, int& iVoid2,
              bool periodicX, bool periodicY, bool periodicZ) {

  float rdist, dist[3];

  dist[0] = fabs(catalog1.voids[iVoid1].barycenter[0] - 
                 catalog2.voids[iVoid2].barycenter[0]);
  dist[1] = fabs(catalog1.voids[iVoid1].barycenter[1] - 
                 catalog2.voids[iVoid2].barycenter[1]);
  dist[2] = fabs(catalog1.voids[iVoid1].barycenter[2] - 
                 catalog2.voids[iVoid2].barycenter[2]);

  if (periodicX) dist[0] = fmin(dist[0], 1.0 - dist[0]);
  if (periodicY) dist[1] = fmin(dist[1], 1.0 - dist[1]);
  if (periodicZ) dist[2] = fmin(dist[2], 1.0 - dist[2]);

  rdist = sqrt(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2]);

  return rdist;

} // end getDist


// ----------------------------------------------------------------------------
void sortMatches(std::vector<MATCHPROPS>& matches, 
                 CATALOG& catalog2, int mode) {

  //printf("SORTING %d\n", matches.size());
  MATCHPROPS tempMatch;
  bool swapped;
  int matchID, matchID2;

  if (matches.size() <= 1) return;

  swapped = true;
  while (swapped) {
    swapped = false;
    for (int iMatch = 0; iMatch < matches.size() - 1; iMatch++) {
      matchID = matches[iMatch].matchID;
      matchID2 = matches[iMatch+1].matchID;
      if ((mode == 0 && matches[iMatch].merit < matches[iMatch+1].merit) ||
          (mode == 1 && matches[iMatch].dist > matches[iMatch+1].dist)) {
        tempMatch = matches[iMatch+1];
        matches[iMatch+1] = matches[iMatch];
        matches[iMatch] = tempMatch;
        swapped = true;
      }
    }
  }

  return;
} // end sortMatches
