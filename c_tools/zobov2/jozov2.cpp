#include <boost/format.hpp>
#include <exception>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <assert.h>
/* jovoz.c by Mark Neyrinck */
#include "jozov2.hpp"
#include "zobov.hpp"

using namespace std;
using boost::format;



int main(int argc,char **argv)
{
  ofstream txt;
  PARTICLE *p;
  pid_t np;

  ZONE *z;
  int nzones;

  string adjfile, volfile, zonfile, zonfile2, txtfile;
  int *jumped;
  int *zonenum;
  float voltol, prob;

  float maxvol, minvol;
  double *sorter;
  int *iord;

  int mockIndex;

  if (argc != 8) {
    printf("Wrong number of arguments.\n");
    printf("arg1: adjacency file\n");
    printf("arg2: volume file\n");
    printf("arg3: output file containing particles in each zone\n");
    printf("arg4: output file containing zones in each void\n");
    printf("arg5: output text file\n");
    printf("arg6: Density threshold (0 for no threshold)\n");
    printf("arg7: Beginning index of mock galaxies\n\n");
    exit(0);
  }
  adjfile = argv[1];
  volfile = argv[2];
  zonfile = argv[3];
  zonfile2 = argv[4];
  txtfile = argv[5];
  if (sscanf(argv[6],"%f",&voltol) == 0) {
    printf("Bad density threshold.\n");
    exit(0);
  }
  if (sscanf(argv[7],"%d",&mockIndex) == 0) {
    printf("Bad mock galaxy index.\n");
    exit(0);
  }
  printf("TOLERANCE: %f\n", voltol);
  if (voltol <= 0.) {
    printf("Proceeding without a density threshold.\n");
    voltol = 1e30;
  }

  try
    {
      readAdjacencyFile(adjfile, p, np);
    }
  catch(const FileError& e)
    {
      return 1;
    }
  if (mockIndex < 0)
    mockIndex = np;

  /* Check that we got all the pairs */
  for (int i = 0; i < np; i++)
    {
      if (p[i].ncnt != p[i].nadj && i < mockIndex) {
	cout
	  << format("We didn't get all of %d's adj's; %d != %d.") 
	  % i % p[i].ncnt % p[i].nadj
	  << endl;
	p[i].nadj = p[i].ncnt;
      }
    }

  /* Volumes */
  try
    {
      readVolumeFile(volfile, p, np, mockIndex);
    }
  catch (const FileError& e)
    {
      return 2;
    }

  buildZones(p, np, jumped, z, nzones, zonenum);
  writeZoneFile(zonfile, p, np, z, nzones, zonenum, jumped);

  maxvol = 0.;
  minvol = BIGFLT;
  for(pid_t i = 0; i < np; i++){
    if (p[i].dens > maxvol)
      maxvol = p[i].dens;
    if (p[i].dens < minvol) 
      minvol = p[i].dens;
  }
  cout << format("Densities range from %e to %e.") % minvol % maxvol << endl;
  FF;

  doWatershed(p, np, z, nzones, maxvol, voltol);
  writeVoidFile(zonfile2, z, nzones);

  sorter = new double[nzones+1];
  /* Assign sorter by probability (could use volume instead) */
  for (int h = 0; h < nzones; h++)
    sorter[h] = (double)z[h].denscontrast;
    
  /* Text output file */

  printf("about to sort ...\n");FF;

  iord = new int[nzones];

  findrtop(sorter, nzones, iord, nzones);
  delete[] sorter;

  txt.open(txtfile.c_str());
  txt << format("%d particles, %d voids.") % np % nzones << endl;
  txt << "Void# FileVoid# CoreParticle CoreDens ZoneVol Zone#Part Void#Zones VoidVol Void#Part VoidDensContrast VoidProb" << endl;
  for (int h=0; h<nzones; h++)
    {
      int i = iord[h];
      prob = exp(-5.12*(z[i].denscontrast-1.) - 0.8*pow(z[i].denscontrast-1.,2.8));
      if (z[i].np == 1) 
	continue;
      
      txt << format("%d %d %d %e %e %d %d %e %d %f %6.2e")
	   % (h+1) % i % z[i].core % p[z[i].core].dens % z[i].vol
           % z[i].np % z[i].nhl % z[i].voljoin % z[i].npjoin
	   % z[i].denscontrast % prob << endl;

    } /* h+1 to start from 1, not zero */
  txt.close();

  delete[] iord;
  delete[] z;
  delete[] p;

  cout << "Done!" << endl;
  return(0);
}
