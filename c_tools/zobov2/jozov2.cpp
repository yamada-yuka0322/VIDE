#include <boost/format.hpp>
#include <exception>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <assert.h>
/* jovoz.c by Mark Neyrinck */

#include "zobov.hpp"

using namespace std;
using boost::format;

#define BIGFLT 1e30 /* Biggest possible floating-point number */
#define NLINKS 1000 /* Number of possible links with the same rho_sl */
#define FF cout.flush()

extern "C" void findrtop(double *a, int na, int *iord, int nb);

class FileError: virtual std::exception
{
};

void readAdjacencyFile(const string& adjfile, PARTICLE*& p, pid_t& np)
 throw(FileError)
{
  ifstream adj(adjfile.c_str());

  if (!adj) {
    cout << format("Unable to open %s") % adjfile << endl;
    throw FileError();
  }
  adj.read((char *)&np, sizeof(pid_t));
  
  cout << format("adj: %d particles") % np << endl;
  FF;

  p = new PARTICLE[np];
  
  /* Adjacencies*/
  for (pid_t i=0; i < np; i++) {
    adj.read((char*)&p[i].nadj, sizeof(pid_t));
    if (!adj)
      throw FileError();

    /* The number of adjacencies per particle */
    if (p[i].nadj > 0)
      p[i].adj = new pid_t[p[i].nadj];
    else 
      p[i].adj = 0;
    p[i].ncnt = 0; /* Temporarily, it's an adj counter */
  }

  for (pid_t i = 0; i < np; i++)
    {
      pid_t nin;
      adj.read((char*)&nin, sizeof(pid_t));
      if (!adj)
        throw FileError();

      if (nin > 0)
        {
          for (int k=0; k < nin; k++) {
            pid_t j;
            
            adj.read((char *)&j,sizeof(pid_t));
            if (j < np)
              {
                /* Set both halves of the pair */
                assert(i < j);
                if (p[i].ncnt == p[i].nadj)
                  {
		    cout << format("OVERFLOW for particle %d (pending %d). List of accepted:") % i % j << endl;
                    for (int q = 0; q < p[i].nadj; q++)
                      cout << format("  %d") % p[i].adj[q] << endl;
                    throw FileError();
                  }
                
                if (p[j].ncnt == p[j].nadj)
                  {
		    cout << format("OVERFLOW for particle %d (pending %d). List of accepted:") % j % i << endl;
                    for (int q = 0; q < p[j].nadj; q++)
                      cout << format("  %d\n") % p[j].adj[q] << endl;
                    throw FileError();
                  }
                
                p[i].adj[p[i].ncnt] = j;
                p[j].adj[p[j].ncnt] = i;
                p[i].ncnt++; 
                p[j].ncnt++;
              }
            else
              {
                cout << format("%d: adj = %d") % i % j << endl;
              }
          }
        }
    }
}

void readVolumeFile(const std::string& volfile, PARTICLE *p, pid_t np, 
                    pid_t mockIndex)
  throw(FileError)
{
  ifstream vol(volfile.c_str());
  pid_t np2;
  
  if (!vol)
    {
      cout << "Unable to open volume file " << volfile << endl;
      throw FileError();
    }
  vol.read((char *)&np2, sizeof(pid_t));
  if (np != np2) {
    cout << format("Number of particles doesn't match! %d != %d") % np %np2 << endl;
    throw FileError();
  }
  cout << format("%d particles") % np << endl;
  FF;

  for (pid_t i = 0; i < np; i++) {
    vol.read((char*)&p[i].dens, sizeof(float));
    if (((p[i].dens < 1e-30) || (p[i].dens > 1e30)) && (i < mockIndex)) {
      cout << format("Whacked-out volume found, of particle %d: %f") % i % p[i].dens << endl;
      p[i].dens = 1.;
    }
    p[i].dens = 1./p[i].dens; /* Get density from volume */
  }
}

void buildInitialZones(PARTICLE *p, pid_t np, pid_t* jumped, 
                      pid_t *numinh, pid_t& numZones)
{  
  pid_t *jumper = new pid_t[np];
  float minvol;

  /* find jumper */
  for (pid_t i = 0; i < np; i++) {
    PARTICLE& cur_p = p[i];

    minvol = cur_p.dens;
    jumper[i] = -1;
    for (int j = 0; j < cur_p.nadj; j++) {
      if (p[cur_p.adj[j]].dens < minvol) {
	jumper[i] = cur_p.adj[j];
	minvol = p[jumper[i]].dens;
      }
    }
    numinh[i] = 0;
  }

  (cout << "About to jump ... " << endl).flush();
  
  /* Jump */
  for (pid_t i = 0; i < np; i++) {
    jumped[i] = i;
    while (jumper[jumped[i]] > -1)
      jumped[i] = jumper[jumped[i]];
    numinh[jumped[i]]++;
  }
  (cout << "Post-jump ..." << endl).flush();
  
  numZones = 0;
  for (pid_t i = 0; i < np; i++)
    if (numinh[i] > 0) 
      numZones++;
  cout << format("%d initial zones found") % numZones << endl;

  delete[] jumper;
}

void buildZoneAdjacencies(PARTICLE *p, pid_t np,
                          ZONE *z, ZONET *zt,
                          int numZones,
                          pid_t *jumped,
                          int *zonenum,
                          int *numinh)
{
  /* Count border particles */
  for (pid_t i = 0; i < np; i++)
    for (int j = 0; j < p[i].nadj; j++) {
      pid_t testpart = p[i].adj[j];
      if (jumped[i] != jumped[testpart])
	zt[zonenum[jumped[i]]].nadj++;
    }

  try
    {
      for (int h = 0; h < numZones; h++) {
        zt[h].adj = new pid_t[zt[h].nadj];
        zt[h].slv = new float[zt[h].nadj];
        zt[h].nadj = 0;
      }
    }
  catch(const std::bad_alloc& a)
    {
      cout << "Could not allocate memory for zone adjacencies." << endl;
      throw a;
    }

  /* Find "weakest links" */
  for (pid_t i = 0; i < np; i++)
    {
      int h = zonenum[jumped[i]];
      ZONET& zt_h = zt[h];
    
      for (int j = 0; j < p[i].nadj; j++)
        {
          pid_t testpart = p[i].adj[j];
          bool already;

          if (h != zonenum[jumped[testpart]]) {
            if (p[testpart].dens > p[i].dens) {
              /* there could be a weakest link through testpart */
              already = false;
              for (int za = 0; za < zt_h.nadj; za++)
                if (zt_h.adj[za] == zonenum[jumped[testpart]]) {
                  already = true;
                  if (p[testpart].dens < zt_h.slv[za]) {
                    zt_h.slv[za] = p[testpart].dens;
                  }
                }
              if (!already) {
                zt_h.adj[zt_h.nadj] = zonenum[jumped[testpart]];
                zt_h.slv[zt_h.nadj] = p[testpart].dens;
                zt_h.nadj++;
              }
            } else { /* There could be a weakest link through i */
              already = false;
              for (int za = 0; za < zt_h.nadj; za++)
                if (zt_h.adj[za] == zonenum[jumped[testpart]]) {
                  already = true;
                  if (p[i].dens < zt_h.slv[za]) {
                    zt_h.slv[za] = p[i].dens;
                  }
                }
              if (!already) {
                zt_h.adj[zt_h.nadj] = zonenum[jumped[testpart]];
                zt_h.slv[zt_h.nadj] = p[i].dens;
                zt_h.nadj++;
              }
            }
          }
        }
    }
  (cout <<"Found zone adjacencies" << endl).flush();

  /* Free particle adjacencies */
  for (pid_t i = 0; i < np; i++)
    {
      if (p[i].adj != 0)
        delete[] p[i].adj;
    }

  /* Use z instead of zt */
  for (int h = 0; h < numZones; h++) {
    z[h].nadj = zt[h].nadj;
    z[h].adj = new pid_t[zt[h].nadj];
    z[h].slv = new float[zt[h].nadj];
    for (int za = 0; za < zt[h].nadj; za++) {
      z[h].adj[za] = zt[h].adj[za];
      z[h].slv[za] = zt[h].slv[za];
    }
    delete[] zt[h].adj;
    delete[] zt[h].slv;
    z[h].np = numinh[z[h].core];
  }
}

void buildZones(PARTICLE *p, pid_t np, pid_t *&jumped, 
                ZONE*& z, int& nzones,
                int*& zonenum)
{
  pid_t *numinh;
  ZONET *zt;

  jumped = new pid_t[np];
  numinh = new pid_t[np];

  buildInitialZones(p, np, jumped, numinh, nzones);

  try
    {
      z = new ZONE[nzones];
      zt = new ZONET[nzones];
      zonenum = new int[np];
    }
  catch (const std::bad_alloc& e)
    {
      cout << "Unable to do zone allocations" << endl;
      throw e;
    }

  int h = 0;
  for (pid_t i = 0; i < np; i++)
    {
      if (numinh[i] > 0) {
        z[h].core = i;
        zonenum[i] = h;
        h++;
      } else {
        zonenum[i] = -1;
      }
    }

  buildZoneAdjacencies(p, np, z, zt,
                       nzones, jumped, zonenum, numinh);

  delete[] zt;
  delete[] numinh;

  for (pid_t i=0; i < np; i++) {
    int h = zonenum[i];
    z[h].vol += 1.0/(double)p[i].dens;
    z[h].numzones = 0;
    z[h].zonelist = 0;
  }
}


void writeZoneFile(const std::string& zonfile, PARTICLE* p, pid_t np,
                   ZONE *z, int numZones, int* zonenum, int *jumped)
{
  pid_t **m = new pid_t *[numZones];
  pid_t *nm = new pid_t[numZones];

  /* Not in the zone struct since it'll be freed up (contiguously, we hope)
     soon */
  for (int h=0; h < numZones; h++)
    {
      m[h] = new pid_t[z[h].np];
      nm[h] = 0;
      z[h].vol = 0.;
    }

  for (pid_t i = 0; i < np; i++)
    {
      int h = zonenum[jumped[i]];
      if (z[h].core == i)
        {
          m[h][nm[h]] = m[h][0];
          m[h][0] = i; /* Put the core particle at the top of the list */
        } 
      else
        {
          m[h][nm[h]] = i;
        }
      nm[h] ++;
    }
  delete[] nm;

  ofstream zon(zonfile.c_str());
  if (!zon) 
    {
      cout << format("Problem opening zonefile %s.") % zonfile << endl;
      throw FileError();
    }
  
  zon.write((char *)&np, sizeof(pid_t));
  zon.write((char *)&numZones, sizeof(int));
  for (int h = 0; h < numZones; h++)
    {
      zon.write((char *)&(z[h].np), sizeof(pid_t));
      zon.write((char *)m[h],z[h].np * sizeof(pid_t));
      delete[] m[h];
    }
  delete[] m;
}

void doWatershed(const std::string& zonfile2, PARTICLE *p, pid_t np, ZONE *z, int numZones, float maxvol, float voltol)
{
  char *inyet, *inyet2;
  int *zonelist, *zonelist2;
  int nhl;
  int *links = new int[NLINKS];
  int *iord;
  float maxdenscontrast = 0;
  bool *done_zones;

  ofstream zon2(zonfile2.c_str());
  if (!zon2)
    {
      cout << format("Problem opening zonefile %s.)") % zonfile2 << endl;
      throw FileError();
    }
  zon2.write((char *)&numZones,sizeof(int));

  inyet = new char[numZones];
  inyet2 = new char[numZones];
  zonelist = new int[numZones];
  zonelist2 = new int[numZones];
  done_zones = new bool[numZones];

  fill(inyet, inyet + numZones, 0);
  fill(inyet2, inyet2 + numZones, 0);
  fill(done_zones, done_zones + numZones, false);

  double *sorter = new double[numZones+1];
  /* Assign sorter by probability (could use volume instead) */
  for (int h = 0; h < numZones; h++)
    sorter[h] = (double)z[h].core;
    
  /* Text output file */

  printf("about to sort (pre-watershed sort) ...\n");FF;

  iord = new int[numZones];
 
  findrtop(sorter, numZones, iord, numZones);
  delete[] sorter;

  nhl = 0;
  for (int ii = 0; ii < numZones; ii++)
    {
      int nhlcount = 0;
      int h = iord[ii];
      float lowvol;
      bool beaten;
      
      for (int hl = 0; hl < nhl; hl++)
        inyet[zonelist[hl]] = 0;

      zonelist[0] = h;
      inyet[h] = 1;
      nhl = 1;
      z[h].npjoin = z[h].np;
      do {
        /* Find the lowest-volume (highest-density) adjacency */
	int nl = 0;
	
	beaten = false;
	lowvol = BIGFLT;

        for (int hl = 0; hl < nhl; hl++)
          {
            int h2 = zonelist[hl];
            if (inyet[h2] == 1) { /* If it's not already identified as 
                                     an interior zone, with inyet=2 */
              bool interior = true;
              for (int za = 0; za < z[h2].nadj; za++)
                {
                  if (inyet[z[h2].adj[za]] == 0)
                    {
                      interior = false;
                      if (z[h2].slv[za] == lowvol)
                        {
                          links[nl] = z[h2].adj[za];
                          nl ++;
                          if (nl == NLINKS) 
                            {
                              printf("Too many links with the same rho_sl!  Increase NLINKS from %d\n",nl);
                              exit(0);
                            }
                        }
                      if (z[h2].slv[za] < lowvol)
                        {
                          lowvol = z[h2].slv[za];
                          links[0] = z[h2].adj[za];
                          nl = 1;
                        }
                    }
                }
              if (interior)
		inyet[h2] = 2; /* No bordering exter. zones */
            }
          }
        
        if (nl == 0)
          {
            beaten = true;
            z[h].leak = maxvol;
            continue;
          }
	
        if (lowvol > voltol)
          {
            beaten = true;
            z[h].leak = lowvol;
            continue;
          }
        
        for (int l = 0; l < nl; l++)
          {
            if (!done_zones[links[l]]) /* Equivalent to p[z[links[l]].core].dens < p[z[h].core].dens) as zones are sorted. */
              beaten = true;
          }

        if (beaten)
          {
            z[h].leak = lowvol;
            continue;
          }

        /* Add everything linked to the link(s) */
        int nhl2 = 0;
        for (int l = 0; l < nl; l++)
          {
            if (inyet2[links[l]] == 0)
              {
                zonelist2[nhl2] = links[l];
                inyet2[links[l]] = 1;
                nhl2 ++;
                bool added = true;
                while (added && !beaten)
                  {
                    added = false;
                    for (int hl = 0; (hl < nhl2) && (!beaten); hl++)
                      {
                        int h2 = zonelist2[hl];

                        if (inyet2[h2] == 1) {
                          bool interior = true; /* Guilty until proven innocent */
                          for (int za = 0; za < z[h2].nadj; za ++) {
                            int link2 = z[h2].adj[za];
                            if ((inyet[link2]+inyet2[link2]) == 0) {
                              interior = false;
                              if (z[h2].slv[za] <= lowvol) {
                                if (!done_zones[link2]) { // Equivalent to p[z[link2].core].dens < p[z[h].core].dens)
                                  beaten = true;
                                  break;
                                }
                                zonelist2[nhl2] = link2;
                                inyet2[link2] = 1;
                                nhl2++;
                                added = true;
                              }
                            }
                          }
                          if (interior)
			    inyet2[h2] = 2;
                        }
                      }
                  }
              }
          }
        for (int hl = 0; hl < nhl2; hl++)
          inyet2[zonelist2[hl]] = 0;
        
        /* See if there's a beater */
        if (beaten) {
          z[h].leak = lowvol;
        } else {
          for (int h2 = 0; h2 < nhl2; h2++) {
            zonelist[nhl] = zonelist2[h2];
            inyet[zonelist2[h2]] = 1;
            nhl++;
            z[h].npjoin += z[zonelist2[h2]].np;
          }
        }
        if (nhl/10000 > nhlcount) {
          nhlcount = nhl/10000;
          printf(" %d",nhl); FF;
        }
      } while((lowvol < BIGFLT) && (!beaten));
      
      z[h].denscontrast = z[h].leak/p[z[h].core].dens;
      if (z[h].denscontrast < 1.) z[h].denscontrast = 1.;
      
      /* find biggest denscontrast */
      if (z[h].denscontrast > maxdenscontrast) {
        maxdenscontrast = (double)z[h].denscontrast;
      }
      
      /* Don't sort; want the core zone to be first */
      
      if (nhlcount > 0) { /* Outputs the number of zones in large voids */
        printf(" h%d:%d\n",h,nhl);
        FF;
      }
      /* Calculate volume */
      z[h].voljoin = 0.;
      z[h].zonelist = new int[nhl];
      z[h].numzones = nhl;
      for (int q = 0; q < nhl; q++) {
        z[h].voljoin += z[zonelist[q]].vol;
	z[h].zonelist[q] = zonelist[q];
      }
      
      done_zones[h] = true;
      z[h].nhl = nhl;
    }
  delete[] zonelist;
  delete[] zonelist2;
  delete[] links;
  delete[] iord;

  cout << "Writing void/zone relations..." << endl;
  for (int h = 0; h < numZones; h++)
    {
      zon2.write((char *)&z[h].nhl, sizeof(int));
      zon2.write((char *)z[h].zonelist, z[h].nhl*sizeof(int));
    }

  cout << format("Maxdenscontrast = %f.") % maxdenscontrast << endl;
}

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

  doWatershed(zonfile2, p, np, z, nzones, maxvol, voltol);

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
