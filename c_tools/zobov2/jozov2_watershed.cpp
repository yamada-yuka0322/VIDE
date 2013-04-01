#ifdef OPENMP
#include <omp.h>
#endif

#include <set>
#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include <string>
#include "jozov2.hpp"
#include "zobov.hpp"

using namespace std;
using boost::format;

void doWatershed(PARTICLE *p, pid_t np, ZONE *z, int numZones, float maxvol, float voltol)
{
  int *iord;

  double *sorter = new double[numZones+1];
  /* Assign sorter by probability (could use volume instead) */
  for (int h = 0; h < numZones; h++)
    sorter[h] = (double)z[h].core;
    
  /* Text output file */

  printf("about to sort (pre-watershed sort) ...\n");FF;

  iord = new int[numZones];
 
  findrtop(sorter, numZones, iord, numZones);
  delete[] sorter;

#pragma omp parallel
  {
    char *inyet, *inyet2;
    int *zonelist, *zonelist2;
    int nhl;
    int *links = new int[NLINKS];
    bool *done_zones;
    int prev_ii = -1;

    inyet = new char[numZones];
    inyet2 = new char[numZones];
    zonelist = new int[numZones];
    zonelist2 = new int[numZones];
    done_zones = new bool[numZones];

    fill(inyet, inyet + numZones, 0);
    fill(inyet2, inyet2 + numZones, 0);
    fill(done_zones, done_zones + numZones, false);

    nhl = 0;
#pragma omp for schedule(dynamic,1)
    for (int h = 0; h < numZones; h++)
      {
        int nhlcount = 0;
        float lowvol;
        bool beaten;
        set<int> to_process;

        for (int hl = 0; hl < nhl; hl++)
          inyet[zonelist[hl]] = 0;
        
        zonelist[0] = h;
        inyet[h] = 1;
        nhl = 1;
        to_process.insert(h);
        z[h].npjoin = z[h].np;
        do {
          /* Find the lowest-volume (highest-density) adjacency */
          int nl = 0;
          
          beaten = false;
          lowvol = BIGFLT;
          
          set<int>::iterator iter = to_process.begin();
          while (iter != to_process.end())
            {
              int h2 = *iter;
              ZONE& z_cur = z[h2];
              bool interior = true, touched = false;
              assert(inyet[h2] == 1);

              for (int za = 0; za < z_cur.nadj; za++)
                {
                 if (inyet[z_cur.adj[za]] == 0)
                    {
                      interior = false;
                      if (z_cur.slv[za] == lowvol)
                        {
                          links[nl] = z_cur.adj[za];
                          touched = true;
                          nl ++;
                          if (nl == NLINKS) 
                            {
                              printf("Too many links with the same rho_sl!  Increase NLINKS from %d\n",nl);
                              exit(0);
                            }
                        }
                      else
                        if (z_cur.slv[za] < lowvol)
                          {
                            lowvol = z_cur.slv[za];
                            links[0] = z_cur.adj[za];
                            nl = 1;
                            touched = true;
                          }
                    }
                }
              ++iter;
              if (interior || !touched)
                to_process.erase(h2);
              if (interior)
                inyet[h2] = 2; /* No bordering exter. zones */
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
              //if (!done_zones[links[l]]) /* Equivalent to p[z[links[l]].core].dens < p[z[h].core].dens) as zones are sorted. */
              if (p[z[links[l]].core].dens < p[z[h].core].dens)
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
                                  //if (!done_zones[link2]) { // Equivalent to p[z[link2].core].dens < p[z[h].core].dens)
                                  if (p[z[link2].core].dens < p[z[h].core].dens) {
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
          
          /* See if there's a beater */
          if (beaten) {
            z[h].leak = lowvol;
          } else {
            for (int h2 = 0; h2 < nhl2; h2++) {
              int new_h = zonelist2[h2];

              zonelist[nhl] = new_h;
              assert(inyet[new_h] == 0);
              inyet[new_h] = 1;
              if (inyet2[new_h] != 2)
                to_process.insert(new_h);
              nhl++;
              z[h].npjoin += z[new_h].np;
            }
          }

          for (int hl = 0; hl < nhl2; hl++)
            inyet2[zonelist2[hl]] = 0;
          if (nhl/10000 > nhlcount) {
            if (nhlcount == 0)
              (cout << format("Zone %d: %d") % h % nhl).flush(); 
            else
              (cout << format(" %d [%d]") % nhl % to_process.size()).flush();
            nhlcount = nhl/10000;
          }
        } while((lowvol < BIGFLT) && (!beaten));
        
        z[h].denscontrast = z[h].leak/p[z[h].core].dens;
        if (z[h].denscontrast < 1.) z[h].denscontrast = 1.;
        
        
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
    delete[] inyet;
    delete[] inyet2;
    delete[] done_zones;

  }
  delete[] iord;

  double maxdenscontrast = 0;
#pragma omp parallel shared(maxdenscontrast)
  {
    double maxdenscontrast_local = 0;

#pragma omp for schedule(static)
    for (int h = 0; h < numZones; h++)
      {
        /* find biggest denscontrast */
        if (z[h].denscontrast > maxdenscontrast_local) {
          maxdenscontrast_local = (double)z[h].denscontrast;
        }
      }

#pragma omp critical
    {
      if (maxdenscontrast_local > maxdenscontrast)
        maxdenscontrast = maxdenscontrast_local;
    }
  }

  cout << format("Maxdenscontrast = %f.") % maxdenscontrast << endl;
}
