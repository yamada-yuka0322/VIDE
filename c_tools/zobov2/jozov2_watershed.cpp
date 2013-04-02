#ifdef OPENMP
#include <omp.h>
#endif

#include <queue>
#include <set>
#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include <string>
#include "jozov2.hpp"
#include "zobov.hpp"

using namespace std;
using boost::format;

struct ZoneDensityPair
{
  int h;
  double density;

  bool operator<(const ZoneDensityPair& p2) const
  {
    return density > p2.density;
  }
};

typedef priority_queue<ZoneDensityPair> ZoneQueue;

static void build_process_queue(ZoneQueue& q, ZONE *z, char *inyet, int h)
{
  ZoneDensityPair zdp;
  ZONE& z_h = z[h];
  bool interior = true;

  assert(inyet[h] == 1);
  for (int za = 0; za < z_h.nadj; za++)
   {
     zdp.h = z_h.adj[za];
     zdp.density = z_h.slv[za];
     if (inyet[zdp.h] == 0)
       {
	 q.push(zdp);
	 interior = false;
       }
   }
  if (interior)
    inyet[h] = 2;
}

void doWatershed(PARTICLE *p, pid_t np, ZONE *z, int numZones, float maxvol, float voltol)
{
  /* Text output file */

#pragma omp parallel
  {
    char *inyet, *inyet2;
    int *zonelist, *zonelist2;
    int nhl;
    int prev_ii = -1;

    inyet = new char[numZones];
    inyet2 = new char[numZones];
    zonelist = new int[numZones];
    zonelist2 = new int[numZones];

    fill(inyet, inyet + numZones, 0);
    fill(inyet2, inyet2 + numZones, 0);

    nhl = 0;
#pragma omp for schedule(dynamic,1)
    for (int h = 0; h < numZones; h++)
      {
        int nhlcount = 0;
        float lowvol, z_cur_core_dens;
        bool beaten;
        priority_queue<ZoneDensityPair> to_process;
	int link0;

        for (int hl = 0; hl < nhl; hl++)
          inyet[zonelist[hl]] = 0;
        
        zonelist[0] = h;
        inyet[h] = 1;
        nhl = 1;
        z[h].npjoin = z[h].np;
        z_cur_core_dens = p[z[h].core].dens;

        build_process_queue(to_process, z, inyet, h);
	
        do {
          /* Find the lowest-volume (highest-density) adjacency */
          int nl = 0;
          
          beaten = false;
          
          if (to_process.empty())
            {
              beaten = true;
              z[h].leak = maxvol;
              continue;
            }

	  do
	    {
	      lowvol = to_process.top().density;
	      link0 = to_process.top().h; 
	      to_process.pop();
	    }
	  while ((inyet[link0] != 0) && (!to_process.empty()));

	  if (to_process.empty())
	    {
	      beaten = true;
	      z[h].leak = maxvol;
	      continue;
	    }

          nl++;
          
	  if (lowvol > voltol)
            {
              beaten = true;
              z[h].leak = lowvol;
              continue;
            }
          
          if (p[z[link0].core].dens < z_cur_core_dens)
            {
	      beaten = true;
              z[h].leak = lowvol;
              continue;
            }
          
          /* Add everything linked to the link(s) */
          int nhl2 = 0;
	  zonelist2[0] = link0;
	  inyet2[link0] = 1;
	  nhl2=1;
	  
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
			  if (p[z[link2].core].dens < z_cur_core_dens) {
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
		build_process_queue(to_process, z, inyet, new_h);
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
	}
	while((lowvol < BIGFLT) && (!beaten));
	
	z[h].denscontrast = z[h].leak/p[z[h].core].dens;
        if (z[h].denscontrast < 1.) 
	  z[h].denscontrast = 1.;
        
        
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
        
        z[h].nhl = nhl;
      }
    delete[] zonelist;
    delete[] zonelist2;
    delete[] inyet;
    delete[] inyet2;
    
  }
  
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
