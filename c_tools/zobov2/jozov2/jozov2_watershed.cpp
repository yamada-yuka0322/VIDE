/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/zobov2/jozov2/jozov2_watershed.cpp
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
  double core;

  bool operator<(const ZoneDensityPair& p2) const
  {
    return (density > p2.density) || (density==p2.density && core > p2.core);
  }
};

typedef priority_queue<ZoneDensityPair> ZoneQueue;

static void build_process_queue(ZoneQueue& q, ZONE *z, PARTICLE *p, char *inyet, int h)
{
  ZoneDensityPair zdp;
  ZONE& z_h = z[h];
  bool interior = true;

  assert(inyet[h] == 1);
  for (int za = 0; za < z_h.nadj; za++)
   {
     zdp.h = z_h.adj[za];
     zdp.density = z_h.slv[za];
     zdp.core = p[z[zdp.h].core].dens;
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
        float previous_lowvol = BIGFLT, lowvol, z_cur_core_dens;
        bool beaten;
        priority_queue<ZoneDensityPair> to_process;
	int link0;
        int save_nhl;
        pid_t save_npjoin;

        for (int hl = 0; hl < nhl; hl++)
          inyet[zonelist[hl]] = 0;
        
        zonelist[0] = h;
        inyet[h] = 1;
        save_nhl = nhl = 1;
        save_npjoin = z[h].npjoin = z[h].np;
        z_cur_core_dens = p[z[h].core].dens;

        build_process_queue(to_process, z, p, inyet, h);
	
        do {
          /* Find the lowest-volume (highest-density) adjacency */
          beaten = false;
          
          if (to_process.empty())
            {
              beaten = true;
              z[h].leak = maxvol;
              save_npjoin = z[h].npjoin;
              save_nhl = nhl;
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
              save_npjoin = z[h].npjoin;
              save_nhl = nhl;
	      beaten = true;
	      z[h].leak = maxvol;
	      continue;
	    }

          /* See if there's a beater */
          if (previous_lowvol != lowvol)
            {
              save_npjoin = z[h].npjoin;
              save_nhl = nhl;
              previous_lowvol = lowvol;
            }

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
              z[h].npjoin += z[new_h].np;
              inyet[new_h] = 1;
              if (inyet2[new_h] != 2)
                build_process_queue(to_process, z, p, inyet, new_h);
              nhl++;
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

        if (!beaten)
          { 
            save_npjoin = z[h].npjoin;
            save_nhl = nhl;
          }
	
	z[h].denscontrast = z[h].leak/p[z[h].core].dens;
        if (z[h].denscontrast < 1.) 
	  z[h].denscontrast = 1.;
        
        
        /* Don't sort; want the core zone to be first */
        
        if (nhlcount > 0) { /* Outputs the number of zones in large voids */
          printf(" h%d:%d\n",h,nhl);
          FF;
        }
        /* Calculate volume */
        z[h].npjoin = save_npjoin;
        z[h].voljoin = 0.;
        z[h].zonelist = new int[save_nhl];
        z[h].numzones = save_nhl;
        for (int q = 0; q < save_nhl; q++) {
          z[h].voljoin += z[zonelist[q]].volume; /* voljoin should not be weighted */
          z[h].zonelist[q] = zonelist[q];
        }
        
        z[h].nhl = save_nhl;
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
