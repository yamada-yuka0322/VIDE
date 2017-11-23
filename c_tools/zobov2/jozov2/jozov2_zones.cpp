/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/zobov2/jozov2/jozov2_zones.cpp
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
#include <iostream>
#include <fstream>
#include "jozov2.hpp"
#include "zobov.hpp"
#include <boost/format.hpp>

using namespace std;
using boost::format;

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
  
  pid_t loc_NumZones = 0;
#pragma omp parallel for schedule(static) reduction(+:loc_NumZones)
  for (pid_t i = 0; i < np; i++)
    if (numinh[i] > 0) 
      loc_NumZones++;

  numZones = loc_NumZones;
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

  size_t nadjAlloced = 0;
  cout << "Total numZones = " << numZones << endl;
  try
    {
      for (int h = 0; h < numZones; h++) {
        cout << "Zone " << h << " nadj = " << zt[h].nadj << endl;
        if (zt[h].nadj > 0) {
          zt[h].adj = new pid_t[zt[h].nadj];
          zt[h].slv = new float[zt[h].nadj];
          nadjAlloced += zt[h].nadj;
          zt[h].nadj = 0;
        } else {
          zt[h].adj = 0;
          zt[h].slv = 0;
        }
      }
    }
  catch(const std::bad_alloc& a)
    {
      cout << "Could not allocate memory for zone adjacencies (nadj so far: " << nadjAlloced << ", memory needed: " << (nadjAlloced*(sizeof(pid_t)+sizeof(float))) << ")" << endl;
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
    if (zt[h].nadj > 0) {
      z[h].adj = new pid_t[zt[h].nadj];
      z[h].slv = new float[zt[h].nadj];
      for (int za = 0; za < zt[h].nadj; za++) {
        z[h].adj[za] = zt[h].adj[za];
        z[h].slv[za] = zt[h].slv[za];
      }
      delete[] zt[h].adj;
      delete[] zt[h].slv;
    } else {
      z[h].adj = 0;
      z[h].slv = 0;
    }
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
        z[h].vol = 0;
        z[h].np = z[h].npjoin = z[h].nadj = z[h].nhl = 0;
        z[h].leak = 0;
        z[h].adj = 0;
        z[h].slv = 0;
        z[h].denscontrast = 0;
        z[h].vol = z[h].voljoin = 0;
        zonenum[i] = h;
        h++;
      } else {
        zonenum[i] = -1;
      }
    }
  for (size_t z = 0; z < nzones; z++) {
    zt[z].nadj = 0;
    zt[z].adj = 0;
    zt[z].slv = 0;
  }

  buildZoneAdjacencies(p, np, z, zt,
                       h, jumped, zonenum, numinh);

  delete[] zt;
  delete[] numinh;

  for (pid_t i=0; i < np; i++) {
    int h = zonenum[jumped[i]];
    z[h].vol += 1.0/(double)p[i].dens;
    z[h].numzones = 0;
    z[h].zonelist = 0;
  }
}
