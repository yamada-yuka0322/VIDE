/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/zobov2/jozov2/jozov2_io.cpp
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
#include <string>
#include <boost/format.hpp>
#include "jozov2.hpp"
#include "zobov.hpp"

using namespace std;
using boost::format;

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
//                    throw FileError();
                  }
                else               
                if (p[j].ncnt == p[j].nadj)
                  {
		    cout << format("OVERFLOW for particle %d (pending %d). List of accepted:") % j % i << endl;
                    for (int q = 0; q < p[j].nadj; q++)
                      cout << format("  %d\n") % p[j].adj[q] << endl;
//                    throw FileError();
                  }
                else{
                p[i].adj[p[i].ncnt] = j;
                p[j].adj[p[j].ncnt] = i;
                p[i].ncnt++; 
                p[j].ncnt++; }
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

void writeVoidFile(const string& zonfile2, ZONE *z, int numZones)
{
  ofstream zon2(zonfile2.c_str());
  if (!zon2)
    {
      cout << format("Problem opening zonefile %s.)") % zonfile2 << endl;
      throw FileError();
    }
  zon2.write((char *)&numZones,sizeof(int));
  cout << "Writing void/zone relations..." << endl;
  for (int h = 0; h < numZones; h++)
    {
      zon2.write((char *)&z[h].nhl, sizeof(int));
      zon2.write((char *)z[h].zonelist, z[h].nhl*sizeof(int));
    }
}
