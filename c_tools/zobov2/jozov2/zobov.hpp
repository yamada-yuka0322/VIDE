/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/zobov2/jozov2/zobov.hpp
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
#ifndef __ZOBOV_HPP
#define __ZOBOV_HPP

typedef struct Particle {
  float dens;
  float weight;
  int nadj;
  int ncnt;
  int *adj;
} PARTICLE;

typedef struct Zone {
  int core; /* Identity of peak particle */
  int np; /* Number of particles in zone */
  int npjoin; /* Number of particles in the joined void */
  int nadj; /* Number of adjacent zones */
  int nhl; /* Number of zones in final joined void */
  float leak; /* Volume of last leak zone*/
  int *adj; /* Each adjacent zone, with ... */
  float *slv; /* Smallest Linking Volume */
  float denscontrast; /* density contrast */
  double vol; /* Total weighted volume of all particles in the zone */
  double volume; /* Total volume of all particles in the zone */
  double voljoin; /* Total volume of all particles in the joined void */

  int *zonelist; /* Zones bound to the void. */
  int numzones; /* Number of zones bound. */
} ZONE;

typedef struct ZoneT {
  int nadj; /* Number of zones on border */
  int *adj; /* Each adjacent zone, with ... */
  float *slv; /* Smallest Linking Volume */
} ZONET;

#endif
