/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/zobov2/jozov2/jozov2.hpp
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
#ifndef __JOZOV2_HPP
#define __JOZOV2_HPP

#include <string>
#include <exception>
#include "zobov.hpp"

#define BIGFLT 1e30 /* Biggest possible floating-point number */
#define NLINKS 1000 /* Number of possible links with the same rho_sl */
#define FF cout.flush()

class FileError: virtual std::exception
{
};

void readAdjacencyFile(const std::string& adjfile, PARTICLE*& p, pid_t& np)
  throw(FileError);

void readVolumeFile(const std::string& volfile, PARTICLE *p, pid_t np, 
                    pid_t mockIndex)
  throw(FileError);

void buildInitialZones(PARTICLE *p, pid_t np, pid_t* jumped, 
                       pid_t *numinh, pid_t& numZones);

void buildZoneAdjacencies(PARTICLE *p, pid_t np,
                          ZONE *z, ZONET *zt,
                          int numZones,
                          pid_t *jumped,
                          int *zonenum,
                          int *numinh);

void buildZones(PARTICLE *p, pid_t np, pid_t *&jumped,
                ZONE*& z, int& nzones,
                int*& zonenum);

void doWatershed(PARTICLE *p, pid_t np, ZONE *z, int numZones, float maxvol, float voltol);

void writeZoneFile(const std::string& zonfile, PARTICLE* p, pid_t np,
                   ZONE *z, int numZones, int* zonenum, int *jumped);

void writeVoidFile(const std::string& zonfile2, ZONE *z, int numZones);


extern "C" void findrtop(double *a, int na, int *iord, int nb);

#endif
