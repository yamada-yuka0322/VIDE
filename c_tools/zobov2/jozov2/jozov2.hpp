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
