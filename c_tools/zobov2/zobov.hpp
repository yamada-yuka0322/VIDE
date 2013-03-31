#ifndef __ZOBOV_HPP
#define __ZOBOV_HPP

typedef struct Particle {
  float dens;
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
  double vol; /* Total volume of all particles in the zone */
  double voljoin; /* Total volume of all particles in the joined void */
} ZONE;

typedef struct ZoneT {
  int nadj; /* Number of zones on border */
  int *adj; /* Each adjacent zone, with ... */
  float *slv; /* Smallest Linking Volume */
} ZONET;

#endif
