#ifndef  __CONTOUR_PIXELS_HPP
#define  __CONTOUR_PIXELS_HPP

#include <vector>
#include <healpix_map.h>

void computeContourPixels(Healpix_Map<float>& map, std::vector<int> contour);

#endif
