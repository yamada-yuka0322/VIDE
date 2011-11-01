#include <vector>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include "contour_pixels.hpp"

using namespace std;

static const bool DEBUG = true;

void computeContourPixels(Healpix_Map<float>& m, vector<int> contour)
{
  for (int p = 0; p < m.Npix(); p++)
    {
      fix_arr<int, 8> result;

      m.neighbors(p, result);
      for (int q = 0; q < 8; q++)
	{
	  if (result[q] < 0)
	    continue;

	  delta = (m[p]-0.5)*(m[result[q]]-0.5);
	  if (delta < 0)
	    {
	      contour.push_back(p);
	      // This is boundary go to next pixel
	      break;
	    }
	}      
    }

  if (DEBUG)
    {
      Healpix_Map<int> contour_map;

      contour_map.SetNside(RING, m.Nside());
      contour_map.fill(0);
      for (int p = 0; p < contour.size(); p++)
	{
	  contour_map[contour[p]]=1;
	}

      fitshandle h;
      h.create("!contour_map.fits");
      write_Healpix_map_to_fits(h, contour_map, planckType<int>());
    }
}
