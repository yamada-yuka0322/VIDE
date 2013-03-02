/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/libzobov/contour_pixels.cpp
    Copyright (C) 2010-2013 Guilhem Lavaux
    Copyright (C) 2011-2013 Paul M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
+*/

#include <vector>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include "contour_pixels.hpp"

using namespace std;

static const bool DEBUG = true;

void computeFilledPixels(Healpix_Map<float>& m, vector<int>& filled)
{
  filled.clear();
  for (int p = 0; p < m.Npix(); p++)
    if (m[p] > 0)
      filled.push_back(p);
}

void computeContourPixels(Healpix_Map<float>& m, vector<int>& contour)
{
  contour.clear();
  for (int p = 0; p < m.Npix(); p++)
    {
      fix_arr<int, 8> result;

      m.neighbors(p, result);
      for (int q = 0; q < 8; q++)
	{
	  if (result[q] < 0)
	    continue;

	  float delta = (m[p]-0.5)*(m[result[q]]-0.5);
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

      contour_map.SetNside(m.Nside(), RING);
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

void computeMaskPixels(Healpix_Map<float>& m, vector<int>& contour)
{
  for (int p = 0; p < m.Npix(); p++)
    {

    if (m[p]>0)
      {
        contour.push_back(p);
        // This is boundary go to next pixel
      }
    }

  if (DEBUG)
    {
      Healpix_Map<int> contour_map;

      contour_map.SetNside(m.Nside(), RING);
      contour_map.fill(0);
      for (int p = 0; p < contour.size(); p++)
  {
    contour_map[contour[p]]=1;
  }

      fitshandle h;
      h.create("!mask_map.fits");
      write_Healpix_map_to_fits(h, contour_map, planckType<int>());
    }
}

