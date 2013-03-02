/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/libzobov/contour_pixels.hpp
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

#ifndef  __CONTOUR_PIXELS_HPP
#define  __CONTOUR_PIXELS_HPP

#include <vector>
#include <healpix_map.h>


void computeContourPixels(Healpix_Map<float>& map, std::vector<int>& contour);
void computeFilledPixels(Healpix_Map<float>& map, std::vector<int>& contour);
void computeMaskPixels(Healpix_Map<float>& map, std::vector<int>& contour);

#endif
