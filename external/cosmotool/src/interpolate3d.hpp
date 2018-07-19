/*+
This is CosmoTool (./src/interpolate3d.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

guilhem.lavaux@gmail.com

This software is a computer program whose purpose is to provide a toolbox for cosmological
data analysis (e.g. filters, generalized Fourier transforms, power spectra, ...)

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
+*/

#ifndef __COSMO_INTERPOLATE3D_HPP
#define __COSMO_INTERPOLATE3D_HPP

#include "config.hpp"
#include "field.hpp"
#include <cmath>

namespace CosmoTool
{

  template<typename IType>
  class GridSampler
  {
  public:
    typedef IType result_type;

    GridSampler(IType *array_, int Nx_, int Ny_, int Nz_, int stride_)
      : array(array_), Nx(Nx_), Ny(Ny_), Nz(Nz_), stride(stride_)
    {
    }

    ~GridSampler()
    {
    }

    IType& operator()(int x, int y, int z)
      throw ()
    {
      while (x < 0)
	x += Nx;
      x %= Nx;
      while (y < 0)
	y += Ny;
      y %= Ny;
      while (z < 0)
	z += Nz;
      z %= Nz;
      
      uint32_t idx = x + Nx * (y + Ny * z);

      return array[idx*stride];
    }

  private:
    IType *array;
    int Nx, Ny, Nz, stride;
  };


  // IType is the quantity to interpolate, 
  template<typename SampledFunction, typename PosType = float>
  class Interpolate3D
  {
  public:
    typedef typename SampledFunction::result_type IType;

    Interpolate3D(SampledFunction& f)
      : sampler(f)
    {
    };

    ~Interpolate3D()
    {
    };
    
    IType get(PosType x, PosType y, PosType z)
      throw (InvalidArgumentException)
    {
      int ix = (int)std::floor(x);
      int iy = (int)std::floor(y);
      int iz = (int)std::floor(z);

      PosType rx = x-ix;
      PosType ry = y-iy;
      PosType rz = z-iz;

      IType v000 = sampler(ix,iy,iz);
      IType v001 = sampler(ix,iy,iz+1);
      IType v010 = sampler(ix,iy+1,iz);
      IType v011 = sampler(ix,iy+1,iz+1);

      IType v100 = sampler(ix+1,iy,iz);
      IType v101 = sampler(ix+1,iy,iz+1);
      IType v110 = sampler(ix+1,iy+1,iz);
      IType v111 = sampler(ix+1,iy+1,iz+1);

      return
	((1-rx) * (1-ry) * (1-rz)) * v000 +
	((1-rx) * (1-ry) *    rz)  * v001 +
	((1-rx) *    ry  * (1-rz)) * v010 +
	((1-rx) *    ry  *    rz)  * v011 +
	(   rx  * (1-ry) * (1-rz)) * v100 +
	(   rx  * (1-ry) *    rz)  * v101 +
	(   rx  *    ry  * (1-rz)) * v110 +
	(   rx  *    ry  *    rz)  * v111;
    }

  private:
    SampledFunction& sampler;    
    int Nx, Ny, Nz;
  };


  template<typename IType, typename ArrayType>
  void singleInterpolation(IType *input_array, ArrayType *x, ArrayType *y, ArrayType *z, ArrayType *scalers)
  {
  }

};

#endif
