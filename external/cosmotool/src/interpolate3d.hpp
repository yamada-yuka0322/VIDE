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

};

#endif
