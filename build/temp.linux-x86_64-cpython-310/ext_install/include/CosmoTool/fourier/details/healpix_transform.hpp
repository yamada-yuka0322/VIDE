/*+
This is CosmoTool (./src/fourier/details/healpix_transform.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __COSMOTOOL_FOURIER_HEALPIX_DETAILS_TRANSFORM_HPP
#define __COSMOTOOL_FOURIER_HEALPIX_DETAILS_TRANSFORM_HPP

#include <valarray>

namespace CosmoTool
{

  template<typename T> struct HealpixJobHelper__ {};

  template<> struct HealpixJobHelper__<double>
    { enum {val=SHARP_DP}; };

  template<> struct HealpixJobHelper__<float>
    { enum {val=0}; };


  template<typename T>
  class HealpixFourierTransform: public FourierTransform<T>
  {
  private:
    sharp_alm_info *ainfo;
    sharp_geom_info *ginfo;
    HealpixFourierMap<T> realMap;
    HealpixFourierALM<T> fourierMap;
    int m_iterate;
  public:
    HealpixFourierTransform(long nSide, long Lmax, long Mmax, int iterate = 0)
      : realMap(nSide), fourierMap(Lmax, Mmax), ainfo(0), ginfo(0), m_iterate(iterate)
    {
      sharp_make_healpix_geom_info (nSide, 1, &ginfo);
      sharp_make_triangular_alm_info (Lmax, Mmax, 1, &ainfo);
    }

    HealpixFourierTransform(long nSide, long Lmax, long Mmax, int iterate, const std::valarray<double>& weights )
      : realMap(nSide), fourierMap(Lmax, Mmax), ainfo(0), ginfo(0), m_iterate(iterate)
    {
      sharp_make_weighted_healpix_geom_info (nSide, 1, &weights[0], &ginfo);
      sharp_make_triangular_alm_info (Lmax, Mmax, 1, &ainfo);
    }

    virtual ~HealpixFourierTransform()
    {
      sharp_destroy_geom_info(ginfo);
      sharp_destroy_alm_info(ainfo);
    }

    virtual const FourierMap<std::complex<T> >& fourierSpace() const { return fourierMap; }

    virtual FourierMap<std::complex<T> >& fourierSpace() { return fourierMap; }

    virtual const FourierMap<T>& realSpace() const { return realMap; }

    virtual FourierMap<T>& realSpace() { return realMap; }
    
    virtual FourierTransform<T> *mimick() const 
    {
      return new HealpixFourierTransform<T>(realMap.Nside(), fourierMap.Lmax(), fourierMap.Mmax());
    }

    virtual void analysis()
    {
      void *aptr=reinterpret_cast<void *>(fourierMap.data()), *mptr=reinterpret_cast<void *>(realMap.data());

      sharp_execute (SHARP_MAP2ALM, 0, &aptr, &mptr, ginfo, ainfo, 1,
        HealpixJobHelper__<T>::val,0,0);
      for (int i = 0; i < m_iterate; i++)
        {
          HealpixFourierMap<T> tmp_map(realMap.Nside());          
          void *tmp_ptr=reinterpret_cast<void *>(tmp_map.data());
          typename HealpixFourierMap<T>::MapType m0 = tmp_map.eigen();
          typename HealpixFourierMap<T>::MapType m1 = realMap.eigen();

          sharp_execute (SHARP_ALM2MAP, 0, &aptr, &tmp_ptr, ginfo, ainfo, 1,
            HealpixJobHelper__<T>::val,0,0);
          m0 = m1 - m0;
          sharp_execute (SHARP_MAP2ALM, 0, &aptr, &tmp_ptr, ginfo, ainfo, 1,
            HealpixJobHelper__<T>::val | SHARP_ADD,0,0);
        }
    }

    virtual void synthesis()
    {
      void *aptr=reinterpret_cast<void *>(fourierMap.data()), *mptr=reinterpret_cast<void *>(realMap.data());

      sharp_execute (SHARP_ALM2MAP, 0, &aptr, &mptr, ginfo, ainfo, 1,
        HealpixJobHelper__<T>::val,0,0);
    }

    virtual void analysis_conjugate()
    {
      synthesis();
      realMap.scale(4*M_PI/realMap.size());
    }

    virtual void synthesis_conjugate()
    {
      analysis();
      fourierMap.scale(realMap.size()/(4*M_PI));
    }
    
  };

};

#endif
