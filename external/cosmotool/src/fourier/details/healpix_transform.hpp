#ifndef __COSMOTOOL_FOURIER_HEALPIX_DETAILS_TRANSFORM_HPP
#define __COSMOTOOL_FOURIER_HEALPIX_DETAILS_TRANSFORM_HPP

namespace CosmoTool
{

  template<typename T> struct HealpixJobHelper__ {};

  template<> struct HealpixJobHelper__<double>
    { enum {val=1}; };

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

      sharp_execute (SHARP_MAP2ALM, 0, 0, &aptr, &mptr, ginfo, ainfo, 1,
        HealpixJobHelper__<T>::val,0,0,0);
      for (int i = 0; i < m_iterate; i++)
        {
          HealpixFourierMap<T> tmp_map(realMap.Nside());          
          void *tmp_ptr=reinterpret_cast<void *>(tmp_map.data());
          typename HealpixFourierMap<T>::MapType m0 = tmp_map.eigen();
          typename HealpixFourierMap<T>::MapType m1 = realMap.eigen();

          sharp_execute (SHARP_ALM2MAP, 0, 0, &aptr, &tmp_ptr, ginfo, ainfo, 1,
            HealpixJobHelper__<T>::val,0,0,0);
          m0 = m1 - m0;
          sharp_execute (SHARP_MAP2ALM, 0, 1, &aptr, &tmp_ptr, ginfo, ainfo, 1,
            HealpixJobHelper__<T>::val,0,0,0);
        }
    }

    virtual void synthesis()
    {
      void *aptr=reinterpret_cast<void *>(fourierMap.data()), *mptr=reinterpret_cast<void *>(realMap.data());

      sharp_execute (SHARP_ALM2MAP, 0, 0, &aptr, &mptr, ginfo, ainfo, 1,
        HealpixJobHelper__<T>::val,0,0,0);
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
