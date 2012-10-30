#ifndef __COSMOTOOL_FOURIER_HEALPIX_HPP
#define __COSMOTOOL_FOURIER_HEALPIX_HPP

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <alm.h>
#include <healpix_base.h>
#include <psht_cxx.h>

namespace CosmoTool
{
    
  template<typename T>
  class HealpixFourierMap: public FourierMap<T>, public Healpix_Base
  {
  private:
    T *m_data;
    Eigen::aligned_allocator<T> alloc;
  public:
    HealpixFourierMap(long nSide)
      : Healpix_Base(RING, nSide, SET_NSIDE)
    {
      m_data = alloc.allocate(Npix);
    }

    virtual ~HealpixFourierMap()
    {
      alloc.deallocate(m_data);
    }

    virtual const T* data() const { return m_data; }
    virtual T *data() { return m_data; }
    virtual long size() const { return Npix(); }
    
    virtual FourierMap<T> *mimick() const
    {
      return new HealpixFourierMap<T>(Nside());
    }
  };


  template<typename T>
  class HealpixFourierALM: public FourierMap<std::complex<T> >, public Alm_Base
  {
  private:
    std::complex<T> *alms;
    long size;
    Eigen::aligned_allocator<std::complex<T> > alloc;
  public:
    HealpixFourierALM(long Lmax, long Mmax)
      : Alm_Base(Lmax, Mmax)
    {
      size = Num_Alms(Lmax, Mmax);
      alms = alloc.allocate(size);
    }

    virtual ~HealpixFourierALM()
    {
      alloc.deallocate(alms);
    }

    virtual const T* data() const { return alms; }
    virtual T * data() { return alms;} 
    virtual long size() const { return size; }

    virtual FourierMap<T> *mimick() const
    {
      return new HealpixFourierALM<T>(Lmax(), Mmax());
    }
  };

  template<typename T>
  class HealpixFourierTransform: public FourierTransform<T>
  {
  private:
    HealpixFourierMap<T> realMap;
    HealpixFourierALM<T> fourierMap;
    psht_joblist<T> jobs;
  public:
    HealpixFourierTransform(long nSide, long Lmax, long Mmax)
      : realMap(nSide), fourierMap(Lmax, Mmax)
    {
      jobs.set_Healpix_geometry(nSide);
      jobs.set_triangular_alm_info(Lmax, Mmax);
    }

    virtual ~HealpixFourierTransform() {}

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
      jobs.add_map2alm(realMap.data(), 
		       reinterpret_cast<xcomplex<T> *>(fourierMap.data()), 
		       false);
      jobs.execute();
      jobs.clear_jobs();
    }

    virtual void synthesis()
    {
      jobs.add_alm2map(reinterpret_cast<xcomplex<T> *>(fourierMap.data()), 
		       realMap.data(),
		       false);
      jobs.execute();
      jobs.clear_jobs();
    }

    virtual void analysis_conjugate()
    {
      synthesis();
      realMap.scale(4*M_PI/realMap.Npix());
    }

    virtual void synthesis_conjugate()
    {
      analysis();
      fourierMap.scale(realMap.Npix()/(4*M_PI));
    }
    
  };
};

#endif
