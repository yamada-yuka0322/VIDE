#ifndef __COSMOTOOL_FOURIER_EUCLIDIAN_HPP
#define __COSMOTOOL_FOURIER_EUCLIDIAN_HPP

#include <vector>
#include <boost/shared_ptr.hpp>
#include "base_types.hpp"
#include "fft/fftw_calls.hpp"

namespace CosmoTool
{
  
  template<typename T>
  class EuclidianFourierMap: public FourierMap<T>
  {
  private:
    boost::shared_ptr<T> m_data;
    std::vector<int> m_dims;
    long m_size;
  public:
    EuclidianFourierMap(boost::shared_ptr<T> indata, std::vector<int> indims)
    {
      m_data = indata;
      m_dims = indims;
      m_size = 1;
      for (int i = 0; i < m_dims.size(); i++)
	m_size *= m_dims[i];
    }

    virtual ~EuclidianFourierMap()
    {
    }

    virtual const T *data() const { return m_data.get(); }
    virtual T *data() { return m_data.get(); }
    virtual long size() const  { return m_size; } 
    
    virtual FourierMap<T> *copy() const
    {
      FourierMap<T> *m = this->mimick();
      m->eigen() = this->eigen();
      return m;
    }

    virtual FourierMap<T> *mimick() const
    {
      return new EuclidianFourierMap<T>(
		      boost::shared_ptr<T>((T *)fftw_malloc(sizeof(T)*m_size), 
					   std::ptr_fun(fftw_free)),
		       m_dims);
    }
  };

  
  template<typename T>
  class EuclidianFourierTransform: public FourierTransform<T>
  {
  private:
    typedef FFTW_Calls<T> calls;
    EuclidianFourierMap<T> *realMap;
    EuclidianFourierMap<std::complex<T> > *fourierMap;
    typename calls::plan_type m_analysis, m_synthesis;
    double volume;
    long N;
    std::vector<int> m_dims;
    std::vector<double> m_L;
  public:
    EuclidianFourierTransform(const std::vector<int>& dims, const std::vector<double>& L)
    {
      assert(L.size() == dims.size());

      m_dims = dims;
      m_L = L;

      N = 1;
      volume = 1;
      for (int i = 0; i < dims.size(); i++)
	{
	  N *= dims[i];
	  volume *= L[i];
	}

      realMap = new EuclidianFourierMap<T>(
		   boost::shared_ptr<T>(calls::alloc_real(N),
					std::ptr_fun(calls::free)),
		   dims);
      fourierMap = new EuclidianFourierMap<std::complex<T> >(
		   boost::shared_ptr<std::complex<T> >((std::complex<T>*)calls::alloc_complex(N),
						       std::ptr_fun(calls::free)), 
		   dims);
      m_analysis = calls::plan_dft_r2c(dims.size(), &dims[0],
				       realMap->data(), (typename calls::complex_type *)fourierMap->data(),
				       FFTW_MEASURE);
      m_synthesis = calls::plan_dft_c2r(dims.size(), &dims[0],
					(typename calls::complex_type *)fourierMap->data(), realMap->data(),
					FFTW_MEASURE);
    }
    
    virtual ~EuclidianFourierTransform()
    {
      delete realMap;
      delete fourierMap;
      calls::destroy_plan(m_synthesis);
      calls::destroy_plan(m_analysis);
    }

    void synthesis()
    {
      calls::execute(m_synthesis);
      realMap->scale(1/volume);
    }
    
    void analysis()
    {
      calls::execute(m_analysis);
      fourierMap->scale(volume/N);
    }
    
    void synthesis_conjugate()
    {
      calls::execute(m_analysis);
      fourierMap->scale(1/volume);
    }
    
    void analysis_conjugate()
    {
      calls::execute(m_synthesis);
      realMap->scale(volume/N);
    }

    const FourierMap<std::complex<T> >& fourierSpace() const
    {
      return *fourierMap;
    }

    FourierMap<std::complex<T> >& fourierSpace()
    {
      return *fourierMap;
    }

    const FourierMap<T>& realSpace() const
    {
      return *realMap;
    }

    FourierMap<T>& realSpace()
    {
      return *realMap;
    }

    FourierTransform<T> *mimick() const
    {
      return new EuclidianFourierTransform(m_dims, m_L);
    }
  };

  template<typename T>
  class EuclidianFourierTransform_2d: public EuclidianFourierTransform<T>
  {
  private:
    template<typename T2>
    static std::vector<T2> make_2d_vector(T2 a, T2 b)
    {
      T2 arr[2] = { a, b};
      return std::vector<T2>(&arr[0], &arr[2]);
    }
  public:
    EuclidianFourierTransform_2d(int Nx, int Ny, double Lx, double Ly)
      : EuclidianFourierTransform<T>(make_2d_vector<int>(Nx, Ny), make_2d_vector<double>(Lx, Ly))
    {
    }

    virtual ~EuclidianFourierTransform_2d() {}

  };

  template<typename T>
  class EuclidianFourierTransform_3d: public EuclidianFourierTransform<T>
  {
  private:
    template<typename T2>
    static std::vector<T2> make_3d_vector(T2 a, T2 b, T2 c)
    {
      T2 arr[2] = { a, b, c};
      return std::vector<T2>(&arr[0], &arr[3]);
    }

  public:
    EuclidianFourierTransform_3d(int Nx, int Ny, int Nz, double Lx, double Ly, double Lz)
      : EuclidianFourierTransform<T>(make_3d_vector<int>(Nx, Ny, Nz), make_3d_vector<double>(Lx, Ly, Lz))
    {
    }

    virtual ~EuclidianFourierTransform_3d() {}
  };



};

#endif
