/*+
This is CosmoTool (./src/fourier/details/euclidian_transform.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __DETAILS_EUCLIDIAN_TRANSFORM
#define __DETAILS_EUCLIDIAN_TRANSFORM

namespace CosmoTool
{
  
  template<typename T>
  class EuclidianFourierTransform: public FourierTransform<T>
  {
  public:
    typedef typename EuclidianFourierMapBase<T>::DimArray DimArray;
  private:
    typedef FFTW_Calls<T> calls;
    EuclidianFourierMapReal<T> *realMap;
    EuclidianFourierMapComplex<T> *fourierMap;
    typename calls::plan_type m_analysis, m_synthesis;
    double volume;
    long N, Nc;
    DimArray m_dims, m_dims_hc;
    std::vector<double> m_L;
  public:
    EuclidianFourierTransform(const DimArray& dims, const std::vector<double>& L)
    {
      realMap = 0;
      create_plan(dims, L);
    }

    void create_plan(const DimArray& dims, const std::vector<double>& L)
    {
      assert(L.size() == dims.size());
      std::vector<double> dk(L.size());
      std::vector<int> swapped_dims(dims.size());
 
      if (realMap != 0)
        {
           delete realMap;
           delete fourierMap;
           calls::destroy_plan(m_synthesis);
           calls::destroy_plan(m_analysis);
        }
 
      m_dims = dims;
      m_dims_hc = dims;
      m_dims_hc[0] = dims[0]/2+1;
      m_L = L;

      N = 1;
      Nc = 1;
      volume = 1;
      for (int i = 0; i < dims.size(); i++)
        {
          N *= dims[i];
          Nc *= m_dims_hc[i];
          volume *= L[i];
          dk[i] = 2*M_PI/L[i];
          swapped_dims[dims.size()-1-i] = dims[i];
        }

      realMap = new EuclidianFourierMapReal<T>(
                boost::shared_ptr<T>(calls::alloc_real(N),
                                    std::ptr_fun(calls::free)),
                                    m_dims);
      fourierMap = new EuclidianFourierMapComplex<T>(
                boost::shared_ptr<std::complex<T> >(
                  (std::complex<T>*)calls::alloc_complex(Nc),
                  std::ptr_fun(calls::free)), 
                dims[0], m_dims_hc, dk);
      {
        m_analysis = calls::plan_dft_r2c(
              dims.size(), &swapped_dims[0],
              realMap->data(), 
              (typename calls::complex_type *)fourierMap->data(),
              FFTW_DESTROY_INPUT|FFTW_MEASURE);
        m_synthesis = calls::plan_dft_c2r(
              dims.size(), &swapped_dims[0],
              (typename calls::complex_type *)fourierMap->data(), realMap->data(),
              FFTW_DESTROY_INPUT|FFTW_MEASURE);
      }
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

    void synthesis_unnormed()
    {
      calls::execute(m_synthesis);
    }
    
    void analysis()
    {
      calls::execute(m_analysis);
      fourierMap->scale(volume/N);
    }

    void analysis_unnormed()
    {
      calls::execute(m_analysis);
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
  class EuclidianFourierTransform_1d: public EuclidianFourierTransform<T>
  {
  private:
    template<typename T2>
    static std::vector<T2> make_1d_vector(T2 a)
    {
      T2 arr[2] = { a};
      return std::vector<T2>(&arr[0],&arr[1]);
    }
  public:
    EuclidianFourierTransform_1d(int Nx, double Lx)
      : EuclidianFourierTransform<T>(make_1d_vector<int>(Nx), make_1d_vector<double>(Lx))
    {
    }

    virtual ~EuclidianFourierTransform_1d() {}

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
      T2 arr[3] = { a, b, c};
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
