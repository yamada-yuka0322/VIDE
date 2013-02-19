#ifndef __DETAILS_EUCLIDIAN_SPECTRUM_1D_BIN
#define __DETAILS_EUCLIDIAN_SPECTRUM_1D_BIN

#include <boost/bind.hpp>
#include <cmath>

namespace CosmoTool
{
  template<typename T>
  class EuclidianSpectrum_1D_Binned: public EuclidianSpectrum_1D<T>
  {
  protected:
    T *m_data;
    long m_size;
    T m_kmin, m_kmax, m_logdeltak;
    std::vector<T> m_dof;
  public:
    typedef typename SpectrumFunction<T>::FourierMapType FourierMapType;
    typedef typename SpectrumFunction<T>::SpectrumFunctionPtr SpectrumFunctionPtr;
    typedef boost::shared_ptr<FourierMapType> ptr_map;

    T interpolate_spectrum(T k)
    {
      T q = std::log(k/m_kmin)/m_logdeltak;
      long ik = std::floor(q);

      if (ik >= m_size-1)
	return m_data[m_size-1];
      else if (ik < 0)
	return k/m_kmin*m_data[0];

      return std::exp((q-ik)*m_data[ik+1] + (1-(q-ik))*m_data[ik]);
    }

    EuclidianSpectrum_1D_Binned(int numBin, T kmin, T kmax)
      : EuclidianSpectrum_1D<T>(boost::bind(&EuclidianSpectrum_1D_Binned::interpolate_spectrum, this, _1))
    {
      m_data = new T[numBin];
      m_kmin = kmin;
      m_kmax = kmax;
      m_size = numBin;
      m_logdeltak = std::log(m_kmax/m_kmin);
    }

    SpectrumFunctionPtr copy() const
    {
      EuclidianSpectrum_1D_Binned *s = new EuclidianSpectrum_1D_Binned(m_size, m_kmin, m_kmax);
      std::copy(m_data, m_data+m_size, s->m_data);
      return SpectrumFunctionPtr(s);
    }

    void set_dof(std::vector<T>& dof_array) 
    {
      assert(m_size == dof_array.size());
      m_dof = dof_array;
    }

    const T *data() const { return m_data; }
    long size() const { return m_size; }
    const T *dof() const { return &m_dof[0]; }
  };

};

#endif
