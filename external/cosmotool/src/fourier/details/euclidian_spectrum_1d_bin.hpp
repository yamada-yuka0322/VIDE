/*+
This is CosmoTool (./src/fourier/details/euclidian_spectrum_1d_bin.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
