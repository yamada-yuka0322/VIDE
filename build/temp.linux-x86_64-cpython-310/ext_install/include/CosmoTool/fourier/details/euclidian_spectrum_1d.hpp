/*+
This is CosmoTool (./src/fourier/details/euclidian_spectrum_1d.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __DETAILS_EUCLIDIAN_SPECTRUM_1D
#define __DETAILS_EUCLIDIAN_SPECTRUM_1D

#include  <iostream>
#include <boost/function.hpp>


namespace CosmoTool
{
  template<typename T>
  class EuclidianOperator
  {
  public:
    typedef boost::function1<T, T> Function;
    
    Function base, op;
    T operator()(T k) {
      return op(base(k));
    }
  };

  template<typename T>
  class EuclidianSpectrum_1D: public SpectrumFunction<T>
  {
  public:
    typedef boost::function1<T, T> Function;
  protected:
    Function f;

    static T msqrt(T a) { return std::sqrt(a); }
  public:
    typedef typename SpectrumFunction<T>::FourierMapType FourierMapType;
    typedef typename SpectrumFunction<T>::SpectrumFunctionPtr SpectrumFunctionPtr;
    typedef boost::shared_ptr<FourierMapType> ptr_map;

    EuclidianSpectrum_1D(Function P)
      : f(P)
    {
    }

    void newRandomFourier(gsl_rng *rng, FourierMapType& out_map) const;

    SpectrumFunctionPtr copy() const {
      return SpectrumFunctionPtr(new EuclidianSpectrum_1D(f));
    }

    void sqrt() {
      EuclidianOperator<T> o;
      o.base = f;
      o.op = &EuclidianSpectrum_1D<T>::msqrt;
      f = (Function(o));
    }

    void mul(FourierMapType& m) const;
    void mul_sqrt(FourierMapType& m) const;
    void mul_inv(FourierMapType& m) const;
    void mul_inv_sqrt(FourierMapType& m) const;
  };


  template<typename T>
  void EuclidianSpectrum_1D<T>::newRandomFourier(gsl_rng *rng, FourierMapType& out_map) const
  {
    typedef EuclidianFourierMapComplex<T> MapT;
    typedef typename EuclidianSpectrum_1D<T>::ptr_map ptr_map;
    typedef typename MapT::DimArray DimArray;

    MapT& rand_map = dynamic_cast<MapT&>(out_map);

    std::complex<T> *d = rand_map.data();
    long idx;
    const DimArray& dims = rand_map.getDims();
    const std::vector<double>& delta_k = rand_map.get_delta_k();
    long plane_size;
    bool alleven = rand_map.allDimensionsEven();
    double V = 1;

    for (int p = 0; p < delta_k.size(); p++)
      V *= (2*M_PI/delta_k[p]);

    for (long p = 1; p < rand_map.size(); p++)
      {
	double A_k = std::sqrt(0.5*V*f(rand_map.get_K_p(p)));
	d[p] = std::complex<T>(gsl_ran_gaussian(rng, A_k),
			       gsl_ran_gaussian(rng, A_k));
      }
    // Generate the mean value
    d[0] = std::complex<T>(gsl_ran_gaussian(rng, std::sqrt(V*f(0))), 0);

    if (!rand_map.firstDimensionEven())
      return;

    // Correct the Nyquist plane
    idx = dims[0]-1; // Stick to the last element of the first dimension
    d[idx] = std::complex<T>(d[idx].real() + d[idx].imag(), 0);
    // 1D is special case
    if (dims.size() == 1)
      return;
    
    plane_size = 1;
    for (int q = 1; q < dims.size(); q++)
      {
	plane_size *= dims[q];
      }

    for (long p = 1; p < plane_size/2+1; p++)
      {
	long q = (p+1)*dims[0]-1;
	long q2 = (plane_size-p+1)*dims[0]-1;
	assert(q < plane_size*dims[0]);
	assert(q2 < plane_size*dims[0]);
	d[q] = conj(d[q2]);
      }

    for (long p = 1; p < plane_size/2+1; p++)
      {
	long q = (p)*dims[0];
	long q2 = (plane_size-p)*dims[0];
	d[q] = conj(d[q2]);
      }

    if (alleven)
      {
	long q = 0;
	for (int i = dims.size()-1; i >= 1; i--)
	  q = dims[i]*q + dims[i]/2;
	q += dims[0]-1;
	d[q] = std::complex<T>(d[q].real()+d[q].imag(),0);
      }
  }

  template<typename T>
  void EuclidianSpectrum_1D<T>::mul(FourierMapType& m) const
  {
    EuclidianFourierMapComplex<T>& m_c = dynamic_cast<EuclidianFourierMapComplex<T>&>(m);
    std::complex<T> *d = m.data();

    for (long p = 0; p < m_c.size(); p++)
      d[p] *= f(m_c.get_K_p(p));
  }

  template<typename T>
  void EuclidianSpectrum_1D<T>::mul_sqrt(FourierMapType& m) const
  {
    EuclidianFourierMapComplex<T>& m_c = dynamic_cast<EuclidianFourierMapComplex<T>&>(m);
    std::complex<T> *d = m.data();

    for (long p = 0; p < m_c.size(); p++)
      d[p] *= std::sqrt(f(m_c.get_K_p(p)));
  }

  template<typename T>
  void EuclidianSpectrum_1D<T>::mul_inv(FourierMapType& m) const
  {
    EuclidianFourierMapComplex<T>& m_c = dynamic_cast<EuclidianFourierMapComplex<T>&>(m);
    std::complex<T> *d = m.data();

    for (long p = 0; p < m_c.size(); p++)
     {
        T A = f(m_c.get_K_p(p));
        if (A==0)
          d[p] = 0;
        else
          d[p] /= A;
     }
  }

  template<typename T>
  void EuclidianSpectrum_1D<T>::mul_inv_sqrt(FourierMapType& m) const
  {
    EuclidianFourierMapComplex<T>& m_c = dynamic_cast<EuclidianFourierMapComplex<T>&>(m);
    std::complex<T> *d = m.data();

    for (long p = 0; p < m_c.size(); p++)
      {
        T A = std::sqrt(f(m_c.get_K_p(p)));
        if (A == 0)
          d[p] = 0;
        else
          d[p] /= A;
      }
  }

};


#endif
