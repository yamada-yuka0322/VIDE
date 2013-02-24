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
	double A_k = std::sqrt(0.5*V*f(rand_map.get_K(p)));
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

    for (long p = 1; p < plane_size/2; p++)
      {
	long q = (p+1)*dims[0]-1;
	long q2 = (plane_size-p+1)*dims[0]-1;
	assert(q < plane_size*dims[0]);
	assert(q2 < plane_size*dims[0]);
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
      d[p] *= f(m_c.get_K(p));
  }

  template<typename T>
  void EuclidianSpectrum_1D<T>::mul_sqrt(FourierMapType& m) const
  {
    EuclidianFourierMapComplex<T>& m_c = dynamic_cast<EuclidianFourierMapComplex<T>&>(m);
    std::complex<T> *d = m.data();

    for (long p = 0; p < m_c.size(); p++)
      d[p] *= std::sqrt(f(m_c.get_K(p)));
  }

  template<typename T>
  void EuclidianSpectrum_1D<T>::mul_inv(FourierMapType& m) const
  {
    EuclidianFourierMapComplex<T>& m_c = dynamic_cast<EuclidianFourierMapComplex<T>&>(m);
    std::complex<T> *d = m.data();

    for (long p = 0; p < m_c.size(); p++)
     {
        T A = f(m_c.get_K(p));
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
        T A = std::sqrt(f(m_c.get_K(p)));
        if (A == 0)
          d[p] = 0;
        else
          d[p] /= A;
      }
  }

};


#endif
