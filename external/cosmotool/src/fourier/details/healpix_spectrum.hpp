#ifndef __COSMOTOOL_FOURIER_HEALPIX_DETAILS_SPECTRUM_HPP
#define __COSMOTOOL_FOURIER_HEALPIX_DETAILS_SPECTRUM_HPP

namespace CosmoTool
{
  template<typename T>
  class HealpixSpectrum: public SpectrumFunction<T>
  {
  protected:
    std::vector<T> cls;
    int *m_dof;
  public:
    typedef typename SpectrumFunction<T>::FourierMapType FourierMapType;
    typedef boost::shared_ptr<FourierMapType> ptr_map;
    typedef typename SpectrumFunction<T>::SpectrumFunctionPtr SpectrumFunctionPtr;
    typedef unsigned long LType;

    HealpixSpectrum(LType Lmax)
      : cls(Lmax+1), m_dof(new int[Lmax+1]) 
    {
      for (int l = 0; l <= Lmax; l++)
	m_dof[l] = l + 1;
    }

    T *data() { return &cls[0]; }
    const T *data() const { return &cls[0]; }
    const int *dof() const { return m_dof; } 
    
    void set_dof(LType l, int dof) { m_dof[l] = dof; }

    LType Lmax() const { return cls.size()-1; }
    long size() const { return cls.size(); }

    void newRandomFourier(gsl_rng *rng, FourierMapType& like_map) const; 

    SpectrumFunctionPtr copy() const {
      HealpixSpectrum<T> *s = new HealpixSpectrum<T>(Lmax());
      s->cls = cls;
      return SpectrumFunctionPtr(s);
    }

    void sqrt() {
      std::transform(cls.begin(), cls.end(), cls.begin(), std::ptr_fun<T,T>(std::sqrt));
    }

    void mul(FourierMapType& m) const;
    void mul_sqrt(FourierMapType& m) const;
    void mul_inv(FourierMapType& m) const;
    void mul_inv_sqrt(FourierMapType& m) const;
  };

  template<typename T>
  void HealpixSpectrum<T>::newRandomFourier(gsl_rng *rng, FourierMapType& out_map) const
  { 
    HealpixFourierALM<T>& alms = dynamic_cast<HealpixFourierALM<T>&>(out_map);
    long lmaxGen = std::min(cls.size()-1, alms.Lmax());
    std::complex<T> *new_data = alms.data();

    for (LType l = 0; l <= lmaxGen; l++)
      {
        double Al = std::sqrt(cls[l]);

        new_data[alms.index(l,0)] = gsl_ran_gaussian(rng, Al);
        Al *= M_SQRT1_2;
        for (LType m = 1; m <= std::min(l,alms.Mmax()); m++)
          {
            std::complex<T>& c = new_data[alms.index(l,m)];
            c.real() = gsl_ran_gaussian(rng, Al); 
            c.imag() = gsl_ran_gaussian(rng, Al); 
          }
      }
  } 

 template<typename T>
 void HealpixSpectrum<T>::mul(FourierMapType& like_map) const
  {
    HealpixFourierALM<T>& alms = dynamic_cast<HealpixFourierALM<T>&>(like_map);
    std::complex<T> *data = alms.data();

    for (LType l = 0; l <= alms.Lmax(); l++)
      {
        double Al = cls[l];

        for (LType m = 0; m <= std::min(l,alms.Mmax()); m++)
          {
            data[alms.index(l,m)] *= Al;
          }
      }
  }

 template<typename T>
 void HealpixSpectrum<T>::mul_sqrt(FourierMapType& like_map) const
  {
    HealpixFourierALM<T>& alms = dynamic_cast<HealpixFourierALM<T>&>(like_map);
    std::complex<T> *data = alms.data();

    for (LType l = 0; l <= alms.Lmax(); l++)
      {
        double Al = std::sqrt(cls[l]);

        for (LType m = 0; m <= std::min(l,alms.Mmax()); m++)
          {
            data[alms.index(l,m)] *= Al;
          }
      }
  }

 template<typename T>
 void HealpixSpectrum<T>::mul_inv(FourierMapType& like_map) const
  {
    HealpixFourierALM<T>& alms = dynamic_cast<HealpixFourierALM<T>&>(like_map);
    std::complex<T> *data = alms.data();

    for (LType l = 0; l <= alms.Lmax(); l++)
      {
        double Al = (cls[l] <= 0) ? 0 : (1/cls[l]);

        for (LType m = 0; m <= std::min(l,alms.Mmax()); m++)
          {
            data[alms.index(l,m)] *= Al;
          }
      }
  }

 template<typename T>
 void HealpixSpectrum<T>::mul_inv_sqrt(FourierMapType& like_map) const
  {
    HealpixFourierALM<T>& alms = dynamic_cast<HealpixFourierALM<T>&>(like_map);
    std::complex<T> *data = alms.data();

    for (LType l = 0; l <= alms.Lmax(); l++)
      {
        double Al = (cls[l] <= 0) ? 0 : std::sqrt(1/cls[l]);

        for (LType m = 0; m <= std::min(l,alms.Mmax()); m++)
          {
            data[alms.index(l,m)] *= Al;
          }
      }
  }

};

#endif
