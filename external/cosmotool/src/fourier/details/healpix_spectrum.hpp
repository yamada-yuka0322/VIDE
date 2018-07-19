/*+
This is CosmoTool (./src/fourier/details/healpix_spectrum.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
            c.real(gsl_ran_gaussian(rng, Al)); 
            c.imag(gsl_ran_gaussian(rng, Al)); 
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
