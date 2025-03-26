/*+
This is CosmoTool (./src/fourier/details/healpix_alms.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __COSMOTOOL_FOURIER_HEALPIX_DETAILS_ALM_HPP
#define __COSMOTOOL_FOURIER_HEALPIX_DETAILS_ALM_HPP

namespace CosmoTool
{
  template<typename T>
  class HealpixFourierALM: public FourierMap<std::complex<T> >
  {
  private:
    std::complex<T> *alms;
    long m_size;
    long Lmax_, Mmax_, TVal_;
    Eigen::aligned_allocator<std::complex<T> > alloc;
  public:
    typedef unsigned long LType;

    LType Lmax() const { return Lmax_; }
    LType Mmax() const { return Mmax_; }

    LType Num_Alms() const
    {
      return ((Mmax_+1)*(Mmax_+2))/2 + (Mmax_+1)*(Lmax_-Mmax_);
    }

    LType index_l0(LType m) const
    {
      return ((m*(TVal_-m))/2);
    }

    LType index(LType l, LType m) const
    {
      return index_l0(m) + l;
    }
    
    HealpixFourierALM(LType lmax, LType mmax)
      : Lmax_(lmax), Mmax_(mmax), TVal_(2*lmax+1)
    {
      m_size = Num_Alms();
      alms = alloc.allocate(m_size);
    }

    virtual ~HealpixFourierALM()
    {
      alloc.deallocate(alms, m_size);
    }

    virtual const std::complex<T>* data() const { return alms; }
    virtual std::complex<T> * data() { return alms;} 
    virtual long size() const { return m_size; }

    virtual FourierMap<std::complex<T> > *mimick() const
    {
      return new HealpixFourierALM<T>(Lmax_, Mmax_);
    }

    virtual std::complex<T> dot_product(const FourierMap<std::complex<T> >& other) const
      throw(std::bad_cast)
    {
      const HealpixFourierALM<T>& mfm = dynamic_cast<const HealpixFourierALM<T>&>(other);
      typedef typename FourierMap<std::complex<T> >::MapType MapType;
      std::complex<T> S;

      if (m_size != mfm.m_size)
        throw std::bad_cast();

      MapType m1(alms, m_size);
      MapType m2(mfm.alms, mfm.m_size);
      
      S = (m1.block(0,0,1,Lmax_+1).conjugate() * m2.block(0,0,1,Lmax_+1)).sum();
      S += std::complex<T>(2,0)*(m1.block(0,1+Lmax_,1,m_size-1-Lmax_).conjugate() * m2.block(0,1+Lmax_,1,m_size-1-Lmax_)).sum();
      return S;
    }
  };
};

#endif
