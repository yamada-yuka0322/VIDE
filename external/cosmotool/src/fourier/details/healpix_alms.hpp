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
