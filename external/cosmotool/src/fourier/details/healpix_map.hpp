#ifndef __COSMOTOOL_FOURIER_HEALPIX_DETAILS_MAP_HPP
#define __COSMOTOOL_FOURIER_HEALPIX_DETAILS_MAP_HPP

namespace CosmoTool
{

  template<typename T>
  class HealpixFourierMap: public FourierMap<T>
  {
  private:
    T *m_data;
    long Npix, m_Nside;
    Eigen::aligned_allocator<T> alloc;
  public:
    HealpixFourierMap(long nSide)
      : Npix(12*nSide*nSide), m_Nside(nSide)
    {
      m_data = alloc.allocate(Npix);
    }

    virtual ~HealpixFourierMap()
    {
      alloc.deallocate(m_data, Npix);
    }

    long Nside() const { return m_Nside; }
    virtual const T* data() const { return m_data; }
    virtual T *data() { return m_data; }
    virtual long size() const { return Npix; }

    virtual T dot_product(const FourierMap<T>& other) const
      throw(std::bad_cast)
    {
      typedef typename FourierMap<T>::MapType MapType;

      const HealpixFourierMap<T>& mfm = dynamic_cast<const HealpixFourierMap<T>&>(other);
      if (Npix != mfm.size())
        throw std::bad_cast();

      MapType m1(m_data, Npix);
      MapType m2(mfm.m_data, mfm.Npix);
      
      return (m1*m2).sum();
    }
    
    virtual FourierMap<T> *mimick() const
    {
      return new HealpixFourierMap<T>(m_Nside);
    }
  };

};

#endif
