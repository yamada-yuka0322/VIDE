#ifndef __BASE_FOURIER_TYPES_HPP
#define __BASE_FOURIER_TYPES_HPP

#include <string>
#include <Eigen/Dense>
#include <complex>

namespace CosmoTool
{
  template<typename T>
  class FourierMap
  {
  protected:
    FourierMap() {}
    
  public:
    typedef Eigen::Matrix<T, 1, Eigen::Dynamic> VecType;
    typedef Eigen::Map<VecType, Eigen::Aligned> MapType;
    typedef Eigen::Map<VecType const, Eigen::Aligned> ConstMapType;

    virtual ~FourierMap() {}

    virtual const T* data() const = 0;
    virtual T* data()  = 0;

    virtual long size() const = 0;

    MapType eigen()
    {
      return MapType(data(), size());
    }

    ConstMapType eigen() const
    {
      return ConstMapType(data(), size());
    }
 
    void scale(const T& factor)
    {
      MapType m(data(), size());
      m *= factor;
    }

    void scale(const FourierMap<T> *map2)
    {
      assert(size() == map2->size());
      MapType m(data(), size());
      MapType m2(map2->data(), map2->size());
      m *= m2;
    }

    void add(const T&  factor)
    {
      eigen() += factor;
    }

    void add(const FourierMap<T> *map2)
    {
      assert(size() == map2->size());
      MapType m(data(), size());
      MapType m2(map2->data(), map2->size());

      eigen() += map2->eigen();
    }

    virtual FourierMap<T> *copy() const
    {
      FourierMap<T> *m = this->mimick();
      
      m->eigen() = this->eigen();
      return m;
    }

    virtual FourierMap<T> *mimick() const = 0;
  };

  template<typename T>
  class FourierTransform
  {
  protected:
    FourierTransform() {}
  public:
    virtual ~FourierTransform() { }
    
    virtual const FourierMap<std::complex<T> >& fourierSpace() const = 0;
    virtual FourierMap<std::complex<T> >& fourierSpace() = 0;

    virtual const FourierMap<T>& realSpace() const = 0;
    virtual FourierMap<T>& realSpace() = 0;

    virtual FourierTransform<T> *mimick() const = 0;

    virtual void analysis() = 0;
    virtual void synthesis() = 0;
    virtual void analysis_conjugate() = 0;
    virtual void synthesis_conjugate() = 0;
  };


};

#endif
