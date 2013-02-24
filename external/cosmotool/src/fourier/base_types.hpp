#ifndef __BASE_FOURIER_TYPES_HPP
#define __BASE_FOURIER_TYPES_HPP

#include <gsl/gsl_rng.h>
#include <boost/shared_ptr.hpp>
#include <string>
#include <Eigen/Dense>
#include <complex>
#include <exception>

namespace CosmoTool
{
  class IncompatibleMap: virtual std::exception {};

  template<typename T> class FourierMap;

  template<typename T>
  class SpectrumFunction
  {
  protected:
    SpectrumFunction() {}
  public:
    typedef Eigen::Array<T, 1, Eigen::Dynamic> VecType;
    typedef Eigen::Map<VecType, Eigen::Aligned> MapType;
    typedef Eigen::Map<VecType const, Eigen::Aligned> ConstMapType;
    typedef FourierMap<std::complex<T> > FourierMapType;
    typedef boost::shared_ptr<FourierMapType> FourierMapPtr;
    typedef boost::shared_ptr<SpectrumFunction<T> > SpectrumFunctionPtr;

    virtual ~SpectrumFunction() {}

    virtual void 
       newRandomFourier(gsl_rng *rng, FourierMapType& out_map) const = 0;

    virtual SpectrumFunctionPtr copy() const = 0;

    virtual const T *data() const { return 0; }
    virtual const int *dof() const { return 0; }
    virtual long size() const { return -1; }

    virtual void sqrt() = 0;

    virtual void mul(FourierMapType& m) const = 0;
    virtual void mul_sqrt(FourierMapType& m) const = 0;
    virtual void mul_inv(FourierMapType& m) const = 0;
    virtual void mul_inv_sqrt(FourierMapType& m) const = 0;
  };

  template<typename T>
  class FourierMap
  {
  protected:
    FourierMap() {}
    
  public:
    typedef Eigen::Array<T, 1, Eigen::Dynamic> VecType;
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

    FourierMap<T>& operator=(const FourierMap<T>& m)
    {
      eigen() = m.eigen();
      return *this;
    }

    void sqrt()
    {
      MapType m = eigen();
      m = m.sqrt();
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
      ConstMapType m2(map2->data(), map2->size());
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

    virtual T dot_product(const FourierMap<T>& second) const
       throw(std::bad_cast) = 0;

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

  template<typename T>
  class MapUtilityFunction
  {
  public:
    typedef SpectrumFunction<T> Spectrum;
    typedef boost::shared_ptr<Spectrum> Spectrum_ptr;
    typedef FourierMap<std::complex<T> > FMap;

    virtual Spectrum_ptr estimateSpectrumFromMap(const FMap& m) const = 0;
    virtual Spectrum_ptr newSpectrumFromRaw(T *data, long size,
                                            const Spectrum& like_spec) const = 0;
  };

};

#endif
