/*+
This is CosmoTool (./src/fourier/base_types.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
    typedef T type;
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
    typedef T type;
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
    FourierTransform(const FourierTransform<T>& a) { abort(); }
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
    typedef T type;
    typedef SpectrumFunction<T> Spectrum;
    typedef boost::shared_ptr<Spectrum> Spectrum_ptr;
    typedef FourierMap<std::complex<T> > FMap;

    virtual Spectrum_ptr estimateSpectrumFromMap(const FMap& m) const = 0;
    virtual Spectrum_ptr newSpectrumFromRaw(T *data, long size,
                                            const Spectrum& like_spec) const = 0;
  };

};

#endif
