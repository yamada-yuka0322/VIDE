/*+
This is CosmoTool (./src/fourier/details/euclidian_maps.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __DETAILS_EUCLIDIAN_MAPS
#define __DETAILS_EUCLIDIAN_MAPS

#include <cmath>
#include <boost/multi_array.hpp>

namespace CosmoTool
{
  
  namespace details {
    static void no_free_euclidian_map(void *) {}
  }
  
  template<typename T>
  class EuclidianFourierMapBase: public FourierMap<T>
  {
  public:
    typedef std::vector<int> DimArray;
  private:
    boost::shared_ptr<T> m_data;
    DimArray m_dims;
    long m_size;
  public:
    
    EuclidianFourierMapBase(boost::shared_ptr<T> indata, const DimArray& indims)
    {
      m_data = indata;
      m_dims = indims;
      m_size = 1;
      for (int i = 0; i < m_dims.size(); i++)
        m_size *= m_dims[i];
    }

    template<typename ArrayType>
    EuclidianFourierMapBase(ArrayType& indata)
    {
      m_data = boost::shared_ptr<T>(
            indata.origin(), 
            std::ptr_fun(details::no_free_euclidian_map));
      m_dims = DimArray(indata.num_dimensions());
      m_size = indata.num_elements();
      for (int i = 0; i < m_dims.size(); i++)
        m_dims[i] = indata.shape()[i];
    }


    virtual ~EuclidianFourierMapBase()
    {
    }

    const DimArray& getDims() const { return m_dims; }

    virtual const T *data() const { return m_data.get(); }
    virtual T *data() { return m_data.get(); }
    virtual long size() const  { return m_size; } 
    
    boost::multi_array_ref<T, 1>& array() { 
      return boost::multi_array_ref<T, 1>(m_data.get(), boost::extents[size]);
    }

    boost::const_multi_array_ref<T, 1>& array() const { 
      return boost::const_multi_array_ref<T, 1>(m_data.get(), boost::extents[size]);
    }
    
    virtual FourierMap<T> *copy() const
    {
      FourierMap<T> *m = this->mimick();
      m->eigen() = this->eigen();
      return m;
    }

  };

  template<typename T>
  class EuclidianFourierMapReal: public EuclidianFourierMapBase<T>
  {
  public:
    typedef typename EuclidianFourierMapBase<T>::DimArray DimArray;

    EuclidianFourierMapReal(boost::shared_ptr<T> indata, const DimArray& indims)
      : EuclidianFourierMapBase<T>(indata, indims) 
    {}

    template<typename ArrayType>
    EuclidianFourierMapReal(ArrayType& indata) 
      : EuclidianFourierMapBase<T>(indata) 
    {}


    virtual FourierMap<T> *mimick() const
    {
      return new EuclidianFourierMapReal<T>(
		      boost::shared_ptr<T>((T *)fftw_malloc(sizeof(T)*this->size()), 
					   std::ptr_fun(fftw_free)),
		      this->getDims());
    }

    virtual T dot_product(const FourierMap<T>& other) const
      throw(std::bad_cast)
    {
      const EuclidianFourierMapReal<T>& m2 = dynamic_cast<const EuclidianFourierMapReal<T>&>(other);
      if (this->size() != m2.size())
        throw std::bad_cast();

      return (this->eigen()*m2.eigen()).sum();
    }
  };

  template<typename T>
  class EuclidianFourierMapComplex: public EuclidianFourierMapBase<std::complex<T> >
  {
  protected:
    typedef boost::shared_ptr<std::complex<T> > ptr_t;
    std::vector<double> delta_k;
    int m_dim0;
    bool even0, alleven;
    long plane_size;
  public:
    typedef typename EuclidianFourierMapBase<std::complex<T> >::DimArray DimArray;

    EuclidianFourierMapComplex(ptr_t indata,
			       int dim0,
			       const DimArray& indims, 
			       const std::vector<double>& dk)
      : EuclidianFourierMapBase<std::complex<T> >(indata, indims), 
        delta_k(dk), m_dim0(dim0), even0((dim0 % 2)==0)
    {
      assert(dk.size() == indims.size()); 
      plane_size = 1;
      alleven = true;
      for (int q = 1; q < indims.size(); q++) {
        plane_size *= indims[q];
        alleven = alleven && ((indims[q]%2)==0);
      }
    }

    template<typename ArrayType>
    EuclidianFourierMapComplex(ArrayType& indata, 
                      	       int dim0,
                               const std::vector<double>& dk) 
      : EuclidianFourierMapBase<std::complex<T> >(indata),
        delta_k(dk), m_dim0(dim0), even0((dim0 % 2)==0)
    {
      assert(dk.size() == indata.num_dimensions()); 
      plane_size = 1;
      alleven = true;
      for (int q = 1; q < this->m_dims.size(); q++) {
        plane_size *= this->m_dims[q];
        alleven = alleven && ((this->m_dims[q]%2)==0);
      }
    }

    virtual FourierMap<std::complex<T> > *mimick() const
    {
      return 
	new EuclidianFourierMapComplex<T>(
		   ptr_t((std::complex<T> *)
			 fftw_malloc(sizeof(std::complex<T>)*this->size()), 
			 std::ptr_fun(fftw_free)),
		   m_dim0,
		   this->getDims(),
		   this->delta_k);
    }

    const std::vector<double>& get_delta_k() const
    {
      return this->delta_k;
    }

    template<typename Array, typename Array2>
    void get_Kvec(const Array& ik, Array2& kvec) const
    {
      const DimArray& dims = this->getDims();
      assert(ik.size() == dims.size());
      assert(kvec.size() == dims.size());

      kvec[0] = ik[0] * delta_k[0];
      for (int q = 1; q < ik.size(); q++) {
        int dk = ik[q];
        if (dk > dims[q]/2)
          dk = dk - dims[q];

        kvec[q] = dk*delta_k[q];
      }
    }

    template<typename Array2>
    void get_Kvec_p(long p, Array2& kvec) const
    {
      const DimArray& dims = this->getDims();
      DimArray d(delta_k.size());
      get_IKvec(p, d);
      get_Kvec(d, kvec);
    }

    void get_IKvec(long p, DimArray& ikvec) const
    {
      const DimArray& dims = this->getDims();
      assert(dims.size()==ikvec.size());
      for (int q = 0; q < ikvec.size(); q++)
       {
         ikvec[q] = p%dims[q];
         p = (p-ikvec[q])/dims[q];
       }
    }


    template<typename Array>
    double get_K(const Array& ik) const
    {
      const DimArray& dims = this->getDims();
      assert(ik.size() == dims.size());
      double k2 = 0;
      k2 += CosmoTool::square(ik[0]*delta_k[0]);

      for (int q = 1; q < ik.size(); q++) {
        int dk = ik[q];

        if (dk > dims[q]/2)
          dk = dk - dims[q];
        
        k2 += CosmoTool::square(delta_k[q]*dk);
      }
      return std::sqrt(k2);
    }

    double get_K_p(long p) const
    {
      const DimArray& dims = this->getDims();
      DimArray d(delta_k.size());
      for (int q = 0; q < d.size(); q++)
       {
         d[q] = p%dims[q];
         p = (p-d[q])/dims[q];
       }
      return get_K(d);
    }

    bool allDimensionsEven() const { return alleven; }
    bool firstDimensionEven() const { return even0; }

    virtual std::complex<T> dot_product(const FourierMap<std::complex<T> >& other) const
      throw(std::bad_cast)
    {
      const EuclidianFourierMapComplex<T>& m2 = 
        dynamic_cast<const EuclidianFourierMapComplex<T>&>(other);
      if (this->size() != m2.size())
        throw std::bad_cast();
      
      const std::complex<T> *d1 = this->data();
      const std::complex<T> *d2 = m2.data();
      const DimArray& dims = this->getDims();
      int N0 = dims[0] + (even0 ? 0 : 1);
      std::complex<T> result = 0;

      for (long q0 = 1; q0 < N0-1; q0++) {
        for (long p = 0; p < plane_size; p++) {
          long idx = q0+dims[0]*p;
          assert(idx < this->size());
          result += T(2)*(std::conj(d1[idx]) * d2[idx]).real();
        }
      }
      if (even0) {
        for (long p = 0; p < plane_size; p++)
          {
            long q0 = N0*p, q1 = (p+1)*N0-1;
            result += T(2)*std::conj(d1[q0]) * d2[q0];
            result += T(2)*std::conj(d1[q1]) * d2[q1];
          }
        }
      return result;
    }

  };
};

#endif
