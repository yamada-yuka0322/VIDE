#ifndef __DETAILS_EUCLIDIAN_MAPS
#define __DETAILS_EUCLIDIAN_MAPS


namespace CosmoTool
{
  
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

    virtual ~EuclidianFourierMapBase()
    {
    }

    const DimArray& getDims() const { return m_dims; }

    virtual const T *data() const { return m_data.get(); }
    virtual T *data() { return m_data.get(); }
    virtual long size() const  { return m_size; } 
    
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
      : EuclidianFourierMapBase<std::complex<T> >(indata, indims), delta_k(dk), m_dim0(dim0), even0((dim0 % 2)==0)
    {
      assert(dk.size() == indims.size()); 
      plane_size = 1;
      alleven = true;
      for (int q = 1; q < indims.size(); q++)
	{
	  plane_size *= indims[q];
	  alleven = alleven && ((indims[q]%2)==0);
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
      for (int q = 1; q < ik.size(); q++)
	{
	  int dk = ik[q];
	  if (dk > dims[q]/2)
	    dk = dk - dims[q];

          kvec[q] = dk*delta_k[q];
        }
    }

    template<typename Array2>
    void get_Kvec(long p, Array2& kvec) const
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

      for (int q = 1; q < ik.size(); q++)
	{
	  int dk = ik[q];

	  if (dk > dims[q]/2)
	    dk = dk - dims[q];
	  
	  k2 += CosmoTool::square(delta_k[q]*dk);
	}
      return std::sqrt(k2);
    }

    double get_K(long p) const
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
      const EuclidianFourierMapComplex<T>& m2 = dynamic_cast<const EuclidianFourierMapComplex<T>&>(other);
      if (this->size() != m2.size())
        throw std::bad_cast();
      
      const std::complex<T> *d1 = this->data();
      const std::complex<T> *d2 = m2.data();
      const DimArray& dims = this->getDims();
      int N0 = dims[0] + (even0 ? 0 : 1);
      std::complex<T> result = 0;

      for (long q0 = 1; q0 < N0-1; q0++)
	{
	  for (long p = 0; p < plane_size; p++)
	    {
	      long idx = q0+dims[0]*p;
	      assert(idx < this->size());
	      result += 2*(conj(d1[idx]) * d2[idx]).real();
	    }
	}
      if (even0)
	{
	  for (long p = 0; p < plane_size; p++)
	    {
	      long q0 = N0*p, q1 = (p+1)*N0-1;
	      result += conj(d1[q0]) * d2[q0];
	      result += conj(d1[q1]) * d2[q1];
	    }
	}
      return result;
    }

  };
};

#endif
