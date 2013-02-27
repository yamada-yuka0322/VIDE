#ifndef __COSMOTOOL_FOURIER_HEALPIX_DETAILS_SPECTRUM_HP
#define __COSMOTOOL_FOURIER_HEALPIX_DETAILS_SPECTRUM_HP

namespace CosmoTool
{

  template<typename T>
  class HealpixUtilityFunction: public MapUtilityFunction<T>
  {
  public:
    typedef typename MapUtilityFunction<T>::Spectrum Spectrum;
    typedef typename MapUtilityFunction<T>::Spectrum_ptr Spectrum_ptr;
    typedef typename MapUtilityFunction<T>::FMap FMap;
    typedef typename HealpixSpectrum<T>::LType LType; 
    
    Spectrum_ptr estimateSpectrumFromMap(const FMap& m) const
    {
      const HealpixFourierALM<T>& alms = dynamic_cast<const HealpixFourierALM<T>&>(m);
      LType Lmax = alms.Lmax(), Mmax = alms.Mmax();
      HealpixSpectrum<T> *spectrum = new HealpixSpectrum<T>(Lmax);
      T *data_cls = spectrum->data();
      const std::complex<T> *data_alms = alms.data();

      // make an estimate of the cls
      std::fill(data_cls, data_cls+Lmax+1, 0);
      LType q = 0;
      for (LType m = 0; m <= Mmax; ++m )
	{
	  for (LType l = m; l <= Lmax; ++l)
	    {
	      // Triangular storage
	      data_cls[l] += std::norm(data_alms[l+q]);
	    }
	  q += Lmax+1-m;
	}
      assert(q==alms.size());
      
      for (LType l = 0; l <= Lmax; ++l)
	{
	  int dof = 1 + std::min(l, Mmax);
	  spectrum->set_dof(l, dof);
	  data_cls[l] /= dof;
	}

      return Spectrum_ptr(spectrum);
    }

    Spectrum_ptr newSpectrumFromRaw(T *data, long size,
				    const Spectrum& like_spec) const
    {
      const HealpixSpectrum<T>& in_spec = dynamic_cast<const HealpixSpectrum<T>&>(like_spec);
      HealpixSpectrum<T> *new_spectrum = new HealpixSpectrum<T>(in_spec.Lmax());
      T *out_d = new_spectrum->data();

      std::copy(data, data + std::min(size,new_spectrum->size()), out_d);

      return Spectrum_ptr(new_spectrum);
    }
  };

};

#endif
