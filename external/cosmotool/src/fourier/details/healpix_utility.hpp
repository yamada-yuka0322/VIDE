/*+
This is CosmoTool (./src/fourier/details/healpix_utility.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
