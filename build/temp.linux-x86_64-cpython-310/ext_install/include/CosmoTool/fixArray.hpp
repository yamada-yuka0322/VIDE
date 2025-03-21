/*+
This is CosmoTool (./src/fixArray.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __FIX_ARRAY_HPP
#define __FIX_ARRAY_HPP

namespace CosmoTool
{

  template <typename T, unsigned int sz> class fixArray
  {
  private:
    T d[sz];
    
  public:
    /*! Returns the size of the array. */
    long size() const { return sz; }

    /*! Returns a reference to element \a #n */
    template<typename T2> T &operator[] (T2 n) {
      return d[n];
    }

    /*! Returns a constant reference to element \a #n */
    template<typename T2> const T &operator[] (T2 n) const {
      return d[n];
    }

    template<typename T2> void importArray(T2 *indata) {
      for (int i = 0; i < sz; i++)
	d[i] = indata[i];
    }

    template<typename T2> void exportArray(T2 *outdata) {
      for (int i = 0; i < sz; i++)
	outdata[i] = d[i];
    }

  };

};

#endif
