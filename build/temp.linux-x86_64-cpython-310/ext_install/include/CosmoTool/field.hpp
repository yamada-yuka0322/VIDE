/*+
This is CosmoTool (./src/field.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __COSMOTOOL_FIELD
#define __COSMOTOOL_FIELD

#include "config.hpp"
#include <iostream>
#include <cassert>

namespace CosmoTool {

  template<typename BaseType>
  struct ScalarField
  {
    BaseType value;    
  };

  template<typename BaseType, int N>
  struct VectorField
  {
    BaseType vec[N];
    
    VectorField& operator=(const VectorField& a)
    {
      for (int i = 0; i < N; i++)
	vec[i] = a.vec[i];
      return *this;
    }

    VectorField()
    {
      for (int i = 0; i < N; i++)
	vec[i] = 0;      
    }

    VectorField(double a)
    {
      assert(a == 0);
      for (int i = 0; i < N; i++)
	vec[i] = 0;      
    }
  };

  template<typename BaseType, int N>
  VectorField<BaseType,N> operator*(BaseType s, const VectorField<BaseType,N>& a)
  {
    VectorField<BaseType,N> v;

    for (int i = 0; i < N; i++)
      v.vec[i] = a.vec[i]*s;
    
    return v;
  }  
  
  template<typename BaseType, int N>
  VectorField<BaseType,N> operator+(const VectorField<BaseType,N>& a, const VectorField<BaseType,N>& b)
  {
    VectorField<BaseType,N> v;

    for (int i = 0; i < N; i++)
      v.vec[i] = a.vec[i]+b.vec[i];
    
    return v;
  }  

  template<typename BaseType, int N>
  VectorField<BaseType,N>& operator+=(VectorField<BaseType,N>& a, const VectorField<BaseType,N>& b)
  {
    for (int i = 0; i < N; i++)
      a.vec[i]+=b.vec[i];
    
    return a;
  }  

};

template<typename BaseType, int N>
std::ostream& operator<<(std::ostream& s, const CosmoTool::VectorField<BaseType,N>& a)
{
  for (int i = 0; i < N; i++)
    s << a.vec[i] << " " ;
  return s;
}

#endif
