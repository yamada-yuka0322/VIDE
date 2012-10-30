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
