#ifndef __COSMOTOOL_ALGO_HPP
#define __COSMOTOOL_ALGO_HPP

#include "config.hpp"

namespace  CosmoTool
{

  template<typename T>
  T cube(T a)
  {
    return a*a*a;
  }

  template<typename T>
  T square(T a)
  {
    return a*a;
  }

  template<int N>
  class SPowerBase
  {
  public:
    template<typename T>
    static T spower(T a)
    {
      if (N<0)
        {
          return 1/SPowerBase<-N>::spower(a);
        }
      if ((N%2)==0)
	{
	  T b = SPowerBase<N/2>::spower(a);
	  return b*b;
	}
      T b = SPowerBase<(N-1)/2>::spower(a);
      return a*b*b;
    }
  };

  template<>
  class SPowerBase<0>
  {
  public:
    template<typename T>
    static T spower(T a)
    {
      return T(1);
    }
  };

  template<int N, typename T>
  T spower(T a)
  {
    return SPowerBase<N>::spower(a);
  }



};

#endif
