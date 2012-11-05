#ifndef __COSMO_MACHINE_TEST_HPP
#define __COSMO_MACHINE_TEST_HPP

#include <iostream>

template<typename T>
T mach_epsilon()
{
  T eps = (T)1;

  do
    {
      eps /= 2;
    }
  while ((T)(1 + (eps/2)) != (T)1);

  return eps;
}

#endif
