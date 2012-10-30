#ifndef _COSMOTOOL_NEWTON_HPP
#define _COSMOTOOL_NEWTON_HPP

#include <cmath>

namespace CosmoTool
{
  template<typename T, typename FunT>
  T newtonSolver(T x0, FunT function, double residual = 1e-3)
  {
    T x, xold = x0;
    T f_x = function.eval(x0);
    T df_x = function.derivative(x0);
    
    x = xold - f_x/df_x;

    while (std::abs(xold-x) > residual)
      {	
	xold = x;
	f_x = function.eval(x);
	df_x = function.derivative(x);
	x = xold - f_x/df_x;
      }
	   
    return x;
  }

};

#endif
