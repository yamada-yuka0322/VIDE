/*+
This is CosmoTool (./src/newton.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
