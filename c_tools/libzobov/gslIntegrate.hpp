/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/libzobov/gslIntegrate.hpp
    Copyright (C) 2010-2014 Guilhem Lavaux
    Copyright (C) 2011-2014 P. M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
+*/



#ifndef __MYGSL_INTEGRATE_HPP
#define __MYGSL_INTEGRATE_HPP

#include <gsl/gsl_integration.h>

template<typename FunT>
double gslSpecialFunction(double x, void *param)
{
  FunT *f = (FunT *)param;

  return (*f)(x);
}

template<typename FunT>
double gslIntegrate(FunT& v, double a, double b, double prec, int NPTS = 1024)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(NPTS);
  gsl_function f;
  double result;
  double abserr;

  f.function = &gslSpecialFunction<FunT>;
  f.params = &v;

  gsl_integration_qag(&f, a, b, prec, 0, NPTS, GSL_INTEG_GAUSS61,
		      w, &result, &abserr);
  
  gsl_integration_workspace_free(w);

  return result;
}

#endif
