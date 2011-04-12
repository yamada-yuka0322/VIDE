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
