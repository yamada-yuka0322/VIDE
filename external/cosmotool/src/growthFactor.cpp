/*+
This is CosmoTool (./src/growthFactor.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#include <cmath>
#include <gsl/gsl_integration.h>
#include "interpolate.hpp"
#include "growthFactor.hpp"

using namespace CosmoTool;

#define AMIN 1e-5
#define AMAX 1.0
#define NUM_WORK 5000
#define TOLERANCE 1e-6

typedef struct {
  double OmegaLambda;
  double OmegaMatter;
  double Hubble;
} Cosmology;

static double computeOmegaMatter(Cosmology *cosmo, double a)
{
  return cosmo->OmegaMatter / (cosmo->OmegaMatter + a*a*a * cosmo->OmegaLambda);
}

static double computeHdotH(Cosmology *cosmo, double a)
{
  return -1.5 * cosmo->OmegaMatter / (a * (cosmo->OmegaMatter + a*a*a*cosmo->OmegaLambda));
}

static double computeE(double OmegaMatter, double OmegaLambda, double a)
{
  double H2;
  double OmegaK = (1 - OmegaMatter - OmegaLambda);
  
  H2 = OmegaMatter/(a*a*a) + OmegaLambda + OmegaK/(a*a);

  return sqrt(H2);
}

static double computeEprime(Cosmology *cosmo, double a)
{
  double H2;
  double OmegaK = (1 - cosmo->OmegaMatter - cosmo->OmegaLambda);
  
  H2 = -3*cosmo->OmegaMatter/(a*a*a*a) - 2*OmegaK/(a*a*a);

  return 0.5*H2/computeE(cosmo->OmegaMatter, cosmo->OmegaLambda, a);  
}

static inline double cube(double x)
{
  return x*x*x;
}

static double integrandGrowthFactor(double a, void *params)
{
  Cosmology *cosmo = (Cosmology *)params;

  return 1/cube(computeE(cosmo->OmegaMatter, cosmo->OmegaLambda, a)*a);
}

Interpolate CosmoTool::buildLinearGrowth(double OmegaLambda, double OmegaMatter, double Hubble, int numPoints)
{
  Cosmology cosmology;
  gsl_integration_workspace *work = gsl_integration_workspace_alloc(NUM_WORK);
  gsl_function f;
  double *a_input, *D_output;

  cosmology.OmegaLambda = OmegaLambda;
  cosmology.OmegaMatter = OmegaMatter;
  cosmology.Hubble = Hubble;

  a_input = new double[numPoints];
  D_output = new double[numPoints];

  f.params = &cosmology;
  f.function = integrandGrowthFactor;
  
  a_input[0] = 0;
  D_output[0] = 0;

  for (int i = 1; i < numPoints; i++)
    {
      double a_dest = 0 + 1.0*i/(numPoints-1);
      double result, abserr;
      double E = computeE(cosmology.OmegaMatter, cosmology.OmegaLambda, a_dest);
      double Eprime = computeEprime(&cosmology, a_dest);

      gsl_integration_qag(&f, 0, a_dest, 0, TOLERANCE, NUM_WORK, 
			  GSL_INTEG_GAUSS61, work, &result, &abserr);


      result *= 2.5 * computeE(cosmology.OmegaMatter, cosmology.OmegaLambda, a_dest) * OmegaMatter;

      D_output[i] = result;
      a_input[i] = a_dest;

    }
  gsl_integration_workspace_free(work);

  for (int i = 0; i < numPoints; i++)
    {
      D_output[i] /= D_output[numPoints-1];
    }

  Interpolate p(a_input, D_output, numPoints, true, false, true);

  delete[] a_input;
  delete[] D_output;

  return p;
}
