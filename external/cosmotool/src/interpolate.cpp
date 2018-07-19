/*+
This is CosmoTool (./src/interpolate.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#include <cassert>
#include <cmath>
#include <algorithm>
#include "interpolate.hpp"
#include <vector>
#include <fstream>
#include <inttypes.h>
#include <sstream>

using namespace std;
using namespace CosmoTool;
  
Interpolate::Interpolate(double *xs, double *values, uint32_t N, bool autofree,
			 bool logx, bool logy)
{
  spline = gsl_spline_alloc (gsl_interp_linear, N);
  gsl_spline_init (spline, xs, values, N);
  accel_interp = gsl_interp_accel_alloc();

  this->logx = logx;
  this->logy = logy;
  this->autoFree = autofree;
}

Interpolate::~Interpolate()
{
  if (spline != 0)
    gsl_spline_free (spline);
  if (accel_interp != 0)
    gsl_interp_accel_free (accel_interp);
}

double Interpolate::compute(double x)
{
  double y;
  
  if (logx)
    x = log(x);

  int err = gsl_spline_eval_e(spline, x, accel_interp, &y);

  if (err)
    throw InvalidRangeException("Interpolate argument outside range");

  if (logy)
    return exp(y);
  else
    return y;
}

double Interpolate::compute(double x) const
{
  double y;

  if (logx)
    x = log(x);

  int err = gsl_spline_eval_e(spline, x, 0, &y);

  if (err)
    throw InvalidRangeException("Interpolate argument outside range");

  if (logy)
    return exp(y);
  else
    return y;
}




double Interpolate::derivative(double x)
{
  double y, dy, x0 = x;
  
  if (logx)
    x0 = log(x0);

  int err = gsl_spline_eval_deriv_e(spline, x0, accel_interp, &dy);

  if (err)
    throw InvalidRangeException("Interpolate argument outside range");

  if (logy)
    {
      int err = gsl_spline_eval_e(spline, x0, accel_interp, &y);
      
      assert(err == 0);
      
      return dy*exp(y)/x0;
    }
  else
    return dy;
}

uint32_t Interpolate::getNumPoints() const
{
  return spline->size;
}

void Interpolate::fillWithXY(double *x, double *y) const
{
  if (x != 0)
    memcpy(x, spline->x, sizeof(double)*spline->size);
  if (y != 0)
    memcpy(y, spline->y, sizeof(double)*spline->size);
}

const Interpolate& Interpolate::operator=(const Interpolate& a)
{
  double *x, *y;

  if (spline != NULL)
    {
      gsl_spline_free (spline);
      gsl_interp_accel_free (accel_interp);
    }

  autoFree = true;
  spline = gsl_spline_alloc (gsl_interp_linear, a.spline->size);
  accel_interp = gsl_interp_accel_alloc ();
  gsl_spline_init(spline, a.spline->x, a.spline->y, a.spline->size);
  logx = a.logx;
  logy = a.logy;
}

double Interpolate::getMaxX() const
{
  if (logx)
    return exp(spline->x[spline->size-1]);
  else
    return spline->x[spline->size-1];
}

typedef struct {
  double a, b;
} MyPair;

bool operator<(const MyPair& a, const MyPair& b)
{
  return a.a < b.a;
}

Interpolate CosmoTool::buildFromVector(const InterpolatePairs& v)
{
  double *x = new double[v.size()];
  double *y = new double[v.size()];
  
  for (uint32_t i = 0; i < v.size(); i++)
    {
      x[i] = v[i].first;
      y[i] = v[i].second;
    }

  Interpolate inter = Interpolate(x, y, v.size(), true);

  delete[] x;
  delete[] y;

  return inter;  
}

Interpolate CosmoTool::buildInterpolateFromFile(const char *fname)
{
  vector<MyPair> allData;
  ifstream f(fname);

  if (!f)
    throw NoSuchFileException(fname);

  do 
    {
      MyPair m;

      if (!(f >> m.a >> m.b))
	break;

      allData.push_back(m);
    }
  while (1);

  sort(allData.begin(), allData.end());

  double *x = new double[allData.size()];
  double *y = new double[allData.size()];
  
  for (uint32_t i = 0; i < allData.size(); i++)
    {
      x[i] = allData[i].a;
      y[i] = allData[i].b;
    }

  Interpolate inter = Interpolate(x, y, allData.size(), true);

  delete[] x;
  delete[] y;

  return inter;
}

Interpolate CosmoTool::buildInterpolateFromColumns(const char *fname, uint32_t col1, uint32_t col2, bool logx,
						   bool logy)
{
  vector<MyPair> allData;
  ifstream f(fname);

  if (!f)
    throw NoSuchFileException(fname);
  
  bool swapped = (col1 > col2);
  uint32_t colMin = min(col1, col2);
  uint32_t colMax = max(col1, col2);

  do 
    {
      MyPair m;
      string line;

      if (getline(f, line).eof())
	break;

      istringstream iss(line);
      double dummy;
      double val1, val2;

      for (uint32_t i = 0; i < colMin; i++)
	iss >> dummy;
      if (!(iss >> val1))
	throw InvalidRangeException("Invalid first column");

      if (col2 != col1)
	{
	  for (uint32_t i = 0; i < (colMax-colMin-1); i++)
	    iss >> dummy;
	  if (!(iss >> val2))
	    throw InvalidRangeException("Invalid second column");
	}
      else
	val2 = val2;

      if (!swapped)
	{
	  m.a = val1;
	  m.b = val2;
	}
      else
	{
	  m.a = val2;
	  m.b = val1;
	}

      allData.push_back(m);
    }
  while (1);

  sort(allData.begin(), allData.end());

  double *x = new double[allData.size()];
  double *y = new double[allData.size()];
  
  for (uint32_t i = 0; i < allData.size(); i++)
    {
      x[i] = logx ? log(allData[i].a) : allData[i].a;
      y[i] = logy ? log(allData[i].b) : allData[i].b;
    }

  Interpolate inter = Interpolate(x, y, allData.size(), true, logx, logy);

  delete[] x;
  delete[] y;
  
  return inter;
}
