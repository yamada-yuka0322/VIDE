/*+
This is CosmoTool (./src/interpolate.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __CTOOL_INTERPOLATE_HPP
#define __CTOOL_INTERPOLATE_HPP

#include "config.hpp"
#include <inttypes.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <vector>
#include <utility>

namespace CosmoTool
{
  
  class Interpolate
  {
  public:
    Interpolate() : spline(0), accel_interp(0) { }
    Interpolate(double *xs, double *values, uint32_t N, bool autofree = false,
		bool logx = false, bool logy = false);
    ~Interpolate();

    double compute(double x);
    double compute(double x) const;
    double derivative(double x);

    const Interpolate& operator=(const Interpolate& a);

    uint32_t getNumPoints() const;
    void fillWithXY(double *x, double *y) const;
    double getMaxX() const;
    double getXi(int i) const { return spline->x[i]; }
    double getYi(int i) const { return spline->y[i]; }
  protected:
    gsl_interp_accel *accel_interp;
    gsl_spline *spline;
    bool autoFree;
    bool logx, logy;
  };

  typedef std::vector< std::pair<double,double> > InterpolatePairs;

  Interpolate buildInterpolateFromFile(const char *fname);
  Interpolate buildInterpolateFromColumns(const char *fname, uint32_t col1, uint32_t col2, bool logx = false, bool logy = false);
  Interpolate buildFromVector(const InterpolatePairs& v);


  class FunctorInterpolate
  {
  public:
    FunctorInterpolate(Interpolate& i) : m_i(i) {}    

    double eval(double x) { return m_i.compute(x); }
    double derivative(double x) { return m_i.derivative(x); }
  private:
    Interpolate& m_i;
  };
  
};

#endif
