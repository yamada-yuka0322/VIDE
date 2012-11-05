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

    double compute(double x)
      throw (InvalidRangeException);
    double compute(double x) const
      throw (InvalidRangeException);
    double derivative(double x)
      throw (InvalidRangeException);

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

  Interpolate buildInterpolateFromFile(const char *fname)
    throw (NoSuchFileException);
  Interpolate buildInterpolateFromColumns(const char *fname, uint32_t col1, uint32_t col2, bool logx = false, bool logy = false)
    throw (NoSuchFileException,InvalidRangeException);
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
