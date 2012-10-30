#ifndef __FFTW_UNIFIED_CALLS_HPP
#define __FFTW_UNIFIED_CALLS_HPP

#include <fftw3.h>

namespace CosmoTool
{

static inline void init_fftw_wisdom()
{
  fftw_import_system_wisdom();
  fftw_import_wisdom_from_filename("fft_wisdom");
}

static inline void save_fftw_wisdom()
{
  fftw_export_wisdom_to_filename("fft_wisdom");
}

template<typename T> class FFTW_Calls {};


#define FFTW_CALLS_BASE(rtype, prefix) \
  template<>				\
class FFTW_Calls<rtype> {		\
public: \
  typedef rtype real_type; \
  typedef prefix ## _complex complex_type; \
  typedef prefix ## _plan plan_type; \
  \
  static complex_type *alloc_complex(int N) { return prefix ## _alloc_complex(N); } \
  static real_type *alloc_real(int N) { return prefix ## _alloc_real(N); } \
  static void free(void *p) { fftw_free(p); } \
\
  static void execute(plan_type p) { prefix ## _execute(p); } \
  static plan_type plan_dft_r2c_2d(int Nx, int Ny,  \
		       real_type *in, complex_type *out, \
		       unsigned flags) \
  { \
    return prefix ## _plan_dft_r2c_2d(Nx, Ny, in, out, \
				flags); \
  } \
  static plan_type plan_dft_c2r_2d(int Nx, int Ny,  \
                       complex_type *in, real_type *out, \
                       unsigned flags) \
  { \
    return prefix ## _plan_dft_c2r_2d(Nx, Ny, in, out, \
                                flags); \
  } \
  static plan_type plan_dft_r2c_3d(int Nx, int Ny, int Nz, \
                                   real_type *in, complex_type *out, \
                                   unsigned flags) \
  { \
    return prefix ## _plan_dft_r2c_3d(Nx, Ny, Nz, in, out, flags); \
  } \
  static plan_type plan_dft_r2c(int rank, const int *n, real_type *in, \
                                complex_type *out, unsigned flags) \
  { \
    return prefix ## _plan_dft_r2c(rank, n, in, out, flags); \
  } \
  static plan_type plan_dft_c2r(int rank, const int *n, complex_type *in, \
                                real_type *out, unsigned flags) \
  { \
    return prefix ## _plan_dft_c2r(rank, n, in, out, flags); \
  } \
  static void destroy_plan(plan_type plan) { prefix ## _destroy_plan(plan); } \
}


FFTW_CALLS_BASE(double, fftw);
FFTW_CALLS_BASE(float, fftwf);

};

#endif
