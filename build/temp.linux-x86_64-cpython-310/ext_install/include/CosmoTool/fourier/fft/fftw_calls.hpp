/*+
This is CosmoTool (./src/fourier/fft/fftw_calls.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __FFTW_UNIFIED_CALLS_HPP
#define __FFTW_UNIFIED_CALLS_HPP

#include <fftw3.h>
#include <complex>

namespace CosmoTool {

  static inline void init_fftw_wisdom() {
    fftw_import_system_wisdom();
    fftw_import_wisdom_from_filename("fft_wisdom");
  }

  static inline void save_fftw_wisdom() {
    fftw_export_wisdom_to_filename("fft_wisdom");
  }

  template <typename T>
  class FFTW_Calls {};

#define FFTW_CALLS_BASE(rtype, prefix)                                         \
  template <>                                                                  \
  class FFTW_Calls<rtype> {                                                    \
  public:                                                                      \
    typedef rtype real_type;                                                   \
    typedef prefix##_complex complex_type;                                     \
    typedef prefix##_plan plan_type;                                           \
                                                                               \
    static complex_type *alloc_complex(size_t N) {                             \
      return prefix##_alloc_complex(N);                                        \
    }                                                                          \
    static real_type *alloc_real(size_t N) { return prefix##_alloc_real(N); }  \
    static void free(void *p) { fftw_free(p); }                                \
                                                                               \
    static void execute(plan_type p) { prefix##_execute(p); }                  \
    static void execute_r2c(plan_type p, real_type *in, complex_type *out) {   \
      prefix##_execute_dft_r2c(p, in, out);                                    \
    }                                                                          \
    static void execute_c2r(plan_type p, complex_type *in, real_type *out) {   \
      prefix##_execute_dft_c2r(p, in, out);                                    \
    }                                                                          \
    static void                                                                \
    execute_r2c(plan_type p, real_type *in, std::complex<real_type> *out) {    \
      prefix##_execute_dft_r2c(p, in, (complex_type *)out);                    \
    }                                                                          \
    static void                                                                \
    execute_c2r(plan_type p, std::complex<real_type> *in, real_type *out) {    \
      prefix##_execute_dft_c2r(p, (complex_type *)in, out);                    \
    }                                                                          \
    static void execute_c2c(                                                   \
        plan_type p, std::complex<real_type> *in,                              \
        std::complex<real_type> *out) {                                        \
      prefix##_execute_dft(p, (complex_type *)in, (complex_type *)out);        \
    }                                                                          \
    static plan_type plan_dft_r2c_1d(                                          \
        int Nx, real_type *in, complex_type *out, unsigned flags) {            \
      return prefix##_plan_dft_r2c_1d(Nx, in, out, flags);                     \
    }                                                                          \
    static plan_type plan_dft_c2r_1d(                                          \
        int Nx, complex_type *in, real_type *out, unsigned flags) {            \
      return prefix##_plan_dft_c2r_1d(Nx, in, out, flags);                     \
    }                                                                          \
    static plan_type plan_dft_r2c_2d(                                          \
        int Nx, int Ny, real_type *in, complex_type *out, unsigned flags) {    \
      return prefix##_plan_dft_r2c_2d(Nx, Ny, in, out, flags);                 \
    }                                                                          \
    static plan_type plan_dft_c2r_2d(                                          \
        int Nx, int Ny, complex_type *in, real_type *out, unsigned flags) {    \
      return prefix##_plan_dft_c2r_2d(Nx, Ny, in, out, flags);                 \
    }                                                                          \
    static plan_type plan_dft_r2c_3d(                                          \
        int Nx, int Ny, int Nz, real_type *in, complex_type *out,              \
        unsigned flags) {                                                      \
      return prefix##_plan_dft_r2c_3d(Nx, Ny, Nz, in, out, flags);             \
    }                                                                          \
    static plan_type plan_dft_c2r_3d(                                          \
        int Nx, int Ny, int Nz, complex_type *in, real_type *out,              \
        unsigned flags) {                                                      \
      return prefix##_plan_dft_c2r_3d(Nx, Ny, Nz, in, out, flags);             \
    }                                                                          \
                                                                               \
    static plan_type plan_dft_r2c(                                             \
        int rank, const int *n, real_type *in, complex_type *out,              \
        unsigned flags) {                                                      \
      return prefix##_plan_dft_r2c(rank, n, in, out, flags);                   \
    }                                                                          \
    static plan_type plan_dft_c2r(                                             \
        int rank, const int *n, complex_type *in, real_type *out,              \
        unsigned flags) {                                                      \
      return prefix##_plan_dft_c2r(rank, n, in, out, flags);                   \
    }                                                                          \
    static plan_type plan_dft_3d(                                              \
        int Nx, int Ny, int Nz, complex_type *in, complex_type *out, int sign, \
        unsigned flags) {                                                      \
      return prefix##_plan_dft_3d(Nx, Ny, Nz, in, out, sign, flags);           \
    }                                                                          \
    static plan_type plan_dft_2d(                                              \
        int Nx, int Ny, complex_type *in, complex_type *out, int sign,         \
        unsigned flags) {                                                      \
      return prefix##_plan_dft_2d(Nx, Ny, in, out, sign, flags);               \
    }                                                                          \
    static plan_type plan_dft_1d(                                              \
        int Nx, complex_type *in, complex_type *out, int sign,                 \
        unsigned flags) {                                                      \
      return prefix##_plan_dft_1d(Nx, in, out, sign, flags);                   \
    }                                                                          \
    static void destroy_plan(plan_type plan) { prefix##_destroy_plan(plan); }  \
  }

  FFTW_CALLS_BASE(double, fftw);
  FFTW_CALLS_BASE(float, fftwf);
#undef FFTW_CALLS_BASE
}; // namespace CosmoTool

#endif
