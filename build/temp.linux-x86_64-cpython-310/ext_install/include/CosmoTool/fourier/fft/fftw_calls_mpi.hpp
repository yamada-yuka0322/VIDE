#ifndef __MPI_FFTW_UNIFIED_CALLS_HPP
#define __MPI_FFTW_UNIFIED_CALLS_HPP

#include <complex>
#include <mpi.h>
#include <fftw3-mpi.h>

namespace CosmoTool {

  static inline void init_fftw_mpi() { fftw_mpi_init(); }

  static inline void done_fftw_mpi() { fftw_mpi_cleanup(); }

  template <typename T>
  class FFTW_MPI_Calls {};

#define FFTW_MPI_CALLS_BASE(rtype, prefix)                                     \
  template <>                                                                  \
  class FFTW_MPI_Calls<rtype> {                                                \
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
    template <size_t Nd>                                                       \
    static ptrdiff_t local_size(                                               \
        std::array<ptrdiff_t, Nd> const &N, MPI_Comm comm,                     \
        ptrdiff_t *local_n0, ptrdiff_t *local_0_start) {                       \
      return prefix##_mpi_local_size(                                          \
          Nd, N.data(), comm, local_n0, local_0_start);                        \
    }                                                                          \
    static ptrdiff_t local_size_2d(                                            \
        ptrdiff_t N0, ptrdiff_t N1, MPI_Comm comm, ptrdiff_t *local_n0,        \
        ptrdiff_t *local_0_start) {                                            \
      return prefix##_mpi_local_size_2d(                                       \
          N0, N1, comm, local_n0, local_0_start);                              \
    }                                                                          \
                                                                               \
    static ptrdiff_t local_size_3d(                                            \
        ptrdiff_t N0, ptrdiff_t N1, ptrdiff_t N2, MPI_Comm comm,               \
        ptrdiff_t *local_n0, ptrdiff_t *local_0_start) {                       \
      return prefix##_mpi_local_size_3d(                                       \
          N0, N1, N2, comm, local_n0, local_0_start);                          \
    }                                                                          \
                                                                               \
    static void execute(plan_type p) { prefix##_execute(p); }                  \
    static void                                                                \
    execute_c2c(plan_type p, complex_type *in, complex_type *out) {            \
      prefix##_mpi_execute_dft(p, in, out);                                    \
    }                                                                          \
    static void execute_c2c(                                                   \
        plan_type p, std::complex<real_type> *in,                              \
        std::complex<real_type> *out) {                                        \
      prefix##_mpi_execute_dft(p, (complex_type *)in, (complex_type *)out);    \
    }                                                                          \
    static void execute_r2c(plan_type p, real_type *in, complex_type *out) {   \
      prefix##_mpi_execute_dft_r2c(p, in, out);                                \
    }                                                                          \
    static void                                                                \
    execute_c2r(plan_type p, std::complex<real_type> *in, real_type *out) {    \
      prefix##_mpi_execute_dft_c2r(p, (complex_type *)in, out);                \
    }                                                                          \
    static void execute_c2r(plan_type p, complex_type *in, real_type *out) {   \
      prefix##_mpi_execute_dft_c2r(p, in, out);                                \
    }                                                                          \
    static void                                                                \
    execute_r2c(plan_type p, real_type *in, std::complex<real_type> *out) {    \
      prefix##_mpi_execute_dft_r2c(p, in, (complex_type *)out);                \
    }                                                                          \
                                                                               \
    static plan_type plan_dft_r2c_1d(                                          \
        int n, real_type *in, complex_type *out, MPI_Comm, unsigned flags) {   \
      return prefix##_plan_dft_r2c_1d(n, in, out, flags);                      \
    }                                                                          \
                                                                               \
    static plan_type plan_dft_r2c_2d(                                          \
        int Nx, int Ny, real_type *in, complex_type *out, MPI_Comm comm,       \
        unsigned flags) {                                                      \
      return prefix##_mpi_plan_dft_r2c_2d(Nx, Ny, in, out, comm, flags);       \
    }                                                                          \
                                                                               \
    static plan_type plan_dft_c2r_1d(                                          \
        int n, complex_type *in, real_type *out, MPI_Comm, unsigned flags) {   \
      return prefix##_plan_dft_c2r_1d(n, in, out, flags);                      \
    }                                                                          \
    static plan_type plan_dft_c2r_2d(                                          \
        int Nx, int Ny, complex_type *in, real_type *out, MPI_Comm comm,       \
        unsigned flags) {                                                      \
      return prefix##_mpi_plan_dft_c2r_2d(Nx, Ny, in, out, comm, flags);       \
    }                                                                          \
                                                                               \
    static plan_type plan_dft_r2c_3d(                                          \
        int Nx, int Ny, int Nz, real_type *in, complex_type *out,              \
        MPI_Comm comm, unsigned flags) {                                       \
      return prefix##_mpi_plan_dft_r2c_3d(Nx, Ny, Nz, in, out, comm, flags);   \
    }                                                                          \
    static plan_type plan_dft_c2r_3d(                                          \
        int Nx, int Ny, int Nz, complex_type *in, real_type *out,              \
        MPI_Comm comm, unsigned flags) {                                       \
      return prefix##_mpi_plan_dft_c2r_3d(Nx, Ny, Nz, in, out, comm, flags);   \
    }                                                                          \
                                                                               \
    static plan_type plan_dft_r2c(                                             \
        int rank, const ptrdiff_t *n, real_type *in, complex_type *out,        \
        MPI_Comm comm, unsigned flags) {                                       \
      return prefix##_mpi_plan_dft_r2c(rank, n, in, out, comm, flags);         \
    }                                                                          \
    static plan_type plan_dft_c2r(                                             \
        int rank, const ptrdiff_t *n, complex_type *in, real_type *out,        \
        MPI_Comm comm, unsigned flags) {                                       \
      return prefix##_mpi_plan_dft_c2r(rank, n, in, out, comm, flags);         \
    }                                                                          \
    static plan_type plan_dft_3d(                                              \
        int Nx, int Ny, int Nz, complex_type *in, complex_type *out,           \
        MPI_Comm comm, int sign, unsigned flags) {                             \
      return prefix##_mpi_plan_dft_3d(Nx, Ny, Nz, in, out, comm, sign, flags); \
    }                                                                          \
    static plan_type plan_dft_2d(                                              \
        int Nx, int Ny, complex_type *in, complex_type *out, MPI_Comm comm,    \
        int sign, unsigned flags) {                                            \
      return prefix##_mpi_plan_dft_2d(Nx, Ny, in, out, comm, sign, flags);     \
    }                                                                          \
    static plan_type plan_dft_1d(                                              \
        int Nx, complex_type *in, complex_type *out, MPI_Comm comm, int sign,  \
        unsigned flags) {                                                      \
      return prefix##_plan_dft_1d(Nx, in, out, sign, flags);                   \
    }                                                                          \
    static void destroy_plan(plan_type plan) { prefix##_destroy_plan(plan); }  \
  }

  FFTW_MPI_CALLS_BASE(double, fftw);
  FFTW_MPI_CALLS_BASE(float, fftwf);

#undef FFTW_MPI_CALLS_BASE

}; // namespace CosmoTool

#endif
