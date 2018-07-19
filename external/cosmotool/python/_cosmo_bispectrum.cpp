#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>
#include <boost/format.hpp>
#include <boost/multi_array.hpp>
#include <sys/types.h>
#include <cmath>
#include "symbol_visible.hpp"
#include "algo.hpp"

using std::cout;
using std::endl;
using boost::format;
using CosmoTool::square;


struct ModeSet
{
    ssize_t N1, N2, N3;
    bool half_copy;
    
    struct TriangleIterator
    {
        ssize_t i1, i2, i3;
        ssize_t N1, N2, N3;
        ssize_t first_iteration;
        
        
        TriangleIterator& operator++() {
            i3++;
            if (i3==(N3/2+1)) { i3 = first_iteration; i2++; }
            if (i2==(N2/2+1)) { i2 = -N2/2; i1++; }
            return *this;
        }
    
        bool operator!=(const TriangleIterator& t) const {
            return i1!=t.i1 || i2!=t.i2 || i3 != t.i3;
        }
        
        bool in_box() const {
            ssize_t hN1 = N1/2, hN2 = N2/2, hN3 = N3/2;
            
            return (i1 >= -hN1) && (i1 <= hN1) && (i2 >= -hN2) && (i2 <= hN2) && (i3 >= -hN3) && (i3 <= hN3);
        }
        
        TriangleIterator operator+(const TriangleIterator& other_t) const {
            TriangleIterator t = *this;
            
            t.i1 = (t.i1+other_t.i1);
            t.i2 = (t.i2+other_t.i2);
            t.i3 = (t.i3+other_t.i3);
            return t;
        }

        TriangleIterator& inp_array() {
            if (i1 < 0)
                i1 += N1;  
            if (i2 < 0)
                i2 += N2;  
            if (i3 < 0)
                i3 += N3;  
            return *this;
        }

        TriangleIterator array() const {
          TriangleIterator t = *this;
          
          t.inp_array();
          return t;
        }
    
        TriangleIterator real() const {
            TriangleIterator t = *this;
            if (t.i1 >= N1/2)
                t.i1 -= N1;  
            if (t.i2 >= N2/2)
                t.i2 -= N2;  
            if (t.i3 >= N3/2)
                t.i3 -= N3;  
            return t;
        }
        
        double norm() const {
            double r1 = i1, r2 = i2, r3 = i3;
            return std::sqrt(r1*r1 + r2*r2 + r3*r3);
        }
        void reverse() { i1=-i1; i2=-i2; i3=-i3; }
  
        TriangleIterator& operator*() { return *this; }
    };

    ModeSet(size_t N1_, size_t N2_, size_t N3_, bool _half_copy = false)
        : N1(N1_), N2(N2_), N3(N3_),half_copy(_half_copy) {
    } 

    TriangleIterator begin() const {
        TriangleIterator t;
        t.i1 = -N1/2;
        t.i2 = -N2/2;
        if (half_copy)
          t.first_iteration = t.i3 = 0;
        else
          t.first_iteration = t.i3 = -N3/2;
        t.N1 = N1;
        t.N2 = N2;
        t.N3 = N3;
        return t;
    }

    TriangleIterator end() const {
        TriangleIterator t;
        t.first_iteration = (half_copy ? 0 : (-N3/2));
        t.i3 = t.first_iteration;
        t.i2 = -N2/2;
        t.i1 = N1/2+1;
        t.N1 = N1;
        t.N2 = N2;
        t.N3 = N3;
        return t;
    }

};

std::ostream& operator<<(std::ostream& o, const ModeSet::TriangleIterator& t)
{
  o << t.i1 << "," << t.i2 << "," << t.i3;
  return o;
}

template<typename T>
static T no_conj(const T& a) { return a; }

template<typename SubArrayB,typename SubArrayCnt,typename Delta>
static inline void accum_bispec(const Delta& delta_mirror, SubArrayB b_Nt, SubArrayCnt b_B,
    const typename Delta::element& v1, const typename Delta::element& v2,
    const ModeSet::TriangleIterator& rm1, 
    const ModeSet::TriangleIterator& rm2, 
    const ModeSet::TriangleIterator& rm3,
    double delta_k,
    size_t Nk)
{
  typedef std::complex<double> CType;

  size_t q1 = std::floor(rm1.norm()/delta_k);
  if (q1 >= Nk)
    return;
    
  size_t q2 = std::floor(rm2.norm()/delta_k);
  if (q2 >= Nk)
    return;

  size_t q3 = std::floor(rm3.norm()/delta_k);
  if (q3 >= Nk)
    return;

  CType prod = v1*v2;
  ModeSet::TriangleIterator m3 = rm3;
  // We use hermitic symmetry to get -m3, it is just the mode in m3 but conjugated.
  m3.reverse();
  m3.inp_array();
  prod *= delta_mirror[m3.i1][m3.i2][m3.i3];

  b_Nt[q1][q2][q3] ++;
  b_B[q1][q2][q3] += prod;
}

extern "C" CTOOL_DLL_PUBLIC
void CosmoTool_compute_bispectrum(
    double *delta_hat, size_t Nx, size_t Ny, size_t Nz,
    size_t *Ntriangles,
    double* B, double delta_k, size_t Nk ) 
{
    // First remap to multi_array for easy access
    size_t kNz = Nz/2+1;
#ifdef _OPENMP
    int Ntasks = omp_get_max_threads();
#else
    int Ntasks = 1;
#endif
    boost::multi_array_ref<std::complex<double>, 3> a_delta(reinterpret_cast<std::complex<double>*>(delta_hat), boost::extents[Nx][Ny][kNz]); 
    boost::multi_array_ref<size_t, 3> a_Nt(Ntriangles, boost::extents[Nk][Nk][Nk]); 
    boost::multi_array_ref<std::complex<double>, 3> a_B(reinterpret_cast<std::complex<double>*>(B), boost::extents[Nk][Nk][Nk]); 
    boost::multi_array<std::complex<double>, 4> b_B(boost::extents[Ntasks][Nk][Nk][Nk]);
    boost::multi_array<size_t, 4> b_Nt(boost::extents[Ntasks][Nk][Nk][Nk]);
    typedef std::complex<double> CType;
    boost::multi_array<std::complex<double>, 3> delta_mirror(boost::extents[Nx][Ny][Nz]); 
    
    // Add hermiticity
    for (auto m : ModeSet(Nx, Ny, Nz, true)) {
      auto n1 = m;
      auto n2 = m.array();
      n1.reverse();
      n1.inp_array();
      delta_mirror[n2.i1][n2.i2][n2.i3] = (a_delta[n2.i1][n2.i2][n2.i3]);
      delta_mirror[n1.i1][n1.i2][n1.i3] = std::conj(delta_mirror[n2.i1][n2.i2][n2.i3]);
    }

#ifdef _OPENMP
    // First loop over m1
#pragma omp parallel
    {
#pragma omp single
      {
        for (auto m1 : ModeSet(Nx, Ny, Nz)) {
          auto am1 = m1.array();
          CType v1 = delta_mirror[am1.i1][am1.i2][am1.i3];
          int tid = omp_get_thread_num();
          
#pragma omp task
          {
            auto rm1 = m1.real();
            // Second mode m2
            for (auto m2 : ModeSet(Nx, Ny, Nz)) {
                // Now derive m3
                auto am2 = m2.array();
                auto m3 = (m1+m2);
                CType v2 = delta_mirror[am2.i1][am2.i2][am2.i3];

                // Not in Fourier box, stop here     
                if (!m3.in_box())
                    continue;

                accum_bispec(delta_mirror, b_Nt[tid], b_B[tid], v1, v2, m1, m2, m3, delta_k, Nk);

            }
          }          

        }
      }
    }

#pragma omp taskwait
  for (int tid = 0; tid < Ntasks; tid++) {
    size_t *b_p = b_Nt[tid].origin();
    size_t *a_p = a_Nt.data();
    std::complex<double> *b_B_p = b_B[tid].origin();
    std::complex<double> *a_B_p = a_B.origin();
//#pragma omp simd
#pragma omp parallel for
    for (size_t q = 0; q < Nk*Nk*Nk; q++) {
      a_p[q] += b_p[q];
      a_B_p[q] += b_B_p[q];
    }
  }
#else
#warning Serial version not implemented
#endif
}


extern "C" CTOOL_DLL_PUBLIC
void CosmoTool_compute_powerspectrum(
    double *delta_hat, size_t Nx, size_t Ny, size_t Nz,
    size_t *Ncounts,
    double* P, double delta_k, size_t Nk ) 
{
    // First remap to multi_array for easy access
    size_t kNz = Nz/2+1;
    boost::multi_array_ref<std::complex<double>, 3> a_delta(reinterpret_cast<std::complex<double>*>(delta_hat), boost::extents[Nx][Ny][kNz]); 
    boost::multi_array_ref<size_t, 1> a_Nc(Ncounts, boost::extents[Nk]); 
    boost::multi_array_ref<double, 1> a_P(reinterpret_cast<double*>(P), boost::extents[Nk]); 
    typedef std::complex<double> CType;

    // First loop over m1
    for (auto m : ModeSet(Nx, Ny, kNz)) {
        auto m1 = m.array();
        CType& v1 = a_delta[m1.i1][m1.i2][m1.i3];

        size_t q1 = std::floor(m.norm()/delta_k);

        if (q1 >= Nk)
            continue;

        a_Nc[q1] ++;
        a_P[q1] += std::norm(v1);
    }
}
