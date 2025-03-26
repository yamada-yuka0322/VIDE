#ifndef __COSMOTOOL_FFT_COMPLEX_HPP
#define __COSMOTOOL_FFT_COMPLEX_HPP

#include <complex>
#include <fftw3.h>

namespace CosmoTool
{
    template<typename T>
    struct adapt_complex {
    };

    template<> struct adapt_complex<fftw_complex> {
        typedef fftw_complex f_type;
        typedef std::complex<double> cpp_complex;
        
        static inline cpp_complex *adapt(f_type *a) {
            return reinterpret_cast<cpp_complex *>(a);
        }
    };
    
    template<> struct adapt_complex<fftwf_complex> {
        typedef fftwf_complex f_type;
        typedef std::complex<float> cpp_complex;
        
        static inline cpp_complex *adapt(f_type *a) {
            return reinterpret_cast<cpp_complex *>(a);
        }
    };

    template<> struct adapt_complex<fftwl_complex> {
        typedef fftwl_complex f_type;
        typedef std::complex<long double> cpp_complex;
        
        static inline cpp_complex *adapt(f_type *a) {
            return reinterpret_cast<cpp_complex *>(a);
        }
    };

}

#endif