#ifndef __CTOOL_OPENMP_HPP
#define __CTOOL_OPENMP_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

namespace CosmoTool {

    static int smp_get_max_threads() {
#ifdef _OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    }
    
    static int smp_get_thread_id() {
#ifdef _OPENMP
        return omp_get_thread_num();
#else
        return 0;
#endif
    }

    static int smp_get_num_threads() {
#ifdef _OPENMP
        return omp_get_num_threads();
#else
        return 1;
#endif

    }
    
    static void smp_set_nested(bool n) {
#ifdef _OPENMP
        omp_set_nested(n ? 1 : 0);
#endif
    }

    
};

#endif
