#ifndef __COSMOTOOL_NUMPY_ADAPTOR_HPP
#define __COSMOTOOL_NUMPY_ADAPTOR_HPP

namespace CosmoTool {

  template<typename T, typename IT>
  void parallel_ufunc_dd_d(char **args, IT* dimensions, IT* steps, void *func) {
    IT i;
    IT n = dimensions[0];
    char *in = args[0], *in2 = args[1], *out = args[2];
    IT in_step = steps[0], in2_step = steps[1], out_step = steps[2];

    double tmp;
    typedef double (*F_t)(double,double);

    F_t f = (F_t)func;

#pragma omp parallel for schedule(static)
    for (i = 0; i < n; i++) {
        T *out_t = (T *)(out + i * out_step);
        T *in_t = (T *)(in + i * in_step);
        T *in2_t = (T *)(in2 + i * in2_step);
        *out_t = f(*in_t, *in2_t);
    }
  }


}

#endif
