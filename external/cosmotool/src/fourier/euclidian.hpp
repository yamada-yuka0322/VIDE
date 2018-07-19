/*+
This is CosmoTool (./src/fourier/euclidian.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __COSMOTOOL_FOURIER_EUCLIDIAN_HPP
#define __COSMOTOOL_FOURIER_EUCLIDIAN_HPP

#include <cmath>
#include <boost/function.hpp>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <gsl/gsl_randist.h>
#include "base_types.hpp"
#include "fft/fftw_calls.hpp"
#include "../algo.hpp"
#include "details/euclidian_maps.hpp"
#include "details/euclidian_spectrum_1d.hpp"
#include "details/euclidian_spectrum_1d_bin.hpp"
#include "details/euclidian_transform.hpp"

#endif
