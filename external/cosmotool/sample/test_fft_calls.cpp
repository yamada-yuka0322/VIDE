/*+
This is CosmoTool (./sample/test_fft_calls.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#include <boost/format.hpp>
#include "yorick.hpp"
#include <gsl/gsl_rng.h>
#include <iostream>
#include <cmath>
#include "fourier/euclidian.hpp"

using namespace CosmoTool;
using boost::format;
using boost::str;
using namespace std;

double spectrum_generator(double k)
{
  if (k==0)
    return 0;
  return 1/(0.01+pow(k, 3.0));
}

template<typename T>
void test_2d(int Nx, int Ny)
{
  EuclidianFourierTransform_2d<T> dft(Nx,Ny,1.0,1.0);
  EuclidianSpectrum_1D<T> spectrum(spectrum_generator);
  double volume = Nx*Ny;
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

  dft.realSpace().eigen().setRandom();
  dft.fourierSpace().scale(std::complex<T>(0,0));
  dft.analysis();
  cout << format("Testing (%d,%d)") % Nx % Ny << endl;
  cout << "Map dot-product = " << dft.realSpace().dot_product(dft.realSpace()) << endl;
  cout << "Fourier dot-product = " << dft.fourierSpace().dot_product(dft.fourierSpace()).real()*volume << endl;
  dft.synthesis();
  cout << "Resynthesis dot-product = " << dft.realSpace().dot_product(dft.realSpace()) << endl;

  dft.realSpace().scale(2.0);
  dft.fourierSpace().scale(0.2);

  spectrum.newRandomFourier(rng, dft.fourierSpace()); 
  dft.synthesis();

  uint32_t dims[2] = { Ny, Nx };
  CosmoTool::saveArray(str(format("generated_map_%d_%d.nc") %Nx % Ny) , dft.realSpace().data(), dims, 2);

  gsl_rng_free(rng);
}

int main(int argc, char **argv)
{
  test_2d<double>(128,128);
  test_2d<double>(131,128);  
  test_2d<double>(130,128);  
  test_2d<float>(128,128);
  test_2d<float>(128,131);  
  test_2d<float>(128,130);  
  test_2d<float>(131,128);  
  return 0;
}
