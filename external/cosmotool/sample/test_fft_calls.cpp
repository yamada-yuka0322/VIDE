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
  return 0;
}
