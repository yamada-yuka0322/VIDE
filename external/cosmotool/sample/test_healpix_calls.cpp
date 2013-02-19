#include <iostream>
#include "fourier/healpix.hpp"

using namespace CosmoTool;
using namespace std;

int main()
{
  HealpixFourierTransform<double> dft(128,2*128,2*128, 40);
  HealpixUtilityFunction<double> utils;
  long Npix = dft.realSpace().size();

  dft.realSpace().eigen().setRandom();
  dft.analysis();
  cout << "Map dot-product = " << dft.realSpace().dot_product(dft.realSpace()) << endl;
  cout << "Fourier dot-product = " << dft.fourierSpace().dot_product(dft.fourierSpace()).real()*Npix/(4*M_PI) << endl;
  cout << "Spectrum analysis" << endl;
  HealpixUtilityFunction<double>::Spectrum_ptr s = utils.estimateSpectrumFromMap(dft.fourierSpace());

  dft.synthesis();
  cout << "Resynthesis dot-product = " << dft.realSpace().dot_product(dft.realSpace()) << endl;

  return 0;

}
