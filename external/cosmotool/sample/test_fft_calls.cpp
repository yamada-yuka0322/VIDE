#include "fourier/euclidian.hpp"

using namespace CosmoTool;

int main()
{
  EuclidianFourierTransform_2d<double> dft(128,128,1.0,1.0);

  dft.realSpace().eigen().setRandom();
  dft.analysis();
  return 0;
}
