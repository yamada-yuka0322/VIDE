#include <iostream>
#include <cmath>

using namespace std;

#include "newton.hpp"

using namespace CosmoTool;

struct SimplePolynom
{
  double eval(double x) { return x*x*x+2*x-1; }
  double derivative(double x) { return 3*x*x+2; }
};

int main()
{
  SimplePolynom poly;
  double r;

  double solution = -2*pow(2./(3*(9+sqrt(177.))), 1./3) + (pow(0.5*(9+sqrt(177.))/9, 1./3));

  r = newtonSolver(3.0, poly);

  cout << "Result = " << r << " delta = " << r-solution << endl;
  return 0;
}
