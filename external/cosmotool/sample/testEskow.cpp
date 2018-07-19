/*+
This is CosmoTool (./sample/testEskow.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#include <fstream>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
#include "eskow.hpp"

using namespace std;

double Hartmann_Matrix[6][6] = {
  { 14.8253,   -6.4243,   7.8746,  -1.2498,  10.2733,  10.2733 },
  { -6.4243,   15.1024,  -1.1155,  -0.2761,  -8.2117,  -8.2117 },
  {  7.8746,   -1.1155,  51.8519, -23.3482,  12.5902,  12.5902 },
  { -1.2498,   -0.2761, -23.3482,  22.7962,  -9.8958,  -9.8958 },
  { 10.2733,   -8.2117,  12.5902,  -9.8958,  21.0656,  21.0656 },
  { 10.2733,   -8.2117,  12.5902,  -9.8958,  21.0656,  21.0656 }
};

struct MatrixOp
{
  vector<double> M;
  int N;

  double& operator()(int i, int j)
  {
    return M[i*N + j];
  }
};

int main(int argc, char **argv)
{
  MatrixOp M;
  double norm_E;
  ifstream fi(argv[1]);
  ofstream f("eskowed.txt");
  CholeskyEskow<double,MatrixOp> chol;

  fi >> M.N;
  M.M.resize(M.N*M.N);

  for (int i = 0; i < M.N; i++)
    {
      for (int j = 0; j < M.N; j++)
	{
	  fi >> M(i,j);
	  if (j > i)
	    M(i,j) =0;
	}
    }

  chol.cholesky(M, M.N, norm_E);

  cout << "norm_E = " << norm_E << endl;

  for (int i = 0; i < M.N; i++)
    {
      for (int j = 0; j < M.N; j++)
	{
	  if (j > i) 
            f << "0 ";
          else
 	    f << setprecision(25) << M(i,j) << " ";
	}
      f << endl;
    }
  return 0;
}
