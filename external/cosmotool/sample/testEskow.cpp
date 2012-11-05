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
