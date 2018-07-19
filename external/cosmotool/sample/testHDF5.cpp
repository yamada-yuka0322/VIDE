/*+
This is CosmoTool (./sample/testHDF5.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
#include <iostream>
#include "hdf5_array.hpp"
#include <H5Cpp.h>

using namespace std;

struct MyStruct
{
  int a;
  double b;
  char c;
};

struct MyStruct2
{
   MyStruct base;
   int d;
};

enum MyColors
{
  RED, GREEN, BLUE
};

CTOOL_STRUCT_TYPE(MyStruct, hdf5t_MyStruct,
    ((int, a))
    ((double, b))
    ((char, c))
)

CTOOL_STRUCT_TYPE(MyStruct2, hdf5t_MyStruct2,
    ((MyStruct, base))
    ((int, d))
)

CTOOL_ENUM_TYPE(MyColors, hdf5t_MyColors,
   (RED) (GREEN) (BLUE)
)

int main()
{
  typedef  boost::multi_array<float, 2> array_type;
  typedef  boost::multi_array<float, 3> array3_type;
  typedef  boost::multi_array<MyStruct, 1> array_mys_type;
  typedef  boost::multi_array<MyColors, 1> array_mys_color;
  typedef  boost::multi_array<bool, 1> array_mys_bool;
  typedef  boost::multi_array<MyStruct2, 1> array_mys2_type;
  typedef  boost::multi_array<std::complex<double>, 2> arrayc_type;
  typedef array_type::index index;

  H5::H5File f("test.h5", H5F_ACC_TRUNC);
  
  H5::Group g = f.createGroup("test_group");
  
  array_type A(boost::extents[2][3]);
  array_type B, Bprime(boost::extents[1][2]);
  array3_type C(boost::extents[2][3][4]);
  arrayc_type D, E;
  array_mys_type F(boost::extents[10]), G;
  array_mys2_type H(boost::extents[10]);
  array_mys_color I(boost::extents[2]);
  array_mys_bool J(boost::extents[2]);
  
  I[0] = RED;
  I[1] = BLUE;
  J[0] = false;
  J[1] = true;
  
  int values = 0;
  for (index i = 0; i != 2; i++)
    for (index j = 0; j != 3; j++)
      A[i][j] = values++;
      
  for (index i = 0; i != 10; i++)
    {
        F[i].a = i;
        F[i].b = double(i)/4.;
        F[i].c = 'r'+i;
        H[i].base = F[i];
        H[i].d = 2*i;
    } 
  std::cout << " c = " << ((char *)&F[1])[offsetof(MyStruct, c)] << endl;
    
  CosmoTool::hdf5_write_array(g, "test_data", A);
  CosmoTool::hdf5_write_array(g, "test_struct", F);
  CosmoTool::hdf5_write_array(g, "test_struct2", H);
  CosmoTool::hdf5_write_array(g, "colors", I);
  CosmoTool::hdf5_write_array(g, "bools", J);
  CosmoTool::hdf5_read_array(g, "test_data", B);

  int verify = 0;
  for (index i = 0; i != 2; i++)
    for (index j = 0; j != 3; j++)
      if (B[i][j] != verify++) {
        std::cout << "Invalid array content" << endl;
        abort();
      }
  
  std::cout << "Testing C " << std::endl;
  try
    {
      CosmoTool::hdf5_read_array(g, "test_data", C);
      std::cout << "Did not throw InvalidDimensions" << endl;
      abort();
    }
  catch (const CosmoTool::InvalidDimensions&)
    {}

  std::cout << "Testing Bprime " << std::endl;
  try
    {
      CosmoTool::hdf5_read_array(g, "test_data", Bprime, false, true);
      for (index i = 0; i != 1; i++)
        for (index j = 0; j != 2; j++)
          if (B[i][j] != Bprime[i][j]) {
             std::cout << "Invalid array content in Bprime" << endl;
             abort();
          }
    }
  catch (const CosmoTool::InvalidDimensions&)
    {
      std::cout << "Bad! Dimensions should be accepted" << std::endl;
      abort();
    }
  
  D.resize(boost::extents[2][3]);
  D = A;
  
  CosmoTool::hdf5_write_array(g, "test_data_c", D);

  CosmoTool::hdf5_read_array(g, "test_data_c", E);

  verify = 0;
  for (index i = 0; i != 2; i++)
    for (index j = 0; j != 3; j++)
      if (E[i][j].real() != verify++) {
        std::cout << "Invalid array content" << endl;
        abort();
      }

  CosmoTool::hdf5_read_array(g, "test_struct", G);
  for (index i = 0; i != 10; i++)
    if (G[i].a != F[i].a || G[i].b != F[i].b || G[i].c != F[i].c) {
        std::cout << "Invalid struct content" << endl;
        abort();
    }

  return 0;
}
