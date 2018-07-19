/*+
This is CosmoTool (./src/yorick_nc3.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#include "ctool_netcdf_ver.hpp"
#include "config.hpp"
#ifdef NETCDFCPP4
#include <netcdf>
using namespace netCDF;
#else
#include <netcdfcpp.h>
#endif
#include <fstream>
#include "yorick.hpp"
#include <assert.h>

using namespace CosmoTool;
using namespace std;

class NetCDF_handle
{
public:
  NcFile *outFile;
  NcVar *curVar;
  long *curPos;
  long *counts;
  long *dimList;
  uint32_t rank;
  
  NetCDF_handle(NcFile *f, NcVar *v, long *dimList, uint32_t rank);
  virtual ~NetCDF_handle();
};

NetCDF_handle::NetCDF_handle(NcFile *f, NcVar *v, long *dimList, uint32_t rank)
{
  this->outFile = f;
  this->curVar = v;
  this->dimList = dimList;
  this->rank = rank;
  this->counts = new long[rank];
  this->curPos = new long[rank];

  for (long i = 0; i < rank; i++)
    this->curPos[i] = 0;

  for (long i = 0; i < rank; i++)
    this->counts[i] = 1;
}

NetCDF_handle::~NetCDF_handle()
{
  delete[] dimList;
  delete outFile;
}

template<typename T>
class InputGenCDF: public NetCDF_handle, public ProgressiveInputImpl<T>
{
public:
  InputGenCDF(NcFile *f, NcVar *v, long *dimList, uint32_t rank)
    : NetCDF_handle(f,v,dimList,rank)
  {}
  virtual ~InputGenCDF() {} 

  virtual T read()
  {
    T a;

    curVar->set_cur(curPos);
    curVar->get(&a, counts);
    
    curPos[rank-1]++;
    for (long i = rank-1; i >= 1; i--)
      {
	if (curPos[i] == dimList[i])
	  {
	    curPos[i-1]++;
	  curPos[i] = 0;
	}
      }    
    return a;
  }

  virtual void seek(uint32_t *pos)
  {
    for (long i = rank-1; i >= 0; i--)
      curPos[i] = pos[rank-1-i];
  }
};

template<typename T>
class OutputGenCDF: public NetCDF_handle, public ProgressiveOutputImpl<T>
{
public:
  OutputGenCDF(NcFile *f, NcVar *v, long *dimList, uint32_t rank)
    : NetCDF_handle(f,v,dimList,rank)
  {}
  virtual ~OutputGenCDF() {} 

  virtual void put(T a)
  {
    curVar->set_cur(curPos);
    curVar->put(&a, counts);
    
    curPos[rank-1]++;
    for (long i = rank-1; i >= 1; i--)
      {
	if (curPos[i] == dimList[i])
	  {
	    curPos[i-1]++;
	  curPos[i] = 0;
	}
      }    
  }
};

template<typename T>
class NetCDF_type
{
public:
  static const NcType t = (NcType)-1;
};

template<>
class NetCDF_type<int>
{
public:
  static const NcType t = ncInt;
};

template<>
class NetCDF_type<unsigned int>
{
public:
  static const NcType t = ncInt;
};

template<>
class NetCDF_type<float>
{
public:
  static const NcType t = ncFloat;
};

template<>
class NetCDF_type<double>
{
public:
  static const NcType t = ncDouble;
};

namespace CosmoTool {
  template<typename T>
  ProgressiveOutput<T>
  ProgressiveOutput<T>::saveArrayProgressive(const std::string& fname, uint32_t *dimList,
					     uint32_t rank)
  {
    NcFile *f = new NcFile(fname.c_str(), NcFile::Replace, 0, 0, NcFile::Netcdf4);
    
    assert(f->is_valid());
    
    const NcDim **dimArray = new const NcDim *[rank];
    for (uint32_t i = 0; i < rank; i++)
      {
	char dimName[255];
	
	sprintf(dimName, "dim%d", i);
	dimArray[i] = f->add_dim(dimName, dimList[rank-1-i]);
      }
    
    NcVar *v = f->add_var("array", NetCDF_type<T>::t, rank, dimArray);  

    long *ldimList = new long[rank];

    for (uint32_t i = 0; i < rank; i++)
      ldimList[rank-1-i] = dimList[i];

    OutputGenCDF<T> *impl = new OutputGenCDF<T>(f, v, ldimList, rank);
    return ProgressiveOutput<T>(impl);   
  }

  template<typename T>
  ProgressiveInput<T>
  ProgressiveInput<T>::loadArrayProgressive(const std::string& fname, uint32_t *&dimList,
					      uint32_t& rank)
  {
    NcFile *f = new NcFile(fname.c_str(), NcFile::ReadOnly);

    assert(f->is_valid());

    NcVar *v = f->get_var("array");

    rank = v->num_dims();
    long *ldimList = v->edges();
    dimList = new uint32_t[rank];
    for (uint32_t i = 0; i < rank; i++)
      {
	dimList[rank-i-1] = ldimList[i];
      }
    InputGenCDF<T> *impl = new InputGenCDF<T>(f, v, ldimList, rank);

    return ProgressiveInput<T>(impl);
  }

  template<typename T>
  void saveArray(const std::string& fname,
		 const T *array, uint32_t *dimList, uint32_t rank)
  {
    NcFile f(fname.c_str(), NcFile::Replace, 0, 0, NcFile::Netcdf4);
    
    assert(f.is_valid());
    
    const NcDim **dimArray = new const NcDim *[rank];
    for (uint32_t i = 0; i < rank; i++)
      {
	char dimName[255];
	
	sprintf(dimName, "dim%d", i);
	dimArray[i] = f.add_dim(dimName, dimList[i]);
      }

    NcVar *v = f.add_var("array", NetCDF_type<T>::t, rank, dimArray);
    
    long *edge = v->edges();
    v->put(array, edge);
    delete[] edge;  
  }

  template<typename T>
  void loadArray(const std::string& fname,
		 T*&array, uint32_t *&dimList, uint32_t& rank)
	throw (NoSuchFileException)
  {
    NcFile f(fname.c_str(), NcFile::ReadOnly);
    
    if (!f.is_valid())
      throw NoSuchFileException(fname);
    
    NcVar *v = f.get_var("array"); 
    rank = v->num_dims();
    long *edge = v->edges();
    uint32_t fullSize = 1;
    dimList = new uint32_t[rank];
    for (int i = 0; i < rank; i++)
      {
	dimList[i] = edge[i];
	fullSize *= edge[i];
      }
    if (fullSize != 0) {
      array = new T[fullSize];
      v->get(array, edge);
    }
    delete[] edge;   
  }

  template class ProgressiveInput<int>;
  template class ProgressiveInput<float>;
  template class ProgressiveInput<double>;

  template class ProgressiveOutput<int>;
  template class ProgressiveOutput<float>;
  template class ProgressiveOutput<double>;

  template void loadArray<int>(const std::string& fname,
			       int*& array, uint32_t *&dimList, uint32_t& rank);
  template void loadArray<float>(const std::string& fname,
				 float*& array, uint32_t *&dimList, uint32_t& rank);
  template void loadArray<double>(const std::string& fname,
				  double*& array, uint32_t *&dimList, uint32_t& rank);

  template void saveArray<int>(const std::string& fname,
			       const int *array, uint32_t *dimList, uint32_t rank);
  template void saveArray<float>(const std::string& fname,
				 const float *array, uint32_t *dimList, uint32_t rank);
  template void saveArray<double>(const std::string& fname,
				  const double *array, uint32_t *dimList, uint32_t rank);
  
}
