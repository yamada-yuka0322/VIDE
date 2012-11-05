#include "ctool_netcdf_ver.hpp"
#include "config.hpp"
#ifdef NETCDFCPP4
#include <netcdf>
using namespace netCDF
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
    NcFile *f = new NcFile(fname.c_str(), NcFile::Replace);
    
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
		 T *array, uint32_t *dimList, uint32_t rank)
  {
    NcFile f(fname.c_str(), NcFile::Replace);
    
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
			       int *array, uint32_t *dimList, uint32_t rank);
  template void saveArray<float>(const std::string& fname,
				 float *array, uint32_t *dimList, uint32_t rank);
  template void saveArray<double>(const std::string& fname,
				  double *array, uint32_t *dimList, uint32_t rank);
  
}
