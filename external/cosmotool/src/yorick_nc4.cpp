#include "ctool_netcdf_ver.hpp"
#include "config.hpp"
#include <netcdf>
using namespace netCDF;
#include <fstream>
#include "yorick.hpp"
#include <assert.h>

using namespace CosmoTool;
using namespace std;

class NetCDF_handle
{
public:
  NcFile *outFile;
  NcVar curVar;
  vector<size_t> curPos;
  vector<size_t> counts;
  vector<NcDim> dimList;
  uint32_t rank;
  
  NetCDF_handle(NcFile *f, NcVar v, vector<NcDim>& dimList, uint32_t rank);
  virtual ~NetCDF_handle();
};

NetCDF_handle::NetCDF_handle(NcFile *f, NcVar v, vector<NcDim>& dimList, uint32_t rank)
{
  this->outFile = f;
  this->curVar = v;
  this->dimList = dimList;
  this->rank = rank;
  this->counts.resize(rank);
  this->curPos.resize(rank);

  for (long i = 0; i < rank; i++)
    this->curPos[i] = 0;

  for (long i = 0; i < rank; i++)
    this->counts[i] = 1;
}

NetCDF_handle::~NetCDF_handle()
{
  delete outFile;
}

template<typename T>
class InputGenCDF: public NetCDF_handle, public ProgressiveInputImpl<T>
{
public:
  InputGenCDF(NcFile *f, NcVar v, vector<NcDim>& dimList, uint32_t rank)
    : NetCDF_handle(f,v,dimList,rank)
  {}
  virtual ~InputGenCDF() {} 

  virtual T read()
  {
    T a;

    curVar.getVar(curPos, counts, &a);
    
    curPos[rank-1]++;
    for (long i = rank-1; i >= 1; i--)
      {
	if (curPos[i] == dimList[i].getSize())
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
  OutputGenCDF(NcFile *f, NcVar v, vector<NcDim>& dimList, uint32_t rank)
    : NetCDF_handle(f,v,dimList,rank)
  {}
  virtual ~OutputGenCDF() {} 

  virtual void put(T a)
  {
    curVar.putVar(curPos, counts, &a);
    
    curPos[rank-1]++;
    for (long i = rank-1; i >= 1; i--)
      {
	if (curPos[i] == dimList[i].getSize())
	  {
	    curPos[i-1]++;
	  curPos[i] = 0;
	}
      }    
  }
};

template<typename T> NcType& get_NetCDF_type();

#define IMPL_TYPE(T,ncT) \
template<> \
NcType& get_NetCDF_type<T>() \
{ \
  return ncT; \
}

IMPL_TYPE(int,ncInt);
IMPL_TYPE(unsigned int,ncInt);
IMPL_TYPE(double,ncDouble);
IMPL_TYPE(float,ncFloat);

namespace CosmoTool {
  template<typename T>
  ProgressiveOutput<T>
  ProgressiveOutput<T>::saveArrayProgressive(const std::string& fname, uint32_t *dimList,
					     uint32_t rank)
  {
    NcFile *f = new NcFile(fname, NcFile::replace);
    
    vector<NcDim> dimArray;
    for (uint32_t i = 0; i < rank; i++)
      {
	char dimName[255];
	
	sprintf(dimName, "dim%d", i);
	dimArray.push_back(f->addDim(dimName, dimList[rank-1-i]));
      }
    
    NcVar v = f->addVar("array", get_NetCDF_type<T>(), dimArray);  

    vector<NcDim> ldimList;

    for (uint32_t i = 0; i < rank; i++)
      ldimList.push_back(dimArray[rank-1-i]);

    OutputGenCDF<T> *impl = new OutputGenCDF<T>(f, v, ldimList, rank);
    return ProgressiveOutput<T>(impl);   
  }

  template<typename T>
  ProgressiveInput<T>
  ProgressiveInput<T>::loadArrayProgressive(const std::string& fname, uint32_t *&dimList,
					      uint32_t& rank)
  {
    NcFile *f = new NcFile(fname, NcFile::read);

    NcVar v = f->getVar("array");

    rank = v.getDimCount();
    vector<NcDim> vdimlist = v.getDims();
    dimList = new uint32_t[rank];
    for (uint32_t i = 0; i < rank; i++)
      {
	dimList[rank-i-1] = vdimlist[i].getSize();
      }
    InputGenCDF<T> *impl = new InputGenCDF<T>(f, v, vdimlist, rank);

    return ProgressiveInput<T>(impl);
  }

  template<typename T>
  void saveArray(const std::string& fname,
		 T *array, uint32_t *dimList, uint32_t rank)
  {
    NcFile f(fname.c_str(), NcFile::replace);
    
    vector<NcDim> dimArray;
    for (uint32_t i = 0; i < rank; i++)
      {
	char dimName[255];
	
	sprintf(dimName, "dim%d", i);
	dimArray.push_back(f.addDim(dimName, dimList[i]));
      }

    NcVar v = f.addVar("array", get_NetCDF_type<T>(), dimArray);
    
    v.putVar(array);
  }

  template<typename T>
  void loadArray(const std::string& fname,
		 T*&array, uint32_t *&dimList, uint32_t& rank)
	throw (NoSuchFileException)
  {
    NcFile f(fname.c_str(), NcFile::read);
    
    //if (!f.is_valid())
    //  throw NoSuchFileException(fname);
    
    NcVar v = f.getVar("array"); 
    vector<NcDim> dims = v.getDims();
    rank = v.getDimCount();
    uint32_t fullSize = 1;
    dimList = new uint32_t[rank];
    for (int i = 0; i < rank; i++)
      {
	dimList[i] = dims[i].getSize();
	fullSize *= dimList[i];
      }
    if (fullSize != 0) {
      array = new T[fullSize];
      v.getVar(array);
    }
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
