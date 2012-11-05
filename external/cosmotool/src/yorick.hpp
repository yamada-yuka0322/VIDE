#ifndef __YORICK_HELPERS_HPP
#define __YORICK_HELPERS_HPP

#include "config.hpp"
#include <stdint.h>
#include <fstream>
#include <string>


namespace CosmoTool 
{

  class ProgressiveDoubleOutputImpl
  {
  public:
    virtual ~ProgressiveDoubleOutputImpl();
    virtual void addDouble(double d) = 0;
  };

  template<class T>
  class ProgressiveInputImpl
  {
  public:
    virtual ~ProgressiveInputImpl() {}
    virtual T read() = 0;
    virtual void seek(uint32_t *pos) = 0;
  };

  template<class T>
  class ProgressiveOutputImpl
  {
  public:
    virtual ~ProgressiveOutputImpl() {}
    virtual void put(T a) = 0;
  };

  class ProgressiveDoubleOutput
  {
  private:
    bool initialized;
    int *ref;
    ProgressiveDoubleOutputImpl *impl;
    
    friend ProgressiveDoubleOutput saveDoubleArrayProgressive(const char *fname, uint32_t *dimList, uint32_t rank);
    void decRef();
  public:
    ProgressiveDoubleOutput();
    ProgressiveDoubleOutput(ProgressiveDoubleOutputImpl *i);
    ProgressiveDoubleOutput(const ProgressiveDoubleOutput& o);
    ~ProgressiveDoubleOutput();

    virtual void addDouble(double a);

    const ProgressiveDoubleOutput& operator=(const ProgressiveDoubleOutput& b); 
  };

  template<class T>
  class ProgressiveInput
  {
  private:
    int *ref;
    ProgressiveInputImpl<T> *impl;
     
    void decRef()
    {
      if (ref == 0)
	return;
      
      (*ref)--;
      if (*ref == 0)
	{
      	  delete ref;
	  delete impl;
	}
      impl = 0;
      ref = 0;
    }

  public:
    static ProgressiveInput<T> 
    loadArrayProgressive(const std::string& fname, uint32_t *&dimList, 
			 uint32_t& rank);

    ProgressiveInput() {
      impl = 0; 
      ref = 0; 
    }
    ProgressiveInput(ProgressiveInputImpl<T> *i) {
      impl = i; 
      ref = new int; 
      *ref = 1;
    }
    ProgressiveInput(const ProgressiveInput<T>& o) {
      ref = o.ref;
      impl = o.impl;
      (*ref)++;
    }
    ~ProgressiveInput() {
      decRef();
    }

    T read()
    {
      return impl->read();
    }

    void seek(uint32_t *pos)
    {
      impl->seek(pos);
    }

    const ProgressiveInput<T>& operator=(const ProgressiveInput<T>& b)
    {
      decRef();
      ref = b.ref;
      impl = b.impl;
      if (ref != 0)
	(*ref)++;
      return *this;
    }
  };

  template<class T>
  class ProgressiveOutput
  {
  private:
    int *ref;
    ProgressiveOutputImpl<T> *impl;
     
    void decRef()
    {
      if (ref == 0)
	return;
      
      (*ref)--;
      if (*ref == 0)
	{
      	  delete ref;
	  delete impl;
	}
      impl = 0;
      ref = 0;
    }

  public:
    static ProgressiveOutput<T> 
    saveArrayProgressive(const std::string& fname, uint32_t *dimList, 
			 uint32_t rank);

    ProgressiveOutput() {
      impl = 0; 
      ref = 0; 
    }
    ProgressiveOutput(ProgressiveOutputImpl<T> *i) {
      impl = i; 
      ref = new int; 
      *ref = 1;
    }
    ProgressiveOutput(const ProgressiveOutput<T>& o) {
      ref = o.ref;
      impl = o.impl;
      (*ref)++;
    }
    ~ProgressiveOutput() {
      decRef();
    }

    void put(T a)
    {
      impl->put(a);
    }

    const ProgressiveOutput<T>& operator=(const ProgressiveOutput<T>& b)
    {
      decRef();
      ref = b.ref;
      impl = b.impl;
      if (ref != 0)
	(*ref)++;
      return *this;
    }
  };

  template<typename T>
  void saveArray(const std::string& fname,
		 T *array, uint32_t *dimList, uint32_t rank);

  template<typename T>
  void loadArray(const std::string& fname,
		 T*& array, uint32_t *& dimList, uint32_t& rank)
    throw (NoSuchFileException);

  ProgressiveDoubleOutput saveDoubleArrayProgressive(const char *fname, uint32_t *dimList, uint32_t rank);  
};

#endif
