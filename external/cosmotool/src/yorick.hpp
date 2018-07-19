/*+
This is CosmoTool (./src/yorick.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
		 const T *array, uint32_t *dimList, uint32_t rank);

  template<typename T>
  void loadArray(const std::string& fname,
		 T*& array, uint32_t *& dimList, uint32_t& rank)
    throw (NoSuchFileException);

  ProgressiveDoubleOutput saveDoubleArrayProgressive(const char *fname, uint32_t *dimList, uint32_t rank);  
};

#endif
