#ifndef __FIX_ARRAY_HPP
#define __FIX_ARRAY_HPP

namespace CosmoTool
{

  template <typename T, unsigned int sz> class fixArray
  {
  private:
    T d[sz];
    
  public:
    /*! Returns the size of the array. */
    long size() const { return sz; }

    /*! Returns a reference to element \a #n */
    template<typename T2> T &operator[] (T2 n) {
      return d[n];
    }

    /*! Returns a constant reference to element \a #n */
    template<typename T2> const T &operator[] (T2 n) const {
      return d[n];
    }

    template<typename T2> void importArray(T2 *indata) {
      for (int i = 0; i < sz; i++)
	d[i] = indata[i];
    }

    template<typename T2> void exportArray(T2 *outdata) {
      for (int i = 0; i < sz; i++)
	outdata[i] = d[i];
    }

  };

};

#endif
