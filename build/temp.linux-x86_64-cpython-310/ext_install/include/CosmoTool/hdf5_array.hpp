/*+
This is CosmoTool (./src/hdf5_array.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
sthat may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
+*/
#ifndef __COSMO_HDF5_ARRAY_HPP
#define __COSMO_HDF5_ARRAY_HPP


#include <string>
#include <stdint.h>
#include <boost/static_assert.hpp>
#include <boost/utility.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/multi_array.hpp>
#include <boost/type_traits/is_same.hpp> 
#include <boost/type_traits/is_floating_point.hpp> 
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <vector>
#include <H5Cpp.h>

namespace CosmoTool {
#if H5_VERSION_GE(1,8,20)
#if H5_VERSION_GE(1,10,1)
  typedef H5::H5Location H5_CommonFileGroup;
#else
  typedef H5::Group H5_CommonFileGroup;
#endif
#else
  typedef H5::CommonFG H5_CommonFileGroup;
#endif

  //!_______________________________________________________________________________________
  //!     
  //!     map types to HDF5 types
  //!         
  //!     
  //!     Leo Goodstadt (04 March 2013), improved with enable_if by Guilhem Lavaux (May 2014)
  //!_______________________________________________________________________________________ 

  template<typename T, class Enable = void> struct get_hdf5_data_type
  {
    static H5::DataType type()  
      {   
          BOOST_MPL_ASSERT_MSG(0, Unknown_HDF5_data_type, ()); 
          return H5::PredType::NATIVE_DOUBLE; 
      }
  };
  
  #define HDF5_TYPE(tl, thdf5) \
  template<typename T> struct get_hdf5_data_type<T, typename boost::enable_if<boost::is_same<T,tl> >::type > \
    { static H5::DataType type() { return H5::PredType::thdf5; }; }

  #define HDF5_SAFE_TYPE(tl, othertl, thdf5) \
  template<typename T> struct get_hdf5_data_type<T, \
          typename boost::enable_if< \
               boost::integral_constant<bool, \
                      boost::is_same<T, tl>::value  \
                  && !boost::is_same<T,othertl>::value >  \
               >::type \
          >  \
    { static H5::DataType type() { return H5::PredType::thdf5; }; }
  

  HDF5_SAFE_TYPE(long, int                        , NATIVE_LONG);
  HDF5_SAFE_TYPE(unsigned long, unsigned int      , NATIVE_ULONG);
  HDF5_SAFE_TYPE(long long, long                  , NATIVE_LLONG);
  HDF5_SAFE_TYPE(unsigned long long, unsigned long, NATIVE_ULLONG);
  HDF5_TYPE(char              , NATIVE_CHAR);
  HDF5_TYPE(unsigned char     , NATIVE_UCHAR);
  HDF5_TYPE(int               , NATIVE_INT);
  HDF5_TYPE(unsigned int      , NATIVE_UINT);
  HDF5_TYPE(float             , NATIVE_FLOAT);
  HDF5_TYPE(double            , NATIVE_DOUBLE);
  
  #undef HDF5_TYPE
  #undef HDF5_SAFE_TYPE

  // Extent generator
  template<std::size_t r>
  struct hdf5_extent_gen {
    typedef typename boost::detail::multi_array::extent_gen<r> type;
    
    static inline type build(hsize_t *d)
    {
      return (hdf5_extent_gen<r-1>::build(d))[d[r-1]];
    }
  };

  template<>
  struct hdf5_extent_gen<0> {    
    static inline boost::multi_array_types::extent_gen build(hsize_t *)
    {
      return boost::extents;
    }
  };
  

//!_______________________________________________________________________________________
//!     
//!     write_hdf5 multi_array
//!         
//!     \author Guilhem Lavaux (2014-2015)
//!     \author leo Goodstadt (04 March 2013)
//!     
//!_______________________________________________________________________________________
  template<typename ArrayType,  typename hdf5_data_type>
  void hdf5_write_array(H5_CommonFileGroup& fg, const std::string& data_set_name, 
                        const ArrayType& data, 
                        const hdf5_data_type& datatype,
                        const std::vector<hsize_t>& dimensions,
                        bool doCreate = true,
                        bool useBases = false)
  {
      std::vector<hsize_t> memdims(data.shape(), data.shape() + data.num_dimensions());
      H5::DataSpace dataspace(int(dimensions.size()), dimensions.data());
      H5::DataSpace memspace(int(memdims.size()), memdims.data());

      if (useBases) {
        std::vector<hsize_t> offsets(data.index_bases(), data.index_bases() + data.num_dimensions());
        dataspace.selectHyperslab(H5S_SELECT_SET, memdims.data(), offsets.data());
      }

      H5::DataSet dataset;
      if (doCreate)
        dataset = fg.createDataSet(data_set_name, datatype, dataspace);
      else
        dataset = fg.openDataSet(data_set_name);

      dataset.write(data.data(), datatype, memspace, dataspace);
  }


  template<typename ArrayType,  typename hdf5_data_type>
  void hdf5_write_array(H5_CommonFileGroup& fg, const std::string& data_set_name, 
                        const ArrayType& data, 
                        const hdf5_data_type& datatype,
                        bool doCreate = true,
                        bool useBases = false)
  {
      std::vector<hsize_t> dimensions(data.shape(), data.shape() + data.num_dimensions());
      hdf5_write_array(fg, data_set_name, data, datatype, dimensions, doCreate, useBases);
  }

  /* HDF5 complex type */
  template<typename T>
  class hdf5_ComplexType
  {
  public:
    H5::CompType type;
    
    hdf5_ComplexType()
      : type(sizeof(std::complex<T>))
    {
      get_hdf5_data_type<T> hdf_data_type;
      type.insertMember("r", 0, hdf_data_type.type());
      type.insertMember("i", sizeof(T), hdf_data_type.type());
      type.pack();
    }
    
    static const hdf5_ComplexType<T> *ctype()
    {
      static hdf5_ComplexType<T> singleton;
      
      return &singleton;
    }
  };

  template<> struct get_hdf5_data_type<std::complex<float> > { 
    static H5::DataType type() { 
      return hdf5_ComplexType<float>::ctype()->type; 
    }
  };

  template<> struct get_hdf5_data_type<std::complex<double> > { 
    static H5::DataType type() { 
      return hdf5_ComplexType<double>::ctype()->type; 
    }
  };

  class hdf5_StringType
  {
  public:
    H5::StrType type;
    
    hdf5_StringType()
        : type(0, H5T_VARIABLE)
    {
    }
    
    static const hdf5_StringType *ctype()
    {
        static hdf5_StringType singleton;
        return &singleton;
    }
  };

  template<> struct get_hdf5_data_type<std::string> { 
    static H5::DataType type() {
      return hdf5_StringType::ctype()->type;
    }
  };

  class hdf5_BoolType
  {
  public:
    H5::EnumType type;

    hdf5_BoolType()
      : type(sizeof(bool))
    {
      bool v;
      
      v = true;
      type.insert("TRUE", &v);
      v = false;
      type.insert("FALSE", &v);
    }
    static const hdf5_BoolType *ctype()
    {
      static hdf5_BoolType singleton;
      
      return &singleton;
    }
  };

  template<> struct get_hdf5_data_type<bool> { 
    static H5::DataType type() {
      return hdf5_BoolType::ctype()->type;
    }
  };

  struct CosmoString { 
     char const * s;
     operator char const *() { return s; }
     char const *operator=(char const *s0) {
       s = s0;
       return s0;
     }
  };

  class hdf5_CosmoStringType {
  public:
    H5::StrType type;

    hdf5_CosmoStringType() : type(H5::PredType::C_S1, H5T_VARIABLE) {
    }

    static const hdf5_CosmoStringType *ctype()
    {
      static hdf5_CosmoStringType singleton;
      return &singleton;
    }
  };

  template<> struct get_hdf5_data_type<CosmoString> {
    static H5::DataType type() {
      return hdf5_CosmoStringType::ctype()->type;
    }
  };

  template<typename ArrayType>
  void hdf5_write_array(H5_CommonFileGroup& fg, const std::string& data_set_name, const ArrayType& data )
  {
      typedef typename ArrayType::element T;
      get_hdf5_data_type<T> hdf_data_type;
      
      hdf5_write_array(fg, data_set_name, data, hdf_data_type.type());
  }

  // HDF5 array reader
  //
  // Author Guilhem Lavaux (May 2014)
  
  class InvalidDimensions: virtual std::exception {
  };
  
  
  // ----------------------------------------------------------------------
  // Conditional resize support
  // If the Array type support resize then it is called. Otherwise
  // the dimensions are checked and lead to a failure if they are different
  
  template<typename Array> class array_has_resize {
    struct Fallback { int resize; };
    struct Derived: Array, Fallback {};

    typedef char yes[1];
    typedef char no[2];

    template<typename U, U> struct Check;

    template<typename U> 
    static yes& func(Check<int Fallback::*, &U::resize> *);

    template<typename U>
    static no& func(...);
  public:
    typedef array_has_resize type;
    enum { value = sizeof(func<Derived>(0)) == sizeof(no) };
  };

  

  template<typename ArrayType>
    typename boost::enable_if<
      array_has_resize<ArrayType> 
  >::type 
    hdf5_resize_array(ArrayType& data, std::vector<hsize_t>& dims) {
    data.resize(
        hdf5_extent_gen<ArrayType::dimensionality>::build(dims.data())
    );
  }

  template<typename ArrayType>
  void hdf5_check_array(ArrayType& data, std::vector<hsize_t>& dims) {
    for (size_t i = 0; i < data.num_dimensions(); i++) {
        if (data.shape()[i] != dims[i]) {
          throw InvalidDimensions();
        }
    }
  }

  template<typename ArrayType>
  void hdf5_weak_check_array(ArrayType& data, std::vector<hsize_t>& dims) {
    for (size_t i = 0; i < data.num_dimensions(); i++) {
        if (data.index_bases()[i] < 0) {
          // Negative indexes are not supported right now.
          throw InvalidDimensions();
        }
        if (data.index_bases()[i]+data.shape()[i] > dims[i]) {
          throw InvalidDimensions();
        }
    }
  }


  template<typename ArrayType>
    typename boost::disable_if<
      array_has_resize<ArrayType> 
  >::type 
    hdf5_resize_array(ArrayType& data, std::vector<hsize_t>& dims) {
    hdf5_check_array(data, dims);
  }

  // ----------------------------------------------------------------------
  
  template<typename ArrayType, typename hdf5_data_type>
  void hdf5_read_array_typed(H5_CommonFileGroup& fg, const std::string& data_set_name, 
                             ArrayType& data, 
                             const hdf5_data_type& datatype, bool auto_resize = true, bool useBases = false)
  {
      H5::DataSet dataset = fg.openDataSet(data_set_name);
      H5::DataSpace dataspace = dataset.getSpace();
      std::vector<hsize_t> dimensions(data.num_dimensions());
      
      if ((size_t)dataspace.getSimpleExtentNdims() != (size_t)data.num_dimensions())
        {
          throw InvalidDimensions();
        }
            
      dataspace.getSimpleExtentDims(dimensions.data());
      if (auto_resize)
          hdf5_resize_array(data, dimensions);
      else {
          if (useBases) {
            hdf5_weak_check_array(data, dimensions);

            std::vector<hsize_t> memdims(data.shape(), data.shape() + data.num_dimensions());
            H5::DataSpace memspace(int(memdims.size()), memdims.data());
            
            std::vector<hsize_t> offsets(data.index_bases(), data.index_bases() + data.num_dimensions());            
            dataspace.selectHyperslab(H5S_SELECT_SET, memdims.data(), offsets.data());
            
            dataset.read(data.data(), datatype, memspace, dataspace);
            return;
          } else {
            hdf5_check_array(data, dimensions);
          }
      }
      dataset.read(data.data(), datatype);
  }

  template<typename ArrayType>
  void hdf5_read_array(H5_CommonFileGroup& fg, const std::string& data_set_name, ArrayType& data, bool auto_resize = true, 
                      bool useBases = false )
  {
      typedef typename ArrayType::element T;
      
      hdf5_read_array_typed(fg, data_set_name, data, get_hdf5_data_type<T>::type(), auto_resize, useBases);
  }


#define CTOOL_HDF5_NAME(STRUCT) BOOST_PP_CAT(hdf5_,STRUCT)
#define CTOOL_HDF5_INSERT_ELEMENT(r, STRUCT, element) \
     { \
       ::CosmoTool::get_hdf5_data_type<BOOST_PP_TUPLE_ELEM(2, 0, element)> t; \
       long position = HOFFSET(STRUCT, BOOST_PP_TUPLE_ELEM(2, 1, element)); \
       const char *field_name = BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2, 1, element)); \
       type.insertMember(field_name,  position, t.type()); \
     } 
  
#define CTOOL_STRUCT_TYPE(STRUCT, TNAME, ATTRIBUTES) \
namespace CosmoTool { \
    class TNAME {  \
    public:  \
      H5::CompType type;  \
  \
      TNAME() : type(sizeof(STRUCT)) \
      {  \
        BOOST_PP_SEQ_FOR_EACH(CTOOL_HDF5_INSERT_ELEMENT, STRUCT, ATTRIBUTES) \
      }  \
      \
      static const TNAME *ctype() \
      { \
        static TNAME singleton; \
        return &singleton; \
      } \
    }; \
    template<> struct get_hdf5_data_type<STRUCT> { \
        static H5::DataType type() { return TNAME::ctype()->type; }; \
    }; \
};


#define CTOOL_HDF5_INSERT_ENUM_ELEMENT(r, STRUCT, element) \
     { \
       const char *field_name = BOOST_PP_STRINGIZE(element); \
       STRUCT a = element; \
       type.insert(field_name, &a); \
     } 


#define CTOOL_ENUM_TYPE(STRUCT, TNAME, ATTRIBUTES) \
namespace CosmoTool { \
    class TNAME {  \
    public:  \
      H5::EnumType type;  \
  \
      TNAME() : type(sizeof(STRUCT)) \
      {  \
        BOOST_PP_SEQ_FOR_EACH(CTOOL_HDF5_INSERT_ENUM_ELEMENT, STRUCT, ATTRIBUTES) \
      }  \
      \
      static const TNAME *ctype() \
      { \
        static TNAME singleton; \
        return &singleton; \
      } \
    }; \
    template<> struct get_hdf5_data_type<STRUCT> { \
        static H5::DataType type() { return TNAME::ctype()->type; }; \
    }; \
};

#define CTOOL_ARRAY_TYPE(ARRAY_TYPE, DIM, TNAME) \
namespace CosmoTool { \
    class TNAME { \
    public: \
      H5::ArrayType *type; \
\
      TNAME() \
      { \
        hsize_t dims[1] = { DIM }; \
        type = new H5::ArrayType(get_hdf5_data_type<ARRAY_TYPE>::type(), 1, dims); \
      } \
      ~TNAME() { delete type; } \
\
      static const TNAME *ctype() \
      { \
        static TNAME singleton; \
        return &singleton; \
      } \
    }; \
\
    template<> struct get_hdf5_data_type< ARRAY_TYPE[DIM] > { \
      static H5::DataType type() { return *(TNAME::ctype()->type); }; \
    }; \
};

}

#endif


