#ifndef __COSMOTOOL_SYMBOL_VISIBLE_HPP
#define __COSMOTOOL_SYMBOL_VISIBLE_HPP


#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_DLL
    #ifdef __GNUC__
      #define CTOOL_DLL_PUBLIC __attribute__ ((dllexport))
    #else
      #define CTOOL_DLL_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #else
    #ifdef __GNUC__
      #define CTOOL_DLL_PUBLIC __attribute__ ((dllimport))
    #else
      #define CTOOL_DLL_PUBLIC __declspec(dllimport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #endif
  #define CTOOL_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define CTOOL_DLL_PUBLIC __attribute__ ((visibility ("default")))
    #define CTOOL_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define CTOOL_DLL_PUBLIC
    #define CTOOL_DLL_LOCAL
  #endif
#endif



#endif
