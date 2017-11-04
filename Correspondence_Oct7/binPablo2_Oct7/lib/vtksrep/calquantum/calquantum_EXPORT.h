#ifndef __calquantum_EXPORT_h_INCLUDED__
#define __calquantum_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols calquantum_EXPORT and calquantum_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (calquantum_BUILD_SHARED)
  #ifdef calquantum_EXPORT_SYMBOLS
    #define calquantum_EXPORT __declspec( dllexport )
  #else
    #define calquantum_EXPORT __declspec( dllimport )
  #endif
  #define calquantum_CDECL __cdecl
#else
  #define calquantum_EXPORT
  #define calquantum_CDECL
#endif // defined(_WIN32)

#endif

