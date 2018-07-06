#ifndef __zlib_EXPORT_h_INCLUDED__
#define __zlib_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols zlib_EXPORT and zlib_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (zlib_BUILD_SHARED)
  #ifdef zlib_EXPORT_SYMBOLS
    #define zlib_EXPORT __declspec( dllexport )
  #else
    #define zlib_EXPORT __declspec( dllimport )
  #endif
  #define zlib_CDECL __cdecl
#else
  #define zlib_EXPORT
  #define zlib_CDECL
#endif // defined(_WIN32)

#endif

