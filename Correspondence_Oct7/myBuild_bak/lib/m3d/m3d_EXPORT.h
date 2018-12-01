#ifndef __m3d_EXPORT_h_INCLUDED__
#define __m3d_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols m3d_EXPORT and m3d_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (m3d_BUILD_SHARED)
  #ifdef m3d_EXPORT_SYMBOLS
    #define m3d_EXPORT __declspec( dllexport )
  #else
    #define m3d_EXPORT __declspec( dllimport )
  #endif
  #define m3d_CDECL __cdecl
#else
  #define m3d_EXPORT
  #define m3d_CDECL
#endif // defined(_WIN32)

#endif

