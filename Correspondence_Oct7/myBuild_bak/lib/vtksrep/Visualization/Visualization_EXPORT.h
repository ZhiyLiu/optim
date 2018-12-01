#ifndef __Visualization_EXPORT_h_INCLUDED__
#define __Visualization_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols Visualization_EXPORT and Visualization_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (Visualization_BUILD_SHARED)
  #ifdef Visualization_EXPORT_SYMBOLS
    #define Visualization_EXPORT __declspec( dllexport )
  #else
    #define Visualization_EXPORT __declspec( dllimport )
  #endif
  #define Visualization_CDECL __cdecl
#else
  #define Visualization_EXPORT
  #define Visualization_CDECL
#endif // defined(_WIN32)

#endif

