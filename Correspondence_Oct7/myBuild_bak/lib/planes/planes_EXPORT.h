#ifndef __planes_EXPORT_h_INCLUDED__
#define __planes_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols planes_EXPORT and planes_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (planes_BUILD_SHARED)
  #ifdef planes_EXPORT_SYMBOLS
    #define planes_EXPORT __declspec( dllexport )
  #else
    #define planes_EXPORT __declspec( dllimport )
  #endif
  #define planes_CDECL __cdecl
#else
  #define planes_EXPORT
  #define planes_CDECL
#endif // defined(_WIN32)

#endif

