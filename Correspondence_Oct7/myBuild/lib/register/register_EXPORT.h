#ifndef __register_EXPORT_h_INCLUDED__
#define __register_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols register_EXPORT and register_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (register_BUILD_SHARED)
  #ifdef register_EXPORT_SYMBOLS
    #define register_EXPORT __declspec( dllexport )
  #else
    #define register_EXPORT __declspec( dllimport )
  #endif
  #define register_CDECL __cdecl
#else
  #define register_EXPORT
  #define register_CDECL
#endif // defined(_WIN32)

#endif

