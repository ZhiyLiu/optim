#ifndef __paul_code_EXPORT_h_INCLUDED__
#define __paul_code_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols paul_code_EXPORT and paul_code_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (paul_code_BUILD_SHARED)
  #ifdef paul_code_EXPORT_SYMBOLS
    #define paul_code_EXPORT __declspec( dllexport )
  #else
    #define paul_code_EXPORT __declspec( dllimport )
  #endif
  #define paul_code_CDECL __cdecl
#else
  #define paul_code_EXPORT
  #define paul_code_CDECL
#endif // defined(_WIN32)

#endif

