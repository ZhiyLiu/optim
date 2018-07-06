#ifndef __flvw-1.0_EXPORT_h_INCLUDED__
#define __flvw-1.0_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols flvw-1.0_EXPORT and flvw-1.0_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (flvw-1.0_BUILD_SHARED)
  #ifdef flvw-1.0_EXPORT_SYMBOLS
    #define flvw-1.0_EXPORT __declspec( dllexport )
  #else
    #define flvw-1.0_EXPORT __declspec( dllimport )
  #endif
  #define flvw-1.0_CDECL __cdecl
#else
  #define flvw-1.0_EXPORT
  #define flvw-1.0_CDECL
#endif // defined(_WIN32)

#endif

