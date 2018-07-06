#ifndef __match_EXPORT_h_INCLUDED__
#define __match_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols match_EXPORT and match_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (match_BUILD_SHARED)
  #ifdef match_EXPORT_SYMBOLS
    #define match_EXPORT __declspec( dllexport )
  #else
    #define match_EXPORT __declspec( dllimport )
  #endif
  #define match_CDECL __cdecl
#else
  #define match_EXPORT
  #define match_CDECL
#endif // defined(_WIN32)

#endif

