#ifndef __version_EXPORT_h_INCLUDED__
#define __version_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols version_EXPORT and version_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (version_BUILD_SHARED)
  #ifdef version_EXPORT_SYMBOLS
    #define version_EXPORT __declspec( dllexport )
  #else
    #define version_EXPORT __declspec( dllimport )
  #endif
  #define version_CDECL __cdecl
#else
  #define version_EXPORT
  #define version_CDECL
#endif // defined(_WIN32)

#endif

