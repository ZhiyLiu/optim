#ifndef __SRepIO_EXPORT_h_INCLUDED__
#define __SRepIO_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols SRepIO_EXPORT and SRepIO_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (SRepIO_BUILD_SHARED)
  #ifdef SRepIO_EXPORT_SYMBOLS
    #define SRepIO_EXPORT __declspec( dllexport )
  #else
    #define SRepIO_EXPORT __declspec( dllimport )
  #endif
  #define SRepIO_CDECL __cdecl
#else
  #define SRepIO_EXPORT
  #define SRepIO_CDECL
#endif // defined(_WIN32)

#endif

