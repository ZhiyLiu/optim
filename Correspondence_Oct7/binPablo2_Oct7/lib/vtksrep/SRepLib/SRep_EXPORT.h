#ifndef __SRep_EXPORT_h_INCLUDED__
#define __SRep_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols SRep_EXPORT and SRep_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (SRep_BUILD_SHARED)
  #ifdef SRep_EXPORT_SYMBOLS
    #define SRep_EXPORT __declspec( dllexport )
  #else
    #define SRep_EXPORT __declspec( dllimport )
  #endif
  #define SRep_CDECL __cdecl
#else
  #define SRep_EXPORT
  #define SRep_CDECL
#endif // defined(_WIN32)

#endif

