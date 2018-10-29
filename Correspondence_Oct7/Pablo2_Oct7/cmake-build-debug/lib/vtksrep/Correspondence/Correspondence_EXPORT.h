#ifndef __Correspondence_EXPORT_h_INCLUDED__
#define __Correspondence_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols Correspondence_EXPORT and Correspondence_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (Correspondence_BUILD_SHARED)
  #ifdef Correspondence_EXPORT_SYMBOLS
    #define Correspondence_EXPORT __declspec( dllexport )
  #else
    #define Correspondence_EXPORT __declspec( dllimport )
  #endif
  #define Correspondence_CDECL __cdecl
#else
  #define Correspondence_EXPORT
  #define Correspondence_CDECL
#endif // defined(_WIN32)

#endif

