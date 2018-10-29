#ifndef __gui_EXPORT_h_INCLUDED__
#define __gui_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols gui_EXPORT and gui_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (gui_BUILD_SHARED)
  #ifdef gui_EXPORT_SYMBOLS
    #define gui_EXPORT __declspec( dllexport )
  #else
    #define gui_EXPORT __declspec( dllimport )
  #endif
  #define gui_CDECL __cdecl
#else
  #define gui_EXPORT
  #define gui_CDECL
#endif // defined(_WIN32)

#endif

