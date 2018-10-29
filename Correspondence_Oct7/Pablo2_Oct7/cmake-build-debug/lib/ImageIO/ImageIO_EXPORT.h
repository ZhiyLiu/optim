#ifndef __ImageIO_EXPORT_h_INCLUDED__
#define __ImageIO_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols ImageIO_EXPORT and ImageIO_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (ImageIO_BUILD_SHARED)
  #ifdef ImageIO_EXPORT_SYMBOLS
    #define ImageIO_EXPORT __declspec( dllexport )
  #else
    #define ImageIO_EXPORT __declspec( dllimport )
  #endif
  #define ImageIO_CDECL __cdecl
#else
  #define ImageIO_EXPORT
  #define ImageIO_CDECL
#endif // defined(_WIN32)

#endif

