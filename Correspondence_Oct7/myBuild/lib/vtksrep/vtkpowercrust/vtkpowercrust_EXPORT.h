#ifndef __vtkpowercrust_EXPORT_h_INCLUDED__
#define __vtkpowercrust_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols vtkpowercrust_EXPORT and vtkpowercrust_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (vtkpowercrust_BUILD_SHARED)
  #ifdef vtkpowercrust_EXPORT_SYMBOLS
    #define vtkpowercrust_EXPORT __declspec( dllexport )
  #else
    #define vtkpowercrust_EXPORT __declspec( dllimport )
  #endif
  #define vtkpowercrust_CDECL __cdecl
#else
  #define vtkpowercrust_EXPORT
  #define vtkpowercrust_CDECL
#endif // defined(_WIN32)

#endif

