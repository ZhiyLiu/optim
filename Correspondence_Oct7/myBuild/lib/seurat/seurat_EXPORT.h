#ifndef __seurat_EXPORT_h_INCLUDED__
#define __seurat_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols seurat_EXPORT and seurat_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (seurat_BUILD_SHARED)
  #ifdef seurat_EXPORT_SYMBOLS
    #define seurat_EXPORT __declspec( dllexport )
  #else
    #define seurat_EXPORT __declspec( dllimport )
  #endif
  #define seurat_CDECL __cdecl
#else
  #define seurat_EXPORT
  #define seurat_CDECL
#endif // defined(_WIN32)

#endif

