#ifndef __SRepVisualization_EXPORT_h_INCLUDED__
#define __SRepVisualization_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols SRepVisualization_EXPORT and SRepVisualization_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (SRepVisualization_BUILD_SHARED)
  #ifdef SRepVisualization_EXPORT_SYMBOLS
    #define SRepVisualization_EXPORT __declspec( dllexport )
  #else
    #define SRepVisualization_EXPORT __declspec( dllimport )
  #endif
  #define SRepVisualization_CDECL __cdecl
#else
  #define SRepVisualization_EXPORT
  #define SRepVisualization_CDECL
#endif // defined(_WIN32)

#endif

