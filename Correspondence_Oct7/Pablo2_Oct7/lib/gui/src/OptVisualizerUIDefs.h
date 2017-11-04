#ifndef _PRE_DEFS_H_
#define _PRE_DEFS_H_


#include "globalBuildParms.h"
#include <stdio.h>

#ifdef OPTIMIZATION_VISUALIZER
#include "LogManagerChart.h"
#include "Plot3DWindowWrapper.h"
#endif	/* OPTIMIZATION_VISUALIZER */

class P3DControl;
class LogManagerChart;
class Plot3DWindowWrapper;
#include <FL/Fl_Window.H>
#include "movable_Fl_Window.h"
#include "OptVisualizerCallback.h"

#ifdef OPTIMIZATION_VISUALIZER

class P3DUserInterface;
class P3DUserInterfaceCallback;

using namespace std;
using namespace pauly;
#endif	/* OPTIMIZATION_VISUALIZER */




#endif

