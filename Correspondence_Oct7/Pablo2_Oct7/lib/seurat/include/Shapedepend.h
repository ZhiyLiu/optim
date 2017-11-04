/********************************************************************************/
/*																				*/
/*  	File	:  Shapedepend.H												*/
/*																				*/
/*	Description:  class dependency declarations and header includes				*/
/*		for shapemonger library													*/
/*																				*/
/*	Project :  Shapemonger														*/
/*																				*/
/*	Author  :  A. Thall															*/
/*																				*/
/*	Date	:  15. November 1996												*/
/*																				*/
/*	Modifications:  5. September 1998 -- modified all code to use Fltk			*/
/*		interface rather than SMOOCH, and removed all SMOOCH and				*/
/*		Motif dependencies.														*/
/*																				*/
/*		17. May 99 -- replaced Quadfigs by Slicefigs							*/
/*		31. July 00 -- modified for use in Pablo								*/
/*				added DI_ ifdefs and defines									*/
/*				converted all from .H to .h										*/
/*		G Gash mods:  removed <stdio.h>, which was causing problems with		*/
/*			namespaces or some such.  Including separately as necessary.		*/
/*		6. Dec 00 -- added new classes to bring up to date with new				*/
/*			DiatomInterp class													*/
/*		8. Aug 01 -- added Zerofinder class support								*/
/*																				*/
/********************************************************************************/

/********************************************************************************/
/* Establishes #defines for required classes for inclusion by shapemonger.h		*/
/*   includes.																	*/
/********************************************************************************/

// NOTE:  Subdivsurf only needs spatialstate and Diatomgrids, not full Slicefig,
//   so we can avoid problems with circularly including Shaperaster.h
//   Does need to know about the OpenGL stuff, though, so needs FLTKINCLUDES

// THE BELOW IS TO ELIMINATE ERROR MESSAGES ASSOCIATED WITH VC++ TRYING TO
//   GENERATE DEBUG INFO FOR TEMPLATED OBJECTS WITH OVERLY LONG NAMES
#pragma warning(disable:4786)


#ifndef MyFloat
#define MyFloat double
#endif

// Include Xferlist.h but keep it out of the namespace 
//   since Tom's diatoms conflict with ours
#ifdef D_XFERLIST
#ifndef DI_XFERLIST
#define DI_XFERLIST
#include <iostream>
#include "Xferlist.h"
#endif
#endif


#ifdef D_INTERPOLMESHATOM
#define D_DIATOM
#define D_INTERPOLFIG
#define D_BEZIER1D
#endif

#ifdef D_PLIST_SUBDIVCOMP
#define D_POINTLIST_SERVER2
#endif

#ifdef D_POINTLIST_SERVERB
#define D_POINTLIST_SERVER2
#endif

#ifdef D_POINTLIST_SERVER2
#define D_SUBDIVSURF
#define D_PSEUDOSET
#define D_CCSUBDIVSURF
#define D_MESH
#endif

#ifdef D_MESH
#define D_CCSUBDIVSURF
#define D_DIATOMGRID
#define D_TUBEMESH
#define D_QUADMESH
#endif
#ifdef D_CCSUBDIVSURF
#define D_MYLIST
#define D_LINALG
#define D_SUBDIVSURF
#endif

#ifdef D_SUBDIVSURF
#define D_DIATOM
#endif

#ifdef D_INTERPOLFIG
#define D_DIATOMINTERP
#endif

#ifdef D_DIATOMEDGEINTERP
#define D_BEZIER1D
#define D_DIATOM
#endif

#ifdef D_DIATOMQUADINTERP
#define D_BEZIER2D
#define D_SURFMATH
#define D_DIATOM
#endif

#ifdef D_DIATOMINTERP
#define D_BEZIER2D
#define D_BICUBICINTERP
#define D_DIATOMGRID
#define D_SURFMATH
#endif

#ifdef D_BEZIER1D
#define D_LINALG
#endif

#ifdef D_BEZIER2D
#define D_LINALG
#endif

#ifdef D_MESHTESTER
#define D_DIATOMGRID
#endif

#ifdef D_DIATOMGRID
#define D_DIATOM
#endif

#ifdef D_DIATOM
#define D_QUAT
#define D_DIATOMGRID
#endif

#ifdef D_QUAT
#define D_LINALG
#endif

#ifdef D_CONJGRAD
#define D_LINALG
#endif

#ifdef D_ZEROFINDER
#define D_CONJGRAD
#define D_LINALG
#endif

// All dependencies have been checked, so go ahead and begin including files
//   in the correct order.

#ifdef D_SAMPLESTATS
#ifndef DI_SAMPLESTATS
#define DI_SAMPLESTATS
namespace ThallCode {
#include "Samplestats.h"
}
#endif
#endif

#ifdef D_LINALG
#ifndef D_MATHINCLUDED
#include <math.h>
#endif
#ifndef DI_LINALG
#define DI_LINALG
#ifdef WIN32
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#ifndef GLUquadric
#define GLUquadric GLUquadricObj    /* For old versions of OpenGL and Mesa */
#endif
namespace ThallCode {
#include "LinAlg.h"
}
#endif
#endif

#ifdef D_MYLIST
#ifndef DI_MYLIST
#define DI_MYLIST
namespace ThallCode {
#include "MyList.h"
}
#endif
#endif

#ifdef D_PSEUDOSET
#ifndef DI_PSEUDOSET
#define DI_PSEUDOSET
#include <iostream>
#include <limits.h>
namespace ThallCode {
#include "PseudoSet.h"
}
#endif
#endif

#ifdef D_CONJGRAD
#ifndef DI_CONJGRAD
#define DI_CONJGRAD
namespace ThallCode {
#include "Conjgrad2.h"
}
#endif
#endif

#ifdef D_ZEROFINDER
#ifndef DI_ZEROFINDER
#define DI_ZEROFINDER
namespace ThallCode {
#include "Zerofinder.h"
}
#endif
#endif

#ifndef DI_SHAPEHEADER
#define DI_SHAPEHEADER
namespace ThallCode {
#include "Shapeheader.h"
}
#endif

#ifdef D_QUAT
#ifndef DI_QUAT
#define DI_QUAT
//namespace ThallCode {
#include "Quat.h"
//}
#endif
#endif

#ifdef D_DIATOM
#ifndef DI_DIATOM
#define DI_DIATOM
namespace ThallCode {
#include "Diatom.h"
}
#endif
#endif

#ifdef D_DIATOMGRID
#ifndef DI_DIATOMGRID
#define DI_DIATOMGRID
namespace ThallCode {
#include "Diatomgrid.h"
}
#endif
#endif

#ifdef D_MESHTESTER
#ifndef DI_MESHTESTER
#define DI_MESHTESTER
namespace ThallCode {
#include "Meshtester.h"
}
#endif
#endif

#ifdef D_BEZIER1D
#ifndef DI_BEZIER1D
#define DI_BEZIER1D
namespace ThallCode {
#include "Bezier1D.h"
}
#endif
#endif

#ifdef D_BEZIER2D
#ifndef DI_BEZIER2D
#define DI_BEZIER2D
namespace ThallCode {
#include "Bezier2D.h"
}
#endif
#endif

#ifdef D_DIATOMEDGEINTERP
#ifndef DI_DIATOMEDGEINTERP
#define DI_DIATOMEDGEINTERP
namespace ThallCode {
#include "DiatomEdgeInterp.h"
}
#endif
#endif

#ifdef D_DIATOMQUADINTERP
#ifndef DI_DIATOMQUADINTERP
#define DI_DIATOMQUADINTERP
namespace ThallCode {
#include "DiatomQuadInterp.h"
}
#endif
#endif

#ifdef D_SUBDIVSURF
#ifndef DI_SUBDIVSURF
#define DI_SUBDIVSURF
#include <stdio.h>
namespace ThallCode {
#include "Subdivsurf.h"
}
#endif
#endif

#ifdef D_CCSUBDIVSURF
#ifndef DI_CCSUBDIVSURF
#define DI_CCSUBDIVSURF
#include <iostream>
#ifdef WIN32
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
namespace ThallCode {
#include "CCSubdivsurf.h"
}
#endif
#endif

#ifdef D_MESH
#ifndef DI_MESH
#define DI_MESH
namespace ThallCode {
#include "Mesh.h"
}
#endif
#endif

#ifdef D_TUBEMESH
#ifndef DI_TUBEMESH
#define DI_TUBEMESH
namespace ThallCode {
#include "TubeMesh.h"
}
#endif
#endif
#ifdef D_QUADMESH
#ifndef DI_QUADMESH
#define DI_QUADMESH
namespace ThallCode {
#include "QuadMesh.h"
}
#endif
#endif

#ifdef D_POINTLIST_SERVER2
#ifndef DI_POINTLIST_SERVER2
#define DI_POINTLIST_SERVER2
#include <queue>
#include <set>
class ThallCode::Diatomgrid;
class ThallCode::CCSubdivsurf;
class Bpoint;
class ThallCode::Tileneighbors;
class ThallCode::QuadMesh;
namespace ThallCode {
#include "Pointlist_server2.h"
}
#endif
#endif

#ifdef D_POINTLIST_SERVERB
#ifndef DI_POINTLIST_SERVERB
#define DI_POINTLIST_SERVERB
class M3DObject;
class M3DFigureTreeNode;
class ThallCode::CCSubdivsurf;
class Bpoint;
class Bpoint2;
class ThallCode::Tileneighbors;
class ThallCode::Pointlist_server2;
namespace ThallCode {
#include "Pointlist_serverB.h"
}
#endif
#endif

#ifdef D_PLIST_SUBDIVCOMP
#ifndef DI_PLIST_SUBDIVCOMP
#define DI_PLIST_SUBDIVCOMP
namespace ThallCode {
#include "Plist_subdivcomp.h"
}
#endif
#endif

