/* vim: set syn=cpp :
 *
 * Swig wrapper for the m3d library
 * - Rohit Saboo
 */

%{
#include "ControlKeys.h"
#include "control_parms_defaults.h"

// m3d includes.
#include "Hermite.h"
#include "Bezier2D.h"
#include "BYU.h"
#include "ControlParmsAccess.h"
#include "ControlParms.h"
#include "Image3D.h"
#include "GaussianBlur3D.h"
#include "ImageResample3D.h"
#include "InterfiguralConstraints.h"
#include "M3DAtomPredictor.h"
#include "M3DAtomPredictorQuad.h"
#include "M3DAtomPredictorTube.h"
#include "M3DFigure.h"
#include "M3DFigureStats.h"
#include "M3DFigureTreeNode.h"
#include "M3DObjectFile.h"
#include "M3DObject.h"
#include "M3DPGA.h"
#include "M3DPGAPrimitiveStats.h"
#include "M3DPGAStats.h"
#include "M3DPrimitive.h"
#include "M3DEndPrimitive.h"
#include "M3DQuadEndPrimitive.h"
#include "M3DQuadFigure.h"
#include "M3DQuadInterpolator.h"
//#include "M3DQuadPrimitive.h"
#include "M3DTubeEndPrimitive.h"
#include "M3DTubeFigure.h"
//#include "M3DTubePrimitive.h"
#include "Mathdefs.h"
#include "Quat.h"
#include "RAWImageFile.h"
#include "SimilarityTransform3D.h"
#include "SubdivBoundary.h"
#include "TileSet.h"
#include "utility.h"
#include "Vector2D.h"
#include "Vector3D.h"
#include "VectorND.h"
#include "WorldSystem.h"
%}

%include swig-helpers.i

TYPECAST(M3DTubeFigure, M3DFigure)
TYPECAST(M3DQuadFigure, M3DFigure)
TYPECAST(M3DQuadPrimitive, M3DPrimitive)
TYPECAST(M3DTubePrimitive, M3DPrimitive)


// For want of a better place, Control Keys is kept in the m3d library instead of someplace else.
%include ControlKeys.h

%include Hermite.h
%include Bezier2D.h
%include BYU.h
%include ControlParmsAccess.h
%include ControlParms.h
%include Image3D.h
%include GaussianBlur3D.h
%include ImageResample3D.h
%include InterfiguralConstraints.h
%include M3DAtomPredictor.h
%include M3DAtomPredictorQuad.h
%include M3DAtomPredictorTube.h
%include M3DFigure.h
%include M3DFigureStats.h
%include M3DFigureTreeNode.h
%include M3DObjectFile.h
%include M3DObject.h
%include M3DPGA.h
%include M3DPGAPrimitiveStats.h
%include M3DPGAStats.h
%include M3DPrimitive.h
%include M3DEndPrimitive.h
%include M3DQuadEndPrimitive.h
%include M3DQuadFigure.h
%include M3DQuadInterpolator.h
//%include M3DQuadPrimitive.h
%include M3DTubeEndPrimitive.h
%include M3DTubeFigure.h
//%include M3DTubePrimitive.h
%include Mathdefs.h
%include Quat.h
%include RAWImageFile.h
%include SimilarityTransform3D.h
%include SubdivBoundary.h
%include TileSet.h
%include utility.h
%include Vector2D.h
%include Vector3D.h
%include VectorND.h
%include WorldSystem.h

