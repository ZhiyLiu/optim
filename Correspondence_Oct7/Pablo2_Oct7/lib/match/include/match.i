/* vim: set syn=cpp :
 *
 * Swig wrapper for the match library
 * - Rohit Saboo
 */
%{
#include "ControlKeys.h"
#include "control_parms_defaults.h"
#include "ControlParms.h"

// match includes
#include "DistanceMap3D.h"
#include "ImageDistanceMap.h"
#include "Mask.h"
#include "Match.h"
#include "SimpleMaskFile.h"
#include "SurfacePatchEnsemble.h"
#include "TemplateProfiles.h"
#include "Tuning.h"
#include "DistanceVectorList.h"

#ifdef BINARY
#include "MatchUtility.h"
#else
#include "ObjectRelativeSampling.h"
#endif
%}

%include swig-helpers.i

%include DistanceMap3D.h
%include ImageDistanceMap.h
%include Mask.h
%include Match.h
%include SimpleMaskFile.h
%include SurfacePatchEnsemble.h
%include TemplateProfiles.h
%include Tuning.h
%include DistanceVectorList.h

#ifdef BINARY
%include MatchUtility.h
#else
%include ObjectRelativeSampling.h
#endif
