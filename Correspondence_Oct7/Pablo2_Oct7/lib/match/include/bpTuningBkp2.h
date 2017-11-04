#ifndef _BP_TUNING_H
#define _BP_TUNING_H


#include "Tuning.h"


// The enumeration names are, except for TUNE_COUNT, the same as the
// long names in bpTuning.cpp.  However, names that are exclusive to
// Binary Pablo also have the prefix "Bp".
enum bpTune_t
{
	ScriptVersion = 0,	// Pablo version, YYYYMMDD  (MUST BE FIRST IN ENUM)

	LineLength,				// avoid printing more than this many chars/line
	// PROGRAMMER AIDS
	BpDebug1,				// Special debugging flag
	BpDebug2,				// Special debugging flag
	SaveDistanceVectors,	// Whether to save the distance vectors into a file

	// NORMALIZING COORDINATE SYSTEMS
	BpCoordMaxExtent,		// Maximum extent of any dim in any image; 0: calc for the current image
	BpCoordCubeImage,		// Cube the image: enable image padding during -pi
	BpVoxelSpacing,			// New isotropic voxel spacing (world)

	// SPEED OPTIMIZATIONS
	BpIgnoreDistFile,		// Overwrite DDM file with new calc's
	BpDDMCenter,			// Do NOT dilate distance map to match normal-search distances
	BpDDMAngleTolerance,	// Angle threshold that enables normal search if DDMENABLE==2
	BpDDMSearchThreshold,	// "normal search" ceases when boundary is within this voxel distance
	BpDDMEnable,			// Chooses distance map, "normal search", angle-triggered "normal search", "disabled"

	// INITIALIZATION
	BpDoMethodOfMoments,	// Method of moments: bitflags 1=trans 2=scale 4=rotate
	WindowIntensities,		// Use windowed intensities for the optimization; default = true (1)
	BrentLinearSearchBoundFactor,	// Set the brent linear search bound factor. See bpTuning.cpp for details.
	OptimizerRandomSeed,	// Set random number seed

	// THIN OBJECTS
	BpDilationFactor,		// Fit a model to a dilated distance map and then shrink the object.

	// SUBDIVIDE MEDIAL SHEET
	BpSubdivide,			// Subdivide the medial sheet

	WeightMatchAtEnds,		// weight ends in the process of computing image match.

	// FIGURE STAGE (single figure only)
	FigureStartingCoefficients,	// The starting coefficient list.
	FigureSkipFigures,		// List of figures for which the figure stage will be skipped
	BpFigureResetModelBefore,	// reset reference model before figure stage

	BpFigureSurfaceLevel,	// Subdivision surface level: 0 is min			// AGG: This should also apply to Gp
	FigureIterations,		// Figure maximum number of ITERations
	FigureModelPerIteration,	// Printf-style format to save model (per iteration) to file
	FigureModelResult,		// Printf-style format to save last best model to file 

	BpFigureImageMatch,		// Figure Image Match
	BpFigureLandmarkMatch,	// Figure Landmark Match
	BpFigureModelMatch,		// Figure Model Match
	FigureMahalanobisMatch,	// Figure Mahalanobis Match
	FigureMaxPGAModes,		// Figure max # of PGA modes used
	FigureMahalanobisOnly,	// Omit similarity optimization from figure stage; do PGAs only
	BpFigureDoScale,		// Optimize figure's scale
	BpFigureNeighborhood,	// Figure atom neighborhood
	BpFigureDistanceType,	// Figure distance types
    BpFigureRSradPenalty,		// RSrad penalty to enforce legality
	BpFigureCurvinessPenalty,	// Figure stage curviness penalty.
    //Xiaojie
	BpFigureImageNormalMatch,  //Figure Image Normal Match
	BpFigureImagePlaneOrienAtEnd,  //Figure Image Plane Orientation Match at End Points
	BpFigureImageVertexMatch, //Figure Image Vertex 

	// FIGURE STAGE (multi figure) -- some are ATOM stage, too
	BpFigureInterpenetrationMatch,	// Mask figure stage penetration				// AGG: This probably should only be in Gp
	RepulsionMask,			// Bitplanes for other objects, for calculating inter-penetration penalty (atom/fig stages)
	AttractionMask,         // Bitplanes for this object
	BpFigureId,				// Multi-object: selects the figure (atom/fig stages)

	// FIGURE optimizer dimensional factors and bounds
	FigureTranslationFactor,// Translation factor
	FigureRotationFactor,	// Rotation factor
	FigureScaleFactor,		// Scale factor
	FigureElongationFactor,	// Figure elongation factor
	FigureTubePhiFactor,	// Phi alignment factor
	FigurePGAFactor,		// PGA factor
	FigureEpsilon,			// Optimizer epsilons (step size)
	FigureTranslationBound,	// Translation bound
	FigureRotationBound,	// Rotation bound
	FigureScaleBound,		// Scale bound
	FigureElongationBound,	// Figure elongation bound
	FigureTubePhiBound,		// Phi alignment bound
	FigurePGABound,			// PGA bound

	// FIGURE contour based initialization -- see also Contour
	FigureContourMatch,	// Weight on the penalty
	FigureContourLevel,	// Subdivision level
	FigureContourCloseness,	// See DistanceToPointSetFunction::setClosenessLimit()
	FigureContourFigureId,	// Which figure to match to contours
	BpFigureContourPadding,		// Number of axial slices to pad
	BpFigureContourSizeWeight,	// Weight on term that encourages m-rep to take up the same number of axial slices as the contours

	// ATOM STAGE: was (1-alpha)*D + alpha*R + beta*K
	AtomSkipFigures,		// List of figures for which the atom stage will be skipped					Not in use yet
	BpAtomResetModelBefore,	// reset reference model before atom stage

	BpAtomSurfaceLevel,		// Atom subdivision surface level: 0 is min					// AGG: This should also apply to Gp
	AtomIterations,			// Atom maximum number of iterations

	BpAtomAverageNeighborMatch,	// AN (Average of Neighbors: model regularity); was alpha*R
	BpAtomImageMatch,		// IM (Image Match: average distance^2); was (1-alpha)*D
	AtomMahalanobisMatch,	// Atom Mahalanobis Match
	BpAtomLandmarkMatch,	// LM (LandMark Match: average distance^2); was gamma*L
	BpAtomModelMatch,		// MM (Model Match: to reference model); was beta*K
	BpAtomNeighborhood,		// Atom neighborhood
	BpAtomDistanceType,		// Atom distance types
	BpAtomRSradPenalty,		// RSrad penalty to enforce legality
	BpAtomSumMatch,			// Sum over all atoms.
	BpAtomDeltaNeighborMatch,	// AN (difference in Average of Neighbors delta in reference versus candidate); was alpha*R
	BpAtomCurvinessPenalty,	// Atom stage curviness penalty
	//Xiaojie
	BpAtomImageNormalMatch,  //Atom Image Normal Match
	BpAtomImagePlaneOrienAtEnd,  //Atom Image Plane Orientation Match at End Points
	BpAtomImageVertexMatch, //Figure Image Vertex 

	// ATOM STAGE (multi figure)
	BpAtomInterpenetrationMatch,	// Mask atom stage penetration

//	BP_TUNE_AMN,			// Model Normals
//	BP_TUNE_AMS,			// Model sliding
//	BP_TUNE_AMP,			// Model Penetration
//	BP_TUNE_AMC,			// Model Creasing

																// AGG: The set of items below probably can be added to Gp easily
	BpAtomSkipThreshold,	// Atom skip threshold
	BpAtomSkipMultiplier,	// Atom skip multiplier
	BpAtomConjugateGradientIterations,		// Conjugate gradient iterations
	BpAtomSkippedPasses,	// Skipped passes
	BpAtomRandomOrder,		// Randomize the order of optimization of atoms

	// ATOM optimizer dimensional factors and bounds
	AtomTranslationFactor,	// Translation factor
	AtomRotationFactor,		// Rotation factor
	AtomScaleFactor,		// Scale factor
	AtomThetaFactor,		// Theta factor
//	AtomElongationFactor,	// Elongation factor (in atom stage???)
	AtomEtaFactor,			// Eta factor: 0 disables eta optimization 
	AtomPGAFactor,			// Using atom PGA
	AtomEpsilon,			// Optimizer epsilons (step size)

	AtomRestrictMovement,	// Restricts the movement of the atom.

	AtomTranslationBound,	// Translation bound
	AtomRotationBound,		// Rotation bound
	AtomScaleBound,			// Scale bound
	AtomThetaBound,			// Theta bound
	AtomEtaBound,			// Eta bound: 0 disables eta optimization 
	AtomPGABound,			// Atom PGA Bound

	//----------------------

	// SPOKE STAGE:
	SpokeSkipFigures,		// List of figures for which the atom stage will be skipped					Not in use yet
	BpSpokeResetModelBefore,	// reset reference model before atom stage

	BpSpokeSurfaceLevel,		// Spoke subdivision surface level: 0 is min					// AGG: This should also apply to Gp
	SpokeIterations,			// Spoke maximum number of iterations

	BpSpokeAverageNeighborMatch,		// AN (Average of Neighbors: model regularity); was alpha*R
	BpSpokeImageMatch,		// IM (Image Match: average distance^2); was (1-alpha)*D
	SpokeMahalanobisMatch,	// Spoke Mahalanobis Match
	BpSpokeLandmarkMatch,	// LM (LandMark Match: average distance^2); was gamma*L
	BpSpokeModelMatch,		// MM (Model Match: to reference model); was beta*K
	BpSpokeNeighborhood,		// Spoke neighborhood
	BpSpokeDistanceType,		// Spoke distance types
	BpSpokeRSradPenalty,		// RSrad penalty to enforce legality
	BpSpokeSumMatch,			// Sum over all atoms.

	// SPOKE STAGE (multi figure)
	BpSpokeInterpenetrationMatch,	// mask atom stage penetration

	BpSpokeSkipThreshold,	// Spoke skip threshold
	BpSpokeSkipMultiplier,	// Spoke skip multiplier
	BpSpokeConjugateGradientIterations,		// Conjugate gradient iterations
	BpSpokeSkippedPasses,	// Skipped passes
	BpSpokeRandomOrder,		// Randomize the order of optimization of atoms

	// SPOKE optimizer dimensional factors and bounds
	SpokePGAFactor,			// Spoke pga factor
	SpokeScaleFactor,		// Spoke scale factor
	SpokeEpsilon,			// Optimizer epsilons (step size)
	SpokePGABound,			// Spoke PGA bound
	SpokeScaleBound,		// Spoke scale bound

	//----------------------
	RSradPenaltyNorm,		// rSrad penalty norm
	RSradPenaltyThreshold,	// rSrad penalty threshold

	// Tile generation
	TileSurfaceLevel,		// Tile subdivision level (min is 0)
	TileQuads,				// Generate quads (1) or triangles (0)
	TileUnit,				// Unit coordinate system (1) or world (0)

	// Input Files
	Image,					// Image file path
	Model,					// Model file path
	Transform,				// Symmetry transform file path
	BpLandmarkModel,		// Landmark model file path
	OutModel,				// Output model file path
	OutModelImage,			// Output file path for  model image
	OutTile,				// Output tiles file path
	OutImageTile,			// Output file path for image tiles
	Contour,				// Input file containing (x,y,z) of contour vertices in model coordinates

	// Details about export-model-as-image
	OutModelImageX,			// Image's x dimension
	OutModelImageY,			// Image's y dimension
	OutModelImageZ,			// Image's z dimension
	OutModelImageSurfaceLevel,	// Subdivision level of model
	OutModelImageFixedSize,	// Prohibit map size from changing dim's to suit the model

	// Not a tunable parameter; this must always be last in the enum
	TUNE_COUNT				// The length of this enumeration
};


class bpTuning : public Tuning
{

public:

	bpTuning();
	~bpTuning() { }

	bpTune_t indexByLongName(const char * name) const {
		return (bpTune_t) Tuning::indexByLongName(name);
	}

protected:

	// Store initial values
	bool initialize();
};


// This variable is instantiated in P3DControl at the global scope.
extern class bpTuning tuneVals;


// An alias for tuneVals.getEntry(int index)
inline const TuneVal_t tuning(int index) {
	return tuneVals.getEntry(index);
}

// An alias for tuneVals.getEntry(bpTune_t parameterId)
inline const TuneVal_t tuning(bpTune_t parameterId) {
	return tuneVals.getEntry((int) parameterId);
}

// Returns the lock status of the specified parameter
inline bool tuningLocked(int index) {
	return tuneVals.locked(index);
}

// Returns the long name of the specified parameter
inline const char * tuningLongName(int index) {
	return tuneVals.longName(index);
}

// Returns the short name of the specified parameter
inline const char * tuningShortName(int index) {
	return tuneVals.shortName(index);
}

// Returns the weight of the specified parameter.
// The tuning entry is assumed to not be a list.
inline double tuningWt(int index) {
	return tuneVals.weight(index);
}

// Returns a pointer to the list for the specified parameter.
// The tuning entry is assumed to be a list.
inline const int * tuningList(int index, int & listLen) {
	return tuneVals.list(index, listLen);
}

// Returns a pointer to the character string of the specified
// parameter.  The tuning entry is assumed to be a string.
inline const char * tuningStr(int index) {
	return tuneVals.string(index);
}



#endif

