#ifndef _GP_TUNING_H
#define _GP_TUNING_H


#include "Tuning.h"


// The enumeration names are, except for TUNE_COUNT, the same as the
// long names in gpTuning.cpp.  However, names that are exclusive to
// Regular (Greyscale) Pablo also have the prefix "Gp".
enum gpTune_t
{
	ScriptVersion = 0,		// Pablo version, YYYYMMDD  (MUST BE FIRST IN ENUM)

	LineLength,				// avoid printing more than this many chars/line

	GpCmdHistTraining,		// Output a histogram file for the input model/image, exit
	GpCmdMultiPatchTraining,   // Multi patch training on the image and model input.  More flexible than hist.
	GpCmdDQFTraining,		// Testing the bone distance quantile function ideas.
	GpCmdDQFImageMaker,		// Making bone distance quantile function image for a grey image.
	GpCmdDQFModelMaker,		// Making bone distance quantile function for model and image.

	// Tube specific
	IgnoreMatchAtEnds,	// Ignore ends in the process of computing image match.

	// Optimizer
	BrentLinearSearchBoundFactor,	// Set the brent linear search bound factor. See gpTuning.cpp for details.
	OptimizerRandomSeed,	// Set random number seed
	GpOptimizerMatchMethod,	// Optimizer matching method to be used
	GpSurfaceLevel,			// Global surface level for optimization
	GpUsePGA,				// If PGA data is available, use it to optimize; default is true
	FigureStartingCoefficients,		// The starting coefficient list.
	GpHistogramMatchType,	// Type of histogram match to use - simple, Gaussian, difference
	WindowIntensities,		// Use windowed intensities for the optimization; default = true (1)

	// Model Stage
	GpSkipModelStage,		// Skip the model or ensemble stage
	GpModelPenaltyWeight,	// Model stage (registration) penalty weight

	// Figure Stage
	FigureSkipFigures,		// List of figures for which the figure stage will be skipped
/*	GP_TUNE_FRESET,			// Reset reference model before figure stage
	GP_TUNE_FLEV,			// subdivision surface LEVel: 0 is min
*/
	FigureIterations,		// Figure maximum number of ITERations
	FigureModelPerIteration,	// Printf-style format to save model (per iteration) to file
	FigureModelResult,		// Printf-style format to save last best model to file 
	GpFigureHistogramMatch,	// Use histogram match at the figure stage
	GpFigureMultiPatchMatch, //Use the MutliPatch match at the figure stage.
	FigureDQFMatch, 		// Use the bone Distance Quantile Function match at figure stage

	FigureMahalanobisMatch,	// Figure Mahalanobis Match
	FigureMaxPGAModes,		// Figure max # of PGA modes used
	FigureMahalanobisOnly,	// Omit similarity optimization from figure stage; do PGAs only
	GpFigureGeometryWeight,		// Geometry weight for figure stage, [0, 1]
	GpFigureInterObjectWeight,	// Inter-object (constraints) penalty weight for figure stage, [0, 1]
	GpFigureProfileMatch,	// Image match penalty; this parameter cannot be used in a script

	// FIGURE STAGE (multi figure) -- some are ATOM stage, too
	GP_TUNE_FIP,			// Mask figure stage penetration
	GP_TUNE_IPMASK,			// For calculating inter-penetration penalty (atom/fig stages)
	GP_TUNE_FIGURE_ID,		// Multi-object: selects the figure (atom/fig stages)


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
	FigureContourMatch,		// Weight on the penalty
	FigureContourLevel,		// Subdivision level
	FigureContourCloseness,	// See DistanceToPointSetFunction::setClosenessLimit()
	FigureContourFigureId,	// Which figure to match to contours
	GpFigureContourCoords,	// Are contour points in world coordinates?  0 (default) for model space, 1 for world space.

//BpFigureContourPadding, //number of axial slices to pad
//BpFigureContourSizeWeight, // weight on term that encourages m-rep to take up the same number of axial slices as the contours


	// ATOM STAGE
	AtomSkipFigures,		// List of figures for which the atom stage will be skipped
	GpAtomHistogramMatch,	// HIM (Histogram Image Match)
	GpAtomMultiPatchMatch,	// MPF (MultiPatch Image Match)
	AtomMahalanobisMatch,	// Atom Mahalanobis Match
	GpAtomGeometryWeight,	// Geometry weight for atom stage, [0, 1]
	GpAtomNeighborWeight,	// Neighbor penalty weight for atom stage, [0, 1]
	GpAtomInterWeight,		// Inter-atom (constraint) penalty weight for atom stage, [0, 1]
	GpAtomProfileMatch,		// Image match penalty; this parameter cannot be used in a script

    // ATOM STAGE: was (1-alpha)*D + alpha*R + beta*K

		
//	GP_TUNE_ARESET,			// reset reference model before atom stage

//	GP_TUNE_ALEV,			// Atom subdivision surface level: 0 is min

	// The smaller of the two following settings will win
	AtomIterations,			// Atom maximum number of cycles through all atoms in the figure
	GpAtomSubIterations,	// Atom maximum number of iterations (primarily useful for debugging)

/*	GP_TUNE_AAN,			// AN (Average of Neighbors: model regularity); was alpha*R
	GP_TUNE_AIM,			// IM (Image Match: average distance^2); was (1-alpha)*D
	GP_TUNE_ALM,			// LM (LandMark Match: average distance^2); was gamma*L
	GP_TUNE_AMM,			// MM (Model Match: to reference model); was beta*K
	GP_TUNE_ANEI,			// Atom neighborhood
	GP_TUNE_ADIST,			// Atom distance types

	// ATOM STAGE (multi figure)
	GP_TUNE_AIP,			// Mask atom stage penetration
*/
//	GP_TUNE_AMN,			// Model Normals
//	GP_TUNE_AMS,			// Model sliding
//	GP_TUNE_AMP,			// Model Penetration
//	GP_TUNE_AMC,			// Model Creasing

/*	GP_TUNE_AST,			// Atom skip threshold
	GP_TUNE_ASM,			// Atom skip multiplier
	GP_TUNE_ACG,			// Conjugate gradient iterations
	GP_TUNE_ASP,			// Skipped passes
	GP_TUNE_ARAND,			// Randomize atom order
*/
	// ATOM optimizer dimensional factors and bounds
	AtomTranslationFactor,	// Translation factor
	AtomRotationFactor,		// Rotation factor
	AtomScaleFactor,		// Scale factor
	AtomThetaFactor,		// Theta factor
//	AtomElongationFactor,	// Elongation factor (in atom stage???)
	AtomEtaFactor,			// Eta factor: 0 disables eta optimization 
	AtomPGAFactor,	        //Using atom PGA
	AtomEpsilon,			// Optimizer epsilons (step size)

	AtomRestrictMovement,	// Restricts the movement of the atom.

	AtomTranslationBound,	// Translation bound
	AtomRotationBound,		// Rotation bound
	AtomScaleBound,			// Scale bound
	AtomThetaBound,			// Theta bound
	AtomEtaBound,			// Eta bound: 0 disables eta optimization 
	AtomPGABound,	        // Bound for atom PGA

	// SPOKE STAGE:
	SpokeSkipFigures,		// List of figures for which the spoke stage will be skipped
	GpSpokeHistogramMatch,	// HIM (Histogram Image Match)
	SpokeMahalanobisMatch,	// Spoke Mahalanobis Match
	GpSpokeGeometryWeight,	// Geometry weight for spoke stage, [0, 1]
	GpSpokeNeighborWeight,	// Neighbor penalty weight for spoke stage, [0, 1]
	// FIXME: Check what this means
	GpSpokeInterWeight,		// Inter-spoke (constraint) penalty weight for spoke stage, [0, 1]
	GpSpokeProfileMatch,	// Image match penalty; this parameter cannot be used in a script

	SpokeIterations,		// Spoke maximum number of ITERations

	// SPOKE optimizer dimensional factors and bounds
	SpokePGAFactor,			// Spoke pga factor
	SpokeScaleFactor,		// Spoke scale factor
	SpokeEpsilon,			// Optimizer epsilons (step size)
	SpokePGABound,			// Spoke PGA bound
	SpokeScaleBound,		// Spoke scale bound

	// ------------
	RSradPenaltyNorm,		// rSrad penalty norm
	RSradPenaltyThreshold,	// rSrad penalty threshold

	// Tile generation
	TileSurfaceLevel,		// Tile subdivision level (min is 0)
	TileQuads,				// Generate quads (1) or triangles (0)
			// Use only one of the next two parameters
	TileFiguresMask,		// Bit mask indicating for which figures tiles will be output
	TileFiguresList,		// List of figures for which tiles will be output
	TileSeparateFiles,		// Whether or not to use multiple files to receive tile output
	TileUnit,				// Unit coordinate system (1) or world (0)

	// Input Files
	Image,					// Image file path
#ifdef AE2_BUILD
	LoadedImage,			// Image file path
	LoadedModel,			// Use loaded model
#endif
	Model,					// Model file path
	PGAModel,               // PGA model File Path
	Transform,				// Symmetry transform file path
	GpSimpleMaskFile,		// Name of the simple mask (.smsk) file
//	GpMaskFile,				// Name of the old-style mask (.msk) file			// AGG: obsolete
	GpReferenceModel,
	GpTrainingImage,
	GpProfileFile,
	GpHistogramFile,
	GpHistogramMeanFile,
//	LandmarkModel,
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


class gpTuning : public Tuning
{

public:

	gpTuning();
	~gpTuning() { }

	gpTune_t indexByLongName(const char * name) const {
		return (gpTune_t) Tuning::indexByLongName(name);
	}

protected:

	// Store initial values
	bool initialize();
};


// This variable is instantiated in P3DControl at the global scope.
extern class gpTuning tuneVals;


// An alias for tuneVals.getEntry(int index)
inline const TuneVal_t tuning(int index) {
	return tuneVals.getEntry(index);
}

// An alias for tuneVals.getEntry(GP_TUNE_ENUM parameterId)
inline const TuneVal_t tuning(gpTune_t parameterId) {
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

