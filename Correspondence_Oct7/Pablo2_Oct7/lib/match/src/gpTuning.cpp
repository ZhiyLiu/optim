#ifndef BINARY

#include <stdio.h>
#include <string.h>
#include <iostream>
#include "gpTuning.h"

using namespace std;



gpTuning::gpTuning() : Tuning()
{
	v = new TuneVal_t[TUNE_COUNT];
	numEntries = TUNE_COUNT;
	initialize();
}

bool gpTuning::initialize()
{
	// These must be in the same order as in gpTune_t.
	// The longName (second) arguments must begin with a capital letter,
	// because the scriptParms registry is designated as case-independent
	// in pablo.cpp.

	setNext("VER", "ScriptVersion",                     20050117.09,
			"[YYYYMMDD.HH] version of the control script. Value must be <= version contained in program. Identifies feature set, not CVS checkout date",
			true);	// Locked

	setNext("len",	"lineLength", 80, "max number of chars to print, per line; keep under 81 to fit to terminal width");

	setNext("HISTTRAIN", "CmdHistTraining",             0,   "[0,1] Compute training Histogram files, exit, 0=disable, 1=enable");

	setNext("MPTRAIN", "CmdMultiPatchTraining",             0,   "[0,1] Compute training Multi patch files, exit, 0=disable, 1=enable");

	setNext("DQFTRAIN", "CmdDQFTraining",             0,   "[0,1] Debug distance quantile functions for bones, exit, 0=disable, 1=enable");
	setNext("DQFIMAGE", "CmdDQFImageMaker",             0,   "[0,1] Make a distance quantile function image for bones, exit, 0=disable, 1=enable");
	setNext("DQFMODEL", "CmdDQFModelMaker",             0,   "[0,1] Make a distance quantile function file for the model, exit, 0=disable, 1=enable");

	// Tube Specific
	setNext("TUBEIGNOREENDS", "IgnoreMatchAtEnds",		0.0,
			"[0,1] 0 - Use ends while computing image match, 1 - ignore ends while computing image match.");

	// Optimizer
	setNext("BLSBF", "BrentLinearSearchBoundFactor",    2.0,
			"(1.0..) All of the bounds bound the particular component of the vector being optimized and are symmetric about the starting point. Also, they are used by the bracketing method in the brent linear search. So if x was the starting point of the optimization and b is the bound, then an object at x + b/2.0 is constructed and evaluated. The factor 2.0 is there to get one of the search points right in the middle of the starting position and the bound. Also, the factor 2.0 can be changed by this tuning parameter BrentLinearSearchBoundFactor. The usage of bounds cannot be completely disabled, but it's effect can be reduced simply by increasing all the bounds by a large number M and dividing this BrentLinearSearchBoundFactor by the same large number M.");
	setNext("SEED", "OptimizerRandomSeed",              12345,
			"[0..65535] seed for srand(); change this and get slightly different results");
	setNext("METHOD", "OptimizerMatchMethod",           0,	// The following values are hard-coded in P3DControl::runPablo() & P3DUserInterfaceCallback::startInteractive()
			"Optimizer matching method to be used: 0 = Gaussian, 1 = Negative Gaussian, 6 = Simple Mask, 7 = Binary Image, 8 = Histogram, 9 = MultiPatch");
	setNext("LEVEL", "SurfaceLevel",					 0,
			"0 - 3 corresponding to the surface level to be used during optimization; the default is 0");
	setNext("PGA", "UsePGA",						     1,
			"[0,1] Use PGA data, if available, during optimization; default is true");
	setNextString("FIGSTART", "FigureStartingCoefficients", "Space delimited list of starting coefficients to use for the figure stage from where the optimization is to be started.");
	setNext("HMTYPE",  "HistogramMatchType",			 0,		"[0,1,2] initialization stage, Histogram object will be initialied for match type: 0=simple 1=gaussian 2=difference");
	setNext("WININT", "WindowIntensities",				1, "[0, 1] Use windowed intensities for the optimization; default is true");

	// Model Stage
	setNext("SKIPMS", "SkipModelStage",					0, "[0, 1] Skip the model or ensemble stage");
	setNext("MSPENWT", "ModelPenaltyWeight",            0.5, "[0, 1] Model stage: penalty weight");

	// Figure stage
	setNextList("FSKPFIGS", "FigureSkipFigures",		"List of figure numbers for which the figure stage will be skipped");
/*	setNext("FRESET", "FigureResetModelBefore",         1.0,
			"[0,1] figure stage: set the reference model to the model as the stage is entered. The model match is calculated relative to the reference model");
	setNext("FLEV",   "FigureSurfaceLevel",             2,
			"[0,1,2,3] figure stage: 0:no subdivision, 1:each tile is split into 4 tiles. 2:16 tiles, 3:64 tiles");
*/
	setNext("FITER",  "FigureIterations",               40,
			"[0..] figure stage: max number of passes; 0 disables");
	setNextString("OMFPI", "FigureModelPerIteration",
		"Save a model after each figure stage iteration; it's value is a printf-stlye format string with one integer (%.03d) representing the iteration number;");
	setNextString("OMFR", "FigureModelResult",
		"if set, save last best model; it's value is a printf-stlye format string");
	setNext("FHIM",  "FigureHistogramMatch",			0.0,
			"[0..] figure stage: weight of histogram image match.  Use instead of FIM for grey image segmentation");
	setNext("FMPIM",  "FigureMultiPatchMatch",			0.0,
			"[0..] figure stage: weight of histogram image match.  Use instead of FIM for grey image segmentation");
	setNext("FDQFM",  "FigureDQFMatch",			0.0,
			"[0..] figure stage: weight of the bone distance quantile function image match.  As other matches, this can is a complementary match");
	setNext("FMAHAM", "FigureMahalanobisMatch",         0.0,
			"[0..] figure stage: weight of Mahalanobis distance from reference model to optimized model: in units of standard deviations of PGA/Eigenmode Coordinates; #modes is a good (expected) value");
	setNext("FMAHAMAX","FigureMaxPGAModes",             99.0,
			"[0..] figure stage: for Mahalanobis distance: maximum number of PGA/Eigenmode Modes used; 0 disables PGA");
	setNext("FMAHAONLY","FigureMahalanobisOnly",        0,
			"[0,1] figure stage: During optimization allow similarity transform: 0 enable, 1 disable simTrans optimization");
	setNext("FSGEOMWT",  "FigureGeometryWeight",        0.5,
			"[0, 1] figure stage: d^2 geometry penalty weight");
	setNext("FSCONSWT",  "FigureInterObjectWeight",     0.5,
			"[0, 1] figure stage: constraints penalty weight");
	setNext("FPROPEN","FigureProfileMatch",				1,
			"Figure stage: Image match penalty; cannot be changed",
			true);	// Locked; see Match.cpp; this parameter cannot be used in a script: the weight is always 1.0.

	// figure id for multiple single objects
	setNext("FIP",		"FigureInterpenetrationMatch",	0,		"[0..] figure stage: weight of penalty for interpentration");
	setNext("IPM",		"InterpenetrationMask",			0x0,	"[bit flags] atom/figure stage: bitmask of image planes to avoid via penetration penalty; see figureInterpenetrationMatch");
	setNext("FID",		"FigureID",						0,		"[0..] atom/figure stage: specifies the single figure to use in figure/atom stage or to initialize (see momMask); other figures NEVER change");

	setNext("FOT",    "FigureTranslationFactor",        1.0,	"[0..] figure stage: weight of translation component in objective function");
	setNext("FOR",    "FigureRotationFactor",           8.0,	"[0..] figure stage: weight of rotation component in objective function");
	setNext("FOS",    "FigureScaleFactor",              8.0,	"[0..] figure stage: weight of scale component in objective function");
	setNext("FOE",    "FigureElongationFactor",         16.0,	"[0..] figure stage: weight of figure elongation component in objective function; 0 disables");
	setNext("FOPHI",  "FigureTubePhiFactor",                16.0,	"[0..] figure stage: weight of phi component in objective function (used only in case of tubes)");
	setNext("FOPGA",  "FigurePGAFactor",                200.0,	"[0..] figure stage: weight of PGA component in objective function");
	setNext("FEPS",   "FigureEpsilon",                  -1.0,	"(0..] EXPERTS ONLY. figure stage: step size in optimizer; in units of unit cube; omit for built-in default");
	setNext("FBT",    "FigureTranslationBound",        0.25,	"[0..] figure stage: The maximum allowed translation during optimization. See BrentLinearSearchBoundFactor for more details.");
	setNext("FBR",    "FigureRotationBound",           180.0,	"[0..] figure stage: The maximum allowed rotation in degrees. See BrentLinearSearchBoundFactor for more details.");
	setNext("FBS",    "FigureScaleBound",              2.0,	"[0..] figure stage: The maximum allowed change in scale. See BrentLinearSearchBoundFactor for more details.");
	setNext("FBE",    "FigureElongationBound",         2.0,	"[0..] figure stage: The maximum allowed elongation. See BrentLinearSearchBoundFactor for more details.");
	setNext("FBPHI",  "FigureTubePhiBound",            180.0,	"[0..] figure stage: maximum allowed change in phi in degrees. See BrentLinearSearchBoundFactor for more details.");
	setNext("FBPGA",  "FigurePGABound",                3.0,	"[0..] figure stage: maximum amount of change in a PGA direction in units of standard deviations. See BrentLinearSearchBoundFactor for more details.");

	// FIGURE contour based initialization
	setNext("FCM",    "FigureContourMatch",             0,     "[0..] Weight on the contour match penalty; 0 disables");
	setNext("FCL",    "FigureContourLevel",             2,     "[0..] Subdivision level for matching surfaces to contours;  Default is 2");
	setNext("FCC",    "FigureContourCloseness",         0.04,  "[0..1] Constraint on regional correspondence; 1 disables; Default is 0.04");
	setNext("FCF",    "FigureContourFigureId",          0,     "[0..] which figure to match to contours; Default is 0");
	setNext("FCF",    "FigureContourCoords",            0,     "[0 or 1] Contours are in world (vs. model) coordinates; Default is 0");
//	setNext("FCP",    "FigureContourPadding",           1,     "[0..] number of axial slices between the end of the object and the end contour ; Default is 1");
//	setNext("FCS",    "FigureContourSizeWeight",        0,     "[0..] weight of term to encourage model to span the same number of axial slices as the contours ; Default is 0");

	// Atom stage
	setNextList("ASKPFIGS", "AtomSkipFigures",			"List of figure numbers for which the atom stage will be skipped");
	setNext("AHIM",  "AtomHistogramMatch",		        0,
			"[0..] atom stage: weight of histogram image match.  Use instead of FIM for grey image segmentation");
	setNext("AMPFIM",   "AtomMultiPatchMatch",		        0,
			"[0..] atom stage: weight of multipatch image match.  Use instead of FIM for grey image segmentation");
	setNext("AMAHAM", "AtomMahalanobisMatch",			0.0,
		     "[0..] atom stage: weight of Mahalanobis distance from reference model to optimized model: in units of standard deviations of PGA/Eigenmode Coordinates; 1 is a good value");
	setNext("AGEOMWT",  "AtomGeometryWeight",           0.5,
			"[0, 1] atom stage: geometry weight");
	setNext("ANBRWT",  "AtomNeighborWeight",			0.5,
			"[0, 1] atom stage: neighbor penalty weight");
	setNext("ACONSWT",  "AtomInterWeight",              0.5,
			"[0, 1] atom stage: constraints penalty weight");
	setNext("APROPEN","AtomProfileMatch",				1,
			"Atom stage: Image match penalty; cannot be changed",
			true);	// Locked; see Match.cpp; this parameter cannot be used in a script: the weight is always 1.0.

//	setNext("ARESET", "AtomResetModelBefore",           1.0,
//			"[0,1] atom stage: set the reference model to the model as the stage is entered. The model match is calculated relative to the reference model");
//	setNext("ALEV",   "AtomSurfaceLevel",               2,
//			"[0,1,2,3] figure stage: 0:no subdivision, 1:each tile is split into 4 tiles. 2:16 tiles, 3:64 tiles");

	// The smaller of the two following settings will win
	setNext("AITER",  "AtomIterations",                 -1,
			"[0..] atom stage: maximum number of passes through the figure; 0 disables the atom stage; -1 unlimited");
	setNext("ASITER",  "AtomSubIterations",                 -1,
			"[0..] atom stage: maximum number of atom fits to be done; 0 disables the atom stage; -1 unlimited");

/*	setNext("AAN",    "AtomAverageNeighborMatch",       0.5,
			"[0..] atom stage: (model regularity) weight of distance of each atom from it's position as predicted by the atom's neighbors; see atomNeighborhood, atomDistanceType");
	setNext("AIM",    "AtomImageMatch",                 0.5,
			"[0..] atom stage: (distance from truth) weight of distance from model to nearest edge in image");
	setNext("ALM",    "AtomLandmarkMatch",              0.0,
			"[0..] atom stage: (correspondance with truth) weight of distance from model landmarks to corresponding image landmark; 0 disables");
	setNext("AMM",    "AtomModelMatch",                 0.5 * 100000/5,
			"[0..] (model reference) atom stage: weight of distance from reference model to optimized model; currently GEODESIC distance");
	setNext("ANEI",   "AtomNeighborhood",               3,
			"[0,1,2,3] atom stage: neighborhood of each atom: 0:ALL_FIRST_NEIGHBORS, 1:PIN_CORNERS, 2:PIN_CORNERS_AND_EDGES, 3:EDGES_HAVE_2_NEIGHBORS");
	setNext("ADIST",  "AtomDistanceType",               3,
			"[0,1,2,3] atom stage: distances between atoms:0:EUCLIDEAN_DIST, 1:AVE_EUCLIDEAN_DIST (Euclidian divided by number of neighbors); 2:GEODESIC_DIST (radius-scaled pose, position), 3:AVE_GEODESIC_DIST (Geodesic divided by number of neighbors)"),
	setNext("AIP",	   "AtomInterpenetrationMatch",		 0,
			"[0..] atom stage: weight of penalty for interpentration");
*/
//	setNext("AMN",    "ModelNormals",                   0.5,	"[0..] for future use");
//	setNext("AMS",    "ModelSliding",                   0.5,	"[0..] for future use");
//	setNext("AMP",    "ModelPenetration",               1.0,	"[0..] for future use");
//	setNext("AMC",    "ModelCreasing",                  0.0005,	"[0..] for future use");

/*	setNext("AST",    "AtomSkipThreshold",              1.0e-9,
			"(0..] EXPERTS ONLY. atom stage: amount an atom must change to be considered next pass");
	setNext("ASM",    "AtomSkipMultiplier",             2,
			"[1,2] EXPERTS ONLY. atom stage: higher values make it more probable atoms will be skipped");
	setNext("ACG",    "AtomConjugateGradientIterations",2,
			"[1,2] EXPERTS ONLY. atom stage: number of conjugate direction changes within each CG iteration");
	setNext("ASP",    "AtomSkippedPasses",              1,
			"[1,2] EXPERTS ONLY. atom stage: number of times all atoms must be skipped before being considered done");
	setNext("ARAND",  "AtomRandomOrder",                1,
			"[0,1] atom stage: 1:atoms are optimized in random order. 0:atoms are optimized in order found in model");
*/
	setNext("AOT",    "AtomTranslationFactor",          0.375,
			"[0..] atom stage: weight of translation component in objective function");
	setNext("AOR",    "AtomRotationFactor",             11.25,
			"[0..] atom stage: weight of rotation component in objective function");
	setNext("AOS",    "AtomScaleFactor",                6.0,
			"[0..] atom stage: weight of scale component in objective function");
	setNext("AOTH",   "AtomThetaFactor",                10.5,
			"[0..] atom stage: weight of theta (angle between spokes) component in objective function");
//	setNext("AOE",    "AtomElongationFactor",           0.5,
//			"[0..] atom stage: weight of atom elongation component in objective function");
	setNext("AON",    "AtomEtaFactor",                  8.0,
			"[0..] atom stage: weight of eta (B-spoke extension) component in objective function");
	setNext("AOPGA", "AtomPGAFactor",					200.0,
		     "[0..] atom stage: weight of PGA components");
	setNext("AEPS",   "AtomEpsilon",                    0.003,
			"(0..] EXPERTS ONLY. atom stage: unit cube step size in optimizer");

	setNext("ARM",    "AtomRestrictMovement",           0.0,
			"A bit vector indicating restriction of atom motion, bits 0-2 prevent translation "
			"along tangent, normal and binormal directions respectively, bit 4 couples r and "
			"theta changes so that the spoke ends stay in the plane defined by the normal and "
			"the bisector."
			"bit 3 changes the behavior of bits 0-2 to mean translation along x,y and z directions before any rotation is applied.");

	setNext("ABT",    "AtomTranslationBound",        0.10,	"[0..] atom stage: The maximum allowed translation during optimization. See BrentLinearSearchBoundFactor for more details.");
	setNext("ABR",    "AtomRotationBound",           90.0,	"[0..] atom stage: The maximum allowed rotation in degrees. See BrentLinearSearchBoundFactor for more details.");
	setNext("ABS",    "AtomScaleBound",              2.0,	"[0..] atom stage: The maximum allowed change in scale. See BrentLinearSearchBoundFactor for more details.");
	setNext("FBE",    "AtomThetaBound",         90.0,	"[0..] atom stage: The maximum allowed change in theta. See BrentLinearSearchBoundFactor for more details.");
	setNext("FBPHI",  "AtomEtaBound",            2.0,	"[0..] atom stage: The maximum allowed multiplicative change in eta. See BrentLinearSearchBoundFactor for more details.");
	setNext("FBPGA",  "AtomPGABound",                3.0,	"[0..] atom stage: maximum amount of change in a PGA direction in units of standard deviations. See BrentLinearSearchBoundFactor for more details.");

	// SPOKE stage
	setNextList("ASKPFIGS", "SpokeSkipFigures",			"List of figure numbers for which the spoke stage will be skipped");
	setNext("AHIM",  "SpokeHistogramMatch",		        0,
			"[0..] spoke stage: weight of histogram image match.  Use instead of FIM for grey image segmentation");
	setNext("AMAHAM", "SpokeMahalanobisMatch",			0.0,
		     "[0..] spoke stage: weight of Mahalanobis distance from reference model to optimized model: in units of standard deviations of PGA/Eigenmode Coordinates; 1 is a good value");
	setNext("AGEOMWT",  "SpokeGeometryWeight",           0.5,
			"[0, 1] spoke stage: geometry weight");
	setNext("ANBRWT",  "SpokeNeighborWeight",			0.5,
			"[0, 1] spoke stage: neighbor penalty weight");
	setNext("ACONSWT",  "SpokeInterWeight",              0.5,
			"[0, 1] spoke stage: constraints penalty weight");
	setNext("APROPEN","SpokeProfileMatch",				1,
			"Spoke stage: Image match penalty; cannot be changed",
			true);	// Locked; see Match.cpp; this parameter cannot be used in a script: the weight is always 1.0.
	setNext("AITER",  "SpokeIterations",                 0,
			"[0..] spoke stage: max number of passes through the figure; 0 disables; -1 unlimited");

	/*
	setNext("AOT",    "SpokeTranslationFactor",          0.375,
			"[0..] spoke stage: weight of translation component in objective function");
	setNext("AOR",    "SpokeRotationFactor",             11.25,
			"[0..] spoke stage: weight of rotation component in objective function");
	setNext("AOS",    "SpokeScaleFactor",                6.0,
			"[0..] spoke stage: weight of scale component in objective function");
	setNext("AOTH",   "SpokeThetaFactor",                10.5,
			"[0..] spoke stage: weight of theta (angle between spokes) component in objective function");
	setNext("AON",    "SpokeEtaFactor",                  8.0,
			"[0..] spoke stage: weight of eta (B-spoke extension) component in objective function");
	*/
	setNext("SOPGA", "SpokePGAFactor",					200.0,
		     "[0..] spoke stage: weight of PGA components");
	setNext("SOSS",  "SpokeScaleFactor",                    7.0,
			"[0..] spoke stage: weight of spoke scale component in objective function");
	setNext("AEPS",  "SpokeEpsilon",                    0.003,
			"(0..] EXPERTS ONLY. spoke stage: unit cube step size in optimizer");
	setNext("SBPGA",  "SpokePGABound",                3.0,	"[0..] spoke stage: maximum amount of change in a PGA direction in units of standard deviations. See BrentLinearSearchBoundFactor for more details.");
	setNext("SBSS",  "SpokeScaleBound",                2.0,	"[0..] spoke stage: maximum allowed change in the scale of a spoke. See BrentLinearSearchBoundFactor for more details.");

	// rSrad
	setNext("RSPN",  "RSradPenaltyNorm",               1.0,	"(0..] The norm to take when computing the summary of all the rSrad penalty values.");
	setNext("RSPT",  "RSradPenaltyThreshold",          0.0,	"[0..1.0) The threshold above which rSrad penalties count.");

	// Outputs
	setNext("TLEV", "TileSurfaceLevel",					2,
			"[0,1,2,3] Tile output: 0:no subdivision, 1:each tile is split into 4 tiles. 2:16 tiles, 3:64 tiles");
	setNext("TQUAD", "TileQuads",						0,
			"[0,1] tile output: 0:output triangles. 1:output quadrilaterals; there will be twice as many tri's as quad's");
	// Use one of the other of the following two parameters
	setNext("TFM", "TileFiguresMask",					0x0,
			"[bit flags] Bit mask indicating for which figures tiles will be output");
	setNextList("TFL", "TileFiguresList",				"List of figures for which tiles will be output");
	setNext("TSEP", "TileSeparateFiles",				0,
			"[0, 1] Whether or not to use multiple files to receive tile output");
	setNext("TWORLD",  "TileUnit",                      0,
			"[0,1] tile output: 0:output model tile vertices in world coordinates if image is loaded, 1: unit coordinates");

	// Input Files
	// These have no default values
	setNextString("IMG", "Image",                      "Image file path");
#ifdef AE2_BUILD
	setNextString("LDIMG", "LoadedImage",               "Loaded image file path");
	setNext("LDMOD", "LoadedModel",					   0, "Use loaded model");
#endif
	setNextString("MOD", "Model",                      "Model file path");
	setNextString("PGAMM","PGAMeanModel",            "PGA model file path");
	setNextString("SYMT", "Transform",                 "Symmetry transform file path");
	setNextString("SMASK", "SimpleMaskFile",
			"Full path of the simple mask file (*.smsk) to be used during optimization");
//	setNextString("MASK", "MaskFile",
//			"Full path of the mask file (*.msk) to be used during optimization");
	setNextString("REFM", "ReferenceModel",
			"Full path of the reference model to be used during optimization");
	setNextString("TNGI", "TrainingImage",
			"Full path of the training image to be used during optimization");
	setNextString("PROF", "ProfileFile",
			"Full path of the profile file (*.stats) to be used during optimization");
	setNextString("HIST", "HistogramFile",
			"Full path of the surface patch ensemble file to be used during optimization");
	setNextString("HISTMEAN", "HistogramMeanFile",
			"Full path of the surface patch ensemble file that overrides the mean vector in hist");
//	setNextString("LMOD", "LandmarkModel",             "Landmark model file path");
	setNextString("OMOD", "OutModel",                  "Output model file path");
	setNextString("OMI", "OutModelImage",              "Output file path for  model image");
	setNextString("OTIL", "OutTile",                   "Output tiles file path");
	setNextString("OIT", "OutImageTile",               "Output file path for image tiles");
	setNextString("CNT", "Contour",				       "Contour file path");

	setNext("OMIX", "OutModelImageX", 128, "for outModelImage, x dimension of image");
	setNext("OMIY", "OutModelImageY", 128, "for outModelImage, y dimension of image");
	setNext("OMIZ", "OutModelImageZ", 128, "for outModelImage, z dimension of image");
	setNext("OMILEV","OutModelImageSurfaceLevel", 6, "for outModelImage, subdivision level of model");
	setNext("OMIFS", "OutModelImageFixedSize", 0, "for outModelImage, do not allow mapsize to adjust to fit model, eg, image may leak if model's surface tiles are larger than the image the model is being scan-converted onto");

	// TUNE_COUNT is here in gpTuning.h

	resetAllWeights();

	if (index != numEntries) {
		cout << "Construction of gpTuning failed: " << index << " parameters inserted; "
			<< numEntries << " expected" << endl;
		return false;
	}
	return true;
}

#endif
