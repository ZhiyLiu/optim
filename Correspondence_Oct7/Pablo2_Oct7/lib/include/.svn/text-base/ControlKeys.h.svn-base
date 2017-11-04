#ifndef PABLO_CONTROL_KEYS_H
#define PABLO_CONTROL_KEYS_H


/*
    This file specifies the names of user preferences within Pablo.

    See files ControlParms.h, ControlParmsAccess.h and ControlFile.h
    in directory m3d/include for the technical details.

    The pablo_options.txt file distributed with Pablo contains a list
    defining the possible contents of the control file.  It should be
    maintained to correspond to the contents of this file.
*/

#define APPLICATION_NAME    "pablo"
#define NUM_CONTROL_ITEMS	157

/*  The names in the following array must correspond to the members of this
    enumeration, and the number of members must be defined above.
*/
enum control_keys_t {
	// Pablo control parameters.
	// If any of these change, then the corresponding change must be
	// made to the source of ConStruct in RadOnc.
	AtomVectorsLineWidth, AtomVectorsType, BVectorsType, FiguresColorAtoms, ShowAtoms, ShowAtomVectors,
	BackgroundColor_B, BackgroundColor_G, BackgroundColor_R, DrawCutPlaneBoundary,
	CutPlaneBoundaryWidth, CutPlaneMode, MeshConnectorsType, ShowMeshConnectors, MeshConnectorsLineWidth,
	PartialSurfaceRendering, PartialSurfaceColor_B, PartialSurfaceColor_G, PartialSurfaceColor_R,
	PartialSurfaceLevel, PartialSurfaceStyle, ShowConstraints, RockingAngle, RockingIncrement,
	ExtraAtomVectors, ConnectorsColor_R, ConnectorsColor_G, ConnectorsColor_B, SmoothImages, ScaleImages,
	DisplayRangeUnits, DisplayModelUnits, CompressImages, ConvertImages, ConvertImageFormat, ImageFormat,
	AxialSliceDefault, CoronalSliceDefault, SagittalSliceDefault, DrawBoundary, SurfaceLevel,
	SurfaceStyle, SurfaceLineWidth, SurfaceSmoothnessDefault, UndoListLength, LeftHandedMouse,
	ShowLandmarks, LandmarksColor_R, LandmarksColor_G, LandmarksColor_B, LandmarkNarrowWidth,
	LandmarkWideWidth, ReorderModels, BYUOutputType, BYUOutputCoords, DisplayStdAxes, ModelDirectory,
	ImageDirectory, TileSetDirectory, TilesColor_R, TilesColor_G, TilesColor_B, OutputVerbosity,
	ByteOrder, IconifyMode, TwoLights, SimTransformSeparate, SimTransformMatrix,

	// Window position parameters
	MainWindow_X, MainWindow_Y, ModelWindow_X, ModelWindow_Y, DisplayControlWindow_X,
	DisplayControlWindow_Y, VisibilityControlWindow_X, VisibilityControlWindow_Y,
	ConstraintsWindow_X, ConstraintsWindow_Y, AtomEditorWindow_X, AtomEditorWindow_Y,
	CutPlanesControlWindow_X, CutPlanesControlWindow_Y, OptimizerControlWindow_X,
	OptimizerControlWindow_Y, InterpolatedPrimitiveWindow_X, InterpolatedPrimitiveWindow_Y,
	CrestPlaneWindow_X, CrestPlaneWindow_Y, AtomPlaneWindow_X, AtomPlaneWindow_Y,
	BPerpNPlaneWindow_X, BPerpNPlaneWindow_Y, PortPlaneWindow_X, PortPlaneWindow_Y, 
	StarboardPlaneWindow_X, StarboardPlaneWindow_Y, InvolutesPlaneWindow_X,
	InvolutesPlaneWindow_Y, AddQuadFigureWindow_X, AddQuadFigureWindow_Y, AttachSubfigureWindow_X,
	AttachSubfigureWindow_Y, PreferencesEditorWindow_X, PreferencesEditorWindow_Y,
	ImagePrefsSubEditorWindow_X, ImagePrefsSubEditorWindow_Y, EditModelPropsWindow_X,
	EditModelPropsWindow_Y, EditLandmarksWindow_X, EditLandmarksWindow_Y, ElongationWindow_X,
	ElongationWindow_Y, AboutPablo_X, AboutPablo_Y,
	PenaltyWeightsWindow_X, PenaltyWeightsWindow_Y, SlideShowWindow_X, SlideShowWindow_Y,
	RegularizerWindow_X, RegularizerWindow_Y, MatchSurfaces_X, MatchSurfaces_Y, CPNSDeformWindow_X,
	CPNSDeformWindow_Y, PGADeformWindow_X, PGADeformWindow_Y, PCADeformWindow_X, PCADeformWindow_Y, 
	OptVisualizerWindow_X, OptVisualizerWindow_Y,

	// Window open/closed parameters
	MainWindow_Open, ModelWindow_Open, DisplayControlWindow_Open, VisibilityControlWindow_Open,
	ConstraintsWindow_Open, AtomEditorWindow_Open, CutPlanesControlWindow_Open,
	OptimizerControlWindow_Open, InterpolatedPrimitiveWindow_Open, CrestPlaneWindow_Open,
	AtomPlaneWindow_Open, BPerpNPlaneWindow_Open, PortPlaneWindow_Open, StarboardPlaneWindow_Open,
	InvolutesPlaneWindow_Open, AddQuadFigureWindow_Open,  AttachSubfigureWindow_Open,
	PreferencesEditorWindow_Open, ElongationWindow_Open, AboutPablo_Open,
	PenaltyWeightsWindow_Open, SlideShowWindow_Open, RegularizerWindow_Open, MatchSurfaces_Open,
	CPNSDeformWindow_Open, PGADeformWindow_Open, OptVisualizerWindow_Open
};

static const char * CONTROL_NAMES[NUM_CONTROL_ITEMS] = {
	// Pablo control parameters
	"AtomVectorsLineWidth", "AtomVectorsType", "BVectorsType", "FiguresColorAtoms", "ShowAtoms",
	"ShowAtomVectors", "BackgroundColor_B", "BackgroundColor_G", "BackgroundColor_R", "DrawCutPlaneBoundary",
	"CutPlaneBoundaryWidth", "CutPlaneMode", "MeshConnectorsType", "ShowMeshConnectors",
	"MeshConnectorsLineWidth", "PartialSurfaceRendering", "PartialSurfaceColor_B", "PartialSurfaceColor_G",
	"PartialSurfaceColor_R", "PartialSurfaceLevel", "PartialSurfaceStyle", "ShowConstraints",
	"RockingAngle", "RockingIncrement", "ExtraAtomVectors", "ConnectorsColor_R", "ConnectorsColor_G",
	"ConnectorsColor_B", "SmoothImages", "ScaleImages", "DisplayRangeUnits", "DisplayModelUnits",
	"CompressImages", "ConvertImages", "ConvertImageFormat", "ImageFormat", "AxialSliceDefault",
	"CoronalSliceDefault", "SagittalSliceDefault", "DrawBoundary", "SurfaceLevel", "SurfaceStyle",
	"SurfaceLineWidth", "SurfaceSmoothnessDefault", "UndoListLength", "LeftHandedMouse",
	"ShowLandmarks", "LandmarksColor_R", "LandmarksColor_G", "LandmarksColor_B", "LandmarkNarrowWidth",
	"LandmarkWideWidth", "ReorderModels", "BYUOutputType", "BYUOutputCoords", "DisplayStdAxes",
	"ModelDirectory", "ImageDirectory", "TileSetDirectory", "TilesColor_R", "TilesColor_G", "TilesColor_B",
	"OutputVerbosity", "ByteOrder", "IconifyMode", "TwoLights", "SimTransformSeparate", "SimTransformMatrix",

	// Window position parameters
	"MainWindow_X", "MainWindow_Y", "ModelWindow_X", "ModelWindow_Y", "DisplayControlWindow_X",
	"DisplayControlWindow_Y", "VisibilityControlWindow_X", "VisibilityControlWindow_Y",
	"ConstraintsWindow_X", "ConstraintsWindow_Y", "AtomEditorWindow_X", "AtomEditorWindow_Y",
	"CutPlanesControlWindow_X", "CutPlanesControlWindow_Y", "OptimizerControlWindow_X",
	"OptimizerControlWindow_Y", "InterpolatedPrimitiveWindow_X", "InterpolatedPrimitiveWindow_Y",
	"CrestPlaneWindow_X", "CrestPlaneWindow_Y", "AtomPlaneWindow_X", "AtomPlaneWindow_Y",
	"BPerpNPlaneWindow_X", "BPerpNPlaneWindow_Y", "PortPlaneWindow_X", "PortPlaneWindow_Y",
	"StarboardPlaneWindow_X", "StarboardPlaneWindow_Y", "InvolutesPlaneWindow_X",
	"InvolutesPlaneWindow_Y", "AddQuadFigureWindow_X", "AddQuadFigureWindow_Y",
	"AttachSubfigureWindow_X", "AttachSubfigureWindow_Y", "PreferencesEditorWindow_X",
	"PreferencesEditorWindow_Y", "ImagePrefsSubEditorWindow_X", "ImagePrefsSubEditorWindow_Y", 
	"EditModelPropsWindow_X", "EditModelPropsWindow_Y", "EditLandmarksWindow_X",
	"EditLandmarksWindow_Y", "ElongationWindow_X", "ElongationWindow_Y", "AboutPablo_X",
	"AboutPablo_Y", "PenaltyWeightsWindow_X",
	"PenaltyWeightsWindow_Y", "SlideShowWindow_X", "SlideShowWindow_Y", "RegularizerWindow_X",
	"RegularizerWindow_Y", "MatchSurfaces_X", "MatchSurfaces_Y", "CPNSDeformWindow_X", "CPNSDeformWindow_Y",
	"PGADeformWindow_X", "PGADeformWindow_Y", "PCADeformWindow_X", "PCADeformWindow_Y", 
	"OptVisualizerWindow_X", "OptVisualizerWindow_Y",

	// Window open/closed parameters
	"MainWindow_Open", "ModelWindow_Open", "DisplayControlWindow_Open",
	"VisibilityControlWindow_Open", "ConstraintsWindow_Open", "AtomEditorWindow_Open",
	"CutPlanesControlWindow_Open", "OptimizerControlWindow_Open",
	"InterpolatedPrimitiveWindow_Open", "CrestPlaneWindow_Open", "AtomPlaneWindow_Open",
	"BPerpNPlaneWindow_Open", "PortPlaneWindow_Open", "StarboardPlaneWindow_Open",
	"InvolutesPlaneWindow_Open", "AddQuadFigureWindow_Open", "AttachSubfigureWindow_Open",
	"PreferencesEditorWindow_Open", "ElongationWindow_Open", "AboutPablo_Open",
	"PenaltyWeightsWindow_Open", "SlideShowWindow_Open", "RegularizerWindow_Open",
	"MatchSurfaces_Open", "CPNSDeformWindow_Open", "PGADeformWindow_Open", "OptVisualizerWindow_Open"
};




#endif	/* PABLO_CONTROL_KEYS_H */

