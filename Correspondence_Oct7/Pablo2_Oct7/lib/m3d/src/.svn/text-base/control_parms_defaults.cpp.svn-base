
#include <stdlib.h>
#include "ControlParms.h"


// Global function intended to be called from main() at program start-up
void setGlobalControlDefaults()
{
	float default_collision_color[3] = { 0.5, 0.5, 0.5 };
	float default_connector_color[3] = { 0.0, 0.5, 0.0 };
	float default_tiles_color[3] = { 0.0, 0.0, 1.0 };	// Same as in TileSet.cpp

	globalControl->setDefault(IconifyMode, 1);	// Iconify by Main Window
	globalControl->setDefault(PartialSurfaceColor_R, default_collision_color[0]);
	globalControl->setDefault(PartialSurfaceColor_G, default_collision_color[1]);
	globalControl->setDefault(PartialSurfaceColor_B, default_collision_color[2]);
    globalControl->setDefault(ShowConstraints, false);
	globalControl->setDefault(DrawCutPlaneBoundary, true);
	globalControl->setDefault(CutPlaneBoundaryWidth, 1);
	globalControl->setDefault(CutPlaneMode, 0);
#ifdef BINARY
	globalControl->setDefault(SmoothImages, false);	// No smoothing for binary images
#else
	globalControl->setDefault(SmoothImages, true);
#endif
	globalControl->setDefault(ScaleImages, true);
	globalControl->setDefault(DrawBoundary, false);
	globalControl->setDefault(DisplayRangeUnits, false);
	globalControl->setDefault(DisplayModelUnits, false);
	globalControl->setDefault(ByteOrder, 1);	// Native
	globalControl->setDefault(CompressImages, true);
	globalControl->setDefault(ConvertImages, false);
	globalControl->setDefault(ConvertImageFormat, false);
	globalControl->setDefault(ImageFormat, 1);
	globalControl->setDefault(AxialSliceDefault, true);
	globalControl->setDefault(CoronalSliceDefault, false);
	globalControl->setDefault(SagittalSliceDefault, false);
	globalControl->setDefault(MeshConnectorsType, false);
	globalControl->setDefault(MeshConnectorsLineWidth, 1);
	globalControl->setDefault(ShowAtoms, true);
	globalControl->setDefault(FiguresColorAtoms, false);
	globalControl->setDefault(ShowAtomVectors, true);
	globalControl->setDefault(ShowMeshConnectors, true);
	globalControl->setDefault(ExtraAtomVectors, false);
	globalControl->setDefault(AtomVectorsType, false);
	globalControl->setDefault(BVectorsType, 2);	// Show crest B vectors only
	globalControl->setDefault(AtomVectorsLineWidth, 1);
	globalControl->setDefault(ConnectorsColor_R, default_connector_color[0]);
	globalControl->setDefault(ConnectorsColor_G, default_connector_color[1]);
	globalControl->setDefault(ConnectorsColor_B, default_connector_color[2]);
	globalControl->setDefault(UndoListLength, 250);
    globalControl->setDefault(ShowLandmarks, true);
    globalControl->setDefault(ReorderModels, -1);
    globalControl->setDefault(BYUOutputType, 0);
    globalControl->setDefault(BYUOutputCoords, 0);
    globalControl->setDefault(SurfaceSmoothnessDefault, 50);	// Default to very smooth surfaces
    globalControl->setDefault(TwoLights, false);
	globalControl->setDefault(DisplayStdAxes, false);
	globalControl->setDefault(ModelDirectory, ".");
	globalControl->setDefault(ImageDirectory, "");
	globalControl->setDefault(TileSetDirectory, "");
	globalControl->setDefault(TilesColor_R, default_tiles_color[0]);
	globalControl->setDefault(TilesColor_G, default_tiles_color[1]);
	globalControl->setDefault(TilesColor_B, default_tiles_color[2]);
	globalControl->setDefault(SimTransformSeparate, false);
	globalControl->setDefault(SimTransformMatrix, false);
}




