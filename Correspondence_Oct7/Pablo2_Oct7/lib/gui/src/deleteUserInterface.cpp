
#include "P3DUserInterface.h"

/*	This file is generated from P3DUserInterface.cpp by
	grepping for "new Fl_Window" and "new movable_Fl_Window".
	The few Fl_Window's are indicated by comments below.
	The deletions are in reverse order from that of their
	allocation.  Deleting an FLTK window will result in
	all its subordinate widgets being deleted.

	This function is called from pablo.cpp.
*/

void deleteUserInterface(P3DUserInterface * userInterface)
{
#ifndef PRODUCTION_VERSION
    delete userInterface->regularizerWindow;
    userInterface->regularizerWindow = NULL;
    delete userInterface->matchSurfacesDialog;		// Fl_Window
	userInterface->matchSurfacesDialog = NULL;
    delete userInterface->testSeuratDialog;		// Fl_Window
	userInterface->testSeuratDialog = NULL;
    delete userInterface->testConstraintsDialog;		// Fl_Window
	userInterface->testConstraintsDialog = NULL;
    delete userInterface->subfigureTransformationTestDialog;		// Fl_Window
	userInterface->subfigureTransformationTestDialog = NULL;
    delete userInterface->interpolatedPrimitiveWindow;
    userInterface->interpolatedPrimitiveWindow = NULL;
#endif
    delete userInterface->transformationRecordingTestDialog;		// Fl_Window
	userInterface->transformationRecordingTestDialog = NULL;
    delete userInterface->pgaDialog;
    userInterface->pgaDialog = NULL;
    delete userInterface->cpnsDialog;
    userInterface->cpnsDialog = NULL;
    delete userInterface->modelSlideShowDialog;
    userInterface->modelSlideShowDialog = NULL;
    delete userInterface->optimizerSettingsDialog;
    userInterface->optimizerSettingsDialog = NULL;
    delete userInterface->optimizerControlDialog;
    userInterface->optimizerControlDialog = NULL;
    delete userInterface->elongationDialog;
    userInterface->elongationDialog = NULL;
    delete userInterface->attachSubfigureDialog;
    userInterface->attachSubfigureDialog = NULL;
    delete userInterface->editModelWindow;
    userInterface->editModelWindow = NULL;
    delete userInterface->editLandmarksWindow;
    userInterface->editLandmarksWindow = NULL;
    delete userInterface->InvoluteCutPlaneWindow;
    userInterface->InvoluteCutPlaneWindow = NULL;
    delete userInterface->bPerpY1CutPlaneWindow;
    userInterface->bPerpY1CutPlaneWindow = NULL;
    delete userInterface->bPerpY0CutPlaneWindow;
    userInterface->bPerpY0CutPlaneWindow = NULL;
    delete userInterface->bPerpNCutPlaneWindow;
    userInterface->bPerpNCutPlaneWindow = NULL;
    delete userInterface->bBPerpCutPlaneWindow;
    userInterface->bBPerpCutPlaneWindow = NULL;
    delete userInterface->bNCutPlaneWindow;
    userInterface->bNCutPlaneWindow = NULL;
    delete userInterface->cutPlanesControlWindow;
    userInterface->cutPlanesControlWindow = NULL;
    delete userInterface->primitiveEditorWindow;
    userInterface->primitiveEditorWindow = NULL;
    delete userInterface->constraintsWindow;
    userInterface->constraintsWindow = NULL;
    delete userInterface->visibilityControlWindow;
    userInterface->visibilityControlWindow = NULL;
    delete userInterface->addQuadFigureDlg;
    userInterface->addQuadFigureDlg = NULL;
    delete userInterface->reorderPopupWindow;		// Fl_Window
	userInterface->reorderPopupWindow = NULL;
	delete userInterface->prefsImageFilesWindow;
	userInterface->prefsImageFilesWindow = NULL;
    delete userInterface->preferencesEditorWindow;
    userInterface->preferencesEditorWindow = NULL;
#ifdef BINARY
    delete userInterface->aboutBinaryPabloWindow;
    userInterface->aboutBinaryPabloWindow = NULL;
#else
    delete userInterface->aboutPabloWindow;
    userInterface->aboutPabloWindow = NULL;
#endif
    delete userInterface->displayControlWindow;
    userInterface->displayControlWindow = NULL;

	// These are deleted in the destructor
//    delete userInterface->mainWindow;
//    delete userInterface->modelWindow;
}


