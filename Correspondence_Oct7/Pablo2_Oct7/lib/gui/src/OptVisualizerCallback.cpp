#include <stdio.h>
#include "OptVisualizerCallback.h"
#include "globalBuildParms.h"
#ifdef OPTIMIZATION_VISUALIZER

#include <FL/Fl_Choice.H>
#include <FL/Fl_Button.H>
#include "P3DUserInterface.h"
#include "P3DUserInterfaceCallback.h"
#include "LogManagerChart.h"
#include "Plot3DWindowWrapper.h"
#include "OptVisualizerUI.h"


OptVisualizerCallback::~OptVisualizerCallback()
{
	if (attributeBtns != NULL)
		delete [] attributeBtns;
}

void OptVisualizerCallback::redrawChart()
{
	oui->updateChart();
	oui->newChart->damage(1);
	oui->plot3DWrapper->damage(1);
}

void OptVisualizerCallback::init()
{
	if (initialized)
		return;

	if (oui == NULL)
		return;
	p3dUI = oui->p3dUI;
	if (p3dUI->callback == NULL)
		return;
	p3dUICallback = p3dUI->callback;

	if (p3dUICallback == NULL)
		return;

	oui->make_chart(p3dUICallback->getControl());
	oui->chartGrp->add(oui->newChart);
	oui->chartGrp->end();

	globalLogManager.initialize();
	attributeBtns = new Fl_Button*[MAX_ATTRIBUTES];
	if (MAX_ATTRIBUTES != 11) {
		cout <<
			"OptVisualizerCallback::MAX_ATTRIBUTES must be changed\n";
		// See LogManager.h
		return;
	}
	attributeBtns[0] = oui->attributeBtn0;
	attributeBtns[1] = oui->attributeBtn1;
	attributeBtns[2] = oui->attributeBtn2;
	attributeBtns[3] = oui->attributeBtn3;
	attributeBtns[4] = oui->attributeBtn4;
	attributeBtns[5] = oui->attributeBtn5;
	attributeBtns[6] = oui->attributeBtn6;
	attributeBtns[7] = oui->attributeBtn7;
	attributeBtns[8] = oui->attributeBtn8;
	attributeBtns[9] = oui->attributeBtn9;
	attributeBtns[10] = oui->attributeBtn10;

	updateLabels();
	
	oui->make_plot(p3dUICallback->getControl(), p3dUICallback);
	oui->plotGrp->add(oui->plot3DWrapper);
	oui->plotGrp->end();

	redrawChart();

	initialized = true;
}

// This function is called from OptVisualizerUI::show().
void OptVisualizerCallback::show()
{
	if (oui == NULL)
		return;

	oui->plotGrp->show();
	oui->chartGrp->show();
	oui->optVisualizerUIWindow->show();
}

void OptVisualizerCallback::categoryChange()
{
	category = oui->categoryChooser->value();
	oui->newChart->setCategory(category);
	redrawChart();
}

void OptVisualizerCallback::rangeChange()
{
	int min = oui->minIndexInput->value();
	int max = oui->maxIndexInput->value();
	int step = oui->stepSizeInput->value();

	oui->newChart->setDataRange(min, max);
	if (step < 1)
		step = 1;
	oui->newChart->setStepSize(step);
	oui->updateChart();
}

void OptVisualizerCallback::hide()
{
	oui->showModelBtn->value(0);
	showModelBtnPressed();
	oui->optVisualizerUIWindow->hide();
}

void OptVisualizerCallback::toggleAttribute(int attrBtn)
{
	oui->newChart->setAttributeIsVisible(attrBtn,
		(int) attributeBtns[attrBtn]->value());
	oui->updateChart();
}

void OptVisualizerCallback::selectionModeChange(int btnVal)
{
	LogManagerChart::SelectionType selection;

	switch (btnVal) {
		case 0:		selection = LogManagerChart::SelectCenter;
					break;
		case 1:		selection = LogManagerChart::SelectX1;
					break;
		case 2:		selection = LogManagerChart::SelectX2;
					break;
		case 3:		selection = LogManagerChart::SelectDiffX1;
					break;
		case 4:		selection = LogManagerChart::SelectDiffX2;
					break;
		case 5:		selection = LogManagerChart::SelectProjection;
					break;
		default:	selection = LogManagerChart::SelectUnknown;
					break;
	}
	oui->newChart->setSelectionType(selection);
//	std::cout << "Selection is: "  << selection << std::endl;
}

void OptVisualizerCallback::updateBtnPressed()
{
	const char * o = oui->originInput->value();
	double x1M = oui->x1MagnitudeInput->value();
	const char * x1 = oui->x1Input->value();
	double x2M = oui->x2MagnitudeInput->value();
	const char * x2 = oui->x2Input->value();
	const char * p = oui->projectionInput->value();
	oui->plot3DWrapper->callback_update(o, x1, x1M, x2, x2M, p);
}

void OptVisualizerCallback::resetBtnPressed()
{
	oui->plot3DWrapper->callback_reset();
}

void OptVisualizerCallback::scalingBtnPressed()
{
	oui->plot3DWrapper->callback_scaling(oui->scalingBtn->value());
}

void OptVisualizerCallback::showModelBtnPressed()
{
	if (! oui->plot3DWrapper->status()) {
		oui->showModelBtn->value(0);
		return;
	}

	oui->plot3DWrapper->callback_show_model(oui->showModelBtn->value());
}

void OptVisualizerCallback::surfaceChanged()
{
	oui->plot3DWrapper->callback_surface_choice(oui->surfaceChoice->value());
}

void OptVisualizerCallback::commandChanged()
{
	oui->plot3DWrapper->callback_selection_command(
		(Plot3DWindowWrapper::selection_t) oui->commandChoice->value());
	oui->commandChoice->value(0);
}

void OptVisualizerCallback::arrowBtnPressed(int direction)
{
	oui->plot3DWrapper->callback_arrows(direction);
}

void OptVisualizerCallback::projectionChanged()
{
	oui->plot3DWrapper->callback_projection();
}

void OptVisualizerCallback::updateLabels() {
	
	int i = 0;
	for (i = 0; i < MAX_ATTRIBUTES; i++ ){
		if (i < globalLogManager.numAttributes) {
			attributeBtns[i]->label(globalLogManager.attributeNames[i]);
			attributeBtns[i]->labelcolor(i);
		}
		else {
			attributeBtns[i]->label(0);
		}
		attributeBtns[i]->value(0);
		oui->newChart->setAttributeIsVisible(i, 0);
	}
	if (oui && oui->surfaceChoice) {
		oui->surfaceChoice->clear();
		for (i = 0; i < globalLogManager.numAttributes; i++) {
			oui->surfaceChoice->add(globalLogManager.attributeNames[i]);
		}
	}


}

#else	/* OPTIMIZATION_VISUALIZER */


void OptVisualizerCallback::categoryChange()
{
}

void OptVisualizerCallback::rangeChange()
{
}

void OptVisualizerCallback::hide()
{
}

void OptVisualizerCallback::toggleAttribute(int attrBtn)
{
}

void OptVisualizerCallback::selectionModeChange(int btnVal)
{
}

void OptVisualizerCallback::updateBtnPressed()
{
}

void OptVisualizerCallback::resetBtnPressed()
{
}

void OptVisualizerCallback::scalingBtnPressed()
{
}

void OptVisualizerCallback::showModelBtnPressed()
{
}

void OptVisualizerCallback::surfaceChanged()
{
}

void OptVisualizerCallback::commandChanged()
{
}

void OptVisualizerCallback::arrowBtnPressed(int direction)
{
}

void OptVisualizerCallback::projectionChanged()
{
}

void OptVisualizerCallback::updateLabels() {
}


#endif	/* OPTIMIZATION_VISUALIZER */


