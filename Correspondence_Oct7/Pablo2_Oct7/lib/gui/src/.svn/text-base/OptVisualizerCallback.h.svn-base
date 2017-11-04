#ifndef _OBJECTIVE_FUNCTION_UI_CALLBACK_H_
#define _OBJECTIVE_FUNCTION_UI_CALLBACK_H_



class Fl_Button;
class OptVisualizerUI;
class P3DUserInterface;
class P3DUserInterfaceCallback;

class OptVisualizerCallback
{

public:
	OptVisualizerCallback(OptVisualizerUI * ui) {
		oui = ui;
		initialized = false;
		attributeBtns = NULL;
	}
	void init();
	~OptVisualizerCallback();

	void categoryChange();

	void redrawChart();

	void rangeChange();
	void toggleAttribute(int attrBtn);
	void selectionModeChange(int btnVal);

	void updateBtnPressed();
	void resetBtnPressed();
	void scalingBtnPressed();
	void showModelBtnPressed();

	void surfaceChanged();
	void commandChanged();

	void projectionChanged();

	void arrowBtnPressed(int direction);

	void updateLabels();

	void show();
	void hide();

private:
	P3DUserInterface * p3dUI;
	P3DUserInterfaceCallback * p3dUICallback;
	OptVisualizerUI * oui;
	int category;
	bool initialized;
	Fl_Button ** attributeBtns;
};




#endif

