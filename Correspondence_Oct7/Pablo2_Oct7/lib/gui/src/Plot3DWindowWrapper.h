#ifndef _PLOT3D_WINDOW_WRAPPER_H_
#define _PLOT3D_WINDOW_WRAPPER_H_

#include "globalBuildParms.h"
#ifdef OPTIMIZATION_VISUALIZER

/**
 * This file defines classes needed for the objective function
 * visualizer.
 *
 * Plot3DWindowWrapper has the widgets for controlling the 3D view
 * Plot3DWindow is the actual surface plot
 */

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Gl_Window.H>
#include "matrix.h"
#include "Trackball.h"



class M3DObject;
class P3DControl;
class P3DUserInterfaceCallback;
class Plot3DWindowWrapper;
class Fl_Button;
class Fl_Choice;
class Fl_Input;
class Fl_Value_Input;



///////////////////////////////////////////////////////////////////
// Plot3DWindow
///////////////////////////////////////////////////////////////////


class Plot3DWindow : public Fl_Gl_Window
{

public:
	Plot3DWindow(int x, int y, int w, int h, const char * label);
	~Plot3DWindow();

	void setControl(P3DControl * c) { control = c; }
	virtual void draw();
	virtual int handle(int eventID);
	void enableScaling(bool value) { scaling = value; }
	void enableNegating(bool value) { negating = value; }

	void selectionPosition(int direction);
	void setAttribute(int attr);
	void dataChanged(bool forgetProjection = true); //notify the plot that data is new
	void setProjectedVector(const Vector & v);
	void reseedOptimizer();

	// todo: get selected or projected vector accordingly
	// this will let me export to a model or recenter the plot
	Vector * getSelection();
	void axisToCenter();
	// Output data for import into gnuplot or other programs
	void dump();

	void showingModel(bool val);
	void displayNewModel();
	void restoreModel();

	M3DObject * getCurrentModel();
	M3DObject * getTargetObject(Vector & x);

protected:

	void init();
	void drawCube();	// Debugging 
	void drawData();	// Draw the mesh surface
	void drawSelection();	// Draw markers for the two selected points
	void drawAxis();	// Draw the x,y,z axis
/*
	void adjustAttribute(int offset); // change which surface we are drawing, 
																		// offset should be +1 or -1
	void lookAt(int x, int y); // look at (x, y, f(x,y)) -- x,y are indices
	void lookAt(double x, double y, double z); // look at (x, y, z); -- x, y are coordinates
*/
	void calculateProjectionCoefficients(); // Map the projection vector into the space
											// of the current data

private:

	// is openGL initialized
	bool initialized;

	// Which attribute from the event log should I draw
	int attribute;

	// Control the camera angle
	double theta, phi;
	double zoom;

	// Input: are we zooming in or out?
	bool negating;

	// Want to be able to scale from [a,b] to [-1, 1]
	bool scaling;
	double min;
	double max;

	double scale(double val) {
		return scaling ? (2.0 * (val - min) / (max - min)) - 1 : val;
	}

	void colorScale(double val);

	// Be able to select a point
	// There will be one selection point
	// Initially it will correspond to the projection point
	// Moving it will move it along the grid
	// Projecting a new point will reset the selection to the projection point
	enum { SelectedProjection, SelectedGrid } selectionMode;
	int selectionX; // 
	int selectionY;

	// The projection of any point into the plane
	double projectionX;
	double projectionY;
	Vector * toProject;
	Vector * projection;

	// For access to the model slideshow
	P3DControl * control;

	// Control the view point
	double centerX, centerY, centerZ;

	// Trackball for adjusting screen
	Trackball track;
	float transformMatrix[4][4];

	// So that I can swap models with control
	M3DObject * controlModel;
	bool controlIsCompromised;
};


///////////////////////////////////////////////////////////////////
// Plot3DWindowWrapper
///////////////////////////////////////////////////////////////////

class Plot3DWindowWrapper : public Fl_Window
{
public:
	enum selection_t { SelectionNothing, SelectionToCenter, SelectionReseed,
		SelectionExportModel, SelectionExportData, SelectionLoadCenter, SelectionLoadX1, SelectionLoadX2, SelectionLoadProj};

	Plot3DWindowWrapper(int x, int y, int w, int h, P3DControl * control);
	virtual ~Plot3DWindowWrapper();
	bool status() const;

	void setControl(P3DControl * control) { threed->setControl(control); }
	void setWidgets(Fl_Input * center, Fl_Input * x1, Fl_Value_Input * x1Magnitude,
		Fl_Input * x2, Fl_Value_Input * x2Magnitude, Fl_Input * projection);

	virtual void show();
	virtual void damage(uchar c);

	// Related to the vectors that define the frame
	void setCenterVector(Vector * v);

	void setX1Vector(const Vector * v);
	void setX1Vector(const Vector & v) { setX1Vector(&v); }
	void setDiffX1Vector(Vector * v) { if (pCenter && v) setX1Vector((*v) - (*pCenter));}	

	void setX2Vector(const Vector * v);
	void setX2Vector(const Vector & v) { setX2Vector(&v); }
	void setDiffX2Vector(Vector * v) { if (pCenter && v) setX2Vector((*v) - (*pCenter));}

	void setPointVector(Vector * v);

	void displayVector(Fl_Input * o, Vector * v, bool normalize = false);
	void correctX2Vector();

	void setP3DUIC(P3DUserInterfaceCallback * p3duic) { uic = p3duic; }

	void updateP3D();

	void callback_update(const char * origin, const char * x1, double x1Magnitude,
		const char * x2, double x2Magnitude, const char * projection);
	void callback_reset();
	void callback_scaling(int val);
	void callback_show_model(int val);

	void callback_surface_choice(int choice);
	void callback_arrows(int direction);
	void callback_selection_command(selection_t choice);

	void callback_projection();

private:

	Vector * pCenter;
	Vector * pX1;
	Vector * pX2;

	// Pointers to widgets owned by OptVisualizerUI
	Fl_Input * centerInput;
	Fl_Input * x1Input;
	Fl_Value_Input * x1MagnitudeInput;
	Fl_Input * x2Input;
	Fl_Value_Input * x2MagnitudeInput;
	Fl_Input * projectionInput;

	// Actual 3d plotting widget
	Plot3DWindow * threed;

	// P3DUserInterfaceCalback so we can update the display
	// when showing models
	P3DUserInterfaceCallback * uic;

//	M3DObject * getTargetObject(Vector & x);

};


#endif	/* OPTIMIZATION_VISUALIZER */

#endif	/* _PLOT3D_WINDOW_WRAPPER_H_ */

