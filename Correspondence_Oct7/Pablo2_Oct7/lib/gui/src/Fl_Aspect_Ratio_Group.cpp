
#include <iostream>
#include <FL/Fl_Group.H>
#include "Fl_Aspect_Ratio_Group.h"
#include "P3DUserInterface.h"

//#define DEBUG


// These constants must correspond with the definitions in P3DUserInterface.fl.
const int gapWidth = 5;

const int centerButtonWidth = 65;
const int smallButtonWidth = 30;
const int modelWindowButtonGrpHeight = 40;
const int modelWindowButtonGrpWidth = 605;

const int taskDisplayGrpY = 5;
const int taskDisplayGrpWidth = 595;
const int taskDisplayGrpHeight = 50;

extern P3DUserInterface * ui;
extern void triggerModelWindowRedraw();

// This is really the difference between screen height, Fl::h(), and the height
// FLTK uses, when it calls the function below to resize the modelDisplayGrp
// after the user has clicked the full sceen button in the upper right corner
// of the window.
const int windowBannerHeight = 30;

bool Fl_Aspect_Ratio_Group::fullscreen = false;
bool Fl_Aspect_Ratio_Group::last_status = false;


using namespace std;


/*  Make sure window always has the same aspect ratio.  This relies on
	the model display group being resized before the button group is
	resized.  Likewise, in the TASKING version, the taskDisplayGrp is
	expected to be drawn before the model display group.  Thus, the
	order of these groups in P3DUserInterface.fl must be taskDisplayGrp,
	then modelDisplayGrp, and finally modelWindowButtonGrp.

	The heights of the button group and of the task display group are
	held constant here.  The modelView window is kept square and made
	as large as possible, given these constraints and the window size.
*/
void Fl_Aspect_Ratio_Group::resize(int x, int y, int w, int h)
{
	if (ui == NULL) {
		Fl_Group::resize(x, y, w, h);
		return;
	}

#ifdef DEBUG
	cout << "\nFl_Aspect_Ratio_Group::resize():\n";
	cout << "    x = " << x << ", y = " << y << ", w = " << w << ", h = "
		<< h << endl; 
	cout << "    Fl: x = " << Fl::x() << ", y = " << Fl::y() << ", w = "
		<< Fl::w() << ", h = " << Fl::h() << endl;
	cout << "    modelWindow: w = " << ui->modelWindow->w() << ", h = "
		<< ui->modelWindow->h() << endl;
	cout << "    modelWindowButtonGrp: x = " << ui->modelWindowButtonGrp->x() << ", y = "
		<< ui->modelWindowButtonGrp->y() << ", w = " << ui->modelWindowButtonGrp->w()
		<< ", h = " << ui->modelWindowButtonGrp->h() << endl;
#endif

	if (this == ui->taskDisplayGrp)
	{
#ifdef TASKING
		h = taskDisplayGrpHeight;
		y = taskDisplayGrpY;
		if (w > taskDisplayGrpWidth)
			w = taskDisplayGrpWidth;
#ifdef DEBUG
	    cout << "taskDisplayGrp resize(x = " << x << ", y = " << y << ", w = "
		    << w << ", h = " << h << ')' << endl;
#endif
		Fl_Group::resize(x, y, w, h);
#endif	/* TASKING */
	}
	else if (this == ui->modelDisplayGrp)
	{
#ifdef DEBUG
		cout << "Resizing modelDisplayGrp\n";
#endif
		// In Fltk 1.1, there is no other way to tell if the fullscreen button in the
		// upper corner of the model window was selected.
		last_status = fullscreen;
		if (h >= Fl::h() - windowBannerHeight)
			fullscreen = true;
		else
			fullscreen = false;

		int windowHeight;
		if (fullscreen)
			windowHeight = Fl::h() - windowBannerHeight;
		else
			windowHeight = ui->modelWindow->h();

		// Based on the P3DUserInterface.fl file, the modelWindowButtonGrp has a height
		// of 40 pixels and the modelDisplayGrp (or the taskDisplayGrp) begins 5 pixels
		// below the base window.  If the taskDisplayGrp is present, there is a 10 pixel
		// gap below it.  This code adjusts the model display to fit the window, so that
		// the button group maintains its constant height and stays on the window.
		int hnew = windowHeight - modelWindowButtonGrpHeight - gapWidth;
#ifdef TASKING
		hnew -= taskDisplayGrpHeight + gapWidth + gapWidth;
#endif
#ifdef DEBUG
		cout << "Window height = " << windowHeight << ", hnew = " << hnew << ", w = " << w << endl;
#endif

		h = hnew;
		// The model display group must remain square
		if(w > h)
			w =  h;
		else
			h = w;

#ifdef DEBUG
	    cout << "modelDisplayGrp resize(x = " << x << ", y = " << y << ", w = "
		    << w << ", h = " << h << ')' << endl;
#endif
#ifdef TASKING
		// There is a 5 pixel gap above and and 10 pixel gap below the task window
		y = taskDisplayGrpHeight + 3*gapWidth;
#endif
		Fl_Group::resize(x, y, w, h);
	}
	else if (this == ui->modelWindowButtonGrp)
	{
		int ynew;

		//cout << "modelWindowButtonGrp: fullscreen = " << fullscreen << endl;
#ifdef DEBUG
		cout << "Resizing modelWindowButtonGrp\n";
#endif
		// Based on the P3DUserInterface.fl file, the modelWindowButtonGrp has a height
		// of 40 pixels and there is no gap immediately above or below it.  This fixes
		// the origin of the button group at ynew.  It relies on the model display group
		// being drawn first, which will be the case (see above), as long as the window
		// is not resized too rapidly by the user.
		ynew = ui->modelDisplayGrp->y() + ui->modelDisplayGrp->h();
#ifdef DEBUG
		cout << "ynew = " << ynew << endl;
#endif

		if (w > modelWindowButtonGrpWidth)
			w = modelWindowButtonGrpWidth;
#ifdef DEBUG
	    cout << "modelWindowButtonGrp resize(x = " << x << ", y = " << ynew << ", w = "
		    << w << ", h = " << modelWindowButtonGrpHeight << ')' << endl;
#endif
		Fl_Group::resize(x, ynew, w, modelWindowButtonGrpHeight);

		int x_start;
		int x_avail, dx_1, dx_2, dx_3;
		const int coronalViewBtnX = gapWidth + smallButtonWidth*2;

		// Adjust spacing of the buttons in modelWindowButtonGrp.
		// Constants are used for the button sizes, because FLTK tends
		// to shrink them when the parent window is reduced.  This is
		// also why resize() is called, instead of calling position().
		if (w < modelWindowButtonGrpWidth) {
			x_start = w - (centerButtonWidth + gapWidth);
			//cout << "x_start = " << x_start << endl;
			ui->centerBtn->resize(x_start, ui->centerBtn->y(), centerButtonWidth,
				ui->centerBtn->h());

			// Determine inter-button spacings, retaining original ratios.
			// The small buttons are assumed here to be 30 pixels wide.
			x_avail = x_start - coronalViewBtnX - 5*smallButtonWidth;
			//cout << "x_avail = " << x_avail << endl;
			dx_1 = (int) (x_avail*(125.0/320.0));
			dx_2 = (int) (x_avail*(25.0/320.0));
			dx_3 = (int) (x_avail*(20.0/320.0));

			x_start -= (dx_1 + 30);
			//cout << "buttons at: " << x_start;
			ui->rotateRight90->resize(x_start, ui->centerBtn->y(), smallButtonWidth,
				ui->centerBtn->h());
			x_start -= (dx_2 + 30);
			//cout << ", " << x_start;
			ui->rotateVertical180->resize(x_start, ui->centerBtn->y(), smallButtonWidth,
				ui->centerBtn->h());
			x_start -= (dx_3 + 30);
			//cout << ", " << x_start;
			ui->rotateHorizontal180->resize(x_start, ui->centerBtn->y(), smallButtonWidth,
				ui->centerBtn->h());
			x_start -= (dx_2 + 30);
			//cout << ", " << x_start << endl;
			ui->rotateLeft90->resize(x_start, ui->centerBtn->y(), smallButtonWidth,
				ui->centerBtn->h());

			x_start = coronalViewBtnX;
			ui->coronalViewBtn->resize(x_start, ui->centerBtn->y(), smallButtonWidth,
				ui->centerBtn->h());
			x_start -= smallButtonWidth;
			ui->saggitalViewBtn->resize(x_start, ui->centerBtn->y(), smallButtonWidth,
				ui->centerBtn->h());
			x_start = gapWidth;
			ui->axialViewBtn->resize(x_start, ui->centerBtn->y(), smallButtonWidth,
				ui->centerBtn->h());
		}

		// This should not be needed, and thus seems to imply there is a bug in
		// FLTK.  If it is omitted, changing the window's display mode to full
		// screen does not produce a corresponding change in the size of the
		// model display group.
		if (last_status != fullscreen)
				triggerModelWindowRedraw();
	}
#ifndef PRODUCTION_VERSION
	else
		cout << "Resized an unknown Fl_Aspect_Ratio_Group\n";
#endif
}


