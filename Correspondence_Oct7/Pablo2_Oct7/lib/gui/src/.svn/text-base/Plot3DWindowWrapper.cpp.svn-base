#include "LogManager.h"
#ifdef OPTIMIZATION_VISUALIZER

#include <iostream>
#include <strstream>
#include <FL/Fl.H>
#ifdef WIN32
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include <FL/fl_draw.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_File_Chooser.H>
#include "FunctionExplorer.h"
#include "M3DObjectFile.h"
#include "P3DUserInterfaceCallback.h"
#include "Plot3DWindowWrapper.h"



using namespace std;





///////////////////////////////////////////////////////////////////
// iostream utilities
///////////////////////////////////////////////////////////////////

std::ostream & operator << (std::ostream & out, const Vector & v) {
	int n = v.size();
	out << "<" << v(0);
	for (int i = 1; i < n; i++) {
		out << ", " << v(i);
	}
	out << ">";
	return out;
}

std::istream & operator >> (std::istream & in, Vector & v) {
	int n = v.size();
	v.setAll(0);
	int i = 0;
	char peek;
	double temp;
	while ((i < n) && ! in.eof()) {
		peek = in.peek();
		switch(peek) {
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
		case '.':
		case '-': {
			in >> temp;
			v(i) = temp;
			i++;
			break;
							}
		default :
			in.ignore();
		}
	}

	std::cout << "Read vector: " << v << std::endl;
	return in;
}




///////////////////////////////////////////////////////////////////
// Plot3DWindow
///////////////////////////////////////////////////////////////////

Plot3DWindow::Plot3DWindow(int x, int y, int w, int h, const char * label)
	: Fl_Gl_Window(x, y, w, h, label), track(1.0, 1.0)
{
	initialized = false;
	attribute = 2;
	theta = 0;
	phi = M_PI/2.0;
	zoom = 4.0;
	negating = false;
	scaling = false;
	selectionMode = SelectedProjection;
	control = NULL;
	centerX = centerY = centerZ = 0.0;
	projectionX = 0.0;
	projectionY = 0.0;
	projection = NULL;
	toProject = NULL;
	track.setControl(NULL);

	transformMatrix[0][0] = 1.0f;
	transformMatrix[0][1] = 0.0f;
	transformMatrix[0][2] = 0.0f;
	transformMatrix[0][3] = 0.0f;
	transformMatrix[1][0] = 0.0f;
	transformMatrix[1][1] = 1.0f;
	transformMatrix[1][2] = 0.0f;
	transformMatrix[1][3] = 0.0f;
	transformMatrix[2][0] = 0.0f;
	transformMatrix[2][1] = 0.0f;
	transformMatrix[2][2] = 1.0f;
	transformMatrix[2][3] = 0.0f;
	transformMatrix[3][0] = -0.5f;
	transformMatrix[3][1] = -0.5f;
	transformMatrix[3][2] = -0.5f;
	transformMatrix[3][3] = 1.0f;

	track.reset(&(transformMatrix[0][0]));

	controlModel = NULL;
	controlIsCompromised = false;
}

Plot3DWindow::~Plot3DWindow()
{
	// Restore the object in P3DControl
	if (controlIsCompromised)
		restoreModel();
}

M3DObject * Plot3DWindow::getTargetObject(Vector & x)
{
	OptimizerBase * opt;

	opt = globalLogManager.getOptimizer();
	if (opt)
		return opt->createTargetObject(x);
	else
		return NULL;
}

M3DObject * Plot3DWindow::getCurrentModel()
{
	Vector * v;

	v = getSelection();
	return getTargetObject(*v);
}

void Plot3DWindow::init()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	gluPerspective(90, 1, 0.1, 5);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
	glEnable(GL_DEPTH_TEST);

	initialized = true;

/*
	// mock up some data to show in the plot
	globalLogManager.clearCategory(2);
	globalLogManager.setCategory(2);
	for (int x = -10; x<=10; x++) {
		for (int y = -10; y<=10; y++) {
			globalLogManager.beginEvent();
			globalLogManager.setValue(2, -(0.0075 * x * x) + -(0.0075 * y * y));
		}
	}
*/
}

void Plot3DWindow::draw()
{
	bool needsCentering = ! initialized;
	if (! initialized)
		init();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, w(), h());
    glOrtho(-1.0, 1.0, -1.0, 1.0, -500.0, 500.0);

	glClearColor(0.75, 0.75, 0.75, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);

	glPushMatrix();
//		glOrtho(0, 1, 0, 1, 0, 1);
	glLoadIdentity();

	glDisable(GL_LIGHTING);

	track.doViewTransform();

//	drawCube();

	drawAxis();
	drawData();

	glPopMatrix();
	glEnable(GL_LIGHTING);

	if (needsCentering)
		axisToCenter();
}

void Plot3DWindow::drawCube()
{
	glBegin(GL_LINES);
	glColor3f(1, 1, 0);
	glVertex3f(-0.5, -0.5, -0.5);
	glVertex3f(0.5, -0.5, -0.5);
	glVertex3f(0.5, -0.5, -0.5);
	glVertex3f(0.5, 0.5, -0.5);
	glVertex3f(0.5, 0.5, -0.5);
	glVertex3f(-0.5, 0.5, -0.5);
	glVertex3f(-0.5, 0.5, -0.5);
	glVertex3f(-0.5, -0.5, -0.5);
	glColor3f(1, 0, 1);
	glVertex3f(-0.5, -0.5, 0.5);
	glVertex3f(0.5, -0.5, 0.5);
	glVertex3f(0.5, -0.5, 0.5);
	glVertex3f(0.5, 0.5, 0.5);
	glVertex3f(0.5, 0.5, 0.5);
	glVertex3f(-0.5, 0.5, 0.5);
	glVertex3f(-0.5, 0.5, 0.5);
	glVertex3f(-0.5, -0.5, 0.5);
	glColor3f(1, 0, 0);
	glVertex3f(0.5, -0.5, 0.5);
	glVertex3f(0.5, -0.5, -0.5);
	glVertex3f(0.5, 0.5, 0.5);
	glVertex3f(0.5, 0.5, -0.5);
	glVertex3f(-0.5, 0.5, 0.5);
	glVertex3f(-0.5, 0.5, -0.5);
	glVertex3f(-0.5, -0.5, 0.5);
	glVertex3f(-0.5, -0.5, -0.5);
	glEnd();
}

void Plot3DWindow::drawAxis()
{
	glLineWidth(3);
	glBegin(GL_LINES);
	glColor3f(1, 0, 0);
	glVertex2i(0, 0);
	glVertex2i(1, 0);
	glColor3f(1, 0, 1);
	glVertex2i(0, 0);
	glVertex2i(0, 1);
	glColor3f(0, 0, 1);
	glVertex2i(0, 0);
	glVertex3i(0, 0, 1);
	glEnd();
	glLineWidth(1);
}

void Plot3DWindow::drawData()
{
	if (globalLogManager.numEventsPerCategory[2] == 441) {

		glBegin(GL_QUADS);

		int idx = 0;
		for (int x = -10; x < 10; x++) {
			for (int y = -10; y < 10; y++) {
				double dx = 0.1 * x;
				double dy = 0.1 * y;
				colorScale(globalLogManager.eventsByCategory[2][idx][attribute]);

				glVertex3d(dx, dy,
					scale(globalLogManager.eventsByCategory[2][idx][attribute]));
				glVertex3d(dx + 0.1, dy,
					scale(globalLogManager.eventsByCategory[2][idx+21][attribute]));

				glVertex3d(dx + 0.1, dy + 0.1,
					scale(globalLogManager.eventsByCategory[2][idx+22][attribute]));
				glVertex3d(dx, dy + 0.1,
					scale(globalLogManager.eventsByCategory[2][idx+1][attribute]));
				idx++;
			}
			idx++;
		}
		glEnd();

		drawSelection();
	}
}

// Mark off the special positions... where we start, and where the neighbor is
void Plot3DWindow::drawSelection()
{
	float oldWidth; 
	glGetFloatv(GL_LINE_WIDTH, &oldWidth);

	glColor3d(0.39, 0.2, 0.48);

	// Draw a line from the point to be projected onto its projection
	glLineWidth(5);
	if (toProject && projection && globalLogManager.getOptimizer()
		&& globalLogManager.getOptimizer()->getProblem())
	{
		globalLogManager.clearCategory(3);
		globalLogManager.setCategory(3); // don't want to log these values
		globalLogManager.getOptimizer()->getProblem()->evaluate(*toProject);
		globalLogManager.getOptimizer()->getProblem()->evaluate(*projection);


		double trueValue = scale(globalLogManager.eventsByCategory[3][0][attribute]);
		double projectedValue = scale(globalLogManager.eventsByCategory[3][1][attribute]);
		if (fabs(trueValue - projectedValue) > 0.1) {
			glBegin(GL_LINES);
			glVertex3d(projectionX, projectionY, trueValue);
			glVertex3d(projectionX, projectionY, projectedValue);
			glEnd();
		}
		else {
			glPointSize(8);
			glBegin(GL_POINTS);
			glVertex3d(projectionX, projectionY, trueValue);
			glEnd();
		}
	}

	if (selectionMode == SelectedGrid) {
		glColor3d(0.0, 0.5, 0.8);
		glPointSize(8);
		double dx = 0.1 * (selectionX - 10);
		double dy = 0.1 * (selectionY - 10);
		double z = scale(
			globalLogManager.eventsByCategory[2][selectionX * 21 + selectionY][attribute]);
		glBegin(GL_POINTS);
		glVertex3d(dx, dy, z);
		glEnd();
	}

	glLineWidth(oldWidth);
}

void Plot3DWindow::setAttribute(int attr)
{
	attribute = attr % globalLogManager.numAttributes;
	while (attribute < 0) {
		attribute += globalLogManager.numAttributes;
	}
	cout << "Now surface is for [" << globalLogManager.attributeNames[attribute] << "]" << endl;
	damage(0xFF);
}

void Plot3DWindow::dump()
{
	Vector * o = (globalLogManager.allParameters[(2 * MAX_EVENTS) + 21 * 10 + 10]);
	Vector * oX = (globalLogManager.allParameters[(2 * MAX_EVENTS) + 21 * 11 + 10]);
	Vector * oY = (globalLogManager.allParameters[(2 * MAX_EVENTS) + 21 * 10 + 11]);
	if (o && oX && oY) {

		cout << "#Origin: " << *o << endl;
		cout << "#Dir 1: " << 10 * ((*oX) - (*o)) << endl;
		cout << "#Dir 2: " << 10 * ((*oY) - (*o)) << endl;
		cout << "#Surface: " << globalLogManager.attributeNames[attribute] << endl;

		for (int x = -10; x <= 10; x++) {
			for (int y = -10; y <= 10; y++) {
				int idx = (21 * (x + 10)) + (y + 10);
				cout << (1.0 * x)/10 << " " << (1.0 * y)/10 << " " << 
					globalLogManager.eventsByCategory[2][idx][attribute] << endl;
			}
		}
	}		
}

void Plot3DWindow::dataChanged(bool forgetProjection)
{
	if (globalLogManager.numEventsPerCategory[2] == 441) {
		max = globalLogManager.eventsByCategory[2][0][attribute];
		min = max;
		for (int i = 1; i < 441; i++) {
			if (globalLogManager.eventsByCategory[2][i][attribute] < min) {
				 min = globalLogManager.eventsByCategory[2][i][attribute];
			}
			else if (globalLogManager.eventsByCategory[2][i][attribute] > max) {
				 max = globalLogManager.eventsByCategory[2][i][attribute];
			}
		}
	}
	else {
		min = 0;
		max = 0;
	}
	if (forgetProjection) {
		calculateProjectionCoefficients();
	}
}

// Called when a arrow button is pressed
void Plot3DWindow::selectionPosition(int direction)
{	
	if (selectionMode == SelectedProjection) {
		double dx = ((projectionX < -1) ? -1 : ((projectionX > 1) ? 1 : projectionX));
		double dy = ((projectionY < -1) ? -1 : ((projectionY > 1) ? 1 : projectionY));
		selectionX = 10 +(dx * 10);
		selectionY = 10 +(dy * 10);
		selectionMode = SelectedGrid;
	}

	switch (direction) {
		case 0:		// up
			selectionX = (21 + selectionX + 1) % 21;
			break;
		case 1:		// down
			selectionX = (21 + selectionX - 1) % 21;
			break;
		case 2:		// left
			selectionY = (21 + selectionY + 1) % 21;
			break;
		case 3:		// right
			selectionY = (21 + selectionY - 1) % 21;
			break;
		default:
			return;
	}

	damage(0xFF);

	if (controlIsCompromised) {
		M3DObject * o = getCurrentModel();
		control->replaceObjectPtr(o);
	}
}

void Plot3DWindow::colorScale(double val)
{
	if (min == max)
		glColor3d(1.0, 1.0, 0.0);
	else {
		double cumulative = (val - min) / (max - min);
		double red, green, blue;
		if (cumulative < 1.0/3) {
			red = cumulative * 3.0;
			green = 0.0;
			blue = 0.0;
		}
		else {
			red = 1.0;
			if (cumulative < 2.0/3.0) {
				green = 3.0 * cumulative - 1.0;
				blue = 0;
			}
			else {
				green = 1.0;
				blue = 3.0 * cumulative - 2.0;
			}
		}
		glColor3d(red, green, blue);
	}
}

void Plot3DWindow::reseedOptimizer()
{
	Vector * v;

	v = getSelection();
	if (v && globalLogManager.getOptimizer())
		globalLogManager.getOptimizer()->setOptimizerPosition(*v);
}

void Plot3DWindow::setProjectedVector(const Vector & v)
{
	if (toProject)
		delete toProject;
	toProject = new Vector(v);

	if (globalLogManager.numEventsPerCategory[2] == 441)
		calculateProjectionCoefficients();
	if (controlIsCompromised) {
		M3DObject * o = getCurrentModel();
		control->replaceObjectPtr(o);
	}
}

void Plot3DWindow::calculateProjectionCoefficients()
{
	if (toProject) {
		int index = 10 * 21 + 10;
		Vector * origin =
			globalLogManager.allParameters[2 * MAX_EVENTS + 10 * 21 + 10]; // 2 is category
		if (origin == NULL)
			return;
		Vector x1 = *globalLogManager.allParameters[2 * MAX_EVENTS + 20 * 21 + 10] - *origin;
		Vector x2 = *globalLogManager.allParameters[2 * MAX_EVENTS + 10 * 21 + 20] - *origin;
		Vector diff = (*toProject) - *origin;
		projectionX = diff.dotProduct(x1) / x1.dotProduct(x1);
		projectionY = diff.dotProduct(x2) / x2.dotProduct(x2);
		double temp = x1.dotProduct(x2);

		if (projection)
			delete projection;
		projection = new Vector(*origin + projectionX * x1 + projectionY * x2);
		selectionMode = SelectedProjection;
	}
}

Vector * Plot3DWindow::getSelection()
{
	if (selectionMode == SelectedProjection)
		return projection;
	return globalLogManager.allParameters[(2 * MAX_EVENTS) + (selectionX * 21) + selectionY];
}

// Handle mouse events
int Plot3DWindow::handle(int eventId)
{
	int button;
    MouseEventStruct mouse;

    switch(eventId)
    {
        case FL_PUSH:
            button = Fl::event_button();
            mouse.button = button;
            mouse.width = w();
            mouse.height = h();
            mouse.type = MOUSE_PUSH;
            mouse.x = Fl::event_x();
            mouse.y = Fl::event_y();
            mouse.viewType = PERSPECTIVE;

            mouse.modifiers = NO_KEY;
            if(Fl::event_state(FL_CTRL))
                mouse.modifiers |= CTRL_KEY;
            if(Fl::event_state(FL_SHIFT))
                mouse.modifiers |= SHIFT_KEY;
            if(Fl::event_state(FL_ALT))
                mouse.modifiers |= ALT_KEY;

            track.mouseEvent(mouse);
			return 1;

		case FL_DRAG:
            if(track.getMode() != NONE)
            {
                mouse.button = Fl::event_button();
                mouse.width = w();
                mouse.height = h();
                mouse.type = MOUSE_DRAG;
                mouse.x = Fl::event_x();
                mouse.y = Fl::event_y();
                mouse.viewType = PERSPECTIVE;

                mouse.modifiers = NO_KEY;
                if(Fl::event_state(FL_CTRL))
                    mouse.modifiers |= CTRL_KEY;
                if(Fl::event_state(FL_SHIFT))
                    mouse.modifiers |= SHIFT_KEY;
                if(Fl::event_state(FL_ALT))
                    mouse.modifiers |= ALT_KEY;

                track.mouseEvent(mouse);
            }
			damage(0xFF);
            return 1;

        case FL_RELEASE:
            button = Fl::event_button();
			if(track.getMode() != NONE)
            {
                mouse.button = button;
                mouse.width = w();
                mouse.height = h();
                mouse.type = MOUSE_RELEASE;
                mouse.x = Fl::event_x();
                mouse.y = Fl::event_y();
                mouse.viewType = PERSPECTIVE;

                mouse.modifiers = NO_KEY;
                if(Fl::event_state(FL_CTRL))
                    mouse.modifiers |= CTRL_KEY;
                if(Fl::event_state(FL_SHIFT))
                    mouse.modifiers |= SHIFT_KEY;
                if(Fl::event_state(FL_ALT))
                    mouse.modifiers |= ALT_KEY;
                if(Fl::event_state(FL_META))
                    mouse.modifiers |= META_KEY;

                track.mouseEvent(mouse);
	        }
			damage(0xFF);
            break;
		}

	    return Fl_Gl_Window::handle(eventId);
}

void Plot3DWindow::axisToCenter()
{
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	double c[3];

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);

	double centerX = viewport[2]/2.0;
	double centerY = viewport[3]/2.0;
	centerY = viewport[3] - centerY - 1.0;
	gluUnProject(centerX, centerY, 0.5, modelview, projection, viewport, &c[0], &c[1], &c[2]);
	track.translateView(c[0], c[1], c[2]);

	redraw();
}

void Plot3DWindow::showingModel(bool val) {
	if (controlIsCompromised && ! val)
		restoreModel();
	else if (val && ! controlIsCompromised)
		displayNewModel();
	controlIsCompromised = val;
}

void Plot3DWindow::displayNewModel() {
	if (controlModel != NULL) {
		cout << "Invalid call to Plot3DWindow::displayNewModel()\n";
		return;
	}
	M3DObject * o = getCurrentModel();
	controlModel = control->replaceObjectPtr(o);
}

void Plot3DWindow::restoreModel() {
	// If (controlModel == NULL) then the image was discarded by the user.
	// This situation should not happen.
	M3DObject * o = control->replaceObjectPtr(controlModel);
	delete o;
	controlModel = NULL;
}



///////////////////////////////////////////////////////////////////
// Plot3DWindowWrapper
///////////////////////////////////////////////////////////////////

const char msg[] = "Error: initialize origin vector to specify the proper dimension";

Plot3DWindowWrapper::Plot3DWindowWrapper(int x, int y, int w, int h, 
	P3DControl * control) : Fl_Window(x, y, w, h, NULL)
{
	threed = new Plot3DWindow(0, 0, w, h, NULL);
	threed->mode(FL_RGB8);
	threed->setControl(control);
	threed->set_visible_focus();
	threed->box(FL_FLAT_BOX);
	threed->end();

	pCenter = NULL;
	pX1 = NULL;
	pX2 = NULL;		
	end();
}

Plot3DWindowWrapper::~Plot3DWindowWrapper()
{
	delete threed;
	delete pCenter;
	delete pX1;
	delete pX2;
}

void Plot3DWindowWrapper::setWidgets(Fl_Input * center, Fl_Input * x1,
	Fl_Value_Input * x1Magnitude, Fl_Input * x2, Fl_Value_Input * x2Magnitude,
	Fl_Input * projection)
{
	centerInput = center;
	x1Input = x1;
	x1MagnitudeInput = x1Magnitude;
	x2Input = x2;
	x2MagnitudeInput = x2Magnitude;
	projectionInput = projection;
}

void Plot3DWindowWrapper::show()
{
	Fl_Window::show();
	if (threed)
		threed->show();
}

void Plot3DWindowWrapper::damage(uchar c)
{
	Fl_Window::damage(c);
	threed->damage(c);
}

void Plot3DWindowWrapper::setCenterVector(Vector * v)
{
	if (pCenter)
		delete pCenter;
	pCenter = (v ? new Vector(*v) : (Vector *) NULL);
	displayVector(centerInput, pCenter);
}

void Plot3DWindowWrapper::setX1Vector(const Vector * v)
{
	if (pX1)
		delete pX1;
	pX1 = (v ? new Vector(*v) : (Vector *) NULL);
	x1MagnitudeInput->value(pX1 ? pX1->twoNorm() : 1);
	displayVector(x1Input, pX1, true);
}

void Plot3DWindowWrapper::setX2Vector(const Vector * v)
{
	if (pX2)
		delete pX2;
	pX2 = (v ? new Vector(*v) : (Vector *) NULL);
	x2MagnitudeInput->value(pX2 ? pX2->twoNorm() : 1);
	displayVector(x2Input, pX2, true);
}

void Plot3DWindowWrapper::correctX2Vector()
{
	if (pX1 && pX2) {
		Vector unitX1(*pX1);
		unitX1.normalize();
		double temp = unitX1.dotProduct(unitX1);
		double x2ProjX1 = pX2->dotProduct(unitX1);
		Vector inLine = x2ProjX1 * (unitX1);
		Vector ortho = *pX2 - inLine;
		temp = inLine.dotProduct(ortho);
		setX2Vector(ortho);
	}
}

void Plot3DWindowWrapper::setPointVector(Vector * v)
{
	displayVector(projectionInput, v);
	projectionInput->do_callback();
}

void Plot3DWindowWrapper::displayVector(Fl_Input * o, Vector * vv, bool normalize)
{
	if (! o || ! vv)
		return; // Need a valid output
	if (vv) {
		Vector v(*vv);	// Copy the vector
		if (normalize)
			v.normalize();
		char * str = v;
		o->value(str);
	}
	else {
		o->value("");
	}
}

void Plot3DWindowWrapper::callback_scaling(int val)
{
	threed->enableScaling((val == 0) ? false : true);
	damage(0xFF);
}

void Plot3DWindowWrapper::callback_arrows(int direction)
{
	threed->selectionPosition(direction);
	updateP3D();
	damage(0xFF);
}

void Plot3DWindowWrapper::callback_surface_choice(int choice)
{
	threed->setAttribute(choice);	
	threed->dataChanged(false);
	damage(0xFF);
}

void Plot3DWindowWrapper::callback_reset()
{
	threed->axisToCenter();
	damage(0xFF);
}

// Called when the update button is pressed
void Plot3DWindowWrapper::callback_update(const char * origin, const char * x1,
	double x1Magnitude, const char * x2, double x2Magnitude, const char * projection)
{
	if (pCenter == NULL) {
		cout << msg << endl;
		return;
	}
	if (! (globalLogManager.getOptimizer() && globalLogManager.getOptimizer()->getProblem())) {
		cout << "Error: no optimizer problem" << endl;
		return;
	}

	int size = pCenter->size();
	Vector v(size);

	istrstream in0(origin);
	in0 >> v;
	setCenterVector(&v);

	istrstream in1(x1);
	in1 >> v;
	v *= x1Magnitude;
	setX1Vector(v);

	istrstream in2(x2);
	in2 >> v;
	v *= x2Magnitude;
	setX2Vector(v);

	istrstream in3(projection);
	in3 >> v;
	threed->setProjectedVector(v);

	correctX2Vector();

	FunctionExplorer::explorePlane(globalLogManager.getOptimizer()->getProblem(),
		*pCenter, *pX1, *pX2);
	threed->dataChanged();
	damage(0xFF);
}

void Plot3DWindowWrapper::callback_projection()
{
	if (! pCenter) {
		cout << msg << endl;
		return;
	}
	int size = pCenter->size();
	Vector v(size);
	istrstream in(projectionInput->value());
	in >> v;
	threed->setProjectedVector(v);
	updateP3D();
	damage(0xFF);
}

void Plot3DWindowWrapper::callback_selection_command(selection_t choice)
{
	M3DObject * obj;

	switch (choice) {
		case (SelectionNothing):
			break;
		case (SelectionToCenter):
			setCenterVector(threed->getSelection());
			damage(0xFF);
			break;
		case (SelectionExportModel):
			obj = threed->getCurrentModel();
			if (obj) {
				const char * filename = uic->askSingleFilename(ModelDirectory,
					"Save Visualizer Model", "*.m3d");
				if (filename) {
					M3DObjectFile out;
					out.write(filename, *obj);
				}
				delete obj;
			}
			break;
		case (SelectionReseed):
			threed->reseedOptimizer();
			break;
		case (SelectionExportData):
			threed->dump();
			break;
		case (SelectionLoadCenter): {
				const char* filename = uic->askSingleFilename(ModelDirectory, "Load Model As Origin", "*.m3d");
				if (filename) {
					M3DObjectFile in;
					M3DObject* target = in.read(filename);
					if (target && globalLogManager.getOptimizer()) {
						setCenterVector(globalLogManager.getOptimizer()->projectModel(target));
					}
					delete target;
				}
				break;
			}
		case (SelectionLoadX1): {
				const char* filename = uic->askSingleFilename(ModelDirectory, "Load Model As X1", "*.m3d");
				if (filename) {
					M3DObjectFile in;
					M3DObject* target = in.read(filename);
					if (target && globalLogManager.getOptimizer()) {
						setDiffX1Vector(globalLogManager.getOptimizer()->projectModel(target));
					}
					delete target;
				}
				break;
			}				
		case (SelectionLoadX2): {
				const char* filename = uic->askSingleFilename(ModelDirectory, "Load Model As X2", "*.m3d");
				if (filename) {
					M3DObjectFile in;
					M3DObject* target = in.read(filename);
					if (target && globalLogManager.getOptimizer()) {
						setDiffX2Vector(globalLogManager.getOptimizer()->projectModel(target));
					}
					delete target;
				}
				break;
			}				

			case (SelectionLoadProj): {
				const char* filename = uic->askSingleFilename(ModelDirectory, "Load Model As Projection", "*.m3d");
				if (filename) {
					M3DObjectFile in;
					M3DObject* target = in.read(filename);
					if (target && globalLogManager.getOptimizer()) {
						setPointVector(globalLogManager.getOptimizer()->projectModel(target));
					}
					delete target;
				}
				break;
			}		
	}
}

bool Plot3DWindowWrapper::status() const
{
	if (! threed || globalLogManager.numTotalEvents == 0)
		return false;
	return true;
}

void Plot3DWindowWrapper::callback_show_model(int val)
{
	threed->showingModel((val == 0) ? false : true);
	updateP3D();
}

void Plot3DWindowWrapper::updateP3D()
{
	if (uic)
		uic->update();
}




#endif /* OPTIMIZATION_VISUALIZER */


