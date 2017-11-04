
#include <iostream>
#include <math.h>

#include "CutPlanes.h"
#include "M3DObject.h"
#include "ControlParms.h"
#include "P3DCutPlaneView.h"

#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>


//#define DEBUG



#ifdef DEBUG
#include "utility.h"
#endif

using namespace std;


P3DCutPlaneView::P3DCutPlaneView(int x, int y, int w, int h,
	const char * label) : Fl_Gl_Window(x, y, w, h, label)
{
#ifdef DEBUG
	cout << "P3DCutPlaneView: x = " << x << " y = " << y
		<< " w = " << w << " h = " << h << endl;
#endif
	cut = NULL;
	im = NULL;
	vis = NULL;
	texture_name = 0;
	done = false;
	clip = false;
	divider = false;
#ifdef DEBUG
	cout << "P3DCutPlaneView::P3DCutPlaneView(): done set to false" << endl;
#endif
}



P3DCutPlaneView::~P3DCutPlaneView() {
#ifdef DEBUG
	cout << "P3DCutPlaneView::~P3DCutPlaneView()" << endl;
#endif
	delete cut;
	// Vis belongs to the CutPlanes object, and so cannot be deleted
	// Im is a pointer to the pixels in vis.
	if (texture_name)
		glDeleteTextures (1, (GLuint*)&texture_name);
}



void P3DCutPlaneView::init(const M3DObject * object, Image3D * image) {
#ifdef DEBUG
	cout << "P3DCutPlaneView:init(object, image)" << endl;
#endif
	if (cut == NULL)
		cut = new CutPlanes(object, image);
	else
		cut->init(object, image);
	vis = NULL;	// Not deleted, because it belongs to cut.
	im = NULL;	// Im is a pointer to the pixels in vis.
	done = false;
#ifdef DEBUG
	cout << "P3DCutPlaneView::init(): done set to false" << endl;
#endif
}



void P3DCutPlaneView::init(const M3DPrimitive * atom, const M3DObject * object,
						   Image3D * image)
{
#ifdef DEBUG
	cout << "P3DCutPlaneView:init(atom, object, image)" << endl;
#endif
	if (cut == NULL)
		cut = new CutPlanes(atom, object, image);
	else
		cut->init(atom, object, image);
	vis = NULL;	// Not deleted, because it belongs to cut.
	im = NULL;	// Im is a pointer to the pixels in vis.
	done = false;
	clip = false;
	divider = false;
#ifdef DEBUG
	cout << "P3DCutPlaneView::init(): done set to false; atom = "
		<< hex << atom << dec << endl;
#endif
}



CutPlanes * P3DCutPlaneView::getCutViewPtr() {
	return cut;
}



bool P3DCutPlaneView::status() {
	if (cut == NULL) return false;
	return cut->valid();
}



void P3DCutPlaneView::type(plane_t pt = b_n) {
	if (cut == NULL) {
		cout << "P3DCutPlaneView::type() ignored\n";
		return;
	}
	cut->type(pt);
}



void P3DCutPlaneView::boundary(bool toggle) {
	drawing_bdry = toggle;
	if (cut == NULL) {
		cout << "P3DCutPlaneView::boundary() ignored\n";
		return;
	}
	cut->boundary(toggle);
}



void P3DCutPlaneView::gamma(float g) {
	if (cut == NULL) {
		cout << "P3DCutPlaneView::gamma() ignored\n";
		return;
	}
	cut->gamma(g);
}



void P3DCutPlaneView::separate_windowing(bool state) {
	if (cut == NULL) {
		cout << "P3DCutPlaneView::separate_windowing() ignored\n";
		return;
	}
	cut->separate_windowing(state);
}



bool P3DCutPlaneView::run(int width, int height) {

#ifdef DEBUG
	cout << "P3DCutPlaneView::run(" << width << ", " << height << ")\n";
#endif
	if (! status()) {
		cout << "Cannot produce cut-plane visualization\n";
		return false;
	}
	im = NULL;
	done = false;
	if (width == 0) v_width = w();
	else v_width = width;
	if (height == 0) v_height = h();
	else v_height = height;
	if (cut == NULL)
		return false;
	if (cut->run(v_width, v_height) == false)
		return false; 
	vis = cut->visualization();
	if (vis == NULL)
		return false;
	im = vis->getVoxels();
	return true;
}



bool P3DCutPlaneView::rerun(M3DObject * object, int width, int height) {

#ifdef DEBUG
	cout << "P3DCutPlaneView::run(" << width << ", " << height << ")\n";
#endif
	if (! status()) {
		cout << "Cannot produce cut-plane visualization\n";
		return false;
	}
	im = NULL;
	done = false;
	if (width == 0) v_width = w();
	else v_width = width;
	if (height == 0) v_height = h();
	else v_height = height;
	if (cut == NULL)
		return false;
	M3DFigure * figure = (M3DFigure *) cut->getRenderedFigure();
	if (object->getFigureIndex(figure) < 0) {
		delete cut;
		cut = NULL;
		return false;
	}
	if (cut->run(v_width, v_height) == false)
		return false; 
	vis = cut->visualization();
	if (vis == NULL)
		return false;
	im = vis->getVoxels();
	return true;
}



Point P3DCutPlaneView::origin() const {
	if (cut == NULL)
		return Point(-1.0, -1.0, -1.0);
	return cut->view_origin();
}



Point P3DCutPlaneView::corner() const {
	if (cut == NULL)
		return Point(-1.0, -1.0, -1.0);
	return cut->view_corner();
}



Point P3DCutPlaneView::end() const {
	if (cut == NULL)
		return Point(-1.0, -1.0, -1.0);
	return cut->view_end();
}



int P3DCutPlaneView::imageWidth() const {
	if (cut == NULL)
		return 0;
	return cut->width();
}



int P3DCutPlaneView::imageHeight() const {
	if (cut == NULL)
		return 0;
	return cut->height();
}



void P3DCutPlaneView::draw() {

#ifdef DEBUG
	cout << "P3DCutPlaneView::draw():  valid == " << (valid() ? "true" : "false")
		<< endl;
#endif

	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);

	if (! done) {
		if (cut == NULL) return;
		vis = cut->visualization();
		if (vis == NULL) return;

		im = vis->getVoxels();
		v_width = vis->getXDim();
		v_height = vis->getYDim();

		// Find the limiting dimension of the window and determine how much should
		// be visible in each dimension.
		if (w() < h()) {
			// The window's width is the limiting dimension.
			view_size_x = 1.0;
			view_size_y = h()/((double) w());
		}
		else {
			// The window's height is the limiting dimension.
			view_size_x = w()/((double) h());
			view_size_y = 1.0;
		}

		min_x = (1.0 - view_size_x)/2.0;
		min_y = (1.0 - view_size_y)/2.0;
	}

#ifdef DEBUG
	report_intensity_range(*vis, "P3DCutPlaneView::vis");
	cout << "P3DCutPlaneView::draw():  vis = " << hex << vis << dec << endl;
#endif

	// Set up the Gl viewport
    glViewport(0, 0, w(), h()); // left X, bottom Y, etc.
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    glOrtho(min_x, view_size_x + min_x, min_y, view_size_y + min_y, -1, 1);

	paint_image();
	if (drawing_bdry)
		draw_boundary();
	draw_lines();
    glFlush ();
	glPopAttrib();
}



void P3DCutPlaneView::paint_image()
{

#ifdef DEBUG
	cout << "paint_image(): Drawing " << v_width << " x " << v_height << " image." << endl;
#endif
//snapshot(*vis, "snap.raw", true);

	// Clear the rendering window.
    glClear(GL_COLOR_BUFFER_BIT);

#ifdef DEBUG
	cout << "paint_image(): texture name = " << texture_name << endl;
#endif

	const Point c0(0.0, 0.0, 0.0);
	const Point c1(1.0, 0.0, 0.0);
	const Point c2(1.0, 1.0, 0.0);
	const Point c3(0.0, 1.0, 0.0);

	glRasterPos2i(0, 0);

	bool smooth = globalControl->readBool(SmoothImages);
	paint_image((void *) &texture_name, &done, &c0, &c1, &c2, &c3, smooth, false);
}



void P3DCutPlaneView::paint_image(void * texture, bool *painted, const Point * c0,
								  const Point * c1, const Point * c2, const Point * c3,
								  bool smooth, bool half_plane)
{
#ifdef DEBUG
	cout << "paint_image(texture, ...)" << endl;
#endif
	if (im == NULL) return;

	if (*(GLuint *) texture == 0)
		glGenTextures(1, (GLuint *) texture);

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
	glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
	glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);

	glBindTexture(GL_TEXTURE_2D, *(GLuint *) texture);
#ifdef DEBUG
	cout << "paint_image(): Bind texture\n";
	cout << "paint_image(): painted = " << (int) painted << '\n';
#endif

	if (! *painted) {
#ifdef DEBUG
		cout << "paint_image(): Regenerating texture" << endl;
#endif

		if (smooth) {
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		}
		else {
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		}

/*
		glPixelTransferf(GL_RED_SCALE, cut->scale());
		glPixelTransferf(GL_RED_BIAS, cut->bias());
		glPixelTransferf(GL_GREEN_SCALE, cut->scale());
		glPixelTransferf(GL_GREEN_BIAS, cut->bias());
		glPixelTransferf(GL_BLUE_SCALE, cut->scale());
		glPixelTransferf(GL_BLUE_BIAS, cut->bias());
#ifdef DEBUG
		cout << "paint_image(): Scale = " << cut->scale() << "   bias = "
			<< cut->bias() << endl;
#endif
*/

		glTexImage2D(GL_TEXTURE_2D, 0, 1, v_width, v_height,
			0, GL_LUMINANCE, GL_UNSIGNED_SHORT, (GLvoid *) im);
#ifdef DEBUG
		cout << "paint_image(): glTexImage2D() called with im = "
			<< hex << im << dec << "  width = " << v_width << "  height = "
			<< v_height << endl;
#endif
/*
		glPixelTransferf(GL_RED_SCALE, 1.0);
		glPixelTransferf(GL_RED_BIAS, 0.0);
		glPixelTransferf(GL_GREEN_SCALE, 1.0);
		glPixelTransferf(GL_GREEN_BIAS, 0.0);
		glPixelTransferf(GL_BLUE_SCALE, 1.0);
		glPixelTransferf(GL_BLUE_BIAS, 0.0);
*/
		*painted = true;
#ifdef DEBUG
		cout << "paint_image(): done set to true" << endl;
#endif
	}

	glEnable(GL_TEXTURE_2D);

#ifdef DEBUG
	cout << "paint_image(): texture name = " << *(GLuint *) texture << endl;
#endif

    // Draw a quad onto which the texture will be mapped
    glBegin(GL_QUADS);
	glColor3f(1.0f, 1.0f, 1.0f);

	if (half_plane) {
		glTexCoord2f(0.0f, 0.5f);    glVertex3d(c0->getX(), c0->getY(), c0->getZ());
		glTexCoord2f(1.0f, 0.5f);    glVertex3d(c1->getX(), c1->getY(), c1->getZ());
	}
	else {
		glTexCoord2f(0.0f, 0.0f);    glVertex3d(c0->getX(), c0->getY(), c0->getZ());
		glTexCoord2f(1.0f, 0.0f);    glVertex3d(c1->getX(), c1->getY(), c1->getZ());
	}
	glTexCoord2f(1.0f, 1.0f);    glVertex3d(c2->getX(), c2->getY(), c2->getZ());
	glTexCoord2f(0.0f, 1.0f);    glVertex3d(c3->getX(), c3->getY(), c3->getZ());
    glEnd();

	glDisable(GL_TEXTURE_2D);

    return;
}



void P3DCutPlaneView::draw_lines() {

#ifdef DEBUG
	cout << "P3DCutPlaneView::draw_lines()" << endl;
#endif

	Point line_end_0 = cut->end0();
	Point line_end_1 = cut->end1();
	Point line_midpt = cut->midpoint();
	Point line_end_y0 = cut->endY0();
	Point line_end_y1 = cut->endY1();
	double b_ext = cut->b_extension();

	double scale = w()/v_width;		// AGG: what if window is not square?

#ifdef DEBUG
	cout << "window dimensions:  w = " << w() << " h = " << h() << '\n';
	cout << "scale = " << scale << '\n';
	switch (cut->view_type()) {
		case b_n:		cout << "\nB-N View:\n";
						break;
		case b_bperp:	cout << "\nB-Bperp View:\n";
						break;
		case bperp_n_chord:
						cout << "\nBperp-Chord View:\n";
						break;
		case bperp_y0:	cout << "\nBperp-Y0 View:\n";
						break;
		case bperp_y1:	cout << "\nBperp-Y1 View:\n";
						break;
		case bperp_n:	cout << "\nBperp-N View:\n";
						break;
	}
	cout << "\tEnd0: "; line_end_0.print();
	cout << "\tEnd1: "; line_end_1.print();
	cout << "\tMidpt: "; line_midpt.print();
	cout << "\tEndY0: "; line_end_y0.print();
	cout << "\tEndY1: "; line_end_y1.print();
	cout << "\tB-ext: "; cout << b_ext << endl;
#endif

	line_end_0.setX(scale*line_end_0.getX());
	line_end_0.setY(scale*line_end_0.getY());
	line_end_1.setX(scale*line_end_1.getX());
	line_end_1.setY(scale*line_end_1.getY());
	line_midpt.setX(scale*line_midpt.getX());
	line_midpt.setY(scale*line_midpt.getY());
	line_end_y0.setX(scale*line_end_y0.getX());
	line_end_y0.setY(scale*line_end_y0.getY());
	line_end_y1.setX(scale*line_end_y1.getX());
	line_end_y1.setY(scale*line_end_y1.getY());
	b_ext *= scale;

#ifdef DEBUG
	cout << "\tScaled end0: "; line_end_0.print();
	cout << "\tScaled end1: "; line_end_1.print();
	cout << "\tScaled midpt: "; line_midpt.print();
	cout << "\tScaled endY0: "; line_end_y0.print();
	cout << "\tScaled endY1: "; line_end_y1.print();
	cout << "\tScaled B-ext: "; cout << b_ext << endl;
#endif

	// Draw the Y vectors: y0 in magenta, y1 in cyan
	Point a, b;
	const Vector3D u(3.0, 0.0, 0.0);	// A short line segment in the X direction
	const Vector3D v(0.0, 3.0, 0.0);	// A short line segment in the Y direction

	ortho();			// The lines are drawn in pixel coordinates
	glBegin(GL_LINES);
	glColor3d(1.0, 0.0, 1.0);	// magenta
	switch (cut->view_type()) {
		case b_n:		glVertex2d(line_end_y0.getX(), line_end_y0.getY());
						glVertex2d(line_midpt.getX(), line_midpt.getY());

						glColor3d(0.0, 1.0, 1.0);	// cyan
						glVertex2d(line_end_y1.getX(), line_end_y1.getY());
						glVertex2d(line_midpt.getX(), line_midpt.getY());
						break;

		case b_bperp:	a = line_end_y0 + u;
						b = line_end_y0 - u;
						glVertex2d(a.getX(), a.getY());
						glVertex2d(b.getX(), b.getY());
						break;

		case bperp_n:	a = line_end_y0 + u;
						b = line_end_y0 - u;
						glVertex2d(a.getX(), a.getY());
						glVertex2d(b.getX(), b.getY());

						glColor3d(0.0, 1.0, 1.0);	// cyan
						a = line_end_y1 + u;
						b = line_end_y1 - u;
						glVertex2d(a.getX(), a.getY());
						glVertex2d(b.getX(), b.getY());
						a = line_end_y1 + v;
						b = line_end_y1 - v;
						glVertex2d(a.getX(), a.getY());
						glVertex2d(b.getX(), b.getY());
						break;

		case bperp_n_chord:
						if (cut->theta_is_90()) {
							// Draw the Y vectors
							glVertex2d(line_end_y0.getX(), line_end_y0.getY());
							glVertex2d(line_midpt.getX(), line_midpt.getY());

							glColor3d(0.0, 1.0, 1.0);	// cyan
							glVertex2d(line_midpt.getX(), line_midpt.getY());
							glVertex2d(line_end_y1.getX(), line_end_y1.getY());
						}
						else {
							// Draw plus signs to indicate the boundary points
							a = line_end_y0 + u;
							b = line_end_y0 - u;
							glVertex2d(a.getX(), a.getY());
							glVertex2d(b.getX(), b.getY());
							a = line_end_y0 + v;
							b = line_end_y0 - v;
							glVertex2d(a.getX(), a.getY());
							glVertex2d(b.getX(), b.getY());

							glColor3d(0.0, 1.0, 1.0);	// cyan
							a = line_end_y1 + u;
							b = line_end_y1 - u;
							glVertex2d(a.getX(), a.getY());
							glVertex2d(b.getX(), b.getY());
							a = line_end_y1 + v;
							b = line_end_y1 - v;
							glVertex2d(a.getX(), a.getY());
							glVertex2d(b.getX(), b.getY());
						}
						break;

		case bperp_y0:	glColor3d(0.0, 1.0, 1.0);	// cyan
						glVertex2d(line_end_y0.getX(), line_end_y0.getY());
						glVertex2d(line_midpt.getX(), line_midpt.getY());

						glColor3d(1.0, 0.0, 1.0);	// magenta
						if (cut->theta_is_90()) {
							glVertex2d(line_end_y1.getX(), line_end_y1.getY());
							glVertex2d(line_midpt.getX(), line_midpt.getY());
						}
						else {
							a = line_end_y1 + u;
							b = line_end_y1 - u;
							glVertex2d(a.getX(), a.getY());
							glVertex2d(b.getX(), b.getY());

							a = line_end_y1 + v;
							b = line_end_y1 - v;
							glVertex2d(a.getX(), a.getY());
							glVertex2d(b.getX(), b.getY());
						}
						break;

		case bperp_y1:	glColor3d(0.0, 1.0, 1.0);	// cyan
						if (cut->theta_is_90()) {
							glVertex2d(line_end_y0.getX(), line_end_y0.getY());
							glVertex2d(line_midpt.getX(), line_midpt.getY());
						}
						else {
							a = line_end_y0 + u;
							b = line_end_y0 - u;
							glVertex2d(a.getX(), a.getY());
							glVertex2d(b.getX(), b.getY());

							a = line_end_y0 + v;
							b = line_end_y0 - v;
							glVertex2d(a.getX(), a.getY());
							glVertex2d(b.getX(), b.getY());
						}

						glColor3d(1.0, 0.0, 1.0);	// magenta
						glVertex2d(line_end_y1.getX(), line_end_y1.getY());
						glVertex2d(line_midpt.getX(), line_midpt.getY());
						break;
	}

	// Draw the vertical axis
	switch (cut->view_type()) {
		case b_n:
		case b_bperp:	glColor3d(1.0, 0.0, 0.0);		// B is red
						glVertex2d(line_end_0.getX(), line_end_0.getY());
						glVertex2d(line_midpt.getX(), line_midpt.getY());
						break;
		case bperp_n:	glColor3d(1.0, 1.0, 1.0);		// N is white
						glVertex2d(line_end_0.getX(), line_end_0.getY());
						glVertex2d(line_midpt.getX(), line_midpt.getY());
						break;
		case bperp_n_chord:
		case bperp_y0:
		case bperp_y1:
						break;
	}

	if (b_ext > 0.0f) {	// This is zero except for the B-N and B-Bperp views
		// Draw a dashed line with 3-pixel dashes separated by 4-pixel gaps
		double y0 = line_end_0.getY() + 4.0;	// End with a gap
		double y1 = b_ext + line_midpt.getY();
#ifdef DEBUG
		cout << "midpt: "; line_midpt.print();
		cout << "end of B = " << y0 - 4.0;
		cout << "   b_ext = " << b_ext << "\ndashed line from y0 = " << y0 << " to " << y1 << endl;
#endif
		while (y1 > y0) {	// Start at end of extended vector
			glVertex2d(line_end_0.getX(), y1);
			glVertex2d(line_end_0.getX(), y1 - 3.0);
			y1 -= 3.0 + 4.0;
		}
	}

	// Draw the horizontal axis
	if (cut->view_type() == b_n)
		glColor3d(1.0, 1.0, 1.0);		// N is white
	else
		glColor3d(0.0, 0.0, 1.0);		// Bperp is blue

	switch (cut->view_type()) {
		case b_n:
		case b_bperp:
		case bperp_n:	
					glVertex2d(line_end_1.getX(), line_end_1.getY());
					glVertex2d(line_midpt.getX(), line_midpt.getY());
					break;
		case bperp_n_chord:
					if (divider) {
						// Draw a dashed line on the window indicating the half-plane edge
						double x0 = 0.0;	// Start with a dash
						double x1 = line_midpt.getX() + line_midpt.getX();
						while (x0 < x1) {
							glVertex2d(x0, line_midpt.getY());
							glVertex2d(x0 + 3.0, line_midpt.getY());
							x0 += 3.0 + 4.0;
						}
					}
					else {
						// Draw a dashed line with 3-pixel dashes separated by 4-pixel gaps
						double x0 = line_end_1.getX();	// Start with a dash
						double x1 = line_midpt.getX();
						while (x0 < x1) {
							glVertex2d(x0, line_midpt.getY());
							glVertex2d(x0 + 3.0, line_midpt.getY());
							x0 += 3.0 + 4.0;
						}
					}
						break;
		case bperp_y0:
		case bperp_y1:
					if (divider) {
						// Draw a line across the window indicating the half-plane edge
						glVertex2d(0.0, line_midpt.getY());
						glVertex2d(line_midpt.getX() + line_midpt.getX(), line_midpt.getY());
					}
					else {
						glVertex2d(line_end_1.getX(), line_end_1.getY());
						glVertex2d(line_midpt.getX(), line_midpt.getY());
					}
					break;
	}
/*
//	glEnable(GL_LIGHTING);
	cout << "drawing sphere" << endl;
	GLUquadric * q = gluNewQuadric();
	gluDisk(q, 0.0, 1.0, 8, 8);
	gluDeleteQuadric(q);
//	glDisable(GL_LIGHTING);
*/
    glEnd();
}



void P3DCutPlaneView::draw_boundary() {

#ifdef DEBUG
	cout << "P3DCutPlaneView::draw_boundary()" << endl;
#endif

	Point v0;
	Point v1;

	double scale = w()/v_width;		// AGG: what if window is not square?

#ifdef DEBUG
	cout << "window dimensions:  w = " << w() << " h = " << h() << '\n';
	cout << "scale for drawing boundary = " << scale << '\n';
#endif

	ortho();			// The lines are drawn in pixel coordinates  ???????????

	int npoints = cut->boundary_point_count();
	glLineWidth((float) bdry_width);
	glBegin(GL_LINES);
	glColor3fv(cut->boundary_color());		// Set boundary color

#ifdef DEBUG
	cout << "Boundary lists " << npoints << " verticies\n";
#endif
	for (int i = 0; i < npoints; i += 2) {
		v0 = scale*cut->boundary_point(i);
		v1 = scale*cut->boundary_point(i + 1);
		glVertex2d(v0.getX(), v0.getY());
		glVertex2d(v1.getX(), v1.getY());
	}
    glEnd();
	glLineWidth(1.0f);
}


