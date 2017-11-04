#ifndef P3D_CUT_VIEW_H
#define P3D_CUT_VIEW_H


#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>

#include "Image3D.h"
#include "CutPlaneGlobals.h"


/*  A class for generating and rendering a cut-plane visualization
	in a window using OpenGl.

	The constructor takes as its arguments the coordinates of the origin,
	the width and height of the window, and an optional label.  The width
	and height should generally be integer powers of 2.  The visualization
	will always be square, so they should usually be the same.

	To produce a visualization, one of the init() functions must be used.
	These determine the atom to be at the center of the visualization.
	If an object is provided, the first "selected" atom found in it will
	be used.  Typically the second version of init() will be used and
	the currently "marked" atom will be specified.  After calling init(),
	status() will return false, if a problem was encountered.

    The CutPlanes object may be obtained by calling getCutViewPtr().  For
	views other than the B-N plane, it will be needed to specify the
	desired view type.  It may also be need to set the width and height of
	the visualization or for other purposes.

	The run() function, which must be used next, unless CutPlanes::run()
	is called, will generate the visualization by calling CutPlanes::run().
	The width and height used in this call will be that specified to the
	constructor, that is, the window's dimensions, if a different width
	and height is not specified as arguments to run().  If new dimensions
	are given, they will be used to generate the visualizatoin image, and
	it will be zoomed in or out to fit the window when that is shown.

    The draw() function is required by the base class.  It renders the
	image and drawn lines in the window.

	The remaining functions return information stored by the class:
	visualization() returns a pointer to the computed cut plane image;
	origin(), corner() and end() return the model-space coordinates of
	three corners of the visual image, if it has been computed.  If it
	has not yet been computed, they return (-1, -1, -1).

	The imageWidth() and imageHeight() functions return the dimensions
	of the visalization image, if known, but zero if it has not been
	computed.

*/

class M3DPrimitive;
class M3DObject;

class P3DCutPlaneView : public Fl_Gl_Window
{

	public:

		P3DCutPlaneView(int x, int y, int w, int h,
			const char * label = NULL);
		virtual ~P3DCutPlaneView();

		void init(const M3DObject * object, Image3D * image);
		void init(const M3DPrimitive * atom, const M3DObject * object,
			Image3D * image);
		bool status();

		void type(plane_t pt);
		void boundary(bool toggle);
		void boundaryWidth(int w) { bdry_width = w; }
		void gamma(float g);						// Obsolete
		void separate_windowing(bool state);
		void half_planes(bool state) { clip = state; };
		void bperp_line(bool state) { divider = state; };
		CutPlanes * getCutViewPtr();
		bool run(int width = 0, int height = 0);
		bool rerun(M3DObject * object, int width = 0, int height = 0);
		void draw();

		Image3D * visualization() const { return vis; }
		Point origin() const;
		Point corner() const;
		Point end() const;
		int imageWidth() const;
		int imageHeight() const;

		void paint_image(void * texture, bool * done,	// For rendering in other windows
			const Point * c0, const Point * c1,
			const Point * c2, const Point * c3, bool smooth, bool half_plane);


	private:

		CutPlanes * cut;
		Image3D * vis;
		GreyValue * im;
		double view_size_x, view_size_y;	// Gl parameters
		double min_x, min_y;				// Gl parameters
		bool done;
		bool clip;
		bool divider;
		unsigned int texture_name;
		int v_width, v_height;				// Dimensions of visualization
		bool drawing_bdry;		// True to draw the boundary in the cut plane view
		int bdry_width;

		void paint_image();
		void draw_boundary();
		void draw_lines();
		P3DCutPlaneView(const P3DCutPlaneView & cpv);	// Not implemented

};






#endif	/* P3D_CUT_VIEW_H */

