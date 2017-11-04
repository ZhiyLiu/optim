#ifndef CUT_PLANES_H
#define CUT_PLANES_H


#include <vector>
#include "CutPlaneGlobals.h"




/*	Class for generating a cut-plane visualization.

	The constructors and init() functions collectively take as input
	an image from which the visualization will be produced and a single
	atom, which determines the center of the visualization.  If an object
	is provided, the first "selected" atom found in it will be used.
	After construction, valid() will return false, if a problem was
	encountered.

	After construction, the type() function may be called to specify the
	plane to be viewed.  If not, the B-N plane will be assumed.  Function
	size() may optionally be used to determine the area to be visualized.
	It's model-space argument should be a factor of r; a default of 3r is
	used if this function is not called or has no argument.  On some screens
	it will be desirable to call gamma() with a value in the range [0, 1].
	which will propotionally brighten the display.  (This does not really
	implement gamma correction, but has a similar effect).

	Then the visualization may be produced by calling run(), with arguments
	indicating the size of the desired display (Gl viewport) in pixels.  The
	result will be generation of the desired image slice, which can be gotten
	by calling visualization().  The dimensions of the image produced can be
	obtained by calling width() and height().

	The coordinates of lines and points that can be drawn on the image are
	also produced.  Functions end0(), end1() and midpt() return the primary
	axes, with the midpoint always being the center of the viewport.  The
	meaning of these axes depends on the view chosen, as follows:
		View:	b_n	  b_bperp   bperp_n	  bperp_y0	 bperp_y1   bperp_n_chord
		end_0:	B	  B			N		  y0		 y1			N (offset)
		end_1:	N	  B-perp	B-perp	  B-perp	 B-perp		B-perp (offset)
	Function b_extension() indicates how long the B vector should be
	extended for end atoms.  Function theta_is_90() will return true if the
	boundary points are on opposite sides of the atom center.  This information
	may be used to decide whether to draw the "Y" vectors as lines or just to
	show their ends.  Finally, functions scale() and bias() return values to be
	use in the glPixelTransfer() functions (e.g. GL_RED_SCALE and GL_RED_BIAS),
	before the image is displayed using glDrawPixels().  Function view_type()
	simply returns the view specification previously set by the user.

	This class may also be used to produce a set of lines (curve) to outline the
	object boundary in the model window.  To do this, the default constructor
	followed by init2() should be used; then the run() with 6 arguments may be
	called to generate the boundary intersection for a particular figure of the
	model.  The level argument of init2() is used the same way as the lvl
	argument of class M3DObjectSurfaceRenderer, so the boundary produced will
	correspond to that seen when rendering the full surface.

    A copy constructor is provided so vectors of CutPlanes objects can be
	constructed and resized.  However, the copy constructor should only be called
	when the object is in its default state, because no attempt is made to copy
	the images, point lists, etc. that it contains.

*/

class M3DPrimitive;
class M3DFigure;
class M3DObject;
class Image3D;

class CutPlanes {

	public:

		CutPlanes();
		CutPlanes(const M3DObject * object, Image3D * image);
		CutPlanes(const M3DPrimitive * atom, const M3DObject * object,
			Image3D * image);
		CutPlanes(const CutPlanes & cp);
		virtual ~CutPlanes();
		bool init(const M3DObject * object, Image3D * image);
		bool init(const M3DPrimitive * prim, const M3DObject * object,
			Image3D * image);
		void init2(const M3DObject * object, int level = 1);
		bool valid() { return status; }
		void reset();

		void type(plane_t pt = b_n) { plane = pt; }
		void boundary(bool toggle) { drawing_bdry = toggle; }
		void size(double s = 0.0) { min_size = s; }
		void gamma(float g);
		void separate_windowing(bool state) { adjust_window = state; }
		bool run(int width = 0, int height = 0);
		void run(int figureId, Point origin, Point corner, Point end,
			int width = 0, int height = 0);
		void setLevel(int level);

		Image3D * visualization() { return vis_image; }
		int width();
		int height();
		plane_t view_type() { return plane; }
		int boundary_point_count() { return bdry.size(); }
		Point boundary_point(int i) { return bdry[i]; }
		float * boundary_color() { return bdry_color; }

		// Functions for rendering the cut view in display space
		Point end0() { return end_0; }
		Point end1() { return end_1; }
		Point midpoint() { return vertex; }
		Point endY0() { return end_y0; }
		Point endY1() { return end_y1; }
		double b_extension() { return eta; }
		bool theta_is_90() { return opposing_y; }
		float scale() { return windowing_scale; }
		float bias() { return offset; }

		// Functions for rendering the cut view in model space
		Point view_origin() { return origin; }
		Point view_corner() { return corner; }
		Point view_end() { return end; }
		Point view_end0() { return mdl_end_0; }
		Point view_end1() { return mdl_end_1; }
		Point view_midpoint() { return mdl_vertex; }
		Point view_endY0() { return mdl_end_y0; }
		Point view_endY1() { return mdl_end_y1; }
		double view_b_extension() { return mdl_eta; }

		const M3DFigure * getRenderedFigure() const { return figure; }


	private:

		M3DPrimitive * atom;	// Atom being visualized
		const M3DFigure * figure;		// Figure containing the atom
		const M3DObject * obj;		// Object containing the atom
		Image3D *input_image;	// Image to be sampled
		Image3D *vis_image;		// Resultant visualization
		double line_scale;
		double min_size;		// Min. width/height in model space of cut patch
		double msize;
		plane_t plane;			// Visualization type (cut identifier)
		Point origin;			// Model-space corners of the image patch
		Point corner;
		Point end;
		Vector3D normal;		// Normal of display plane in model space
		bool drawing_bdry;		// True to draw the boundary in the cut plane view
		void * plist;

		// Model space values
		Point mdl_end_0;		// Points defining the axial vectors to be drawn
		Point mdl_end_1;
		Point mdl_vertex;
		Point mdl_end_y0;		// Points defining the Y vectors
		Point mdl_end_y1;
		double mdl_eta;			// End-atom extension

		// Display space values
		Point end_0;			// Points defining the axial vectors to be drawn
		Point end_1;
		Point vertex;
		Point end_y0;			// Points defining the Y vectors
		Point end_y1;
		double eta;				// End-atom extension

		bool opposing_y;
		bool status;
		bool adjust_window;
		float windowing_scale;	// Windowing parameters for OpenGl
		float offset;
		float gfact;
		float bdry_color[3];
		int lvl;				// Used for drawing in the model window

		std::vector<Point> bdry;

		void gen_cut_view(double dist);
		void define_atom_cut(int max_width, int max_height);
		void define_lines(int dim);
		void initialize_boundary(bool standard);
		void generate_boundary();
		void convert_boundary(int dim, Point displayCenter);
		// The new parameter tilesetindex coloring the contours
		// see the comment near "ifdef COLOR_CONTOURS" 
		void clip_boundary(int numtiles, void * tileset, int * tilesetindex = 0);
		int intersect_edge_with_plane(Point & pt0, Point & pt1,
			Point & intersection);
		// The new parameters idx0, idx1, idx2 are related to coloring the contours
		// see the comment near "ifdef COLOR_CONTOURS" 
		void clip_triangle(double * vtx0, double * vtx1, double * vtx2, int idx0 = 0, int idx1 = 0, int idx2 = 0);

#ifdef COLOR_CONTOURS
 		// When COLOR_CONTOURS is enabled, I want to use colors from the 
		// ObjectSurfaceVisualizer to color the vertices of the contours
		// keep track of which bpoint contributes to each tile
		// so that the visualizer can recover them

		std::vector<int> bpointIndices;

		public:
			int bpoint_index(int i) { return bpointIndices[i]; }
#endif
};



#endif	/* CUT_PLANES_H */

