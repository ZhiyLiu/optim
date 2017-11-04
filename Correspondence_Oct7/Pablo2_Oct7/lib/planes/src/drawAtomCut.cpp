
#include <iostream>
#include "CutPlanes.h"
#include "M3DObject.h"
#include "Image3D.h"

#ifdef DEBUG
#include "utility.h"
#endif


//#define DEBUG


using namespace std;



/*
	Function to complete production of the cut-plane image.  This is called
	from gen_cut_view().

	Variables origin, corner, and end must contain coordinates of the cut's
	square bounding box in model space; origin and end are opposite points.
	The arguments are the dimensions of the window to contain the image.
*/
void CutPlanes::define_atom_cut(int max_width, int max_height)
{
	int i, j;
	int dim, width;
	float zoom;
	Vector3D u, v;
	Point io, ie, ic;


#ifdef DEBUG
	cout << "\nCutPlanes::define_atom_cut(): initial bounding box:\n";
	cout << "\torigin = ";
	origin.print();
	cout << "\tend = ";
	end.print();
	cout << "\tcorner = ";
	corner.print();
	verify_squareness(origin, corner, end, "bounding box");

	u = corner - origin;
	v = end - corner;
	cout << "|corner - origin| = " << u.norm();
	cout << "   |end - corner| = " << v.norm() << endl;
	if (abs(u.norm() - v.norm()) > 1e-10)
		cout << "These lengths should be equal and are not!\n";

	io = origin;
	input_image->modelToWorldCoordinates(io);		// Image-space origin
	ie = end;
	input_image->modelToWorldCoordinates(ie);		// Image-space end
	ic = corner;
	input_image->modelToWorldCoordinates(ic);		// Image-space corner
	u = ic - io;
	v = ie - ic;
	cout << "|corner - origin| = " << u.norm() << " cm";
	cout << "   |end - corner| = " << v.norm() << " cm" << endl;
	if (abs(u.norm() - v.norm()) > 1e-10)
		cout << "These lengths should be equal and are not!\n";

	cout << "Cut plane is ";
	switch (plane) {
		case b_n:		cout << "B-N\n";
						break;

		case b_bperp:	cout << "B-Bperp\n";
						break;

		case bperp_n:	cout << "Bperp-N\n";
						break;

		case bperp_n_chord:
						cout << "Bperp-Nchord\n";
						break;

		case bperp_y0:	cout << "Bperp-Y0\n";
						break;

		case bperp_y1:	cout << "Bperp-Y1\n";
						break;
	}
#endif

	// Define display space, a mapping of the viewport to pixel space
	io = origin;
	input_image->modelToImageCoordinates(io);		// Image-space origin
	ie = end;
	input_image->modelToImageCoordinates(ie);		// Image-space end
	ic = corner;
	input_image->modelToImageCoordinates(ic);		// Image-space corner

#ifdef DEBUG
	cout << "Cut view image origin is at pixel "; io.print();
	cout << "Cut view image end is at pixel "; ie.print();
	cout << "Cut view image corner is at pixel "; ic.print();
	u = ic - io;
	v = ie - ic;
	cout << "|ic - io| = " << u.norm();
	cout << "   |ie - ic| = " << v.norm() << endl;
#endif

	// Compute the number of pixels per direction in the output image.
	// Note that the sides of the cut usually will not parallel the X
	// and Y axes.
	v = ic - io;
	i = (int) (0.5 + v.norm());
	v = ie - ic;
	j = (int) (0.5 + v.norm());
	if (i < j)
		dim = j;
	else
		dim = i;

	// Zoom in model space to fit the display window
	if (max_width == 0 || max_height == 0) {
		zoom = 1.0;
		width = dim;
	}
	else {
		if (max_width < max_height) {
			zoom = max_height/(float) dim;
			width = max_height;
		}
		else {
			zoom = max_width/(float) dim;
			width = max_width;
		}
	}

#ifdef DEBUG
	cout << "Zoom factor = " << zoom << '\n';
#endif

/*
	if (zoom > 1) {
		cout << "Zoom factor clipped to 1\n";
		zoom = 1;
	}
*/
	dim = (int) (zoom*(float) dim);

	// Force the dimensions to be odd for better drawing of the central lines
	if (2.0*((float) (dim/2)) == dim) {
		if (dim < width)
			dim++;
		else
			dim--;
	}

	// Compute differential vectors in model space, for scan conversion
	// of the cut patch.  Traversal will be from the origin toward the end.
	// Note that the model-space viewport is required to be square.  Model
	// space is a regular space, so dx and dy will have the same length.
	v = corner - origin;
	Vector3D dx = (1.0/(dim - 1))*v;
	v = end - corner;
	Vector3D dy = (1.0/(dim - 1))*v;
	line_scale = dim/v.norm();		// Used for line drawing below

#ifdef DEBUG
	cout << "Viewport is " << dim << " by " << dim << " pixels\n";
	cout << "Line_scale is " << line_scale << '\n';
	cout << "Max width = " << max_width << "   max_height = " << max_height << endl;

	cout << "Differentials for scan conversion: \n\t("
		<< dx.getX() << ", " << dx.getY() << ", " << dx.getZ() << "),\n\t("
		<< dy.getX() << ", " << dy.getY() << ", " << dy.getZ() << ")\n";
#endif

	// Adjust dim to the GL window size, now that the differentials are known
	if (dim < width) dim++;

	int npixels = dim*dim;
#ifdef DEBUG
	cout << "npixels = " << npixels << '\n';
#endif

	// Produce the 2D output image
	vis_image = new Image3D;
	GreyValue * pixels = new GreyValue[npixels];
	vis_image->setVoxels(pixels, dim, dim, 1);

	Point p, q;
	GreyValue intens;

	// Determine the intensities and their range.  Note that the
	// image is inverted here (the Y axis is reversed), because the
	// OpenGL origin is at the bottom of the viewport.
	q = origin;
	q -= dx;
	q -= dy;
	GreyValue min = MAX_GREY_VALUE;
	GreyValue max = MIN_GREY_VALUE;
	if (adjust_window) {
		for (j = dim - 1; j >= 0; j--) {
			q += dy;
			p = q;
			for (i = 0; i < dim; i++) {
				p += dx;
				Vector3D x = p;
				input_image->modelToImageCoordinates(x);	// Convert x
				intens = (GreyValue) (0.5 + input_image->getInterpolatedVoxelValue(
					x.getX(), x.getY(), x.getZ()));
				pixels[i + j*dim] = intens;
				if (min > intens)
					min = intens;
				else if (max < intens)
					max = intens;
			}
		}
	}
	else {
		for (j = dim - 1; j >= 0; j--) {
			q += dy;
			p = q;
			for (i = 0; i < dim; i++) {
				p += dx;
				Vector3D x = p;
				input_image->modelToImageCoordinates(x);	// Convert x
				intens = (GreyValue) (0.5 + input_image->getWindowedInterpolatedVoxelValue(
					x.getX(), x.getY(), x.getZ()));
				pixels[i + j*dim] = intens;
				if (min > intens)
					min = intens;
				else if (max < intens)
					max = intens;
			}
		}
	}
/*
	// For debugging: draw a 3 pixel wide border of medium intensity on two edges
	intens = (min + max)/2;
	for (i = 0; i < dim; i++) {
		pixels[i + 0*dim] = intens;
		pixels[i + 1*dim] = intens;
		pixels[i + 2*dim] = intens;
	}
	for (j = 0; j < dim; j++) {
		pixels[dim - 1 + j*dim] = intens;
		pixels[dim - 2 + j*dim] = intens;
		pixels[dim - 3 + j*dim] = intens;
	}
*/

//if (adjust_window)
//snapshot(*vis_image, "snap.im", true);

	// Set intensity windowing parameters for OpenGl
	if (max != min) {
		windowing_scale = (1.0 + (float) MAX_GREY_VALUE
			- (float) MIN_GREY_VALUE)/(1.0 + max - min);
		offset =  -min/((float) (max - min));
#ifdef DEBUG
		cout << "Defining Gl mapping:\n";
		cout << "\n\tmin = " << min << "    max = " << max;
		cout << "\n\tscale = " << windowing_scale << "    bias = " << offset << endl;
#endif

		// Adjust for "gamma"
		windowing_scale *= (1.0 + gfact/3.0);
//		offset +=  -(windowing_scale - 1.0)*gfact;
		offset *=  (1.0 - 0.5*gfact);
#ifdef DEBUG
		cout << "\tAdjusting for gamma = " << gfact;
		cout << "\n\t\tscale = " << windowing_scale << "    bias = " << offset << endl;
#endif
/*
//		windowing_scale *= (1.0 - gfact);
		offset *=  (1.0 - 0.5*gfact);
*/
#ifdef DEBUG
		cout << "Defining Gl mapping:\n";
		cout << "\tgamma = " << gfact;
		cout << "\n\tmin = " << min << "    max = " << max;
		cout << "\n\tscale = " << windowing_scale << "    bias = " << offset << endl;
#endif
	}
	else {
		windowing_scale = 1.0;
		offset = 0.0;
	}
}



/*
	Function to complete production of the cut-plane image.  This is called
	from gen_cut_view().

	Several class variables are set in gen_cut_view(), before this function is
	called.  Points mdl_end_0, mdl_end_1 and mdl_vertex define the axis lines
	to be drawn, with mdl_vertex being the common point.  Points mdl_end_y0 and
	mdl_end_y1 are the boundary points projected into the drawing plane in
	model space.  All five are in model space initially, and are returned in
	display space.  Mdl_eta, the end-atom extender, is simply converted from
	model space to display space.

  	In addition, line_scale must contain the ratio of the square cut-view
	image's dimensions to the square viewport's dimensions in scale-space.
	This is computed in define_atom_cut.

	The argument is the dimension of the visualization image, set in
	define_atom_cut().
*/
void CutPlanes::define_lines(int dim)
{
	double dx, dy;
	Vector3D u, v;
	Point io, ie, ic;


	// Transform points needed to draw axes to display space
	// Scale is the number of model units per pixel.
	int xcenter, ycenter;
	xcenter = (int) ((dim + 0.5)/2.0);
	ycenter = (int) ((dim + 0.5)/2.0);
#ifdef DEBUG
		cout << "\nCutPlanes::define_lines(): In model space:\n";
		cout << "\tmdl_end_0: "; mdl_end_0.print();
		cout << "\tmdl_end_1: "; mdl_end_1.print();
		cout << "\tmdl_vertex: "; mdl_vertex.print();
		cout << "\tmdl_end_y0: "; mdl_end_y0.print();
		cout << "\tmdl_end_y1: "; mdl_end_y1.print();
#endif

	// Compute ends of the axial vectors in display space
	u = mdl_end_0 - mdl_vertex;
#ifdef DEBUG
	cout << "|mdl_end_0 - mdl_vertex| = " << u.norm() << '\n';
#endif
	end_0.set(xcenter, ycenter + line_scale*u.norm(), 0.0);

	u = mdl_end_1 - mdl_vertex;
#ifdef DEBUG
	cout << "|mdl_end_1 - mdl_vertex| = " << u.norm() << '\n';
#endif
	end_1.set(xcenter - line_scale*u.norm(), ycenter, 0.0);

	// Compute lengths of the y0 and y1 vectors to be drawn.
	// Mdl_end_y0 is the projection of bdry0 on B and mdl_end_y1
	// is the projection of bdry0 on N, in most cases.  In case
	// bperp_y0, mdl_end_y0 will be bdry0 and end_y1 will be the projection
	// of bdry1 onto the bdry0-atom line .  Similarly, for bperp_y1,
	// end_y1 will be bdry1 and mdl_end_y0 will be the projection
	// of bdry0 onto the bdry1-atom line.  See genCutPlanes.cpp.
	u = mdl_end_y0 - mdl_vertex;
#ifdef DEBUG
	cout << "|end_y0 - vertex| = " << u.norm() << '\n';
#endif
	dx = line_scale*u.norm();
	v = mdl_end_y1 - mdl_vertex;
#ifdef DEBUG
	cout << "|end_y1 - vertex| = " << v.norm() << '\n';
#endif
	dy = line_scale*v.norm();

	// Compute ends of the y0 and y1 vectors
	switch (plane) {
		case b_n:
						end_y0.set(xcenter - dy, ycenter + dx, 0.0);
						end_y1.set(xcenter + dy, ycenter + dx, 0.0);
						break;

		case b_bperp:
						end_y0.set(xcenter, ycenter + dx, 0.0);
						end_y1 = end_y0;
						break;

		case bperp_n:
		case bperp_n_chord:		// i == j for bperp_n_chord
						end_y0.set(xcenter, ycenter + dy, 0.0);
						end_y1.set(xcenter, ycenter - dy, 0.0);
						break;

		case bperp_y0:
						end_y0.set(xcenter, ycenter + dx, 0.0);
						end_y1.set(xcenter, ycenter - dy, 0.0);
						break;

		case bperp_y1:
//						end_y0.set(xcenter, ycenter + dx, 0.0);
//						end_y1.set(xcenter, ycenter - dy, 0.0);
						end_y1.set(xcenter, ycenter + dy, 0.0);
						end_y0.set(xcenter, ycenter - dx, 0.0);
						break;
	}

	// Convert eta to display space
	eta = line_scale*mdl_eta;

	// Specify the vertex in display space
	vertex.set(xcenter, ycenter, 0.0);
#ifdef DEBUG
	cout << "In display space:\n";
	cout << "\tcenter: (" <<xcenter << ", " << ycenter << ")\n";
	cout << "\tend_0: "; end_0.print();
	cout << "\tend_1: "; end_1.print();
	cout << "\tvertex: "; vertex.print();
	cout << "\tend_y0: "; end_y0.print();
	cout << "\tend_y1: "; end_y1.print();
	cout << "\teta: " << eta << endl;

    report_intensity_range(*vis_image);	// Diagnostic output
#endif

}


void CutPlanes::convert_boundary(int dim, Point displayCenter)
{
	Vector3D u, v, w;

	int xcenter, ycenter;
	xcenter = (int) ((dim + 0.5)/2.0);
	ycenter = (int) ((dim + 0.5)/2.0);
/*	offset = corner - atom->getX();
	qc.rotateVector(offset);
	offset *= line_scale;
	cout << "Offset = "; offset.print();
*/

#ifdef DEBUG
	cout << "\nCutPlanes::convert_boundary(" << dim << ")\n";
#endif
	u = corner - origin;
	u.normalize();
	v = end - corner;
	v.normalize();
	w = u.cross(v);
	Quat q(u, v, w);	// AGG: This is may now be broken because of the matrixToQuat() change by Rohit
	Quat qc = q.conj();		// Get inverse quaternion

#ifdef DEBUG
	cout << "u = " ; u.print();
	cout << "v = " ; v.print();
	cout << "w = " ; w.print();
	cout << "q = "; q.print();
	cout << "qc = "; qc.print();
#endif

	for (int i = 0; i < bdry.size(); i++) {
#ifdef DEBUG
		cout << "(" << bdry[i].getX() << ", " << bdry[i].getY() << ", " << bdry[i].getZ()
			<< ") -> (";
#endif
		Point p = bdry[i] - displayCenter;	// Placing origin at center of display
		qc.rotateVector(p);
#ifdef DEBUG
 		cout << p.getX() << ", " << p.getY() << ", " << p.getZ() << ")\n";
#endif
		p.setX(p.getX()*line_scale);
		p.setY(-p.getY()*line_scale);
		p.setX(p.getX() + xcenter);
		p.setY(p.getY() + ycenter);
		p.setZ(0.0);
#ifdef DEBUG
		cout << "\t-> (" << p.getX() << ", " << p.getY() << ", " << p.getZ() << ")\n";
#endif
		bdry[i] = p;
	}
#ifdef DEBUG
	cout << endl;
#endif
}

