
#include <vector>
#include <math.h>

#include "CutPlanes.h"
#include "M3DObject.h"
#include "utility.h"
#define D_XFERLIST
#define D_POINTLIST_SERVER2
#include "Shapedepend.h"
#include "renderDefinitions.h"
#include "M3DFigure.h"

//#define DEBUG


using namespace std;


// Default factor of the size of the area around atoms to be viewed
#define ATOMIC_ZOOM    2.0		/* A factor of r */



CutPlanes::~CutPlanes() {
	delete vis_image;
	ThallCode::Pointlist_server2 * pL = (ThallCode::Pointlist_server2 *) plist;
	if (pL)
		delete pL;
	// input_image and atom point to objects not owned by CutPlanes
}


CutPlanes::CutPlanes() {
#ifdef DEBUG
	cout << "CutPlanes::CutPlanes()" << endl;
#endif
	input_image = NULL;
	vis_image = NULL;
	figure = NULL;
	min_size = 0.0;
	atom = NULL;
	plane = b_n;
	opposing_y = false;
	gfact = 0.0;
	adjust_window = false;
	plist = NULL;
	status = false;
}


CutPlanes::CutPlanes(const M3DObject * object, Image3D * image) {
	input_image = NULL;
	vis_image = NULL;
	figure = NULL;
	min_size = 0.0;
	atom = NULL;
	plane = b_n;
	opposing_y = false;
	gfact = 0.0;
	adjust_window = false;
	plist = NULL;
	status = init(object, image);
}


void CutPlanes::reset() {
#ifdef DEBUG
	cout << "CutPlanes::reset()" << endl;
#endif
	delete vis_image;
	vis_image = NULL;

	ThallCode::Pointlist_server2 * pL = (ThallCode::Pointlist_server2 *) plist;
	if (pL) {
		delete pL;
		plist = NULL;
	}

	//	input_image and obj are not reset

	figure = NULL;
	min_size = 0.0;
	atom = NULL;
	plane = b_n;
	opposing_y = false;
	gfact = 0.0;
	adjust_window = false;
	status = false;
}


CutPlanes::CutPlanes(const M3DPrimitive * prim, const M3DObject * object,
					 Image3D * image)
{
	input_image = NULL;
	vis_image = NULL;
	figure = NULL;
	min_size = 0.0;
	atom = NULL;
	plane = b_n;
	opposing_y = false;
	gfact = 0.0;
	adjust_window = false;
	plist = NULL;
	status = init(prim, object, image);
}


/*	This version of init() takes a pointer to the object and finds
	the first selected primitive in it, for which it produces the
	visualization.  It uses the first selected atom, because the
	M3DPrimitive class cannot be queried to find which atom is
	marked.
*/
bool CutPlanes::init(const M3DObject * object, Image3D * image) {
#ifdef DEBUG
	cout << "CutPlanes::init()" << endl;
#endif
	delete vis_image;
	vis_image = NULL;
	obj = object;
	figure = NULL;

	input_image = image;
#ifdef DEBUG
	report_intensity_range(*image);	// Diagnostic output
#endif

	if (input_image == NULL || object == NULL) return false;

	int nprim = object->getPrimitiveCount();
	atom = NULL;
	for (int i = 0; i < nprim; i++) {
		M3DPrimitive * p = object->getPrimitivePtr(i);
		if (p->isSelected()) {							// AGG: is Marked?
			atom = p;
			if (min_size == 0.0)
				msize = ATOMIC_ZOOM*p->getR();
			else
				msize = min_size;
#ifdef DEBUG
			cout << "Minimum size set to " << msize << '\n';
#endif
			// The casts to float effect use of a fuzz
			if ((float) atom->getTheta() == (float) M_PI)
				opposing_y = true;
			else
				opposing_y = false;
			break;			// Work with the first selected atom found
		}
	}
	if (atom == NULL) {
		cout << "No atoms selected: cut planes cannot be visualized" << endl;
		return false;
	}
	return true;
}


bool CutPlanes::init(const M3DPrimitive * prim, const M3DObject * object,
					 Image3D * image)
{
#ifdef DEBUG
	cout << "CutPlanes::init():  atom = " << hex << prim << dec << endl;
#endif
	delete vis_image;
	vis_image = NULL;
	obj = object;
	figure = NULL;

	input_image = image;
	atom = (M3DPrimitive *) prim;

	if (input_image == NULL || atom == NULL) return false;

	if (min_size == 0.0)
		msize = ATOMIC_ZOOM*atom->getR();
	else
		msize = min_size;
#ifdef DEBUG
	cout << "Minimum size set to " << msize << '\n';
#endif
			// The casts to float effect use of a fuzz
	if ((float) atom->getTheta() == (float) M_PI)
		opposing_y = true;
	else
		opposing_y = false;
	return true;
}


/* Initialization function used to generate the curves of figures
    intersecting image planes displayed in the main window.
*/
void CutPlanes::init2(const M3DObject * object, int level)
{
	ThallCode::Pointlist_server2 * pL;

#ifdef DEBUG
	cout << "CutPlanes::init2() called with level = " << level << endl;
#endif
	obj = object;
	lvl = level;
	pL = (ThallCode::Pointlist_server2 *) plist;
        if (pL) {
		delete pL;
                plist = NULL;
        }
}


void CutPlanes::setLevel(int level) { 
#ifdef DEBUG
	cout << "CutPlanes::setLevel() called with level = " << level << endl;
#endif
	lvl = level;
}


/*  The copy constructor should only be used when the object is in its
    default state, since no attempt is made to copy images, point lists,
	etc.
*/
CutPlanes::CutPlanes(const CutPlanes & cp) {
#ifdef DEBUG
	cout << "CutPlanes::CutPlanes(const CutPlanes & cp)" << endl;
#endif
	// This test was added because this code was crashing ConStruct.  In
	// Visual Studio, when a std::vector or CutPlanes objects is resized,
	// the standard library calls the copy constructor, instead of the
	// default constructor, to create the new instances of the class.
	// I am not satisfied that this is entirely the proper solution.  AGG.
        if (cp.status) {
	    input_image = cp.input_image;
	    vis_image = cp.vis_image;
	    figure = cp.figure;
	    min_size = cp.min_size;
	    atom = cp.atom;
	    plane = cp.plane;
	    opposing_y = cp.opposing_y;
	    gfact = cp.gfact;
	    adjust_window = cp.adjust_window;
	    plist = cp.plist;
	    status = cp.status;
        }
        else {
	    input_image = NULL;
	    vis_image = NULL;
	    figure = NULL;
	    min_size = 0.0;
	    atom = NULL;
	    plane = b_n;
	    opposing_y = false;
	    gfact = 0.0;
	    adjust_window = false;
	    plist = NULL;
	    status = false;
        }
}


void CutPlanes::gamma(float g) {
	if (g < 0.0) g = 0.0;
	else if (g > 1.0) g = 1.0;
	else gfact = g;
#ifdef DEBUG
	cout << "CutPlanes::gamma(" << g << ") set gfact to " << gfact << endl;
#endif
}


//  Find the next higher power of 2

unsigned int powerize(unsigned int i)
{
        unsigned int     pwr;

        if (i == 1)
			return i;
        pwr = 1;
        while (i > pwr)
			pwr <<= 1;
        return pwr;
}


bool CutPlanes::run(int width, int height) {

	int dim;

#ifdef DEBUG
	cout << "CutPlanes::run(" << width << ", " << height << "):  atom = "
		<< hex << atom << dec << endl;
#endif
	if (height < width)
		dim = height;
	else
		dim = width;

	// Adjust dim to a power of 2.  The lower possiblity is used, unless
	// the original value of dim is 2/3 the way to the higher possibility.
	int d = (int) powerize((unsigned int) dim);
	if (d > dim) {
		int c = d >> 1;
		int x = c + (int) (2.0*(d - c)/3.0);
#ifdef DEBUG
		cout << "dim = " << dim;
#endif
		if (x < dim)
			dim = d;
		else dim = c;
#ifdef DEBUG
		cout << "  x = " << x << "  c = " << c << "  d = " << d << "  new dim = " << dim << endl;
#endif
	}
	// else d == dim, so dim is a power of 2

	if (atom == NULL || input_image == NULL || input_image->getVoxels() == NULL) {
		cout << "\nCut-plane visualization failed" << endl;
		return false;
	}

	delete vis_image;
	vis_image = NULL;

	// Run the cut plane visualization
#ifdef DEBUG
	cout << "Calling genCutView() with size = " << msize
		<< " and type of ";
	switch (plane) {
		case b_n:		cout << "b_n";
						break;
		case b_bperp:	cout << "b_bperp";
						break;
		case bperp_n:	cout << "bperp_n";
						break;
		case bperp_n_chord:	cout << "bperp_n_chord";
						break;
		case bperp_y0:	cout << "bperp_y0";
						break;
		case bperp_y1:	cout << "bperp_y1";
						break;
		default:		cout << "ERROR";
						break;
	}
	cout << endl;
#endif

	// Note: this forces the visualization to be square for OpenGl
	gen_cut_view(msize);

#ifdef DEBUG
	cout << "\torigin = ";
	origin.print();
	cout << "\tend = ";
	end.print();
	cout << "\tcorner = ";
	corner.print();
#endif

	if (drawing_bdry) {
		initialize_boundary(true);
		generate_boundary();
#ifdef DEBUG
		cout << "Boundary curve size = " << bdry.size() << endl;
#endif
	}

	define_atom_cut(dim, dim);

	if (drawing_bdry)
		convert_boundary(dim, atom->getX());
	define_lines(vis_image->getXDim());

	return true;
}


/*  This function generates the curves of figures intersecting
    image planes displayed in the main window.
*/
void CutPlanes::run(int figureId, Point viewPortOrigin, Point viewPortCorner,
					Point viewPortEnd, int width, int height)
{
	int dim;
	Point center;

#ifdef DEBUG
	cout << "CutPlanes::run(figureId " << figureId << ", Point, Point, Point, "
		<< width << ", " << height << ')' << endl;
#endif
	if (height < width)
		dim = height;
	else
		dim = width;
	if (dim == 0)
		dim = 512;

	origin = viewPortOrigin;
	corner = viewPortCorner;
	end = viewPortEnd;
	center = 0.5*(origin + end);
	figure = (M3DFigure *) obj->getFigurePtr(figureId);
	initialize_boundary(false);
	generate_boundary();
#ifdef DEBUG
	cout << "Produced " << bdry.size() << " verticies" << endl;
#endif
}


int CutPlanes::width() {
        if (vis_image == NULL) return 0;
        return vis_image->getXDim();
}


int CutPlanes::height() {
        if (vis_image == NULL) return 0;
        return vis_image->getYDim();
}


/*  Given the origin, corner and end coordinates of the viewport,
    generate the tiles needed for calculating the boundary curve.
	Standard should be true, if one of the 6 cut plane visualizations
	is to be drawn; and false, if a curve is being drawn in the main
	window.
*/
void CutPlanes::initialize_boundary(bool standard)
{
	ThallCode::Pointlist_server2 * pL;
	int i;

#ifdef DEBUG
	cout << "CutPlanes::initialize_boundary(" << standard << ") called with lvl = " << lvl << endl;
#endif
	if (standard) {
		if (atom == NULL)
			return;
		if (input_image == NULL)
			return;
	}
	if (obj == NULL)
		return;

	bdry.erase(bdry.begin(), bdry.end());
#ifdef COLOR_CONTOURS
	bpointIndices.erase(bpointIndices.begin(), bpointIndices.end());
#endif

	// Compute the normal to the display plane
	Vector3D v = origin - corner;
	Vector3D w = end - corner;
	normal = w.cross(v);
	normal.normalize();		// AGG:  needed?     
#ifdef DEBUG
	cout << "Origin: "; origin.print();
	cout << "Corner: "; corner.print();
	cout << "End: "; end.print();
	cout << "normal: "; normal.print();
#endif

	// Get a pointer to the marked figure
    if (figure == NULL) {
		int num_figures = obj->getFigureCount();
		for (i = 0; i < num_figures; i++) {
			figure = (M3DFigure *) obj->getFigurePtr(i);
			for (int j = 0; j < figure->getPrimitiveCount(); j++) {
				if (atom == figure->getPrimitivePtr(j)) {
					i = num_figures;	// Force break from outer loop
					break;
				}
			}
		}
	}

	const float * color = ((M3DFigure *) figure)->getColor();
#ifdef DEBUG
	cout << "Figure has " << figure->getPrimitiveCount() << " atoms\n";
	cout << "Figure = " << figure << endl;
	cout << "Boundary color: (" << color[0] << ", " << color[1]
		<< ", " << color[2] << ')' << endl;
#endif
	bdry_color[0] = color[0];
	bdry_color[1] = color[1];
	bdry_color[2] = color[2];

	// Generate boundary tiles for the figure
	if (standard)
		lvl = SUBDIV_MEDIAL_SUBDIVISIONS;
    if (! standard && plist != NULL)
	return;

    Xferlist * xferList = convertM3DtoXfer(figure);

	pL = (ThallCode::Pointlist_server2 *) plist;
	if (pL)
		delete pL;  	    // AGG: Is this really necessary?
	    pL = new ThallCode::Pointlist_server2;
	pL->init(xferList);
	    pL->ComputeSubdivBoundaryTiles(lvl);

	    // Added by CL.  Add displacement values to boundary points.
	    pL->ComputeSubdivPointCloud(lvl); 
	    SubdivBoundary * boundary = figure->getBoundaryPtr();
	    if (boundary != NULL)
	    {
		    Displacements * dsp = boundary->getDisplacements(lvl);
		    if (dsp != NULL)
			    pL -> applyDisplacements(dsp->getNumPts(), dsp->getVals());
	    }

	plist = (void *) pL;
}


void CutPlanes::generate_boundary()
{
	Bpoint * tileset;
	int numtiles;

	ThallCode::Pointlist_server2 * pL = (ThallCode::Pointlist_server2 *) plist;
	if (! pL)
		return;

    pL->subdivtileinfo(&numtiles, &tileset);

    int *tileIndices = 0; // this will be used below to pass
                          // indices to clip_boundary
#ifdef COLOR_CONTOURS
// Map the tile vertices back to bpoints so we can recover their
// color information
		int numTilesToIndex = 0;
		pL->subdivtileindexinfo(&numTilesToIndex, &tileIndices);

		// printf("number of tile indices= %d\n", numTilesToIndex);
#endif /* COLOR_CONTOURS */

#ifdef DEBUG
	printf("number of tiles = %d    number of points = %d\n", numtiles, 4*numtiles);
/*
	// Print out vertices
	Point atom_posn = atom->getX();
	for (i = 0; i < 4*numtiles; i++) {
//		if (inter_pt_distance(tileset[i].pnt, atom_posn) < 0.12)
			cout << '(' << tileset[i].pnt[0] << ", " << tileset[i].pnt[1]
				<< ", " << tileset[i].pnt[2] << ")  " << i << endl;
	}
*/
#endif
	clip_boundary(numtiles, (void *) tileset, tileIndices);
}


/*	This function works very much like Cyrus-Beck clipping.  It finds
	the intersection of a line, parameterized by 0 <= t <= 1 between
	the two arguments and the viewport's plane, defined by the normal,
	which must have been caluculated before entry.  If an intersection
	occurs inside the viewport, the function returns 1.  If the line
	intersects outside the viewport or lies in the plane, it returns 0.
	If the line lies in the plane, -1 is returned and the intersection
	argument is not set.  If a single endpoint lies in the plane, 2
	is returned.
*/
int CutPlanes::intersect_edge_with_plane(Point & pt0, Point & pt1,
										   Point & intersection)
{
	Vector3D l, m;
	double denom, numer, t;

	l = pt1 - pt0;
	denom = normal*l;
	if (denom == 0) {
		// The edge lies parallel to the display viewport.
		// Must check that the line is in the cut plane.
		double dp = normal*origin;
		Point q = project_point_to_plane(pt0, normal, dp);
		if (inter_pt_distance(q, pt0) == 0.0)
			return -1;	// Edge is in the plane
		else
			return 0;	// Ignore this edge
	}
	m = pt0 - origin;
	numer = normal*m;
	t = -numer/denom;
	intersection = pt0 + t*l;
	if (t < 0.0) return 0;
	if (t > 1.0) return 0;

	double d = inter_pt_distance(pt0, intersection);
	if (d == 0.0)
		return 2;

	d = inter_pt_distance(pt1, intersection);
	if (d == 0.0)
		return 2;

	return 1;
}


/*  This function processes a single triangle, drawing the
    intersection(s) of it with the view plane.
*/
void CutPlanes::clip_triangle(double * vtx0, double * vtx1, double * vtx2, int idx0, int idx1, int idx2)
{
	int flag0, flag1, flag2, sum, n;
	Point q0, q1, q2, p0, p1, p2;
	double dp;

	p0 = Point(vtx0);
	p1 = Point(vtx1);
	p2 = Point(vtx2);

	flag0 = intersect_edge_with_plane(p0, p1, q0);
	flag1 = intersect_edge_with_plane(p1, p2, q1);
	flag2 = intersect_edge_with_plane(p2, p0, q2);
/*
	cout << "vtx0: "; p0.print();
	cout << "vtx1: "; p1.print();
	cout << "vtx2: "; p2.print();
	cout << "flag0 = " << flag0 << "   q0: "; q0.print();
	cout << "flag1 = " << flag1 << "   q1: "; q1.print();
	cout << "flag2 = " << flag2 << "   q2: "; q2.print();
*/

	if (flag0 != 0 || flag1 != 0 || flag2 != 0) {

		if (flag0 >= 0 && flag1 >= 0 && flag2 >= 0) {
			sum = flag0 + flag1 + flag2;
			switch (sum) {
				case 6:		// 2,2,2
					// This should only occur occassionally in very distorted objects.
					// In theory, the points of intersection define a triangle that
					// should be filled in the display.  For the moment, this triangle
					// is drawn.
//					cout << "Error: 2,2,2 touch" << endl;
					bdry.push_back(q0);
					bdry.push_back(q1);
					bdry.push_back(q1);
					bdry.push_back(q2);
					bdry.push_back(q2);
					bdry.push_back(q0);
#ifdef COLOR_CONTOURS
					bpointIndices.push_back(idx0);
					bpointIndices.push_back(idx1);
					bpointIndices.push_back(idx1);
					bpointIndices.push_back(idx2);
					bpointIndices.push_back(idx2);
					bpointIndices.push_back(idx0);
#endif  /* COLOR_CONTOURS */
					break;

				case 5:		// 1,2,2
					// In this case the two intersections returning 2 are the same
					// point, both resulting in an edge ending in the view plane.
					if (flag0 == 2)
						flag0 = 0;
					else if (flag1 == 2)
						flag1 = 0;
					// Fall through to case 3

				case 3:		// 0,1,2
					// In this case intersections returning 0 can be treated
					// as returning 2; the apparent error is a result of the
					// finite length of variables in intersect_edge_with_plane().
					// Essentially, two intersections are the same point, both
					// resulting in an edge ending in the view plane.
					if (flag0 == 1) {
						bdry.push_back(q0);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx0);
#endif /* COLOR_CONTOURS */
					}
					else if (flag0 == 2) {
						bdry.push_back(q0);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx0);
#endif
					}
					if (flag1 == 1) {
						bdry.push_back(q1);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx1);
#endif
					}
					else if (flag1 == 2) {
						bdry.push_back(q1);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx1);
#endif
					}
					if (flag2 == 1) {
						bdry.push_back(q2);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx2);
#endif
					}
					else if (flag2 == 2) {
						bdry.push_back(q2);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx2);
#endif
					}
					break;

				case 4:		// 0,2,2
					// In this case, two edges meet at the same end point lying in
					// the view plane, returning 2.  The third edge, with flag of 0,
					// links the other endpoints of the first two edges and lies
					// entirely outside the plane.
					dp = normal*origin;
					if (flag0 == 0) {
						bdry.push_back(q1);
						bdry.push_back(q2);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx1);
					  	bpointIndices.push_back(idx2);
#endif /* COLOR_CONTOURS */
					}
					else if (flag1 == 0) {
						bdry.push_back(q0);
						bdry.push_back(q2);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx0);
						bpointIndices.push_back(idx2);					 
#endif /* COLOR_CONTOURS */
					}
					else if (flag2 == 0) {
						bdry.push_back(q0);
						bdry.push_back(q1);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx0);
  						bpointIndices.push_back(idx1);		
#endif /* COLOR_CONTOURS */
					}
					break;

				case 2:		// 0,1,1 or 0,0,2
					// Normal two-sided intersection, 0,1,1
					// (0,0,2 - should never occur)
					n = 0;
					if (flag0 == 1) {
						bdry.push_back(q0);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx0);
#endif /* COLOR_CONTOURS */
						n++;
					}
					if (flag1 == 1) {
						bdry.push_back(q1);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx1);
#endif /* COLOR_CONTOURS */
						n++;
					}
					if (flag2 == 1) {
						bdry.push_back(q2);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx2);
#endif /* COLOR_CONTOURS */
						n++;
					}
					if (n == 1) {
						cout << "Error: 0,0,2 touch" << endl;
						bdry.pop_back();
#ifdef COLOR_CONTOURS
						bpointIndices.pop_back();
#endif /* COLOR_CONTOURS */
					}
					break;

				case 1:		// 0,0,1
					// This should theoretically never occur, however can happen
					// because of arithmetic error.  The fixup here, which may not
					// be correct, averages the two odd intersections, under the
					// assumption that they must be the same point.
//					cout << "Error: 0,0,1 crossing" << endl;
					if (flag0 == 1) {
						Point q = q1 + q2;
						q *= 0.5;
						bdry.push_back(q0);
						bdry.push_back(q);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx0);
						bpointIndices.push_back(idx1);
#endif /* COLOR_CONTOURS */
					}
					else if (flag1 == 1) {
						Point q = q0 + q2;
						q *= 0.5;
						bdry.push_back(q1);
						bdry.push_back(q);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx1);
						bpointIndices.push_back(idx0);
#endif /* COLOR_CONTOURS */
					}
					else if (flag2 == 1) {
						Point q = q0 + q1;
						q *= 0.5;
						bdry.push_back(q2);
						bdry.push_back(q);
#ifdef COLOR_CONTOURS
						bpointIndices.push_back(idx2);
						bpointIndices.push_back(idx0);
#endif /* COLOR_CONTOURS */
					}
					break;

				case 0:		// 0,0,0
					// Triangle is away from the view plane: no action
					break;
			}
		}
		else {
			// One or more edges lie in the view plane (flag < 0)
			n = 0;
			if (flag0 < 0) n++;
			if (flag1 < 0) n++;
			if (flag2 < 0) n++;
			if (n == 1) {
				// One edge lies in the plane
				if (flag0 < 0) {
					bdry.push_back(q1);
					bdry.push_back(q2);
#ifdef COLOR_CONTOURS
					bpointIndices.push_back(idx1);
					bpointIndices.push_back(idx2);
#endif /* COLOR_CONTOURS */	
					if (flag1 != 2 && flag2 != 2)
						cout << "Error: " << flag0 << ',' << flag1 << ',' << flag2 << endl;
				}
				if (flag1 < 0) {
					bdry.push_back(q0);
					bdry.push_back(q2);
#ifdef COLOR_CONTOURS
					bpointIndices.push_back(idx0);
					bpointIndices.push_back(idx2);
#endif /* COLOR_CONTOURS */
					if (flag0 != 2 && flag2 != 2)
						cout << "Error: " << flag0 << ',' << flag1 << ',' << flag2 << endl;
				}
				if (flag2 < 0) {
					bdry.push_back(q0);
					bdry.push_back(q1);
#ifdef COLOR_CONTOURS
					bpointIndices.push_back(idx0);
					bpointIndices.push_back(idx1);
#endif /* COLOR_CONTOURS */
					if (flag0 != 2 && flag1 != 2)
						cout << "Error: " << flag0 << ',' << flag1 << ',' << flag2 << endl;
				}
			}
			else {
				cout << "Error: " << flag0 << ',' << flag1 << ',' << flag2 << endl;
			}
		}
	}
	// else all flags == 0: ignore

}

/*	This function intersects the quardilaterals with the viewport, the
	normal of which must have been set before entry.  The quadrilaterals,
	which are assumed to not be planar, are processed as pairs of triangles.
	Every triangle that crosses the viewport will generate two points.  If
	a triangle touches the viewport at a corner, then the corner is recorded
	twice.  Thus the output, stored in bdry, is a list of line segments.
	The clipping of the output to the viewport is left to OpenGL.
*/
void CutPlanes::clip_boundary(int numtiles, void /*Bpoint*/ * tileset, int * tilesetindices)
{
	int i;
	Bpoint * tset;

	tset = (Bpoint *) tileset;

	// Loop over all quadrilateral tiles
	for (i = 0; i < 4*numtiles; i += 4) {

		if (tilesetindices) {
			// This can only occur if COLOR_CONTOURS has been defined
			// First triangle
			clip_triangle(tset[i].pnt, tset[i + 1].pnt, tset[i + 2].pnt, tilesetindices[i], tilesetindices[i+1], tilesetindices[i+2]);

			// Second triangle
			clip_triangle(tset[i + 2].pnt, tset[i + 3].pnt, tset[i].pnt, tilesetindices[i+2], tilesetindices[i+3], tilesetindices[i]);
		}
		else {
			// First triangle
			clip_triangle(tset[i].pnt, tset[i + 1].pnt, tset[i + 2].pnt);

			// Second triangle
			clip_triangle(tset[i + 2].pnt, tset[i + 3].pnt, tset[i].pnt);
		}
	}
}



