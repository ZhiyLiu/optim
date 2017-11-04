
#include <math.h>
#include "CutPlanes.h"
#include "M3DObject.h"

#ifndef M_PI
#define M_PI	3.1415926535897932384626433832795
#endif


//#define DEBUG


using namespace std;





/* 
	Function to produce a visualization of a cut plane through an atom 
	embedded in the image data.  The rendering of this visualization is
	performed in class P3DCutPlaneView.

	Either one of the six orthogonal planes defined by the frame of the atom
	can be displayed, depending on the value of CutPlanes::plane.  Member
	variable atom must contain a primitive, usually one "marked" by the user.
	The first argument is a minimum distance in model coordinates for the
	width and height of the plane.  This will typically be a multiple of r or
	the inter-atom distance.  Arguments max_width and max_height are the
	maximum dimensions of the display window in pixels.  They may be zero,
	in which case reasonable values will be determined from r.

	The Point variable of the class return points for line drawing in display
	coordinates.  Points mdl_end_0, mdl_vertex, and mdl_end_1 define the main
	axes to be drawn, thusly:
		View:	b_n	  b_bperp   bperp_n	  bperp_y0	 bperp_y1   bperp_n_chord
	mdl_end_0:	B	  B			N		  y0		 y1			N (offset)
	mdl_end_1:	N	  B-perp	B-perp	  B-perp	 B-perp		B-perp (offset)
	Mdl_vertex is the point of origin of the lines (the atom center).  Points
	mdl_end_y0 and mdl_end_y1 similarly return the ends of the Y vectors.  (Note
	that these values are given temporary values here, which are used and
	replaced by define_lines(), called below).  The length of the B-perp vectors
	(for drawing) is set to r; the length of B is set to r*cos(theta), and the
	length of N is set to r*sin(theta).  If the atom is an end primitive,
	mdl_eta will be set to a non-zero value, indicating how long to draw the
	extended B vector.  Mdl_eta will be made zero for views that do not show
	the B vector.

	Class variables origin, corner and end will be set to contain model-space
	coordinates defining the viewport of the visualization.  The mdl_end_0
	vector will point to the origin-corner edge and the mdl_end_1 vector will
	point away from the corner-end edge.  The mdl_vertex (atom) point is always
	located at the center of the viewport.

  	The input image must be stored in CutPlanes::input_image and the output
	image will be allocated and returned in CutPlanes::vis_image.  This image
	will always be square and have size that is a power of 2, because OpenGl and
	many graphics boards prefer or require such dimensions.

*/

void CutPlanes::gen_cut_view(double dist)
{

	Vector3D v;			// Temporary storage
	Vector3D b, y0, y1, n, bperp, long_b;
	double theta, r, d;
	Point p, q;
	//  The origin and end are guaranteed to be diagonally opposed points. 


	p = atom->getX();
	b = atom->getB();
	n = atom->getN();
	theta = atom->getTheta();
	r = atom->getR();
	y0 = atom->getY0();
	y1 = atom->getY1();

	mdl_vertex = p;

#ifdef DEBUG
	cout << "CutPlanes::gen_cut_view()\n";
	cout << "Atom: "; p.print();
	cout << "R: " << atom->getR() << '\n';
	cout << "Unit B: "; b.print(); cout << "|Unit B|: " << b.norm() << '\n';
	cout << "N: "; n.print();
	cout << "Y0: "; y0.print();
	cout << "Y1: "; y1.print();
#endif

	// Generate a square viewport in model space for each cut plane
	switch (plane) {
		case b_n:		if (atom->type() == M3D_END_PRIMITIVE) {
							long_b = atom->getExtendedB();

#ifdef DEBUG
							cout << "Extended B: "; long_b.print();
							cout << "|Ext-B|: " << long_b.norm() << '\n';
#endif
							mdl_eta = long_b.norm();
#ifdef DEBUG
							cout << "Eta = " << mdl_eta << endl;
#endif
						}
						else
							mdl_eta = 0.0;

#ifdef DEBUG
						mdl_end_y0 = r*cos(theta)*b;	// B vector
						cout << "B: "; mdl_end_0.print();
						cout << "|B|: " << mdl_end_0.norm() << '\n';
#endif
						mdl_end_y0 = p + r*cos(theta)*b;	// Projection of bdry0 on B
						mdl_end_y1 = p + r*cos(M_PI/4.0
										- theta)*n;			// Projection of bdry0 on N

						mdl_end_0 = mdl_end_y0;		// mdl_end_0 is the end of the B vector
						mdl_end_1 = mdl_end_y1;		// mdl_end_1 is the end of the N vector
						origin = p + r*b + r*n;		// Initial values
						corner = p + r*b - r*n;

						// Make sure the box has reasonable width and height
						q = p + r*b;		// Projection of origin on B
						v = q - p;			// A vector parallel to B
						d = v.norm();
						if (d < dist) {
							v = ((dist - d)/d)*v;	// Add component in B direction
							origin += v;
							corner += v;
						}
						q = p + r*n;		// Projection of bdry0 on N
						v = q - p;			// A vector parallel to N
						d = v.norm();
						if (d < dist) {
							v = ((dist - d)/d)*v;	// Add component in N direction
							origin += v;
							corner -= v;
						}

						// Complete the bounding box
						v = p - origin;
						end = p + v;
						break;

		case b_bperp:	if (atom->type() == M3D_END_PRIMITIVE) {
							long_b = atom->getExtendedB();
#ifdef DEBUG
							cout << "Extended B: "; long_b.print();
							cout << "|Ext-B|: " << long_b.norm() << '\n';
#endif
							mdl_eta = long_b.norm();
#ifdef DEBUG
							cout << "Eta = " << mdl_eta << endl;
#endif
						}
						else
							mdl_eta = 0.0;

#ifdef DEBUG
						mdl_end_y0 = r*cos(theta)*b;	// B vector
						cout << "B: "; mdl_end_y0.print();
						cout << "|B|: " << mdl_end_y0.norm() << '\n';
#endif
						mdl_end_y0 = p + r*cos(theta)*b;	// Projection of bdry0 on B
						mdl_end_y1 = p + r*cos(M_PI/4.0
										- theta)*n;			// Projection of bdry0 on N

						bperp = atom->getBPerp();	// Length == 1
#ifdef DEBUG
						cout << "Bperp: "; bperp.print();
#endif
						mdl_end_1 = p + r*bperp;	// mdl_end_1 is r along the Bperp vector
						mdl_end_0 = mdl_end_y0;		// mdl_end_0 is the end of the B vector

						q = mdl_end_y0;		// Projection of bdry0 on B

						// Make sure q is a reasonable distance from p
						v = q - p;			// A vector parallel to B
						d = v.norm();
						if (d < dist) {
							q = p + (dist/d)*v;
							d = dist;
						}

						// Determine bounding box
						origin = q + d*bperp;
						corner = q - d*bperp;
						v = p - origin;
						end = p + v;
						break;

		case bperp_n_chord:
		case bperp_n:	mdl_eta = 0.0;
						bperp = atom->getBPerp();	// Length == 1
#ifdef DEBUG
						cout << "Bperp: "; bperp.print();
#endif
						mdl_end_y0 = p + r*cos(theta)*b;	// Projection of bdry0 on B
						mdl_end_y1 = p + r*cos(M_PI/4.0
										- theta)*n;			// Projection of bdry0 on N

						mdl_end_0 = mdl_end_y1;		// mdl_end_0 is the end of the N vector
						mdl_end_1 = p + r*bperp;	// mdl_end_1 is r along the Bperp vector

						q = mdl_end_y1;		// Projection of bdry0 on N

						// Make sure q is a reasonable distance from p
						v = q - p;			// A vector parallel to N
						d = v.norm();
						if (d < dist) {
							q = p + (dist/d)*v;
							d = dist;
						}

						// Determine bounding box
						origin = q + d*bperp;
						corner = q - d*bperp;
						v = p - origin;
						end = p + v;

						if (plane == bperp_n) break;

						// bperp_n_chord code: simply move frame to chord middle
						v = 0.5*(y0 + y1);
						origin += v;
						corner += v;
						end += v;
						vertex += v;
						mdl_end_0 += v;
						mdl_end_1 += v;

						break;

		case bperp_y0:	mdl_eta = 0.0;
						bperp = atom->getBPerp();	// Length == 1
#ifdef DEBUG
						cout << "Bperp: "; bperp.print();
#endif
						mdl_end_y0 = p + y0;	// End of the y0 vector
						// End of the y1 vector projected onto the line through
						// p and bdry0.
						d = y1.norm()*cos(M_PI - 2.0*theta); // Length of projection
						mdl_end_y1 = p - (d/r)*y0;	// Note that ||y0|| = r

						mdl_end_0 = mdl_end_y0;		// mdl_end_0 is the end of the y0 vector
						mdl_end_1 = p + r*bperp;	// mdl_end_1 is r along the Bperp vector

						// Make sure q is a reasonable distance from p
						d = y0.norm();
						if (d < dist) {
							q = p + (dist/d)*y0;
							d = dist;
						}
						else q = end_0;

						// Determine bounding box
						origin = q + d*bperp;
						corner = q - d*bperp;
						v = p - origin;
						end = p + v;
						break;

		case bperp_y1:	mdl_eta = 0.0;
						bperp = atom->getBPerp();	// Length == 1
#ifdef DEBUG
						cout << "Bperp: "; bperp.print();
#endif
						// End of the y0 vector projected onto the line through
						// p and bdry1.
						d = y0.norm()*cos(M_PI - 2.0*theta); // Length of projection
						mdl_end_y0 = p - (d/r)*y1;	// Note that ||y1|| = r
						mdl_end_y1 = p + y1;		// End of the y1 vector

						mdl_end_0 = mdl_end_y1;		// mdl_end_0 is the end of the y1 vector
						mdl_end_1 = p + r*bperp;	// mdl_end_1 is r along the Bperp vector

						// Make sure q is a reasonable distance from p
						d = y1.norm();
						if (d < dist) {
							q = p + (dist/d)*y1;
							d = dist;
						}
						else q = mdl_end_0;

						// Determine bounding box
						origin = q + d*bperp;
						corner = q - d*bperp;
						v = p - origin;
						end = p + v;
						break;
	}

	// Function define_atom_cut() will replace vis_image
}


