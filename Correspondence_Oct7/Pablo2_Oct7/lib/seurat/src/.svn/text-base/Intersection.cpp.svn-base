
#include <iostream>
#include <queue>
#include <set>
#include "M3DObject.h"
#include "M3DDisplayGlobals.h"
#include "M3DQuadFigure.h"
#include "Intersection.h"
#include "utility.h"
#define D_XFERLIST
#define D_SUBDIVSURF
#define D_POINTLIST_SERVER2
#include "Shapedepend.h"

// Uncomment the following to cause the blending regions to completely
// wrap around the hinge.  Otherwise, the two ends are omitted.
//#define CHECK_ENDS	1

// Uncomment the following to cause the blending region end caps to be
// rendered as spheres, when render() is called.
//#define DRAW_ENDS	1


//#define DEBUG

using namespace std;
using namespace ThallCode;

Intersection::Intersection(M3DFigureTreeNode & subNode, M3DObject * subObject,
						   void * parentPList, void * childPList)
{
	int u, v;
	M3DPrimitiveLinkInfo * link;
	M3DFigureTreeNode * parentNode;
	M3DPrimitive * atom0;
	M3DPrimitive * atom1;
	M3DPrimitive * atom2;
	Point p0, p1, p2;
	Vector3D c0, c1;
	Vector3D s0, s1;
	Vector3D parentNormal;
	Pointlist_server2 * parent_pList;
	Pointlist_server2 * child_pList;
	Bpoint bp;
	double dist;
	bool outside;

#ifdef DEBUG
	cout << "Intersection::Intersection()\n";
#endif
	// Determine indexes of the two figures
	subFigureId = subNode.getFigureId();
	parentNode = subNode.getParent();
	mainFigureId = parentNode->getFigureId();	// Needed for .h file
	parent_pList = (Pointlist_server2 *) parentPList;
	child_pList = (Pointlist_server2 *) childPList;

	// Allocate arrays for storage of line segments and corresponding radii
	nlinks = subNode.getLinkCount();
	if (nlinks <= 1) {
		if (nlinks == 1)
			cout << "Figural intersections with a single link atom are not supported" << endl;
		nlinks = 0;
		vertex = NULL;
		radius = NULL;
	}
	vertex = new Point[nlinks + nlinks];
	radius = new double[nlinks + nlinks];

	// For each attached subfigure, determine bounding truncated cone for the blend
	// region.

	M3DFigure * subFigure = subObject->getFigurePtr(subFigureId);

	// Generate first pair of blending region centers and radii
	link = subNode.getLink(0);		// Get link for atom 0
	dynamic_cast<M3DQuadFigure*>( subFigure)->figuralCoordinates(link->primitiveId, u, v);	// Get u and v of atom 0
	atom0 = dynamic_cast<M3DQuadFigure*>( subFigure)->getPrimitivePtr(u, v);	// Get link atom 0
	radius[0] = atom0->getR();	// Initial estimates of radii
	radius[1] = atom0->getR();

	if (! parent_pList->subdivBposition(&bp, link->u, link->v, link->t)) {
		cout << "Plist FAILURE 1!\n";
		destroy();
		return;
	}
	parentNormal = -Vector3D(bp.norm[0], bp.norm[1], bp.norm[2]);

	link = subNode.getLink(1);		// Get link for atom 1
	dynamic_cast<M3DQuadFigure*>( subFigure)->figuralCoordinates(link->primitiveId, u, v);	// Get u and v of atom 1
	atom1 = dynamic_cast<M3DQuadFigure*>( subFigure)->getPrimitivePtr(u, v);	// Get link atom 1
	radius[2] = atom1->getR();	// Initial estimates of radii
	radius[3] = atom1->getR();

	p0 = atom0->getX();
	p1 = atom1->getX();
	c0 = p1 - p0;	// Vector connecting atom 0 with atom 1

	s0 = c0.cross(parentNormal);	// AGG: s0 ~= N; would it be better to use N?
	s0.normalize();
	s0 *= atom0->getR()*sin(atom0->getTheta());		// Multiply by approx distance to surface
	s0 += p0;		// Initial estimate of intersection of surfaces
	if (! child_pList->subdivNearSurface(s0.getX(), s0.getY(), s0.getZ(), &bp, &dist, &outside)) {
		cout << "Plist FAILURE 2!\n";
		destroy();
		return;
	}
	s0 = Point(bp.pnt[0], bp.pnt[1], bp.pnt[2]) - p0;	// Better estimate of vector to intersection

	vertex[0] = p0 + s0;
	vertex[1] = p0 - s0;

	// Refine the radii by use of the parent figure
	if (! parent_pList->subdivNearSurface(vertex[0].getX(), vertex[0].getY(), vertex[0].getZ(),
			&bp, &dist, &outside))
	{
		cout << "Plist FAILURE 3!\n";
		destroy();
		return;
	}
	if (radius[0] < bp.rad)		// Use the larger of the two radii
		radius[0] = bp.rad;
	if (! parent_pList->subdivNearSurface(vertex[1].getX(), vertex[1].getY(), vertex[1].getZ(),
			&bp, &dist, &outside))
	{
		cout << "Plist FAILURE 4!\n";
		destroy();
		return;
	}
	if (radius[1] < bp.rad)		// Use the larger of the two radii
		radius[1] = bp.rad;

	// Intermediate blending region centers and radii
	if (nlinks >= 3)
		for (int n = 1; n < nlinks - 1; n++) {
			if (! parent_pList->subdivBposition(&bp, link->u, link->v, link->t)) {
				cout << "Plist FAILURE 5!\n";
				destroy();
				return;
			}
			parentNormal = -Vector3D(bp.norm[0], bp.norm[1], bp.norm[2]);

			link = subNode.getLink(n + 1);	// Get link
			dynamic_cast<M3DQuadFigure*>( subFigure)->figuralCoordinates(link->primitiveId, u, v);	// Get u and v
			atom2 = dynamic_cast<M3DQuadFigure*>( subFigure)->getPrimitivePtr(u, v);	// Get link atom
			radius[n + n + 2] = atom2->getR();	// Initial estimates of radii
			radius[n + n + 3] = atom2->getR();

			p2 = atom2->getX();
			c1 = p2 - p1;

			double r = atom1->getR()*sin(atom1->getTheta());	// Approx distance to surface
			s0 = c0.cross(parentNormal);
			s0.normalize();
			s0 *= r;		// Multiply by approx distance to surface
			s0 += p1;		// Initial estimate of intersection of surfaces
			if (! child_pList->subdivNearSurface(s0.getX(), s0.getY(), s0.getZ(), &bp, &dist, &outside)) {
				cout << "Plist FAILURE 6!\n";
				destroy();
				return;
			}
			s0 = Point(bp.pnt[0], bp.pnt[1], bp.pnt[2])
				- p1;	// Better estimate of vector to intersection

			s1 = c1.cross(parentNormal);
			s1.normalize();
			s1 *= r;		// Multiply by approx distance to surface
			s1 += p1;		// Initial estimate of intersection of surfaces
			if (! child_pList->subdivNearSurface(s1.getX(), s1.getY(), s1.getZ(), &bp, &dist, &outside)) {
				cout << "Plist FAILURE 7!\n";
				destroy();
				return;
			}
			s1 = Point(bp.pnt[0], bp.pnt[1], bp.pnt[2])
				- p1;	// Better estimate of vector to intersection

			s0 += s1;
			s0 *= 0.5;	// s0 contains average displacement based on two connectors to atom n
			vertex[n + n] = p1 + s0;
			vertex[n + n + 1] = p1 - s0;

			// Refine the radii by use of the parent figure
			if (! parent_pList->subdivNearSurface(vertex[n + n].getX(), vertex[n + n].getY(),
					vertex[n + n].getZ(), &bp, &dist, &outside))
			{
				cout << "Plist FAILURE 8!\n";
				destroy();
				return;
			}
			if (radius[n + n] < bp.rad)		// Use the larger of the two radii
				radius[n + n] = bp.rad;
			if (! parent_pList->subdivNearSurface(vertex[n + n + 1].getX(), vertex[n + n + 1].getY(),
					vertex[n + n + 1].getZ(), &bp, &dist, &outside))
			{
				cout << "Plist FAILURE 9!\n";
				destroy();
				return;
			}
			if (radius[n + n + 1] < bp.rad)		// Use the larger of the two radii
				radius[n + n + 1] = bp.rad;

			atom0 = atom1;
			atom1 = atom2;
			p0 = p1;
			p1 = p2;
			c0 = c1;
		}

	// Last pair of blending region centers and radii
	if (! parent_pList->subdivBposition(&bp, link->u, link->v, link->t)) {
		cout << "Plist FAILURE 10!\n";
		destroy();
		return;
	}
	parentNormal = -Vector3D(bp.norm[0], bp.norm[1], bp.norm[2]);

	s1 = c0.cross(parentNormal);	// AGG: s1 ~= N; would it be better to use N?
	s1.normalize();
	s1 *= atom1->getR()*sin(atom1->getTheta());		// Multiply by approx distance to surface
	s1 += p1;		// Initial estimate of intersection of surfaces
	if (! child_pList->subdivNearSurface(s1.getX(), s1.getY(), s1.getZ(), &bp, &dist, &outside)) {
		cout << "Plist FAILURE 11!\n";
		destroy();
		return;
	}
	s1 = Point(bp.pnt[0], bp.pnt[1], bp.pnt[2]) - p1;	// Better estimate of vector to surface

	vertex[nlinks + nlinks - 2] = p1 + s1;
	vertex[nlinks + nlinks - 1] = p1 - s1;

	// Refine the radii by use of the parent figure
	if (! parent_pList->subdivNearSurface(vertex[nlinks + nlinks - 2].getX(),
			vertex[nlinks + nlinks - 2].getY(), vertex[nlinks + nlinks - 2].getZ(),
			&bp, &dist, &outside))
	{
		cout << "Plist FAILURE 12!\n";
		destroy();
		return;
	}
	if (radius[nlinks + nlinks - 2] < bp.rad)		// Use the larger of the two radii
		radius[nlinks + nlinks - 2] = bp.rad;
	if (! parent_pList->subdivNearSurface(vertex[nlinks + nlinks - 1].getX(),
			vertex[nlinks + nlinks - 1].getY(), vertex[nlinks + nlinks - 1].getZ(),
			&bp, &dist, &outside))
	{
		cout << "Plist FAILURE 13!\n";
		destroy();
		return;
	}
	if (radius[nlinks + nlinks - 1] < bp.rad)		// Use the larger of the two radii
		radius[nlinks + nlinks - 1] = bp.rad;

	factorSq = 1.0;
}

void Intersection::destroy() {
#ifdef DEBUG
	cout << "Intersection::destroy()\n";
#endif
	delete [] vertex;
	delete [] radius;
	vertex = NULL;
	radius = NULL;
	nlinks = 0;
}

Intersection::~Intersection() {
#ifdef DEBUG
	cout << "Intersection::~Intersection()\n";
#endif
	delete [] vertex;
	delete [] radius;
}

Intersection::Intersection(const Intersection & isc) {
#ifdef DEBUG
	cout << "Intersection::Intersection(const Intersection &)\n";
#endif
	nlinks = isc.nlinks;
	vertex = new Point[nlinks + nlinks];
	radius = new double[nlinks + nlinks];
	for (int n = 0; n < nlinks; n++) {
		vertex[n] = isc.vertex[n];
		vertex[n + nlinks] = isc.vertex[n + nlinks];	// Vertex contains 2*nlinks items
		radius[n] = isc.radius[n];
		radius[n] = isc.radius[n + nlinks];		// Radius contains 2*nlinks items
	}
	factorSq = isc.factorSq;
}

#if 0
// This version seems wrong and leads to a warning. -- Foskey
Intersection & Intersection::operator= (const Intersection & isc) {
#ifdef DEBUG
	cout << "Intersection::operator=()\n";
#endif
	delete [] vertex;
	delete [] radius;
    // This next line is going to call operator=, which is the 
    // function we are defining.  Seems wrong.
	*this = isc;  
	return *this;
}
#endif

Intersection & Intersection::operator= (const Intersection & isc) {
#ifdef DEBUG
	cout << "Intersection::operator=()\n";
#endif
	delete [] vertex;
	delete [] radius;
	nlinks = isc.nlinks;
	vertex = new Point[nlinks + nlinks];
	radius = new double[nlinks + nlinks];
	for (int n = 0; n < nlinks; n++) {
		vertex[n] = isc.vertex[n];
		vertex[n + nlinks] = isc.vertex[n + nlinks];	// Vertex contains 2*nlinks items
		radius[n] = isc.radius[n];
		radius[n] = isc.radius[n + nlinks];		// Radius contains 2*nlinks items
	}
	factorSq = isc.factorSq;
	return *this;
}

// Test a point against a single bounding region
inline bool Intersection::test(const Point & p, int m, int n)
{
	double t, r;
	Point q;
	double a, b, c;
	double radsq;

	q = project_point_to_line(p, vertex[n], vertex[m], t);
	// The cone ends are rounded
	if (t < 0.0)
		r = radius[n];
	else if (t > 1.0)
		r = radius[m];
	else
		r = t*radius[m] + (1.0 - t)*radius[n];

	a = p.getX() - q.getX();
	b = p.getY() - q.getY();
	c = p.getZ() - q.getZ();
	radsq = a*a + b*b + c*c;

	if (radsq <= r*r*factorSq)
		return true;	// Inside region
	else
		return false;	// Outside region
}

// Make two tests per atom
bool Intersection::testPoint(const Point & p, int n)
{
	if (test(p, n + n, n + n + 2))
		return true;	// Inside upper region

	if (test(p, n + n + 1, n + n + 3))
		return true;	// Inside lower region
	else
		return false;	// Outside both regions
}

bool Intersection::testPoint(const Point & p)
{
	// Check the inter-atom bounding regions
	for (int n = 0; n < nlinks - 1; n++) {
		if (testPoint(p, n))
			return true;
	}

#ifdef CHECK_ENDS
	// Check the two intra-atom bounding regions on the ends
	if (test(p, 0, 1))
		return true;
	if (test(p, nlinks + nlinks - 2, nlinks + nlinks - 1))
		return true;
#endif
	return false;
}

// Debugging function
void Intersection::print() {
	cout << "Intersection:\n  mainFigureId = " << mainFigureId
		<< ", subFigureId = " << subFigureId << '\n';
	cout << "  Blending boxes:\n";
	for (int n = 0; n < nlinks + nlinks - 2; n++) {
		cout << "  " << n << ". centers: (" << vertex[n].getX() << ", " << vertex[n].getY()
			<< ", " << vertex[n].getZ() << "), (" << vertex[n + 2].getX() << ", "
			<< vertex[n + 2].getY() << ", " << vertex[n + 2].getZ() << ")\n";
		cout << "\tradii = " << radius[n] << " and " << radius[n + 2] << '\n';
	}
	cout << "  multiplier = " << sqrt(factorSq) << endl;
}

// Debugging function
void Intersection::render(double currentScale) {
    double scale, length;
	int n;
	Vector3D arrow;

	float color[4] = { 0.5f, 0.5f, 0.5f, 0.4f };
	float lineColor[3] = { 1.0f, 1.0f, 1.0f };
	const double degrees_per_radian = 90.0/asin(1.0);

	glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LIGHTING);
	scale = LARGE_POINT_SIZE/currentScale;	

	GLUquadric * qCyl = gluNewQuadric();
#ifdef DRAW_ENDS
	GLUquadric * qSph = gluNewQuadric();
#endif

	glColor4fv(color);

	// Draw a semi-transparent surface for the truncated, bounding cones
	// The end caps are omitted

#ifdef CHECK_ENDS
	// Start with the first end region
	arrow = vertex[1] - vertex[0];
	length = arrow.norm();

	// Find the unit vector that bisects arrow and the Z axis
	arrow /= length;
	arrow.setZ(arrow.getZ() + 1.0);
	arrow *= 0.5;

	glPushMatrix();
	glTranslated(vertex[0].getX(), vertex[0].getY(), vertex[0].getZ());
	glRotated(180.0, arrow.getX(), arrow.getY(), arrow.getZ());
	glScaled(scale, scale, scale);
	length /= scale;

	gluCylinder(qCyl, radius[1]*sqrt(factorSq)/scale,
		radius[0]*sqrt(factorSq)/scale, length, 12, 1);
	glPopMatrix();
#endif	/* CHECK_ENDS */

	// Process the side pairs
	for (n = 0; n < nlinks + nlinks - 2; n++) {
		arrow = vertex[n + 2] - vertex[n];
		length = arrow.norm();

		// Find the unit vector that bisects arrow and the Z axis
		arrow /= length;
		arrow.setZ(arrow.getZ() + 1.0);
		arrow *= 0.5;

		glPushMatrix();
		glTranslated(vertex[n].getX(), vertex[n].getY(), vertex[n].getZ());
		glRotated(180.0, arrow.getX(), arrow.getY(), arrow.getZ());
		glScaled(scale, scale, scale);
		length /= scale;

#ifdef DRAW_ENDS
		if (n == 0 || n == 1)
			gluSphere(qSph, radius[n]*sqrt(factorSq)/scale, 12, 12);
#endif
		gluCylinder(qCyl, radius[n + 2]*sqrt(factorSq)/scale,
			radius[n]*sqrt(factorSq)/scale, length, 12, 1);
#ifdef DRAW_ENDS
		if (n == nlinks + nlinks - 4 || n == nlinks + nlinks - 3) {
			glPopMatrix();
			glPushMatrix();
			glTranslated(vertex[n + 2].getX(), vertex[n + 2].getY(), vertex[n + 2].getZ());
			glRotated(180.0, arrow.getX(), arrow.getY(), arrow.getZ());
			glScaled(scale, scale, scale);
			gluSphere(qSph, radius[n + 2]*sqrt(factorSq)/scale, 12, 12);
		}
#endif
		glPopMatrix();
	}

#ifdef CHECK_ENDS
	// Process the remaining end region
	arrow = vertex[nlinks + nlinks - 2] - vertex[nlinks + nlinks - 1];
	length = arrow.norm();

	// Find the unit vector that bisects arrow and the Z axis
	arrow /= length;
	arrow.setZ(arrow.getZ() + 1.0);
	arrow *= 0.5;

	glPushMatrix();
	glTranslated(vertex[nlinks + nlinks - 1].getX(), vertex[nlinks + nlinks - 1].getY(),
		vertex[nlinks + nlinks - 1].getZ());
	glRotated(180.0, arrow.getX(), arrow.getY(), arrow.getZ());
	glScaled(scale, scale, scale);
	length /= scale;


	gluCylinder(qCyl, radius[nlinks + nlinks - 2]*sqrt(factorSq)/scale,
		radius[nlinks + nlinks - 1]*sqrt(factorSq)/scale, length, 12, 1);
	glPopMatrix();
#endif	/* CHECK_ENDS */

	gluDeleteQuadric(qCyl);
#ifdef DRAW_ENDS
	gluDeleteQuadric(qSph);
#endif
    glDisable(GL_LIGHTING);
	glDisable(GL_BLEND);

	// Draw the center axes of the blending regions
	glLineWidth(3.0f);
	glBegin(GL_LINES);
	glColor3fv(lineColor);

	for (n = 0; n < nlinks + nlinks - 2; n++) {
	    glVertex3d(vertex[n].getX(), vertex[n].getY(), vertex[n].getZ());
	    glVertex3d(vertex[n + 2].getX(), vertex[n + 2].getY(), vertex[n + 2].getZ());
	}
	glEnd();
	glLineWidth(1.0f);

}


