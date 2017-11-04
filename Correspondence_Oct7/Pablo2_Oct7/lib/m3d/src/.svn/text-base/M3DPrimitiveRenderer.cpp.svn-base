#include <stdio.h>
#include <math.h>

#include "M3DPrimitiveRenderer.h"
#include "M3DDisplayGlobals.h"
#include "M3DObjectRenderer.h"


const float DEFAULT_POINT_COLOR[3] = { 1.0f, 1.0f, 1.0f };
const float DEFAULT_SELECTED_POINT_COLOR[3] = { 1.0f, 1.0f, 0.0f };
// Same as DEFAULT_LANDMARK_COLOR in P3DView.cpp
const float DEFAULT_TIP_POINT_COLOR[3] = { 1.0f, 165.0f/255.0f, 0.0f  };


void M3DPrimitiveRenderer::drawPoint()
{
    glEnable(GL_LIGHTING);
    GLUquadric * q = gluNewQuadric();
    gluSphere(q, 1.0, 8, 8);
    gluDeleteQuadric(q);
    glDisable(GL_LIGHTING);
}

M3DPrimitiveRenderer::M3DPrimitiveRenderer()
{
    primitive = NULL;

    selectionId = IGNORE_SELECTION_ID;
    pointScaleFactor = 1.0;
    axesOn = M3DObjectRenderer::showExtraAtomVectors();
    fancyAxes = M3DObjectRenderer::atomVectorsType();
	drawBAxes = M3DObjectRenderer::atomB_VectorsType();
    lineWidth = M3DObjectRenderer::atomVectorsLineWidth();
    pointColor[0] = DEFAULT_POINT_COLOR[0];
    pointColor[1] = DEFAULT_POINT_COLOR[1];
    pointColor[2] = DEFAULT_POINT_COLOR[2];

    selectedPointColor[0] = DEFAULT_SELECTED_POINT_COLOR[0];
    selectedPointColor[1] = DEFAULT_SELECTED_POINT_COLOR[1];
    selectedPointColor[2] = DEFAULT_SELECTED_POINT_COLOR[2];

	tipName = NULL;
	tipT = 0;	// default is b-spoke
	tipV = 0;
    tipPointColor[0] = DEFAULT_TIP_POINT_COLOR[0];
    tipPointColor[1] = DEFAULT_TIP_POINT_COLOR[1];
    tipPointColor[2] = DEFAULT_TIP_POINT_COLOR[2];
}

M3DPrimitiveRenderer::~M3DPrimitiveRenderer()
{
}

// dibyendu - modified the render() function on 06/16/2011 to facilitate s-rep display

void M3DPrimitiveRenderer::render()
{
	Vector3D center;
	Vector3D long_b;
	Vector3D b;
	Vector3D n;
    Vector3D bPerp;
	Vector3D y0;
	Vector3D y1;
	Vector3D yn;

	// dibyendu
	Vector3D yEnd ;

    double pointSize;
    double scale;

    if (primitive == NULL)
		return;

	center = primitive->getX();		// Atom position

	if (M3DObjectRenderer::atomVectorVisibility()) {

		// Compute coordinates of vectors/axes
		if (fancyAxes)
			b = primitive->getR()*cos(primitive->getTheta())*primitive->getB();
		// Long_b is R proportional and for end atoms also eta-elongated
		long_b = primitive->getExtendedB();
		n = primitive->getN();
		bPerp = primitive->getBPerp();
		y0 = primitive->getY0() + primitive->getX();
		y1 = primitive->getY1() + primitive->getX();

		// dibyendu - get yEnd for quad primitive
		if( dynamic_cast<M3DQuadEndPrimitive*>(primitive) != NULL )
			yEnd = (dynamic_cast<M3DQuadEndPrimitive*>(primitive))->getYEnd() + primitive->getX() ;

		glLineWidth(lineWidth);
		glBegin(GL_LINES);

		// Draw optional axes
		if (axesOn == true)
		{
			// Draw N axis
			glColor3f(0.0, 1.0, 0.0);
			glVertex3d(center.getX(), center.getY(), center.getZ());
			n = n * primitive->getR() + center;
			glVertex3d(n.getX(), n.getY(), n.getZ());

			// Draw B-Perp axis
			glColor3f(0.0, 0.0, 1.0);
			glVertex3d(center.getX(), center.getY(), center.getZ());
			bPerp = bPerp * primitive->getR() + center;
			glVertex3d(bPerp.getX(), bPerp.getY(), bPerp.getZ());
		}

		// Draw all 8 vectors if tube primitive
		M3DTubePrimitive* tubePrimitive	= dynamic_cast<M3DTubePrimitive*>(primitive);
		if(tubePrimitive != NULL) {
			const int N	= tubePrimitive->getNumberOfSpokes();
			const float hn	= N/2.0f;
			for (int i = 0; i != N; ++i) {
				yn	= tubePrimitive->getYN(i) + tubePrimitive->getX();
				glColor3f( 1.0f - fabsf(i-hn)/hn, fabsf(i-hn)/hn, 1.0f);
				glVertex3d(center.getX(), center.getY(), center.getZ());
				glVertex3d(yn.getX(), yn.getY(), yn.getZ());
			}
		}
		else {
			// Draw Y0 vector
			glColor3f(0.0f, 1.0f, 1.0f);
			glVertex3d(center.getX(), center.getY(), center.getZ());
			glVertex3d(y0.getX(), y0.getY(), y0.getZ());

			// Draw Y1 vector
			glColor3f(1.0f, 0.0f, 1.0f);
			glVertex3d(center.getX(), center.getY(), center.getZ());
			glVertex3d(y1.getX(), y1.getY(), y1.getZ());

			// dibyendu - if it is an end atom, draw yEnd
			bool drawEnd = 1 ;
			if( drawEnd ) {
				if( dynamic_cast<M3DQuadEndPrimitive*>(primitive) != NULL ) {
					glColor3f(1.0, 0.0, 0.0);
					glVertex3d(center.getX(), center.getY(), center.getZ());
					glVertex3d(yEnd.getX(), yEnd.getY(), yEnd.getZ());
				}
			}
		}	

		glEnd() ;
 
		// -------------------------------------------------------------------------------------------------------------
		// dibyendu - commented out the fancyAxes part because it doesnt make sense to show the axes notation in s-reps
		// -------------------------------------------------------------------------------------------------------------
		
		//glLineWidth(lineWidth);
		//glBegin(GL_LINES);

		//if (fancyAxes) {

		//	if (drawBAxes == 1 || (drawBAxes == 2 && primitive->type() == M3D_END_PRIMITIVE)) {
		//		// Draw B axis
		//		glColor3f(1.0, 0.0, 0.0);
		//		glVertex3d(center.getX(), center.getY(), center.getZ());
		//		b = b + center;
		//		glVertex3d(b.getX(), b.getY(), b.getZ());
		//	}
		//	glEnd();

		//	if (drawBAxes == 1 || (drawBAxes > 1 && primitive->type() == M3D_END_PRIMITIVE)) {

		//		// Draw B axis extension in the case of end atoms
		//		glEnable(GL_LINE_STIPPLE);
		//		glLineStipple(1, 0x0F0F);
		//		glBegin(GL_LINES);

		//		long_b = long_b + center;
		//		glColor3f(1.0, 0.0, 0.0);
		//		// Note, the order of the next two statements guarantees that the
		//		// stippling ends at the outer end of the vector.
		//		glVertex3d(long_b.getX(), long_b.getY(), long_b.getZ());
		//		glVertex3d(b.getX(), b.getY(), b.getZ());

		//		glEnd();
		//		glDisable(GL_LINE_STIPPLE);
		//	}
		//}
		//else {

		//	if (drawBAxes == 1 || (drawBAxes > 1 && primitive->type() == M3D_END_PRIMITIVE)) {
		//		// Draw B axis
		//		glColor3f(1.0, 0.0, 0.0);
		//		glVertex3d(center.getX(), center.getY(), center.getZ());
		//		long_b = long_b + center;
		//		glVertex3d(long_b.getX(), long_b.getY(), long_b.getZ());
		//		glEnd();
		//	}
		//	else
		//		glEnd();
		//}

		glLineWidth(1.0f);
	}

	if (M3DObjectRenderer::atomVisibility()) {

		if (primitive->isSelected()) {
			pointSize = LARGE_POINT_SIZE;
			if (! primitive->isSelectedForRegularity())
				pointSize += LARGE_POINT_SIZE;
			glColor3fv(selectedPointColor);
		}
		else {
			pointSize = SMALL_POINT_SIZE;
			if (! primitive->isSelectedForRegularity())
				pointSize += LARGE_POINT_SIZE;
			glColor3fv(pointColor);
		}

		glLoadName(selectionId);

		scale = pointSize * pointScaleFactor;
		glPushMatrix();
		glTranslated(center.getX(), center.getY(), center.getZ());
		glScaled(scale, scale, scale);

		drawPoint();

		glPopMatrix();
		glLoadName(IGNORE_SELECTION_ID);
	}
}

inline Vector3D screenToModel(Vector3D pos)
{
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	double a[3];

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);

	gluUnProject(pos.getX(), pos.getY(), pos.getZ(),
		modelview, projection, viewport, &a[0], &a[1], &a[2]);
	return Vector3D(a[0], a[1], a[2]);
}

inline Vector3D modelToScreen(Vector3D & pos)
{
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	double a[3];

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);

	gluProject(pos.getX(), pos.getY(), pos.getZ(), modelview, projection,
		viewport, &a[0], &a[1], &a[2]);
	return Vector3D(a[0], a[1], a[2]);
}

void M3DPrimitiveRenderer::renderLandmark()
{
    double scale;
    double pointSize;

//	cout << "M3DPrimitiveRenderer::renderLandmark()\n";
    if (primitive == NULL)
		return;

	if (! M3DObjectRenderer::landmark_visibility)
		return;

	// Draw B spoke end & text label if atom's tip is named
	if (tipName && *tipName)
	{
		Vector3D tip;		// A marked and named spoke-end

		Vector3D hub = primitive->getX();		// Atom position
		Vector3D hubscr = modelToScreen(hub);

		pointSize = LARGE_POINT_SIZE;
		scale = pointSize * pointScaleFactor;

		const M3DTubePrimitive * tprim =
			dynamic_cast<const M3DTubePrimitive*>(primitive);
		if (tprim != NULL) {
			if (tipV == -1.0)
				tip = tprim->getExtendedB();
			else
				tip = tprim->getYPhi(tipV * 2 * R_PI);
		}
		else {
			// vectors from center to each tip
			Vector3D y0 = primitive->getY0();
			Vector3D bi = primitive->getExtendedB();
			Vector3D y1 = primitive->getY1();
			if (tipT == 1.0) tip = y1;
			else if (tipT == -1.0) tip = y0;
			else if (tipT == 0.0) tip = bi;
			else if (tipT < 0) tip = bi + -tipT * (y0 - bi); 
			else tip = bi + tipT * (y1 - bi);
		}
		Vector3D spokeEnd = hub + tip;
		
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();

		// Move to the spoke end.
		// Note: tip is a vector, not a position.
		glTranslated(spokeEnd.getX(), spokeEnd.getY(), spokeEnd.getZ());
		glScaled(scale, scale, scale);
		// Now working in raster (screen pixel) coordinates

		glColor3fv(tipPointColor);
		drawPoint();	// Draw spoke-end ball for landmark

		// Offset the label from the surface by 0.25*r
		tip *= 0.25/scale;
		glTranslated(tip.getX(), tip.getY(), tip.getZ());

		// Assumes all characters are the same width
		int sw = M3DObjectRenderer::stringDrawingWidth(tipName);
		glRasterPos2d(0.0, 0.0);	// Set the raster origin at tip
		// Get the screen coordinates of tip
		double start[4];
		glGetDoublev(GL_CURRENT_RASTER_POSITION, start);
		Vector3D scrPos(start[0], start[1], start[2]);
//			cout << "scrPos = "; scrPos.print();

		// Should the label be on the left or right?
		if (scrPos.getX() < hubscr.getX()) {
			// Draw from screen left toward the tip point
			int sw = M3DObjectRenderer::stringDrawingWidth(tipName);
			scrPos.setX(scrPos.getX() - sw);
//				cout << "Moving string by scrPos = "; scrPos.print();
			Vector3D point = screenToModel(scrPos);
			glRasterPos3d(point.getX(), point.getY(), point.getZ());
		}

		M3DObjectRenderer::drawString(tipName);
//			cout << "Drew landmark " << tipName << " at "; tip.print();

		glPopMatrix();
	}
}

