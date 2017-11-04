#include <math.h>
#include <FL/glut.H>

#include "M3DFigure.h"
#include "M3DQuadFigureRenderer.h"
#include "M3DTubeFigureRenderer.h"
#include "ControlParms.h"
#include "renderDefinitions.h"
#include "M3DDisplayGlobals.h"
#include "M3DObjectRenderer.h"

#include <typeinfo>

using namespace std;


// Note: static variables are used for the following controls because
// this permits the values to be accessed directly by the classes
// that do the actual rendering.

// Parameters used to control visibility
bool M3DObjectRenderer::mesh_visibility;
bool M3DObjectRenderer::atom_vector_visibility;
bool M3DObjectRenderer::atom_visibility;

// Parameter used for rendering the medial mesh
bool M3DObjectRenderer::mesh_connectors_type;
float M3DObjectRenderer::mesh_connectors_line_width;
float M3DObjectRenderer::mesh_connectors_color[3];

// Parameters used for rendering primitives
bool M3DObjectRenderer::extra_atom_vectors;
bool M3DObjectRenderer::atom_vectors_type;
int M3DObjectRenderer::b_vectors_type;
float M3DObjectRenderer::atom_vectors_line_width;
bool M3DObjectRenderer::use_figure_color;

// Flag for rendering interfigural constraints
bool M3DObjectRenderer::draw_constraint_arrows;

// Flag for rendering model landmark tips
bool M3DObjectRenderer::landmark_visibility;


M3DObjectRenderer::M3DObjectRenderer(M3DObject * obj)
{
    object = NULL;
    currentScale = 1.0;
    markedPrimitiveId = INVALID_PRIMITIVE_ID;
    markedPrimitiveRenderer = NULL;
    setObjectPtr(obj);

	mesh_connectors_type = globalControl->readBool(MeshConnectorsType);
	mesh_connectors_line_width = (float) globalControl->readInt(MeshConnectorsLineWidth);
	if (mesh_connectors_line_width > MESH_CONNECTORS_WIDTH_MAX)
		mesh_connectors_line_width = MESH_CONNECTORS_WIDTH_MAX;
	mesh_connectors_color[0] = (float) globalControl->readDouble(ConnectorsColor_R);
	mesh_connectors_color[1] = (float) globalControl->readDouble(ConnectorsColor_G);
	mesh_connectors_color[2] = (float) globalControl->readDouble(ConnectorsColor_B);
	mesh_visibility = globalControl->readBool(ShowMeshConnectors);
	extra_atom_vectors = globalControl->readBool(ExtraAtomVectors);
	atom_vectors_type = globalControl->readBool(AtomVectorsType);
	b_vectors_type = globalControl->readBool(BVectorsType);
	atom_vectors_line_width = (float) globalControl->readInt(AtomVectorsLineWidth);
	use_figure_color = globalControl->readBool(FiguresColorAtoms);
	if (atom_vectors_line_width > ATOM_VECTORS_WIDTH_MAX)
		atom_vectors_line_width = ATOM_VECTORS_WIDTH_MAX;
	landmark_visibility = globalControl->readBool(ShowLandmarks);
}

M3DObjectRenderer::~M3DObjectRenderer()
{
    if (figureRenderers.size() > 0) {
        for (int i = 0; i < figureRenderers.size(); i++) {
            delete figureRenderers[i];
            figureRenderers[i] = NULL;
        }
        figureRenderers.clear();
    }
}

void M3DObjectRenderer::setObjectPtr(M3DObject * obj)
{
    int numFigures;
    int selectionId;
    int i;

    M3DFigure * figurePtr;

    M3DFigureRenderer * figureRendererPtr;

    object = obj;

	figureRenderers.erase( figureRenderers.begin(), figureRenderers.end() );

    if(object == NULL)
        return;

    numFigures = object->getFigureCount();
    selectionId = 0;
    for(i = 0; i < numFigures; i++)
    {
        figurePtr = object->getFigurePtr(i);
        if(figurePtr == NULL)
            continue;

		if( typeid(*figurePtr) == typeid(M3DQuadFigure) ) {
            figureRendererPtr = 
                new M3DQuadFigureRenderer(dynamic_cast<M3DQuadFigure*>(figurePtr));
            figureRendererPtr->setStartingSelectionId(selectionId);
            selectionId += figurePtr->getPrimitiveCount();
            figureRendererPtr->setColor(figurePtr->getColor());
            figureRenderers.push_back(figureRendererPtr);
        }
		else if( typeid(*figurePtr) == typeid(M3DTubeFigure) ) {
            figureRendererPtr = 
                new M3DTubeFigureRenderer(dynamic_cast<M3DTubeFigure*>(figurePtr));
            figureRendererPtr->setStartingSelectionId(selectionId);
            selectionId += figurePtr->getPrimitiveCount();
            figureRendererPtr->setColor(figurePtr->getColor());
            figureRenderers.push_back(figureRendererPtr);
		}
		else {
			// FIXME:
			// show an error
			std::cout << typeid(figurePtr).name() << std::endl;
        }
    }
}

void M3DObjectRenderer::setCurrentScale(double scale)
{
    int numFigures;
    int i;

    M3DFigureRenderer * figureRendererPtr;

    double pointScaleFactor;

    if(scale == 0.0)
        return;

    currentScale = scale;
    pointScaleFactor = 1.0 / scale;

    numFigures = figureRenderers.size();
    for(i = 0; i < numFigures; i++)
    {
        figureRendererPtr = figureRenderers[i];
        if(figureRendererPtr != NULL)
            figureRendererPtr->setPointScaleFactor(pointScaleFactor);
    }
}

// A useful routine for testing.  This places a point at the
// current OpenGL model view origin.
void drawPoint()
{
    glEnable(GL_LIGHTING);
	GLUquadric * q = gluNewQuadric();
	gluSphere(q, 1.0, 8, 8);
	gluDeleteQuadric(q);
    glDisable(GL_LIGHTING);
}

void M3DObjectRenderer::renderConstraintArrows() {
	// FIXME: Needs a workover
	M3DPrimitive * atom;
	int i, j;
	int u, v;
    double scale, length;
	const float * color;
	const double degrees_per_radian = 90.0/asin(1.0);

    glEnable(GL_LIGHTING);
	scale = LARGE_POINT_SIZE/currentScale;	

	GLUquadric * qCyl = gluNewQuadric();
	for (i = 0; i < object->getFigureCount(); i++) {
		M3DQuadFigure* figure = dynamic_cast<M3DQuadFigure*>(object->getFigurePtr(i));
		color = figure->getColor();
		InterfiguralConstraints constraints = figure->constraints();

		if (constraints.size() > 0) {
			glColor4f(color[0], color[1], color[2], 0.4f);
//			std::cout << "\nColor of figure " << i << " = (" << color[0] << ", "
//				<< color[1] << ", " << color[2] << ")\n";

			// Determine a point in the figure to serve as the base of the arrow.
			// The center of gravity is averaged with the position of a central
			// atom.  This is done because the cog could be outside the figure
			// and the central atom may be near an end of the figure.
			Vector3D cog = figure->getCOG(true);
			u = figure->getRowCount()/2;
			v = figure->getColumnCount()/2;
			atom = figure->getPrimitivePtr(u, v);
			cog += atom->getX();
			cog *= 0.5;

			for (j = 0; j < constraints.size(); j++) {
				int governedId = constraints.figure(j);

				// Almost the same computation as above, but for the governed figure.
				// The if's are used to provide a slight offset, so arrows pointing
				// in oposite directions aren't colinear.
				M3DQuadFigure* governed = dynamic_cast<M3DQuadFigure*>(object->getFigurePtr(governedId));
				Vector3D governed_cog = governed->getCOG(true);
				u = governed->getRowCount();
				if (u%2 != 0)
					u++;
				u /= 2;
				v = governed->getColumnCount();
				if (v%2 != 0)
					v++;
				v /= 2;
				atom = governed->getPrimitivePtr(u, v);
				governed_cog += atom->getX();
				governed_cog *= 0.5;

				Vector3D arrow = governed_cog - cog;
//				std::cout << "Governor at "; cog.print();
//				std::cout << "Governed figure " << governedId << " at "; governed_cog.print();
//				std::cout << "Arrow = "; arrow.print();
				length = arrow.norm();

				// Find the unit vector that bisects arrow and the Z axis
				arrow /= length;
				arrow.setZ(arrow.getZ() + 1.0);
				arrow *= 0.5;

				glPushMatrix();
				glTranslated(cog.getX(), cog.getY(), cog.getZ());
				glRotated(180.0, arrow.getX(), arrow.getY(), arrow.getZ());
				glScaled(scale, scale, scale);
				length /= scale;

				gluCylinder(qCyl, 2.0, 0.0, length, 12, 1);
				glPopMatrix();
			}
		}
	}
	gluDeleteQuadric(qCyl);
    glDisable(GL_LIGHTING);
}

void M3DObjectRenderer::render()
{
    int numFigures;
    int i;
    M3DFigureRenderer * figureRendererPtr;
    M3DPrimitive * markedPrimitivePtr;
	M3DFigure * markedFigure;

    numFigures = figureRenderers.size();
	if (numFigures == 0)
		return;

    for(i = 0; i < numFigures; i++)
    {
        figureRendererPtr = figureRenderers[i];
        if(figureRendererPtr != NULL)
            figureRendererPtr->render();
    }

    i = object->getFigureAtomIndexes(markedPrimitiveId);

#ifdef BINARY
	if (markedPrimitiveId <= INVALID_PRIMITIVE_ID)
		return;
#endif

	markedFigure = object->getFigurePtr(i);
	if (markedFigure->getVisibility()) {

		markedPrimitivePtr = markedFigure->getPrimitivePtr(markedPrimitiveId);
		if (markedPrimitivePtr != NULL) {
			markedPrimitiveRenderer->setPrimitive(markedPrimitivePtr);

			if(currentScale != 0.0)
				markedPrimitiveRenderer->setPointScaleFactor(1.5 / currentScale);

			markedPrimitiveRenderer->render();
		}
	}

	if (constraintArrows())
		renderConstraintArrows();
}

void M3DObjectRenderer::renderLandmarks()
{
    int numFigures;
    int i;
    M3DFigureRenderer * figureRendererPtr;

    numFigures = figureRenderers.size();
	if (numFigures == 0)
		return;

    for(i = 0; i < numFigures; i++)
    {
        figureRendererPtr = figureRenderers[i];
        if(figureRendererPtr != NULL)
            figureRendererPtr->renderLandmarks();
    }
}

void M3DObjectRenderer::renderFigureNames()
{
    int numFigures;
    int i;
    M3DFigureRenderer * figureRendererPtr;

    numFigures = figureRenderers.size();
	if (numFigures == 0)
		return;

    for(i = 0; i < numFigures; i++)
    {
        figureRendererPtr = figureRenderers[i];
        if(figureRendererPtr != NULL)
            figureRendererPtr->renderFigureName();
    }
}

// Draw a string of characters in an openGL window at the last position
// specified by glRasterPos.
int M3DObjectRenderer::stringDrawingWidth(const char * s)
{
	int w = 0;
	int i = 0;
	while (s[i]) {
		//w += glutBitmapWidth(GLUT_BITMAP_TIMES_ROMAN_24, s[i]);
		i++;
	}
	return w;
}

// Draw a string of characters in an openGL window at the last position
// specified by glRasterPos.
void M3DObjectRenderer::drawString(const char * s)
{
	int i = 0;
	while (s[i]) {
		//glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, s[i]);
		i++;
	}
}

