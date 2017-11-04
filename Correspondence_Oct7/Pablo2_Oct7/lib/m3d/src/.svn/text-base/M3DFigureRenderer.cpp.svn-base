#include "M3DFigureRenderer.h"
#include "M3DDisplayGlobals.h"
#include "M3DObjectRenderer.h"
#include <vector>

#include <typeinfo>

using namespace std;

static
int atomIsLandmark(M3DFigure *fig, int atomIndex)
{
	for (int lm = 0; lm < fig->getLandmarkCount(); lm++)
	{
		if (fig->getLandmarkAtomIndex(lm) == atomIndex)
			return lm;
	}
	return -1;
}

// return a vector of landmark indices that involve this atom.
void findAtomsLandmarkIndices(M3DFigure *fig, int atomIndex, vector<int> &lmIndices)
{
	lmIndices.clear();
	int nLm = fig->getLandmarkCount();
	for (int lm = 0; lm < nLm; lm++)
	{
		if (fig->getLandmarkAtomIndex(lm) == atomIndex)
			lmIndices.push_back(lm);
	}
}

M3DFigureRenderer::M3DFigureRenderer(M3DFigure* _fig)
{
	fig	= _fig;
    primitiveRenderer = new M3DPrimitiveRenderer;
    startingSelectionId = IGNORE_SELECTION_ID;

    meshLineType = M3DObjectRenderer::meshConnectorsType();
    lineWidth = M3DObjectRenderer::meshConnectorsLineWidth();

	// capture tipNames & tipTs - handle fixed size of tipTs.
	// LATER: replace tip* arrays with vector<>
    int numPrimitives = fig->getPrimitiveCount();
	int maxPrim = sizeof(tipTs) / sizeof(tipTs[0]);
	if (numPrimitives > maxPrim)
	{
		static int  warned = 0;
		numPrimitives = maxPrim;
		if (!warned)
			cout << "WARNING: only drawing first " << maxPrim << " atoms" << endl;
		warned = 1;
	}

    for (int atom = 0; atom < numPrimitives; atom++)
	{
		vector<int> lmIndices;
		findAtomsLandmarkIndices(fig, atom, lmIndices);
		int nlm = lmIndices.size();
		if (nlm > 0) {
			vector<int>::iterator lm;
			int tCount = 0;
			for (lm = lmIndices.begin(); lm != lmIndices.end(); lm++)
			{
				tipNames[atom][tCount] = (char *) fig->getLandmarkName (*lm);
				tipTs[atom][tCount++] = fig->getLandmarkAtomT(*lm);
			}
			tipNames[atom][tCount] = NULL;	// null-terminate it
			tipTs[atom][tCount] = 0;
		}
		else {
			tipNames[atom][0] = NULL;
			tipTs[atom][0] = 0;
		}
	}
}

M3DFigureRenderer::~M3DFigureRenderer()
{
    delete primitiveRenderer;
}

void M3DFigureRenderer::setColor(const float * color)
{
	// Compute max and average of color component values
	float avg = 0.0;
	float max = 0.0;
	float diff;
	int i, n;

//	std::cout << "Original color: " << color[0] << ", " << color[1] << ", " << color[2] << '\n';
	for (i = 0; i < 3; i++) {
		rgb[i] = color[i];
		avg += rgb[i];
		if (rgb[i] > max)
			max = rgb[i];
	}
	avg /= 3.0;
	n = 0;	// Number of components greater than avg
	for (i = 0; i < 3; i++)
		if (rgb[i] > avg)
			n++;
	diff = max - avg;
	if (n < 1)
		diff = -diff;

	// Adjust color components
	for (i = 0; i < 3; i++) {
		rgb[i] += diff;
		if (rgb[i] > 1.0)
			rgb[i] = 1.0;
	}
//	std::cout << "Final color: " << rgb[0] << ", " << rgb[1] << ", " << rgb[2] << '\n';
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

void M3DFigureRenderer::renderFigureName()
{
    M3DPrimitive * centerAtom;
    Vector3D point;
    Vector3D hub;
    int numRows, numColumns;
	const char * label;

    if (fig == NULL)
        return;

//	cout << "M3DFigureRenderer::renderFigureName()\n";

	if (fig->getVisibility() == false)
		return;

	label = fig->getName();
	if (label) {

		// Draw the figure name near the center of the figure
		if (typeid(*fig) == typeid(M3DQuadFigure)) {
			const M3DQuadFigure * figure = dynamic_cast<M3DQuadFigure*>(fig);
			numRows = figure->getRowCount();
			numColumns = figure->getColumnCount();
			centerAtom = figure->getPrimitivePtr(numRows/2, numColumns/2);
		}
		else {
			const M3DTubeFigure * figure = dynamic_cast<M3DTubeFigure*>(fig);
			numColumns = figure->getColumnCount();
			centerAtom = figure->getPrimitivePtr(numColumns/2);
		}

		glColor3f(1.0f, 1.0f, 1.0f);

		// Note: drawing figure names is a many-body problem, except in the
		// case of having a single figure.  No effort is made here to
		// consider other figures when rendering the label for a given figure.
		if (centerAtom) {
			hub = centerAtom->getX();
			Vector3D hubscr = modelToScreen(hub);
//			cout << "hub (screen) = "; hubscr.print();
			point = centerAtom->getY0();
			point *= 1.25;	// Offset from the surface by 0.25*r
			point += hub;
//			cout << "point = "; point.print();
			Vector3D scrPos = modelToScreen(point);
//			cout << "scrPos = "; scrPos.print();

			// Heuristically adjust the drawing position for best results.
			// Choosing the spoke that has the smaller Z value did not
			// work as well as doing it this way, although this method
			// sometimes allows the label to go behind the figure.
			if (scrPos.getX() > hubscr.getX()) {
				// Draw from the point toward the screen right
				if (scrPos.getZ() > hubscr.getZ()) {
					// Move label from behind figure
					scrPos.setZ(hubscr.getZ() - (scrPos.getZ() - hubscr.getZ()));
					point = screenToModel(scrPos);
//					cout << "adjusted scrPos = "; scrPos.print();
//					cout << "adjusted point = "; point.print();
				}
			}
			else {
				// Draw from screen left toward the point
				int sw = M3DObjectRenderer::stringDrawingWidth(label);
				scrPos.setX(scrPos.getX() - sw);
				if (scrPos.getZ() > hubscr.getZ()) {
					// Move label from behind figure
					scrPos.setZ(hubscr.getZ() - (scrPos.getZ() - hubscr.getZ()));
					point = screenToModel(scrPos);
//					cout << "offset & adjusted scrPos = "; scrPos.print();
//					cout << "offset & adjusted point = "; point.print();
				}
				else {
					point = screenToModel(scrPos);
//					cout << "offset scrPos = "; scrPos.print();
//					cout << "offset point = "; point.print();
				}
			}

			// Draw the figure's name
			glPushMatrix();
			glTranslated(point.getX(), point.getY(), point.getZ());
			glRasterPos3i(0, 0, 0);
			M3DObjectRenderer::drawString(label);
			glPopMatrix();
		}
	}
}

