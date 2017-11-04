#include "M3DQuadFigureRenderer.h"
#include "M3DObjectRenderer.h"
#include "M3DDisplayGlobals.h"

using namespace std;


M3DQuadFigureRenderer::M3DQuadFigureRenderer(M3DQuadFigure * fig) : M3DFigureRenderer(fig)
{
}

M3DQuadFigureRenderer::~M3DQuadFigureRenderer()
{
}

void M3DQuadFigureRenderer::render()
{
    int numRows,
        numColumns,
        numPrimitives;
    int i, j;

    M3DPrimitive * primPtr1,
                 * primPtr2;

    Vector3D point1,
             point2;

    if (fig == NULL)
        return;

	const M3DQuadFigure* figure	= dynamic_cast<M3DQuadFigure*>(fig);
	if (figure->getVisibility() == false)
		return;

	if (M3DObjectRenderer::meshVisibility()) {
		float * color;

		// Draw Connectors
		if (meshLineType) {
			glEnable(GL_LINE_STIPPLE);
			glLineStipple(1, 0xC0C0);
		}
		glLineWidth(lineWidth);
		glBegin(GL_LINES);
		color = M3DObjectRenderer::meshConnectorsColor();
		glColor3fv(color);

		numRows = figure->getRowCount();
		numColumns = figure->getColumnCount();

		for (i = 0; i < numRows; i++)
		{
			for (j = 0; j < numColumns; j++)
			{
				primPtr1 = figure->getPrimitivePtr(i, j);
				if (primPtr1 == NULL)
					continue;

				if (i < numRows - 1)
				{
					primPtr2 = figure->getPrimitivePtr(i+1, j);

					if (primPtr2 != NULL)
					{
						point1 = primPtr1->getX();
						point2 = primPtr2->getX();

						glVertex3d(point1.getX(), point1.getY(), point1.getZ());
						glVertex3d(point2.getX(), point2.getY(), point2.getZ());
					}
				}

				if (j < numColumns - 1)
				{
					primPtr2 = figure->getPrimitivePtr(i, j+1);

					if (primPtr2 != NULL)
					{
						point1 = primPtr1->getX();
						point2 = primPtr2->getX();

						glVertex3d(point1.getX(), point1.getY(), point1.getZ());
						glVertex3d(point2.getX(), point2.getY(), point2.getZ());
					}
				}
			}
		}

		glEnd();
		glLineWidth(1.0f);
		glDisable(GL_LINE_STIPPLE);
	}

    if (primitiveRenderer == NULL)
        return;

    numPrimitives = figure->getPrimitiveCount();
    primitiveRenderer->setSelectionId(IGNORE_SELECTION_ID);
	if (M3DObjectRenderer::useFigureColor()) {
		primitiveRenderer->setPointColor(rgb[0], rgb[1], rgb[2]);
		primitiveRenderer->setSelectedPointColor(rgb[0], rgb[1], rgb[2]);
	}
	// else use the default colors

    for (i = 0; i < numPrimitives; i++)
    {
        primPtr1 = figure->getPrimitivePtr(i);
        if(primPtr1 == NULL)
            continue;

        primitiveRenderer->setPrimitive(primPtr1);
        if (startingSelectionId != IGNORE_SELECTION_ID)
            primitiveRenderer->setSelectionId(startingSelectionId + i);
        primitiveRenderer->render();
    }
}

void M3DQuadFigureRenderer::renderLandmarks()
{
    int numPrimitives;
    int i;
    M3DPrimitive * primPtr1;

//	cout << "M3DQuadFigureRenderer::renderLandmarks()\n";
#ifdef LM_METHOD_OBJ

    if (fig == NULL)
        return;

	const M3DQuadFigure * figure = dynamic_cast<M3DQuadFigure*>(fig);
	if (figure->getVisibility() == false)
		return;

    if (primitiveRenderer == NULL)
        return;

    numPrimitives = figure->getPrimitiveCount();
    primitiveRenderer->setSelectionId(IGNORE_SELECTION_ID);
	if (M3DObjectRenderer::useFigureColor()) {
		primitiveRenderer->setPointColor(rgb[0], rgb[1], rgb[2]);
		primitiveRenderer->setSelectedPointColor(rgb[0], rgb[1], rgb[2]);
	}
	// else use the default colors

    for (i = 0; i < numPrimitives; i++)
    {
        primPtr1 = figure->getPrimitivePtr(i);
        if (primPtr1 == NULL)
            continue;

        primitiveRenderer->setPrimitive(primPtr1);
        if (startingSelectionId != IGNORE_SELECTION_ID)
            primitiveRenderer->setSelectionId(startingSelectionId + i);

		// Model spoke-end landmarks
		// Each atom can have up to 3 named spoke ends
		if (tipNames[i][0]) {
//			cout << "Landmark on atom " << i << endl;
			//  (this draws the atom up to 3 times - yuck!)
			for (int tCount = 0; tipNames[i][tCount] && tCount < 3; tCount++) {
				primitiveRenderer->setTipName(tipNames[i][tCount]);
				primitiveRenderer->setTipT(tipTs[i][tCount]);
				primitiveRenderer->renderLandmark();
			}
		}
    }
#endif
}

