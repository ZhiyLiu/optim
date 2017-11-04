#include "M3DTubeFigureRenderer.h"
#include "M3DDisplayGlobals.h"
#include "M3DObjectRenderer.h"
#include <vector>

using namespace std;

M3DTubeFigureRenderer::M3DTubeFigureRenderer(M3DTubeFigure * fig) : M3DFigureRenderer(fig)
{
	// capture tipNames & tipVs
    int numPrimitives = fig->getPrimitiveCount();
    for (int atom = 0; atom < numPrimitives; atom++)
	{
		vector<int> lmIndices;
		findAtomsLandmarkIndices(fig, atom, lmIndices);
		int nlm = lmIndices.size();
		if (nlm > 0)
		{
			vector<int>::iterator lm;
			int tCount = 0;
			for (lm = lmIndices.begin(); lm != lmIndices.end(); lm++)
			{
				tipVs   [atom][tCount++] = fig->getLandmarkAtomV(*lm);
			}
			tipVs[atom][tCount] = 0;
		}
		else
		{
			tipVs[atom][0] = 0;
		}
	}
}

M3DTubeFigureRenderer::~M3DTubeFigureRenderer()
{
}

void M3DTubeFigureRenderer::render()
{
    int numColumns,
        numPrimitives;
    int i, j;

    M3DPrimitive * primPtr1,
                 * primPtr2;

    Vector3D point1,
             point2;

    if(fig == NULL)
        return;

	const M3DTubeFigure* figure	= dynamic_cast<M3DTubeFigure*>(fig);
	if (figure->getVisibility() == false)
		return;

    numColumns = figure->getColumnCount();

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

		for (j = 0; j < numColumns; j++)
		{
			primPtr1 = figure->getPrimitivePtr(j);
			if (primPtr1 == NULL)
				continue;

			if (j < numColumns - 1)
			{
				primPtr2 = figure->getPrimitivePtr(j+1);

				if (primPtr2 != NULL)
				{
					point1 = primPtr1->getX();
					point2 = primPtr2->getX();

					glVertex3d(point1.getX(), point1.getY(), point1.getZ());
					glVertex3d(point2.getX(), point2.getY(), point2.getZ());
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
        if (primPtr1 == NULL)
            continue;

        primitiveRenderer->setPrimitive(primPtr1);
        if (startingSelectionId != IGNORE_SELECTION_ID)
            primitiveRenderer->setSelectionId(startingSelectionId + i);
        primitiveRenderer->render();
    }
}

void M3DTubeFigureRenderer::renderLandmarks()
{
    int numPrimitives;
    int i;
    M3DPrimitive * primPtr1;
	
//	cout << "M3DTubeFigureRenderer::renderLandmarks()\n";
#ifdef LM_METHOD_OBJ

    if (fig == NULL)
        return;

	const M3DTubeFigure * figure = dynamic_cast<M3DTubeFigure*>(fig);
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
			for (int tCount=0; tipNames[i][tCount] && tCount < 3; tCount++)
			{
				primitiveRenderer->setTipName(tipNames[i][tCount]);
				primitiveRenderer->setTipT(tipTs[i][tCount]);
				primitiveRenderer->setTipV(tipVs[i][tCount]);
				primitiveRenderer->renderLandmark();
			}
		}
    }
#endif
}

