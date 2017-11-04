
#include <iostream>
#include <queue>
#include <set>

#include "M3DObject.h"

#define D_DIATOM
#define D_XFERLIST
#define D_SAMPLESTATS
#define D_POINTLIST_SERVER2
#include "Shapedepend.h"

#include "SelectedPartialFigures.h"

#include "M3DObjectSurfaceRenderer.h"
#include "M3DObjectSurfaceVisualizer.h"
#include "renderDefinitions.h"


#include <typeinfo>
/********************************************************************************/
/* Modifications:																*/
/*		7. Dec 00 -- changed alphas to etas, as per other Seurat code			*/
/********************************************************************************/

//#define DEBUG
//#define QIONG_DEBUG_LINALG

M3DObjectSurfaceRenderer::M3DObjectSurfaceRenderer(M3DObject * obj, int level)
{
    object = obj;
    surfaceLevel = level;
	list_length = 0;
	list = NULL;
	lineWidth = 1;
	rows = NULL;
	cols = NULL;
	num_figures = 0;
}

void M3DObjectSurfaceRenderer::clean()
{
	if (list_length != 0) {
		for (int i = 0; i < list_length; i++)
			delete list[i];
		delete [] list;
	}
	list = NULL;
	list_length = 0;
}

void M3DObjectSurfaceRenderer::clean(int index)
{
	if (index >= list_length)
		return;
	delete list[index];
	list[index] = NULL;
}

M3DObjectSurfaceRenderer::~M3DObjectSurfaceRenderer()
{
	clean();
	delete [] rows;
	delete [] cols;
	// object does not belong to this class
}

/*	This supports fast updating in case the topology of an individual figure
	changes, even though at the present time, pablo cannot make such changes,
	since it does not support changing the number of rows/columns of figures.
	To detect a topology change, the numbers of rows and columns for each
	figure must simply be compared with those used in the previous rendering.

	Notice that looping through the figures' point lists here must be
	consistent with other similar loops in this class.
*/
void M3DObjectSurfaceRenderer::setObject(M3DObject * obj)
{
	if (obj != NULL) {
		if (num_figures == obj->getFigureCount()) {
			for (int i = 0; i < num_figures; i++) {
				if (typeid(*obj->getFigurePtr(i)) == typeid(M3DQuadFigure) ) {
					M3DQuadFigure * fig = dynamic_cast<M3DQuadFigure*>( obj->getFigurePtr(i));
					if (fig->getRowCount() != rows[i] ||
						fig->getColumnCount() != cols[i]) {
							clean(i);	// Topology of this figure changed
							rows[i] = fig->getRowCount();
							cols[i] = fig->getColumnCount();
					}
				}
				else {
					M3DTubeFigure * fig = dynamic_cast<M3DTubeFigure*>( obj->getFigurePtr(i));
					if (1 != rows[i] ||
						fig->getColumnCount() != cols[i]) {
							clean(i);	// Topology of this figure changed
							rows[i] = 1;
							cols[i] = fig->getColumnCount();
					}
				}
			}
		}
		else {
			clean();	// Number of figures changed
			num_figures = obj->getFigureCount();
			delete [] rows;
			rows = new int[num_figures];
			delete [] cols;
			cols = new int[num_figures];
			for (int i = 0; i < num_figures; i++) {
				if (typeid(*obj->getFigurePtr(i)) == typeid(M3DQuadFigure) ) {
					M3DQuadFigure * fig = dynamic_cast<M3DQuadFigure*>( obj->getFigurePtr(i));
					rows[i] = fig->getRowCount();
					cols[i] = fig->getColumnCount();
				}
				else {
					M3DTubeFigure * fig = dynamic_cast<M3DTubeFigure*>( obj->getFigurePtr(i));
					rows[i] = 1;
					cols[i] = fig->getColumnCount();
				}
			}
		}
	}
	else {
		clean();	// obj is NULL
		num_figures = 0;
		delete [] rows;
		rows = NULL;
		delete [] cols;
		cols = NULL;
	}

    object = obj;
}

ThallCode::Diatomgrid * M3DObjectSurfaceRenderer::convertFigure(M3DFigure * _fig)
{
    ThallCode::Diatomgrid * grid;
    M3DPrimitive * prim;
    ThallCode::Diatom * atom;
    int index,
        i, j;
    Quat rotQ, q;
    Quat q2;

    if(_fig == NULL)
        return NULL;

	M3DQuadFigure* fig = dynamic_cast<M3DQuadFigure*>(_fig);

    if(fig->getRowCount() <= 0 || fig->getColumnCount() <= 0)
        return NULL;

    grid = new ThallCode::Diatomgrid ;

    grid->rows = fig->getRowCount();
    grid->cols = fig->getColumnCount();

    grid->dlist = new ThallCode::Diatom[grid->rows * grid->cols];

    rotQ.setAxisAngle(Vector3D(1.0, 0.0, 0.0), R_HALF_PI);

    index = 0;
    for(i = 0; i < grid->cols; i++)
    {
        for(j = 0; j < grid->rows; j++)
        {
            prim = fig->getPrimitivePtr(j, i);
            if(prim == NULL)
                continue;

            atom = grid->dlist + index;

            q = prim->getQ() * rotQ;
            q2.set(q.getX(), q.getY(), q.getZ(), q.getW());

            atom->set_p(ThallCode::DbVector3(prim->getX().getX(), prim->getX().getY(), prim->getX().getZ()));
            atom->set_q(q2);
            atom->set_phi(R_HALF_PI - prim->getTheta());
            atom->set_r(prim->getR());
			atom->set_r1(dynamic_cast<M3DQuadPrimitive*>(prim)->getR0());
			atom->set_r2(dynamic_cast<M3DQuadPrimitive*>(prim)->getR1());

			std::cout << "M3DObjectSurfaceRenderer" << std::endl;

            if(prim->type() == M3D_END_PRIMITIVE)
                atom->set_eta((dynamic_cast<M3DEndPrimitive*>(prim))->getElongation());
            else
                atom->set_eta(1.0);

            index++;
        }
    }

    return grid;
}

using namespace ThallCode;

void M3DObjectSurfaceRenderer::renderFigure(int figureId)
{
    M3DFigure* figure;

    if(object == NULL)
        return;

    figure = object->getFigurePtr(figureId);

	if(figure == NULL)
        return;

	if (figure->getVisibility() == false)
		return;

#ifdef QIONG_DEBUG_LINALG
	figure->print();
#endif

#ifdef DEBUG
	std::cout << "rendering figure " << figureId << std::endl;
#endif

    glColor3fv(figure->getColor());

	Xferlist * xferList = NULL;
	Pointlist_server2 * thisPList;
	thisPList = list[figureId];
	if (thisPList == NULL)
	{
		xferList = convertM3DtoXfer(figure);
		thisPList = new Pointlist_server2;
		thisPList->init(xferList);
		thisPList->ComputeSubdivBoundaryTiles(surfaceLevel);
		thisPList->ComputeSubdivPointCloud(surfaceLevel);
		thisPList->lineWidth(lineWidth);
		list[figureId] = thisPList;
	}
	else {
		xferList = convertM3DtoXfer(figure);
		thisPList->init(xferList);
		thisPList->UpdateSubdivBoundaryTiles(surfaceLevel, xferList);
		thisPList->UpdateSubdivPointCloud(surfaceLevel, xferList);
	}

#ifdef QIONG_DEBUG_LINALG
	int i, j;
	XferAtom * thisAtom;
	for(i = 0; i < xferList->numrows; i++)
		for(j = 0; j < xferList->numcols; j++)
		{
			thisAtom=&(xferList->atomlist[i*xferList->numcols+j]);
			std::cout << "For the atom[" << i << "," << j << "]" << std::endl;
			std::cout << thisAtom->X_eta << std::endl <<
				thisAtom->X_p[0] << " " << thisAtom->X_p[1] << " " << thisAtom->X_p[2] << std::endl <<
				thisAtom->X_q[0] << " " << thisAtom->X_q[1] << " " << thisAtom->X_q[2] << " " << thisAtom->X_q[3] << std::endl << 
				thisAtom->X_r << " " << thisAtom->X_rho << " " << thisAtom->X_theta << std::endl;
		}
#endif

	delete [] (xferList->atomlist);
	delete xferList;

#ifdef DEBUG
	std::cout << "Computing surface at level " << surfaceLevel << " for full rendering"
		<< std::endl;
#endif
// qiong working 12302002
#ifdef QIONG_DEBUG_LINALG2
	figure->print();
#endif

	// Added by CL.  Add displacements to boundary point list
	SubdivBoundary * boundary = figure -> getBoundaryPtr();

	if (boundary != NULL)
	{
		Displacements * disp = boundary -> getDisplacements(surfaceLevel);
		if (disp != NULL)
			thisPList->applyDisplacements(disp->getNumPts(), disp->getVals());

		switch (surfaceStyle) 
		{
			case M3D_NONE:
				break;
			case M3D_POINTCLOUD:
				thisPList->glRenderSubdivPoints();
				break;
			case M3D_WIREFRAME:
				thisPList->glRenderSubdivWires();
				break;
			case M3D_SOLID:
				glEnable(GL_LIGHTING);
				thisPList->glRenderSubdivTiles(0);
				glDisable(GL_LIGHTING);
				break;
			case M3D_COLORED_POINTCLOUD:
				objectSurfaceVisualizer->lineWidth(lineWidth);
				objectSurfaceVisualizer->renderColoredSubdivPoints(figureId,
					surfaceLevel, surfaceColorsMin, surfaceColorsMax);
				break;
			case M3D_COLORED_WIREFRAME:
				objectSurfaceVisualizer->lineWidth(lineWidth);
				objectSurfaceVisualizer->renderColoredSubdivWires(figureId,
					surfaceLevel, surfaceColorsMin, surfaceColorsMax);
				break;
			case M3D_COLORED_SOLID:
				objectSurfaceVisualizer->renderColoredSubdivTiles(figureId,
					surfaceLevel, surfaceColorsMin, surfaceColorsMax, tile_outlines, 0);
				break;
		}
	}
	else
	{
		switch (surfaceStyle) 
		{
			case M3D_NONE:
				break;
			case M3D_POINTCLOUD: 
				thisPList->glRender(surfaceLevel, VPOINTS);
				break;
			case M3D_WIREFRAME:
				thisPList->glRender(surfaceLevel, WIREFRAME);
				break;
			case M3D_SOLID:
				thisPList->glRender(surfaceLevel, FILLEDPOLY);
				break;
			case M3D_COLORED_POINTCLOUD:
			case M3D_COLORED_WIREFRAME:
			case M3D_COLORED_SOLID:
				break;
		}
	}
}

ThallCode::Pointlist_server2 * M3DObjectSurfaceRenderer::computeFigure(int figureId,
    bool ignoreVisibility, bool ignoreBoundary)
{
    M3DFigure * figure;

    if(object == NULL)
        return NULL;

    figure = object->getFigurePtr(figureId);
    if(figure == NULL)
        return NULL;

	if (ignoreVisibility == false && figure->getVisibility() == false)
		return NULL;

#ifdef DEBUG
	std::cout << "Computing rendering of figure " << figureId << std::endl;
#endif

    glColor3fv(figure->getColor());

	ThallCode::Pointlist_server2 * pList = new ThallCode::Pointlist_server2;
	Xferlist * xferList = convertM3DtoXfer(figure);
	pList->init(xferList);
#ifdef DEBUG
	std::cout << "Computing surface at level " << surfaceLevel << " for full rendering"
		<< std::endl;
#endif
    pList->ComputeSubdivBoundaryTiles(surfaceLevel);
	pList->ComputeSubdivPointCloud(surfaceLevel);
	pList->lineWidth(lineWidth);

	if (! ignoreBoundary) {
		SubdivBoundary * boundary = figure -> getBoundaryPtr();
		if (boundary != NULL)
		{
			Displacements * disp = boundary -> getDisplacements(surfaceLevel);
			if (disp != NULL)
				pList->applyDisplacements(disp->getNumPts(), disp->getVals());
		}
	}

	delete [] (xferList->atomlist);
	delete xferList;
	return pList;
}

ThallCode::Pointlist_server2 * M3DObjectSurfaceRenderer::getPList(int figureId) {
    int numFigures;

	if(object == NULL)
        return NULL;

    numFigures = object->getFigureCount();
	if (figureId >= numFigures || figureId < 0)
		return NULL;

	if (list == NULL)
		return NULL;

	return list[figureId];
}

void M3DObjectSurfaceRenderer::createList(int length)
{
	list = new ThallCode::Pointlist_server2 *[length];
	for (int i = 0; i < length; i++)
		list[i] = NULL;
	list_length = length;
}

void M3DObjectSurfaceRenderer::render(bool feedback, bool ignoreVisibility, bool ignoreBoundary)
{
    int numFigures;
    int i;

#ifdef DEBUG
	std::cout << "M3DObjectSurfaceRenderer::render() called" << std::endl;
	std::cout << "\tsurfaceStyle = " << surfaceStyle << std::endl;
#endif

    if(object == NULL)
        return;

    numFigures = object->getFigureCount();
	// The same list is used for either feedback rendering, or for partial and
	// regular rendering for display.  It cannot be used for both, which is why
	// those two modes cannot be mixed.
	if (list == NULL)
		createList(numFigures);

	if (feedback) {
		for(i = 0; i < numFigures; i++)
			list[i] = computeFigure(i, ignoreVisibility, ignoreBoundary);
	}
	else {
		for(i = 0; i < numFigures; i++)
			renderFigure(i);
	}
}

void M3DObjectSurfaceRenderer::drawDots(double * spot_list, int numSpots,
	double size, float defaultWidth, float markedWidth, const float * color,
	int numBigSpot)
{
#ifdef DEBUG
	std::cout << "Drawing " << numSpots << " dots; number " << numBigSpot << " is marked\n";
#endif
	if (color != NULL)
		glColor3fv(color);
	if (size > 0.0) {
		if (defaultWidth > 1.0f)
			glLineWidth(defaultWidth);
		glBegin(GL_LINES);
		for (int i = 0; i < 3*numSpots; i += 3) {
#ifdef DEBUG
			std::cout << '(' << spot_list[i] << ", " << spot_list[i + 1] << ", "
				<< spot_list[i + 2] << ')' << std::endl;
#endif
			if (i == numBigSpot) {
				double s = 1.5*size;
				glEnd();
				glLineWidth(markedWidth);
				glBegin(GL_LINES);

				// Three crossed lines
				glVertex3d(spot_list[i] - s, spot_list[i + 1], spot_list[i + 2]);
				glVertex3d(spot_list[i] + s, spot_list[i + 1], spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1] - s, spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1] + s, spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1], spot_list[i + 2] - s);
				glVertex3d(spot_list[i], spot_list[i + 1], spot_list[i + 2] + s);

				// Draw lines between the ends of the three crossed lines
				glVertex3d(spot_list[i] - s, spot_list[i + 1], spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1] - s, spot_list[i + 2]);
				glVertex3d(spot_list[i] - s, spot_list[i + 1], spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1] + s, spot_list[i + 2]);
				glVertex3d(spot_list[i] - s, spot_list[i + 1], spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1], spot_list[i + 2] - s);
				glVertex3d(spot_list[i] - s, spot_list[i + 1], spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1], spot_list[i + 2] + s);

				glVertex3d(spot_list[i] + s, spot_list[i + 1], spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1] - s, spot_list[i + 2]);
				glVertex3d(spot_list[i] + s, spot_list[i + 1], spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1] + s, spot_list[i + 2]);
				glVertex3d(spot_list[i] + s, spot_list[i + 1], spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1], spot_list[i + 2] - s);
				glVertex3d(spot_list[i] + s, spot_list[i + 1], spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1], spot_list[i + 2] + s);

				glEnd();
				glLineWidth(defaultWidth);
				glBegin(GL_LINES);
			}
			else {
				// Three crossed lines
				glVertex3d(spot_list[i] - size, spot_list[i + 1], spot_list[i + 2]);
				glVertex3d(spot_list[i] + size, spot_list[i + 1], spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1] - size, spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1] + size, spot_list[i + 2]);
				glVertex3d(spot_list[i], spot_list[i + 1], spot_list[i + 2] - size);
				glVertex3d(spot_list[i], spot_list[i + 1], spot_list[i + 2] + size);
			}
		}
		glEnd();
		if (defaultWidth > 1.0f)
			glLineWidth(1.0f);
	}
	else {
		glBegin(GL_POINTS);
		for (int i = 0; i < 3*numSpots; i += 3)
			glVertex3d(spot_list[i], spot_list[i + 1], spot_list[i + 2]);
		glEnd();
	}
}

void M3DObjectSurfaceRenderer::drawDots(double size, float defaultWidth, bool distinguishMarked,
										float markedWidth, const float * color)
{
    int numFigures;
    int i, mark;
	M3DFigure * figure;
	double * landmarks;
	int numLandmarks;
	float dotColor[3];
	float * c;

#ifdef DEBUG
	std::cout << "M3DObjectSurfaceRenderer::drawDots(size=" << size << ", defaultWidth="
		<< defaultWidth << ", markedWidth=" << markedWidth << ')' << std::endl;
#endif
    if(object == NULL)
        return;

    numFigures = object->getFigureCount();
	for(i = 0; i < numFigures; i++) {
		figure = object->getFigurePtr(i);
		if (color == NULL) {
			const float * figuralColor = figure->getColor();
			if (figuralColor[0] > 0.5f)
				dotColor[0] = 0.0f;
			else
				dotColor[0] = 1.0f;
			if (figuralColor[1] > 0.5f)
				dotColor[1] = 0.0f;
			else
				dotColor[1] = 1.0f;
			if (figuralColor[2] > 0.5f)
				dotColor[2] = 0.0f;
			else
				dotColor[2] = 1.0f;
			c = dotColor;
		}
		else
			c = (float *) color;
		numLandmarks = figure->getLandmarkCount();
		if (numLandmarks > 0) {
			if (distinguishMarked)
				mark = figure->getMarkedLandmark();		// mark = -1, if unmarked
			else
				mark = -1;
			landmarks = new double[3*numLandmarks];
			for (int l = 0; l < numLandmarks; l++) {
				Vector3D vl = figure->getLandmark(l);
				landmarks[3*l] = vl.getX();
				landmarks[3*l + 1] = vl.getY();
				landmarks[3*l + 2] = vl.getZ();
			}
			drawDots(landmarks, numLandmarks, size, defaultWidth, markedWidth, c, 3*mark);
			delete [] landmarks;
		}
	}
}

/*
	// Compute max and average of color component values
	float color[3];
	float avg = 0.0;
	float max = 0.0;
	int i;

	for (i = 0; i < 3; i++) {
		color[i] = (figure->getColor())[i];
		avg += color[i];
		if (color[i] > max)
			max = color[i];
	}
	avg /= 3.0;

	// Adjust color components -
	// AGG: other ideas: use a single unique color
	//	    use alpha blending
	cout << "Original color: " << color[0] << ", " << color[1] << ", " << color[2] << '\n';
	for (i = 0; i < 3; i++) {
		if (color[i] == max)
			color[i] -= avg/2.0;
		else color[i] += avg/2.0;
		if (color[i] > 1.0)
			color[i] = 1.0;
		else if (color[i] < 0.0)
			color[i] = 0,0;
	}
	cout << "Final color: " << color[0] << ", " << color[1] << ", " << color[2] << '\n';
*/

void M3DObjectSurfaceRenderer::partial_render(SelectedPartialFigures * figure_list,
				double cull_distance, bool complete_figures, int level)
{
    int current_figure;
	int num_listed_figures;
    int i;

#ifdef DEBUG
	std::cout << "M3DObjectSurfaceRenderer::partial_render() called\n";
	std::cout << "\tpartialSurfaceStyle = " << partialSurfaceStyle << '\n';
	std::cout << "\tpartialSurfaceLevel = " << level << std::endl;
#endif

	if (object == NULL)
		return;

	num_listed_figures = figure_list->number();
	current_figure = figure_list->displayFigure();

	if (surfaceStyle != M3D_NONE) {
		// Render the marked figure
		int markedFigureId;
		markedFigureId = figure_list->markedFigure();
		if (list == NULL)
			createList(object->getFigureCount());
		renderFigure(markedFigureId);
	}

	ThallCode::Pointlist_server2 ** pLists =
		(ThallCode::Pointlist_server2 **) figure_list->pointLists();

	if (figure_list->pointListsChanged()) {
	    M3DFigure * figure;
#ifdef DEBUG
		std::cout << "Pointlists changed; number of listed figures = "
			<< num_listed_figures << std::endl;
#endif
		for (i = 0; i < num_listed_figures; i++) {
			int figureID = figure_list->figureID(i);
			figure = object->getFigurePtr(figureID);
			//if (figure == NULL)
			//	continue;
			Xferlist * xferList = convertM3DtoXfer(figure);
			pLists[i]->init(xferList);
#ifdef DEBUG
			std::cout << "Computing boundary points at level "
				<< level + CONSTRAINTS_RESOLUTION_BIAS << " for figure "
				<< figure_list->figureID(i) << std::endl;
#endif
			// Generate tile sets and boundary point lists used in the 
			// rendering functions called below.
#ifdef DEBUG
	std::cout << "Computing surface at level " << level + CONSTRAINTS_RESOLUTION_BIAS
		<< " for partial rendering" << std::endl;
#endif
			pLists[i]->ComputeSubdivBoundaryTiles(level + CONSTRAINTS_RESOLUTION_BIAS);
			pLists[i]->ComputeSubdivPointCloud(level + CONSTRAINTS_RESOLUTION_BIAS);

#ifdef DEBUG
			int numBPoints, numQuads;
			Bpoint * bpoints;
			Bpoint * quads;

			pLists[i]->subdivtileinfo(&numQuads, &quads);
			pLists[i]->subdivboundaryinfo(&numBPoints, &bpoints);
			std::cout << "Figure " << figure_list->figureID(i) << " has " << numBPoints
				<< " boundary points and " << numQuads << " tiles." << std::endl;
#endif

			pLists[i]->setColors(figure->getColor(), figure_list->color());

			delete [] (xferList->atomlist);
			delete xferList;
		}
		figure_list->pointListsChangeReset();
	}

	// Render the figures for which constraints are being set
	for (i = 0; i < num_listed_figures; i++) {
		if (current_figure == -1 || current_figure == figure_list->figureID(i)) {
#ifdef DEBUG
			std::cout << "partial-rendering figure " << figure_list->figureID(i) << std::endl;
#endif
			switch (partialSurfaceStyle) {
				case M3D_NONE:
					break;

				case M3D_POINTCLOUD:
					pLists[i]->glPartialRenderSubdivPoints(figure_list->point_distances(i),
						cull_distance, complete_figures);
					break;

				case M3D_WIREFRAME:
					pLists[i]->glPartialRenderSubdivWires(figure_list->quad_distances(i),
						cull_distance, complete_figures);
					break;

				case M3D_SOLID:
					glEnable(GL_LIGHTING);
					pLists[i]->glPartialRenderSubdivTiles(figure_list->quad_distances(i),
						cull_distance, complete_figures);
					glDisable(GL_LIGHTING);
					break;

				case M3D_COLORED_POINTCLOUD:
				case M3D_COLORED_WIREFRAME:
				case M3D_COLORED_SOLID:
					break;
			}
	    }
    }
}

// Used with Pointlist_server2::init(5 args) for rendering a surface directly
// from a list of points
void M3DObjectSurfaceRenderer::renderFigureFromControlPoints
    (int numverts, int numfaces, Bpoint *pntlist, ThallCode::PindexList *faces)
{
   Pointlist_server2 * thisPList;

#ifdef DEBUG
   std::cout << "Rendering object from list of control points" << std::endl;
#endif

   float color[3];
   color[0] = 0.8f;
   color[1] = 0.2f;
   color[2] = 1.0f;
   glColor3fv(color);

   thisPList = new Pointlist_server2;
   thisPList->init(numverts, numfaces, pntlist, faces, surfaceLevel);
   thisPList->ComputeSubdivPointCloud(surfaceLevel);
   thisPList->lineWidth(lineWidth);

#ifdef DEBUG
   std::cout << "Computing surface at level " << surfaceLevel << " for full rendering"
       << std::endl;
#endif

   switch (surfaceStyle)
   {
   case M3D_NONE:
       break;
   case M3D_POINTCLOUD:
       thisPList->glRender(surfaceLevel, VPOINTS);
       break;
   case M3D_WIREFRAME:
       thisPList->glRender(surfaceLevel, WIREFRAME);
       break;
   case M3D_SOLID:
       thisPList->glRender(surfaceLevel, FILLEDPOLY);
       break;
   case M3D_COLORED_POINTCLOUD:
   case M3D_COLORED_WIREFRAME:
   case M3D_COLORED_SOLID:
       break;
   }

   delete thisPList;
}


