
#include <iostream>
#include <queue>
#include <set>

#define D_QUAT
#define D_XFERLIST
#define D_POINTLIST_SERVER2
#define D_POINTLIST_SERVERB
#include "Shapedepend.h"

#include "M3DQuadFigure.h"
#include "M3DObject.h"
#include "M3DFigureTreeNode.h"
#include "M3DObjectSurfaceRenderer.h"
#include "renderDefinitions.h"
#include "SelectedPartialFigures.h"
#include "M3DBlendedRenderer.h"

using namespace ThallCode;
using namespace std;

void (* testInterrupt)(void);
bool M3DBlendedRenderer::halt;


//#define FIGURE_MAP_TO_REMOVE_TESTING

//float * M3DBlendedRenderer::blendExtents;
//float * M3DBlendedRenderer::blendAmounts;
//#define DEBUG
//#define QIONG_DEBUG_LINALG
//#define QIONG_DEBUG_LINALG2

M3DBlendedRenderer::M3DBlendedRenderer(M3DObject * obj, int subLevel)
{
    object = obj;
    level = subLevel;
//	list_length = 0;
	subObject = NULL;

	numTrees = 0;
//	blendExtents = NULL;
//	blendAmounts = NULL;
	blendedPList = NULL;
//    transformMatrix.setToIdentity();

	lineWidth=1.0;
}

M3DBlendedRenderer::~M3DBlendedRenderer()
{
	// FIXME: Why is this commented out?
/*
	if (list_length != 0) {
		for (int i = 0; i < list_length; i++)
			delete list[i];
		delete [] list;
	}
*/
	// object does not belong to this class
}

void M3DBlendedRenderer::init(M3DObject *obj, int subLevel)
{
#ifdef DEBUG
        cout << "M3DBlendedRenderer::init() called" << endl;
#endif
    object = obj;
    level = subLevel;

#if defined(DRAW_ONLY_BLEND) || defined(RECOLOR_BLEND)
#endif

/*
	if(blendExtents!=NULL)
	{
		delete []blendExtents;
		blendExtents=NULL;
	}
	if(blendAmounts!=NULL)
	{
		delete []blendAmounts;
		blendExtents=NULL;
	}
	//delete [] blendExtents;
	//delete [] blendAmounts;
*/
    if (object != NULL) // && object->getFigureTreeCount()>0) 
	{
		numTrees = object->getFigureTreeCount();
		// For simplicity, this allocates more members than the necessary (tree_depth - 1)
//		blendExtents = new float[object->getFigureCount()];
//		blendAmounts = new float[object->getFigureCount()];
	}
	else
	{
		numTrees = 0;
//		blendExtents = NULL;
//		blendAmounts = NULL;
	}

	if(blendedPList!=NULL)
	{
		delete blendedPList;
		blendedPList=NULL;
	}
}

ThallCode::Diatomgrid * M3DBlendedRenderer::convertFigure(M3DQuadFigure * fig)
{
    ThallCode::Diatomgrid * grid;
    M3DPrimitive * prim;
    ThallCode::Diatom * atom;
    int index,
        i, j;
	Quat rotQ, q;
	Quat q2;

    if(fig == NULL)
        return NULL;

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

			std::cout << "M3DBlendedRenderer" << std::endl;

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
void M3DBlendedRenderer::renderFigure(int index)
{
    M3DFigure * figure;

    if(object == NULL)
        return;

    figure = object->getFigurePtr(index);

	if(figure == NULL)
        return;

	if (figure->getVisibility() == false)
		return;

#ifdef QIONG_DEBUG_LINALG
	figure->print();
#endif

#ifdef DEBUG
	std::cout << "rendering figure " << index << std::endl;
#endif

    glColor3fv(figure->getColor());

	ThallCode::Pointlist_server2 pList;
	Xferlist * xferList = convertM3DtoXfer(figure);

#ifdef QIONG_DEBUG_LINALG
	int i, j;
	XferAtom *thisAtom;
	for(i=0; i<xferList->numrows; i++)
		for(j=0; j<xferList->numcols; j++)
		{
			thisAtom=&(xferList->atomlist[i*xferList->numcols+j]);
			std::cout << "For the atom[" << i << "," << j << "]" << std::endl;
			std::cout << thisAtom->X_eta << std::endl <<
				thisAtom->X_p[0] << " " << thisAtom->X_p[1] << " " << thisAtom->X_p[2] << std::endl <<
				thisAtom->X_q[0] << " " << thisAtom->X_q[1] << " " << thisAtom->X_q[2] << " " << thisAtom->X_q[3] << std::endl << 
				thisAtom->X_r << " " << thisAtom->X_rho << " " << thisAtom->X_theta << std::endl;
		}
#endif

	pList.init(xferList);
#ifdef DEBUG
	std::cout << "Computing surface at level " << level << " for full rendering"
		<< std::endl;
#endif

#ifdef QIONG_DEBUG_LINALG2
	figure->print();
#endif

    pList.ComputeSubdivBoundaryTiles(level);
	pList.ComputeSubdivPointCloud(level);
	//pList.InitializeSubdivSurf(level);

	///////////////////////////////////////////////////////////////////
	// Added by CL.  Add displacements to boundary point list

	SubdivBoundary * boundary = figure -> getBoundaryPtr();
	/*if (boundary == NULL)
	{
		std::cout << "Null boundary in M3DobjRenderer::renderFig(" << index <<")." << std::endl;
		return;
	}*/

	if (boundary != NULL)
	{
		Displacements * disp = boundary -> getDisplacements(level);
		if (disp != NULL)
			pList.applyDisplacements(disp->getNumPts(), disp->getVals());

		switch (surfaceStyle) 
		{
		case M3D_NONE:
			break;
		case M3D_POINTCLOUD:
			pList.glRenderSubdivPoints();
			break;
		case M3D_WIREFRAME:
			pList.glRenderSubdivWires();
			break;
		case M3D_SOLID:
			glEnable(GL_LIGHTING);
			pList.glRenderSubdivTiles(0);
			glDisable(GL_LIGHTING);
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
			{
				pList.glRender(level, VPOINTS);
				break;
			}
		case M3D_WIREFRAME:
			{
				pList.glRender(level, WIREFRAME);
				break;
			}
		case M3D_SOLID:
			{
				pList.glRender(level, FILLEDPOLY);
				break;
			}
		}
	}

	delete [] (xferList->atomlist);
	delete xferList;
}

ThallCode::Pointlist_server2 * M3DBlendedRenderer::computeFigure(int index,
    bool ignoreVisibility, bool ignoreBoundary)
{
    M3DFigure * figure;

    if(object == NULL)
        return NULL;

    figure = object->getFigurePtr(index);
    if(figure == NULL)
        return NULL;

	if (ignoreVisibility == false && figure->getVisibility() == false)
		return NULL;

#ifdef DEBUG
	std::cout << "Computing rendering of figure " << index << std::endl;
#endif

    glColor3fv(figure->getColor());

	ThallCode::Pointlist_server2 * pList = new ThallCode::Pointlist_server2;
	Xferlist * xferList = convertM3DtoXfer(figure);
	pList->init(xferList);
#ifdef DEBUG
	std::cout << "Computing surface at level " << level << " for full rendering"
		<< std::endl;
#endif
    pList->ComputeSubdivBoundaryTiles(level);
	pList->ComputeSubdivPointCloud(level);

	if (! ignoreBoundary) {
		SubdivBoundary * boundary = figure -> getBoundaryPtr();
		if (boundary != NULL)
		{
			Displacements * disp = boundary -> getDisplacements(level);
			if (disp != NULL)
				pList->applyDisplacements(disp->getNumPts(), disp->getVals());
		}
	}

	delete [] (xferList->atomlist);
	delete xferList;
	return pList;
}
/*
ThallCode::Pointlist_server2 * M3DBlendedRenderer::getPList(int figureId) {
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
*/
void M3DBlendedRenderer::drawBlendedObject()
{
	blendedPList->lineWidth(lineWidth);

	switch (surfaceStyle) 
	{
	case M3D_NONE:
		break;
	case M3D_POINTCLOUD: 
	{
		//pList.glRender(level, VPOINTS);
		blendedPList->glRender(level, VPOINTS);
		break;
	}
	case M3D_WIREFRAME:
	{
		//pList.glRender(level, WIREFRAME);
		blendedPList->glRender(level, WIREFRAME);
		break;
	}
	case M3D_SOLID:
	{
		//pList.glRender(level, FILLEDPOLY);
		blendedPList->glRender(level, FILLEDPOLY);
		break;
	}
	default:
		break;
	}
}

bool M3DBlendedRenderer::render()
{
	// BUG FIXED
	// object		-> pointer to the entire object complex
	// root			-> root to the current object (single- or multi-figured)
	// subobject	->
	// node			->
	// currentFigureId used to index the current object figure (either in single- or multi-figure
	// object) in the whole complex
	RenderStyle renderMode;		// WIREFRAME, FILLEDPOLY, VPOINTS, NO_RENDER
	switch (surfaceStyle) 
	{
	case M3D_NONE:
		renderMode=NO_RENDER;
		break;
	case M3D_POINTCLOUD: 
		renderMode=VPOINTS;
		break;
	case M3D_WIREFRAME:
		renderMode=WIREFRAME;
		break;
	case M3D_SOLID:
		renderMode=FILLEDPOLY;
		break;
	default:
		renderMode=NO_RENDER;
		break;
	}

	int currentFigureId=0;

    M3DFigureTreeNode * root;
    M3DFigureTreeNode * node;
	bool recompute;
//	double lmin;	// Used to compute cubeSize and bound

#ifdef DEBUG
	cout << "M3DBlendedRenderer::render() called" << endl;
#endif
    if(object == NULL)
        return false;

#ifdef COMPOUND_FIGURES_ONLY
	M3DObjectSurfaceRenderer surfaceRenderer(object, level);
	surfaceRenderer.setSurfaceStyle(surfaceStyle);
#else
	M3DObjectSurfaceRenderer surfaceRenderer(object, level);
	surfaceRenderer.setSurfaceStyle(surfaceStyle);
#endif

	if(blendedPList==NULL)
		recompute=true;
	else
		recompute = !blendedPList->subdivSurfInitialized; //(blendedPList == NULL); //recompute = (vertexLists == NULL);

	if (recompute && numTrees > 0) 
	{
//		vertexLists = new vector<Point>[numTrees];
//		normalsLists = new vector<Vector3D>[numTrees];
#if defined(DRAW_ONLY_BLEND) || defined(RECOLOR_BLEND)
//		inBlendLists = new vector<bool>[numTrees];
#endif
		cout << "Blending ..." << flush;
	}

	int * unused = new int[object->getFigureCount()];

	for (int tree = 0; tree < numTrees; tree++) 
	{
		root = object->getFigureTreeRoot(tree);
		if (root == NULL)
			continue;

#ifdef DEBUG
		//cout << "Tree " << tree << " has " << root->getChildCount() << " children" << endl;
#endif
#ifdef COMPOUND_FIGURES_ONLY
		//if (root->getChildCount() > 0)
#endif
		// need to deal with the case of multi-figure object and single-figure objects separately!
		if(root->getChildCount() > 0)
		{
			glColor3fv(object->getFigurePtr(root->getFigureId())->getColor());
			if(recompute)	// equivalent to (recompute)
			{
				// Create an object containing a copy of the figure tree
				if(subObject!=NULL)
					delete subObject;
				subObject=new M3DObject;
				subObject->addTree(root, object, unused); //subObject->addTree(root, object); changed according to changed in M3DObject

				blendedPList=new Pointlist_serverB();

				// Build a breadth-first list of the nodes in this tree
				node=subObject->getFigureTreeRoot(0);

				if (halt) 
				{
					cout << endl;
					halt = false;
					return false;
				}

				if(node->getChildCount()>0)
				{
					// my plan is to implement for simple case, including 2 figures first
					// then extend the same method to the tree structure

					// Create the tileLists entries for the tree in breadth-first order,
					// and initialize the minimum length of the object in model space

					// The implcit rendering cube size depends on the object size
//					renderFigure(tree);
#ifndef FIGURE_MAP_TO_REMOVE_TESTING
					blendedPList->init(subObject, 1);
					blendedPList->SetRenderMode(renderMode);
					if(!blendedPList->InitializeSubdivSurf(level))
					{
						delete blendedPList;
						blendedPList=NULL;
					}
#else
					blendedPList->init(subObject, level);
					blendedPList->SetRenderMode(renderMode);
					Bpoint2 testPt;
					CCVertex *thisVert;
					int ptsToRemoveNum[2];
					short *vertexMap[2];
					for(int i=0; i<2; i++)
						vertexMap[i]=NULL;
					MyFloat color[3], color_2[3];
					color[0]=1;
					color[1]=0;
					color[2]=0;
					color_2[0]=0;
					color_2[1]=1;
					color_2[2]=0;
					blendedPList->ComputeFiguresMapsToRemove(ptsToRemoveNum, vertexMap, 3);
					for(i=0; i<blendedPList->GetOriginalPList(0)->GetThisSurfPointer()->
						SubMeshPointer(3)->numverts(); i++)
					{
						if(!vertexMap[0][i])
						{
							thisVert=blendedPList->GetOriginalPList(0)->GetThisSurfPointer()->
								SubMeshPointer(3)->VertList(i);
							testPt.setvals(thisVert->limitloc.x(), thisVert->limitloc.y(), thisVert->limitloc.z(),
								0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
							blendedPList->glRenderPoint(testPt, color);
						}
						else
						{
							thisVert=blendedPList->GetOriginalPList(0)->GetThisSurfPointer()->
								SubMeshPointer(3)->VertList(i);
							testPt.setvals(thisVert->limitloc.x(), thisVert->limitloc.y(), thisVert->limitloc.z(),
								0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
							blendedPList->glRenderPoint(testPt, color_2);
						}
					}
					color[2]=1;
					//blendedPList->ComputeFigureVertexMapToRemove(&subPtsToRemoveNum, &vertexMapSub, 4, 1);
					if(vertexMap[1]!=NULL)
					{
						for(i=0; i<blendedPList->GetOriginalPList(1)->GetThisSurfPointer()->
							SubMeshPointer(3)->numverts(); i++)
						{
							if(!vertexMap[1][i])
							{
								thisVert=blendedPList->GetOriginalPList(1)->GetThisSurfPointer()->
									SubMeshPointer(3)->VertList(i);
								testPt.setvals(thisVert->limitloc.x(), thisVert->limitloc.y(), thisVert->limitloc.z(),
									0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
								blendedPList->glRenderPoint(testPt, color);
							}
							else
							{
								thisVert=blendedPList->GetOriginalPList(1)->GetThisSurfPointer()->
									SubMeshPointer(3)->VertList(i);
								testPt.setvals(thisVert->limitloc.x(), thisVert->limitloc.y(), thisVert->limitloc.z(),
									0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
								blendedPList->glRenderPoint(testPt, color_2);
							}
						}
					}
					else
							surfaceRenderer.renderFigure(1);
					for(i=0; i<2; i++)
						if(vertexMap[i]!=NULL)
							delete []vertexMap[i];
#endif

					delete subObject;
					subObject = NULL;
				}
				else
					surfaceRenderer.renderFigure(currentFigureId);

			} // if (recompute)
#ifndef FIGURE_MAP_TO_REMOVE_TESTING
			if(blendedPList!=NULL)
				drawBlendedObject(); //redrawFigure(tree, root->getFigureId());
#endif

//#ifdef QIONG_BLENDING_NO_DEBUG
//			drawBlendedObject(); //redrawFigure(tree, root->getFigureId());
//#endif
		}
		else
			surfaceRenderer.renderFigure(root->getFigureId());
#ifdef COMPOUND_FIGURES_ONLY
//		else
//			surfaceRenderer.renderFigure(root->getFigureId());
#endif

		currentFigureId+=1+root->getChildCount();
	}
	if (recompute && numTrees > 0)
		cout << " done." << endl;

	delete [] unused;

	return true;
}

void M3DBlendedRenderer::partial_render(SelectedPartialFigures * figure_list,
				double cull_distance, bool complete_figures, int lvl)
{
    int current_figure;
	int num_listed_figures;
    int i;

#ifdef DEBUG
	std::cout << "M3DBlendedRenderer::partial_render() called\n";
	std::cout << "\tpartialSurfaceStyle = " << partialSurfaceStyle << '\n';
	std::cout << "\tpartialSurfaceLevel = " << lvl << std::endl;
#endif

	if (object == NULL)
		return;

	num_listed_figures = figure_list->number();
	current_figure = figure_list->displayFigure();

	if (surfaceStyle != M3D_NONE) {
		// Render the marked figure
		int markedFigureId;
		markedFigureId = figure_list->markedFigure();
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
				<< lvl + CONSTRAINTS_RESOLUTION_BIAS << " for figure "
				<< figure_list->figureID(i) << std::endl;
#endif
			// Generate tile sets and boundary point lists used in the 
			// rendering functions called below.
#ifdef DEBUG
	std::cout << "Computing surface at level " << lvl + CONSTRAINTS_RESOLUTION_BIAS
		<< " for partial rendering" << std::endl;
#endif
			pLists[i]->ComputeSubdivBoundaryTiles(lvl + CONSTRAINTS_RESOLUTION_BIAS);
			pLists[i]->ComputeSubdivPointCloud(lvl + CONSTRAINTS_RESOLUTION_BIAS);

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
			}
	    }
    }
}

