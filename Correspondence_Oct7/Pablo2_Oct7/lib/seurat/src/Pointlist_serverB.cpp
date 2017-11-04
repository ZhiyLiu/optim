/**************************************************************************************************************
2003/09/08
	BUG FIXED
		if(iCurveLengthBeforeDilation > iCurveLengthAfterDilation)
		to detect the incorrect orientation of the intersection curve due to the wrong ordering of the
		sub-figure atoms

		if(interPNum>MAX_INTERSECTION_NUMBER || interPNum<=0)
		if(interPNum[0]>MAX_INTERSECTION_NUMBER || interPNum[1]>MAX_INTERSECTION_NUMBER || \
		interPNum[0]<=0 || interPNum[1]<=0)
		if(interPNum[cDex]>MAX_INTERSECTION_NUMBER || interPNum[cDex]<=0)
		to detect the case that too many or none intersection points found in the tri-tri intersections

		if(subInerPNum<=0), if(mCurveNum<=0 || mainInterPNum<=0)						for 2-figure object
		if(subInerPNum[cDex]<=0), if(mCurveNum[cDex]<=0 || mainInterPNum[cDex]<=0)		for Torus
		to detect the case that no CUT or BOUNDARY CURVE can be found on sub- or host-figure at all

2003/07/15 
	POSSIBLE BUG
		(i=mCurveStartDex; i<=mCurveEndDex; i++)
		mCurveStartDex>mCurveEndDex could be true because of the cyclic curve?

2003/07/10
	BUG FIXED
	1. in RegularizeIntersectionCurve

	2. in RegularizeDilatedIntersectionCurve

2003/07/05
	NEW FUNCTION
	CheckIntersectionCurveOrientation
	to check the orientation of the intersection curve directly by the reference of the sub-figure
	by doing this, I can actually do the dilation of the intersection curve chunk by chunk!

2003/07/01
	NEW FUNCTION
	// Re-regularize the intersection curve after the dilation
	// after the dilation step, still need to re-regularize the dilated curve because at high-curvature
	// part, vertices tend to be crunched together
	RegularizeDilatedIntersectionCurve(interCurveD, interCurve, &interPNum);

2003/06/25
	NEW FUNCTION 2003/06/25

	QIONG_DEBUG_20030625

	1. re-write the codes to connect the cut-curve on the subfigure and the dilated-i-curve

2003/06/20
	TO REMOVE:

2003/06/01
	(no memory leak so far found in the blending codes! all the leakage is caused by the old seurat,
	which has been optimized)

	BUG FIXED 2003/06/01

	1. 
		fixed the codes in render(). the old codes treat all the objects in an object complex as 
		multi-figured.

2003/05/25
	Trying to use the new method: calculate the intersection at level 1 instead of level 0
	W/ SUBDIV_LEVEL_TO_USE
		Control the level of meshes used in the intersection calculation!

	W/ QIONG_DEBUG_FILE_TO_WRITE
		open a file called 'debug.dat' to write all the garbage debugging information. if it's turned 
		off, all the garbage will be dumped to the console window

	W/ QIONG_DEBUG_20030525

	BUG FIXED 2003/05/25
	5.
		//int maxU
		double maxU;

	4.
		// BUG FIXED 2003/05/25
		right before connecting all the intersection segments
		// i++;
		i=1;

	3.
		During the calculation of the intersection curve, some of the intersection segments were 
		degenerated	to zero-length. They need to be tagged and counted!
	2.
		EPSILON2 plays a crucial role! , the intersection curve could be accidently determined as open 
		because the numerical problems give me false 'negative' matches between intersection ponits.
		Have to manually reduce EPSILON2 when calculating the intersection between level 1 meshes.

		TO DO?
		make EPSILON2 a dynamic value, if an OPEN intersection curve is found, increase EPSILON2 and
		try again?

	1.
		// int minUPlus, minUPlus2;
		MyFloat minUPlus, minUPlus2;

2003/05/10
	W/ QIONG_DEBUG_20030510

	NEW FUNCTION 2003/05/10
		void RegularizeIntersectionCurveUVT(double uMaxM, double vMaxM, double uMaxS, double vMaxS, \
		Bpoint2 *curv, int pNum);

		bool IsCounterClockwise
		need2ReverseIntersectionCurve
		to check whether I need to reverse the whole intersection curve

	BUG FIXED 2003/05/10

	5.
		TO DO: make it even more robust?
		in Pointlist_serverB::RenderFigure
		Apparently in the codes cutting the sub-figure, the part checking whether to remove a particular 
		vertex is not CORRECT! Can't draw the conclusion just by checking
		if(thisVert->uvt.u-SUB_FIGURE_CUTOFF_U<-1)
		So far so good for now, but how to deal with the case that the intersection curve is on the 
		bottom of the sub-figure? Or the problem is, whether we really need to handle that case?

	4.
		TO DO: make it even more robust!
		in Pointlist_serverB::RenderFigure
		connecting intersection curve to the sub-figure, to figure out which of subVertOldDex[0/1] goes 
		first. the ultimate solution needs to be found from the orientation of the cut-curve on the 
		sub-figure

	3.
		in Pointlist_serverB::RenderFigure
		vertListMain[..] was not updated correctly

	2.
		In Pointlist_serverB::ConvertCoordinates
		It is a numerical problem because the UVT.t of those vertices on top/bottom is supposed to be
		1.0/-1.0, but the real value is 1.0+1.0e-16 or something like that, it's because of the 
		interpolation of the barycentric coordinates of all the intersection vertices from the 2 
		intersected triangles. Meanwhile, I wrote a new function RegularizeIntersectionCurveUVT to 
		regularize those 'irregular' UVTs.
		have to deal with such kind of numerical problem all the time!

		WARNING

		Therefore a new function RegularizeIntersectionCurveUVT is added to regularize the UVTs of all
		the intersection vertices of the intersection curve

	1.
		in RegularizeIntersectionCurve

2003/04/25

	2.
	TO DO
	Update the Pointlist_server2 to accomodate the 2nd set of UVTR coordinates of the intersection curve

	1.
	Calculate the UVT for the intersection points in the tri_tri_intersect_with_isectline_w_UVTRad 
	routine itself

2003/04/09
	simple Delaunay-like metric + Minimal angle
	so far works well, just one more thing: need to deal with illegal connections!

2003/04/08
	Revisited the remeshing codes, now use the minimal angle in the triangle as the creterion
	(still, a ROBUST algorithm will be much more helpful!
	this unfortunately does not work :(

2003/04/01
	Revisited the codes between QIONG_BLENDING_MAIN_FIGURE and QIONG_BLENDING_INTERSECTION_MAIN_FIGURE, the 
	main reason is to get a more robust routine to calcualte the curve on the main figure that should be 
	connected to the dilated intersection curve.

2003/03/25
	Revisited the codes to connect the dilated intersection curve and the main-figure
	used Delaunay-like triangulation instead to connect 2 curves together
	(previously, simply picked the minimal edge of the 2)

***************************************************************************************************************/

/**************************************************************************************************************
2003/05/25
	NEW FUNCTIONS
	int WithinPolygonAtSublvlUVT(CCPolygon *thisPoly, int subLvl, double u, double v, double t);
	int WithinPolyDexAtSublvlUVT(int pDex, int subLvl, double u, double v, double t);

	bool Pointlist_serverB::InitializeBlendedSubdivSurf_1(int vNum, int fNum, Bpoint2 *vList, PindexList *fList, int *vTag, int *fTag, Pointlist_serverB **pListOld)
	bool Pointlist_serverB::InitializeBlendedSubdivSurf_New(int vNum, int fNum, Bpoint2 *vList, PindexList *fList, int *vTag, int *fTag, Pointlist_serverB **pListOld)

	to 'cooperate' with the M3DBlendedRenderer, in order to calculate the intersection curve between
	the level 1 mesh instead of level 0 mesh

2003/04/25
	Bug
		When averaging 2 vertices (u1, v1, t1) and (u2, v2, t2), if t1=+/- 1 but t2 is in between,
		we can not just add them up and divide it by 2
		The right thing to do is to deal with all the cases just like those to dilate the curve
		I decide to make the change in CCSubdivSurf instead of here, and we can suppose that all
		the UVTRs output from the subdivision routine are ALREADY REGULARIZED!

	W/ QIONG_DEBUG_20030425
		debugging the InitializeBlendedSubdivSurf(...)

2003/04/01
	W/ QIONG_BLENDING_PURE_SMOOTHING 
		to show the pure smoothing effect with the constructed mesh with the blended regions,
		more or less a debugging tool to check the validity of the mesh itself

2003/03/25
	in bool Pointlist_serverB::InitializeSubdivSurf()
	called function
	CCSubdivsurf::CompareNormalsFromVerticesAndFromSubdivision(int meshLvl)
	to compare the analytic normals and the numerically computed normals
***************************************************************************************************************/


/*
	111702
			void UpdateSubdivBoundaryTiles(int subLevel, Xferlist *thisList);
			void UpdateSubdivPointCloud(int subLevel, Xferlist *thisList);
			bool UpdateSubdivSurf(int subLevel, Xferlist *thisList);

	111602
			void UpdateSubdivBoundaryTiles(int subLevel, Diatomgrid *thisGrid);
			void UpdateSubdivPointCloud(int subLevel, Diatomgrid *thisGrid);
			bool UpdateSubdivSurf(int subLevel, Diatomgrid *thisGrid);

	111302	Pointlist_serverB::UpdateSubdivSurf(Xferlist *thisList) and 
			Pointlist_serverB::UpdateSubdivSurf(Diatomgrid *thisGrid)
			1. added this new function call, keep the topology untouched!

	111202	Pointlist_serverB::subdivBposition
			1. disabled the codes normalizing the interpolated normals!
			ask Guodong to update his subfigure transformation codes
*/

// Disable warning C4786 on Windows (identifier length too long)
#ifdef WIN32
#pragma warning ( disable : 4786 )
#endif

#include <stdio.h>
#include <math.h>
#include <float.h>

#include <queue>
#include <set>
#include <list>

#include "M3DFigure.h"
#include "M3DObject.h"
#include "M3DFigureTreeNode.h"

#define D_XFERLIST
#define D_POINTLIST_SERVERB
#define D_ZEROFINDER
#include "Shapedepend.h"

#include "MyQueue.h"
#include "M3DObjectSurfaceRenderer.h"
#include "Tritri.h"

using namespace ThallCode;
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SOME CONSTANTS
#define MAX_INTERSECTION_NUMBER 500
#define MAX_MAIN_CURVE_NUMBER 50

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DEBUGGING TOOLS
#define QIONG_BLENDING_NO_DEBUG
#ifndef QIONG_BLENDING_NO_DEBUG
	//#define QIONG_DEBUG_FILE_TO_WRITE

	#define QIONG_DEBUG_20030625
	//#define QIONG_DEBUG_20030525
	//#define QIONG_DEBUG_20030510

	//#define QIONG_BLENDING_INTERSECTION_QUAD
	//#define QIONG_BLENDING_SUB_FIGURE_QUAD_LIST
	//#define QIONG_BLENDING_MAIN_FIGURE_QUAD_LIST
	//#define QIONG_BLENDING_DILATED_ICURVE_RELATED_QUADS
	#define QIONG_BLENDING_MAIN_FIGURE_EXPANDED_QUADS_TO_REMOVE

	#define QIONG_BLENDING_BPOINT2UVT
	#define QIONG_BLENDING_DILATION
	//#define QIONG_BLENDING_DILATION_DEBUG

	#define QIONG_BLENDING_SUB_FIGURE							// display the part of sub-figure excluding blended region
	//#define QIONG_BLENDING_SUB_FIGURE_INTERSECTION_CURVE		// display the cut-curve on the main-figure
	#define QIONG_BLENDING_SUB_FIGURE_INTERSECTION			// display the part of sub-figure connected to the i-curve
	/////////////////////////////////////////////////////
	//#define QIONG_BLENDING_INTERSECTION_CURVE_DEBUG
	#define QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#define QIONG_BLENDING_INTERSECTION_CURVE_REGULARIZED
	#define QIONG_BLENDING_INTERSECTION_CURVE_DILATED
	/////////////////////////////////////////////////////
	#define QIONG_BLENDING_INTERSECTION_MAIN_FIGURE				// display the part of main-figure connected to the i-curve
	//#define QIONG_BLENDING_MAIN_FIGURE_INTERSECTION_CURVE		// display the cut-curve on the main-figure
	#define QIONG_BLENDING_MAIN_FIGURE						// display the part of main-figure excluding blended region
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// define the level of meshes used to calculate the intersection
//#define SUBDIV_LEVEL_TO_USE 1
//#define UVT_UNIT 1.0/(pow(2, SUBDIV_LEVEL_TO_USE))

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// choose which object model to use and 
// define the sets of parameters for various objects, they will be eventually written in the M3D files

// subBlending and mainBlending will be within [0.01, 1.0]

//#define PARAMETERS_FOR_PROTRUSION_OBJECT
//#define PARAMETERS_FOR_INDENTATION_OBJECT
#define PARAMETERS_FOR_LIVER_OBJECT
//#define PARAMETERS_FOR_TORUS_OBJECT

//#define PARAMETERS_FOR_TESTING_OBJECT_1
//#define PARAMETERS_FOR_TESTING_OBJECT_2

#ifdef PARAMETERS_FOR_TORUS_OBJECT
	#define SUB_FIGURE_CUTOFF_U_ORIGINAL 1.01
	#define SUB_FIGURE_CUTOFF_U_RANGE .1511
	#define SUB_FIGURE_CUTOFF_TRANSFORM .05
	#define DILATION_STEPS_ORIGINAL 200
	#define DILATION_STEPS_RANGE 80

	#define SUB_FIGURE_CUTOFF_U_ORIGINAL_2 3.99
	#define SUB_FIGURE_CUTOFF_U_RANGE_2 -.1511
	#define SUB_FIGURE_CUTOFF_TRANSFORM_2 -.05
	#define DILATION_STEPS_ORIGINAL_2 150
	#define DILATION_STEPS_RANGE_2 80

	#define EPSILON2_ORIGINAL 0.000001
	#define DILATION_DT 0.0005
	#define FACTOR_T_2_U .5
	#define FACTOR_T_2_V .3
#endif

#ifdef PARAMETERS_FOR_LIVER_OBJECT
		#define SUB_FIGURE_CUTOFF_U_ORIGINAL 1.01
		#define SUB_FIGURE_CUTOFF_U_RANGE .3511
		#define SUB_FIGURE_CUTOFF_TRANSFORM .05
		#define DILATION_STEPS_ORIGINAL 200
		#define DILATION_STEPS_RANGE 80
	#define EPSILON2_ORIGINAL 0.000001
	#define DILATION_DT 0.0008
	#define FACTOR_T_2_U .5
	#define FACTOR_T_2_V .3
#endif

#ifdef PARAMETERS_FOR_TESTING_OBJECT_1
	#define SUB_FIGURE_CUTOFF_U .35
	#define SUB_FIGURE_CUTOFF_TRANSFORM 0.55
	#define EPSILON2 0.000005
	#define DILATION_DT 0.002
	#define DILATION_STEPS 65
	#define FACTOR_T_2_U .5
	#define FACTOR_T_2_V .5
#endif

#ifdef PARAMETERS_FOR_TESTING_OBJECT_2
	#define SUB_FIGURE_CUTOFF_U 1.01
	#define SUB_FIGURE_CUTOFF_TRANSFORM 0.05
	#define EPSILON2 0.000005
	#define DILATION_DT 0.002
	#define DILATION_STEPS 30
	#define FACTOR_T_2_U .3
	#define FACTOR_T_2_V .3
#endif

#ifdef PARAMETERS_FOR_PROTRUSION_OBJECT
		#define SUB_FIGURE_CUTOFF_U_ORIGINAL .51
		#define SUB_FIGURE_CUTOFF_U_RANGE .5111
		#define SUB_FIGURE_CUTOFF_TRANSFORM 0.01
		#define DILATION_STEPS_ORIGINAL 50
		#define DILATION_STEPS_RANGE 70
	#define EPSILON2_ORIGINAL 0.00001
	#define DILATION_DT 0.001
	#define FACTOR_T_2_U .3
	#define FACTOR_T_2_V .3
#endif

#ifndef PARAMETERS_FOR_INDENTATION_OBJECT
#define IS_IT_AN_INDENTATION_MODEL false
#endif

#ifdef PARAMETERS_FOR_INDENTATION_OBJECT
	#define	IS_IT_AN_INDENTATION_MODEL true
	#define SUB_FIGURE_CUTOFF_U_ORIGINAL .51
	#define SUB_FIGURE_CUTOFF_U_RANGE .7111
	#define SUB_FIGURE_CUTOFF_TRANSFORM 0.1
	#define DILATION_STEPS_ORIGINAL 80
	#define DILATION_STEPS_RANGE 90

	#define EPSILON2_ORIGINAL 0.00001
//	#define MINIMAL_SEGMENT_LENGTH_SQUARE 0.002/(pow(2, SUBDIV_LEVEL_TO_USE+1))
	#define DILATION_DT 0.0005
	#define FACTOR_T_2_U .3
	#define FACTOR_T_2_V .3
#endif

// If the next line is commented out, all figures in the model
// will be implicitely rendered, even those with no subfigures.
// On slow systems, it may be desirable to uncomment this line.
//#define COMPOUND_FIGURES_ONLY

// The #defines below are useful for debugging.
// Normally they should be commented out.

//#define DEBUG		// Uncomment to print debugging output

//#define VERBOSE	// Uncomment to print counts

// Causes lines to be rendered showing the normals
//#define SHOW_NORMALS

// Causes normals to be computed from the triangles generated
// by the implicit render.  This is automatic when the parametric
// surface is being rendered.
//#define CHEAP_NORMALS

// Uncomment to render only the blend region of the surface
//#define DRAW_ONLY_BLEND

// Uncomment to render the cubes containing the surface rather than the surface
//#define DRAW_BOXES

// Uncomment to render the blend region of the surface in a unique color.
// If uncommented, MODULATE_COLOR must not be used.
//#define RECOLOR_BLEND

// Uncomment to render using a color gradient through time.
// If uncommented, RECOLOR_BLEND must not be used.
//#define MODULATE_COLOR

//int interPNum;
//Bpoint interPList[200];
//Bpoint interCurve[100];

//typedef enum { WIREFRAME, FILLEDPOLY, VPOINTS, NO_RENDER } RenderStyle;

#ifdef QIONG_DEBUG_FILE_TO_WRITE
	FILE *fDebug;
#endif


//#include <time.h>
//#include <sys/timeb.h>

// Uncomment to get rainbow coloration of figures
//#define COLOR_DEBUG

#define TOP_BOTTOM 0
#define CREST_U 1
#define CREST_V 2

//#define QIONG_BLENDING_PURE_SMOOTHING
#define QIONG_BLENDING_MAINFIGURE_INTERSECTION_CURVE
#define QIONG_DEBUG_20030425

#define QIONG_DEBUG_02232003
#define QIONG_DEBUG_BLENDING_02192003

// effective zero for testing of double values as differences of near-integers
//   (want a loose fit, here, not mach-eps).
const double TESTEPS = 1e-10;


// destructor()
Pointlist_serverB::~Pointlist_serverB() 
{
	int i;//, figureNum=2;

//	if(!object) delete object;

	if(pListOriginal!=NULL)
	{
		for(i=0; i<figureNum; i++)
			if(pListOriginal[i]!=NULL)
				delete pListOriginal[i];
		delete []pListOriginal;
	}

	ClearPointCloudPointers();
	ClearBoundaryTilesPointers();

	if(thisSurf!=NULL)
		delete thisSurf;

	if(subObject!=NULL)
		delete subObject;

#if defined(DRAW_ONLY_BLEND) || defined(RECOLOR_BLEND)
#endif
#if defined(DRAW_ONLY_BLEND) || defined(RECOLOR_BLEND)
#endif

	if(blendExtents!=NULL)
		delete []blendExtents;
	if(blendAmounts!=NULL)
		delete []blendAmounts;
	if(attachmentsType != NULL)
		delete []attachmentsType;
}

void Pointlist_serverB::DetermineParametersBasedOnSublvlToUse()
{
	UVT_UNIT=1.0/(pow(2.0, SUBDIV_LEVEL_TO_USE));

#ifdef PARAMETERS_FOR_TORUS_OBJECT
	MINIMAL_SEGMENT_LENGTH_SQUARE=0.004/(pow(2.0, SUBDIV_LEVEL_TO_USE+1));
	numPointsOnCutCurve=(int)(2*(pListOriginal[1]->max_v()+2)*pow(2.0, SUBDIV_LEVEL_TO_USE+1)*1.1);	// hard-coded 1.1 here
#endif

#ifdef PARAMETERS_FOR_LIVER_OBJECT
	MINIMAL_SEGMENT_LENGTH_SQUARE=0.008/(pow(2.0, SUBDIV_LEVEL_TO_USE+1));
	numPointsOnCutCurve=(int)(2*(pListOriginal[1]->max_v()+2)*pow(2.0, SUBDIV_LEVEL_TO_USE)*1.1);	// hard-coded 1.1 here
#endif

#ifdef PARAMETERS_FOR_TESTING_OBJECT_1
	MINIMAL_SEGMENT_LENGTH_SQUARE=0.0001;
	numPointsOnCutCurve=(int)(2*(pListOriginal[1]->max_v()+2)*pow(2.0, SUBDIV_LEVEL_TO_USE)*1.1);	// hard-coded 1.1 here
#endif

#ifdef PARAMETERS_FOR_TESTING_OBJECT_2
	MINIMAL_SEGMENT_LENGTH_SQUARE=0.0001;
#endif

// BUG FIXED 2003/05/25
#ifdef PARAMETERS_FOR_PROTRUSION_OBJECT
	MINIMAL_SEGMENT_LENGTH_SQUARE=0.004/6.0/(pow(2.0, SUBDIV_LEVEL_TO_USE+1));
	numPointsOnCutCurve=(int)(2*(pListOriginal[1]->max_v()+2)*pow(2.0, SUBDIV_LEVEL_TO_USE)*1.1);	// hard-coded 1.1 here
#endif

#ifdef PARAMETERS_FOR_INDENTATION_OBJECT
	MINIMAL_SEGMENT_LENGTH_SQUARE=0.002/(pow(2.0, SUBDIV_LEVEL_TO_USE+1));
	numPointsOnCutCurve=(int)(2*(pListOriginal[1]->max_v()+2)*pow(2.0, SUBDIV_LEVEL_TO_USE+1)*1.1);	// hard-coded 1.1 here
#endif
}

/********************************************************************************/
/* init() -- initialize the Pointlist_serverB, passing in a M3DObject			*/
/********************************************************************************/
void Pointlist_serverB::init(const M3DObject *obj, int subLevel)
{
	int i, j, figureId;
	M3DQuadFigure * figure;
	Xferlist * xferList;

	SUBDIV_LEVEL_TO_USE=subLevel;

#ifdef DEBUG
        cout << "Pointlist_serverB::init() called" << endl;
#endif
	if(obj==NULL)
	{
		fprintf(stderr, "NULL pointer of the M3DObject object in Pointlist_serverB::init\n");
		return;
	}

	object = obj;
//    object = obj->assign();
//    level = subLevel;
	figureNum=obj->getFigureCount();

#if defined(DRAW_ONLY_BLEND) || defined(RECOLOR_BLEND)
#endif

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
	if(attachmentsType!=NULL)
	{
		delete []attachmentsType;
		attachmentsType=NULL;
	}

    M3DFigureTreeNode * root;
    M3DFigureTreeNode * node;
	figureId=0;
    if(figureNum>0) 
	{
		// For simplicity, this allocates more members than the necessary (tree_depth - 1)
		blendExtents = new float[figureNum];
		blendAmounts = new float[figureNum];
		attachmentsType = new SubfigureAttachment_t[figureNum];

		numTrees = obj->getFigureTreeCount();
		for(i=0; i<numTrees; i++)
		{
			root=object->getFigureTreeRoot(i);
			blendExtents[figureId] = (float) root->getBlendExtent();
			blendAmounts[figureId] = (float) root->getBlendAmount();
			attachmentsType[figureId++] = (SubfigureAttachment_t)root->getAttachmentMode();
			for(j=0; j<root->getChildCount(); j++)
			{
				node=root->getChild(j);
				blendExtents[figureId] = (float) node->getBlendExtent();
				blendAmounts[figureId] = (float) node->getBlendAmount();
				attachmentsType[figureId++] = (SubfigureAttachment_t)node->getAttachmentMode();
			}
		}
	}
	else
	{
		numTrees = 0;
		blendExtents = NULL;
		blendAmounts = NULL;
		attachmentsType = NULL;
	}

	if(pListOriginal!=NULL)
	{
		for(i=0; i<figureNum; i++)
			if(pListOriginal[i]!=NULL)
			{
				figure=dynamic_cast<M3DQuadFigure*>(object->getFigurePtr(i));
				xferList=convertM3DtoXfer(figure);

				if(i==0 && attachmentsType[1] == INDENT)
				//if(i==0 && IS_IT_AN_INDENTATION_MODEL)
				{
					pListOriginal[i]->UpdateSubdivBoundaryTiles(SUBDIV_LEVEL_TO_USE+1, xferList);
				}
				else
					pListOriginal[i]->UpdateSubdivBoundaryTiles(SUBDIV_LEVEL_TO_USE, xferList);

				delete []xferList->atomlist;
				delete xferList;
			}
	}
	else
	{
		pListOriginal=new Pointlist_server2 *[figureNum];
		for(i=0; i<figureNum; i++)
		{
			figure=dynamic_cast<M3DQuadFigure*>(object->getFigurePtr(i));
			xferList=convertM3DtoXfer(figure);
			pListOriginal[i]=new Pointlist_server2;
			pListOriginal[i]->init(xferList);

			if(i==0 && attachmentsType[1] == INDENT)
			//if(i==0 && IS_IT_AN_INDENTATION_MODEL)
			{
				pListOriginal[i]->ComputeSubdivBoundaryTiles(SUBDIV_LEVEL_TO_USE+1);
			}
			else
			{
				pListOriginal[i]->ComputeSubdivBoundaryTiles(SUBDIV_LEVEL_TO_USE);
			}
			delete []xferList->atomlist;
			delete xferList;
		}
	}

	DetermineParametersBasedOnSublvlToUse();
}

/********************************************************************************/
/* update_mesh() -- update the newgrid with the given medial atom at the given	*/
/*		location, and update the interpolator for this changed location.		*/
/********************************************************************************/
void Pointlist_serverB::update_mesh(int figureId, XferAtom *thisXatom, int Xfer_col, int Xfer_row)
{
	if(ValidFigureID(figureId))
		//newGrid[figureId]->update_mesh(thisXatom, Xfer_col, Xfer_row);
		pListOriginal[figureId]->update_mesh(thisXatom, Xfer_col, Xfer_row);
	else
		fprintf(stderr, "Incorrect 'figureId' in Pointlist_serverB::update_mesh!\n");
}

void Pointlist_serverB::update_mesh(int figureId, Diatom *thisDiatom, int u_col, int v_row)
{
	if(ValidFigureID(figureId))
		//newGrid[figureId]->update_mesh(thisDiatom, u_col, v_row);
		pListOriginal[figureId]->update_mesh(thisDiatom, u_col, v_row);
	else
		fprintf(stderr, "Incorrect 'figureId' in Pointlist_serverB::update_mesh!\n");
}

/********************************************************************************/
/* etafun() -- compute a multiplier for the etavec for a given theta and		*/
/*		a given rotation rotangle from v1 to v2 where angle(v1,v2) = 2*theta.	*/
/********************************************************************************/
inline double Pointlist_serverB::etafun(double eta, double theta, double rotangle)
{
	// try a simple sinusoid---window a cosine between -pi and pi, shift and scale
	// So, rotangle goes from 0 to 2*theta, so subtract theta and multiply by pi/theta

	double sinusend = ((rotangle - theta)*M_PI/theta);
	double multiplier = (cos(sinusend) + 1.0) / 2.0;
	return 1.0 + multiplier*(eta - 1.0);
}


// some functions never called!
// as above, but only lerps rad and (u, v, t)
void lerp_bpoint2_B(Bpoint *Bpnt, double t, Bpoint *ptA, Bpoint *ptB)
{
	Bpnt->rad = (1.0 - t)*ptA->rad + t*ptB->rad;
	Bpnt->u = (1.0 - t)*ptA->u + t*ptB->u;
	Bpnt->v = (1.0 - t)*ptA->v + t*ptB->v;
	Bpnt->t = (1.0 - t)*ptA->t + t*ptB->t;
}

/********************************************************************************
/* interpolate_quad2_B() -- as above, but only interpolates rad					*
/*	 and (u, v, t),	not entire Bpoint.											*
/********************************************************************************/
//void Pointlist_serverB::interpolate_quad2_B(Bpoint *Bpnt, double uval, double vval, 
void interpolate_quad2_B(Bpoint *Bpnt, double uval, double vval, 
		                  Bpoint *ptA, Bpoint *ptB, Bpoint *ptC, Bpoint *ptD)
{
	Bpoint workpointAB, workpointAC, workpointAD;
	double u_t, v_t;

	// if both are zero, just return values from ptA
	if (uval == 0.0 && vval == 0.0) {
		Bpnt->rad = ptA->rad;
		Bpnt->u = ptA->u;
		Bpnt->v = ptA->v;
		Bpnt->t = ptA->t;
	}
	else if (uval > vval) {
		// compute patch coordinates
		u_t = uval;
		v_t = vval/uval;

		lerp_bpoint2_B(&workpointAC, u_t, ptA, ptC);
		lerp_bpoint2_B(&workpointAD, u_t, ptA, ptD);
		lerp_bpoint2_B(Bpnt, v_t, &workpointAD, &workpointAC);
	}
	else {
		// compute patch coordinates
		u_t = uval/vval;
		v_t = vval;

		lerp_bpoint2_B(&workpointAC, v_t, ptA, ptC);
		lerp_bpoint2_B(&workpointAB, v_t, ptA, ptB);
		lerp_bpoint2_B(Bpnt, u_t, &workpointAB, &workpointAC);
	}
}

/********************************************************************************/
/* subdivNearpoint() -- return the nearpoint on the subdivided tileset to the	*/
/*    input (x0, y0, z0)														*/
/********************************************************************************/
double Pointlist_serverB::subdivNearpoint(double x0, double y0, double z0,
										  DbVector3 *nearpointXYZ, DbVector3 *nearpointUVT,
										  int *multiplicity)
{
	double mindistsqr = DBL_MAX;
	double distsqr;
	DbVector3 thispoint(x0, y0, z0), diffvec;
	DbVector3 minpointXYZ;
	DbVector3 minpointUVT;
	DbVector3 testpoint;

	for (int vert = 0; vert < numsubdivpoints; vert++) {
		testpoint = DbVector3(subdivbpoints[vert].pnt[0], subdivbpoints[vert].pnt[1], subdivbpoints[vert].pnt[2]);
		diffvec = testpoint - thispoint;

		if ((distsqr = diffvec.dot(diffvec)) < mindistsqr) {
			mindistsqr = distsqr;
			minpointXYZ = testpoint;
			minpointUVT = DbVector3(subdivbpoints[vert].u, subdivbpoints[vert].v, subdivbpoints[vert].t);
		}
	}
	*nearpointXYZ = minpointXYZ;
	*nearpointUVT = minpointUVT;
	*multiplicity = 1;

	return sqrt(mindistsqr);
}

/********************************************************************************/
/* subdivNearQuads() -- return the nearpoint on the subdivided tileset to the	*/
/*    input (x0, y0, z0)														*/
/********************************************************************************/
double Pointlist_serverB::subdivNearQuads(double x0, double y0, double z0,
										  DbVector3 *nearpointXYZ, DbVector3 *nearpointUVT,
										  int *multiplicity)
{
	double mindistsqr = DBL_MAX;
	double distsqr;
	DbVector3 thispoint(x0, y0, z0), diffvec;
	DbVector3 minpointXYZ;
	DbVector3 minpointUVT;
	DbVector3 testpoint;

	for (int quad = 0; quad < numsubdivtiles; quad++) {
		for (int vert = 0; vert < 4; vert++) {
			testpoint = DbVector3(subdivtileset[4*quad + vert].pnt[0], subdivtileset[4*quad + vert].pnt[1],
				                  subdivtileset[4*quad + vert].pnt[2]);
			diffvec = testpoint - thispoint;

			if ((distsqr = diffvec.dot(diffvec)) < mindistsqr) {
				mindistsqr = distsqr;
				minpointXYZ = testpoint;
				minpointUVT = DbVector3(subdivtileset[4*quad + vert].u, subdivtileset[4*quad + vert].v,
					                    subdivtileset[4*quad + vert].t);
			}
		}
	}
	*nearpointXYZ = minpointXYZ;
	*nearpointUVT = minpointUVT;
	*multiplicity = 1;

	return sqrt(mindistsqr);
} 

// if thispoint is on the positive side of triangle (A-(A+nA)-B or A-(B+nB)-B)
//     and likewise for the other 3 edges BC, CD, and DA, then thispoint
//     has a possible Phong nearpoint on the quad, else return FALSE
bool possible_phongneartile_B(const DbVector3& thispoint,
		const DbVector3& A, const DbVector3& B, const DbVector3& C, const DbVector3& D,
		const DbVector3& nA, const DbVector3& nB, const DbVector3& nC, const DbVector3& nD)
{
	bool retval;
	if (thispoint.outsideface(A, A + nA, B) == -1 && thispoint.outsideface(A, B + nB, B) == -1)
		retval = false;
	else if (thispoint.outsideface(C, C + nC, D) == -1 && thispoint.outsideface(C, D + nD, D) == -1)
		retval = false;
	else if (thispoint.outsideface(B, B + nB, C) == -1 && thispoint.outsideface(B, C + nC, C) == -1)
		retval = false;
	else if (thispoint.outsideface(D, D + nD, A) == -1 && thispoint.outsideface(D, A + nA, A) == -1)
		retval = false;
	else
		retval = true;

	return retval;
}


// 2003/06/10
/********************************************************************************/
/* subdivNearSurface() -- return the nearpoint on the subdivided tileset to the	*/
/*    input (x0, y0, z0) -- first finds the nearest quad-vertex, then computes	*/
/*	  the point on the bilinearly interpolated quad with a Phong normal nearest */
/*	  angularly from the point in space.										*/
/* Returns a bpoint containing XYZ and UVT position on surface, phong normal	*/
/*	  for this position (normalized?  no, I don't think so.) and interpolated	*/
/*	  radius; and return and distance to point in space FROM SURFACE.			*/
/*	  NOTE:  This is not tau, though could be if it	were better so.				*/
/* outsideFLAG:  if NULL, no action												*/
/*				  if non-NULL, set true if point is outside of surface,			*/
/*				  else set false.												*/
/********************************************************************************/
bool Pointlist_serverB::subdivNearSurface(double x0, double y0, double z0, 
										  Bpoint *nearpoint, double *neardist, bool *outsideFLAG)
{
	DbVector3 thispoint(x0, y0, z0), diffvec;
	DbVector3 testpoint;

	IntQueue nearqueue;		// queue container for int values to be searched
	IntQueue checkedqueue;	// second queue to search if no values found during first search of nearqueue
	IntQueue checkedqueue2;	// third queue to search if no values found during first two searches

	if (OldnearsetA.maxsize() != numsubdivtiles) {
		// init all sets for the first time
		OldnearsetA.init(numsubdivtiles);
		nearset.init(numsubdivtiles);
		neighborset.init(numsubdivtiles);
		neighborset2.init(numsubdivtiles);
		neighborset3.init(numsubdivtiles);
	}

	nearset.clear();
	neighborset.clear();
	neighborset2.clear();
	neighborset3.clear();

	// Uncomment to DISABLE CACHING
	//	OldnearsetA.clear();

	int pos;
	DbVector3 A, B, C, D;
	DbVector3 nA, nB, nC, nD;
	int vert, quad, quadpos;
	double mindistsqr = DBL_MAX;
	double distsqr;

	// if Oldnearset is not empty, use its quads to initialize the routine
	// ---if none is a possible quad, do the brute force search
	if (!OldnearsetA.empty()) {

		for (pos = 0; pos < OldnearsetA.size(); pos++) {

			quad = OldnearsetA[pos];
			quadpos = 4*quad;

			A.set(subdivtileset[quadpos].pnt);
			nA.set(subdivtileset[quadpos].norm);
			B.set(subdivtileset[quadpos + 1].pnt);
			nB.set(subdivtileset[quadpos + 1].norm);
			C.set(subdivtileset[quadpos + 2].pnt);
			nC.set(subdivtileset[quadpos + 2].norm);
			D.set(subdivtileset[quadpos + 3].pnt);
			nD.set(subdivtileset[quadpos + 3].norm);

			if (possible_phongneartile_B(thispoint, A, B, C, D, nA, nB, nC, nD)) {

				nearset.insert(quad);
				mindistsqr = MIN(mindistsqr, (A - thispoint).length());
				mindistsqr = MIN(mindistsqr, (B - thispoint).length());
				mindistsqr = MIN(mindistsqr, (C - thispoint).length());
				mindistsqr = MIN(mindistsqr, (D - thispoint).length());
			}
		}

		// if none of the above tiles are now possible, check their neighbors
		if (nearset.empty()) {

			for (pos = 0; pos < OldnearsetA.size(); pos++) {

				int tnum;
				for (int i = 0; i < MAX_TILE_NEIGHBORS; i++) {
					if ((tnum = tneighbors[OldnearsetA[pos]][i]) != -1) {

						if (!OldnearsetA.hasmember(tnum))
							neighborset.insert(tnum);
					}		
				}
			}

			// Now check the neighbors
			for (pos = 0; pos < neighborset.size(); pos++) {

				quad = neighborset[pos];
				quadpos = 4*quad;

				A.set(subdivtileset[quadpos].pnt);
				nA.set(subdivtileset[quadpos].norm);
				B.set(subdivtileset[quadpos + 1].pnt);
				nB.set(subdivtileset[quadpos + 1].norm);
				C.set(subdivtileset[quadpos + 2].pnt);
				nC.set(subdivtileset[quadpos + 2].norm);
				D.set(subdivtileset[quadpos + 3].pnt);
				nD.set(subdivtileset[quadpos + 3].norm);

				if (possible_phongneartile_B(thispoint, A, B, C, D, nA, nB, nC, nD)) {

					nearset.insert(quad);
					mindistsqr = MIN(mindistsqr, (A - thispoint).length());
					mindistsqr = MIN(mindistsqr, (B - thispoint).length());
					mindistsqr = MIN(mindistsqr, (C - thispoint).length());
					mindistsqr = MIN(mindistsqr, (D - thispoint).length());
				}
			}
		}

		// if none of the above tiles are now possible, check their neighbors
		if (nearset.empty()) {

			for (pos = 0; pos < neighborset.size(); pos++) {

				int tnum;
				for (int i = 0; i < MAX_TILE_NEIGHBORS; i++) {
					if ((tnum = tneighbors[neighborset[pos]][i]) != -1) {

						if ((!OldnearsetA.hasmember(tnum)) && (!neighborset.hasmember(tnum)))
							neighborset2.insert(tnum);
					}		
				}
			}

			// Now check the neighbors
			for (pos = 0; pos < neighborset2.size(); pos++) {

				quad = neighborset2[pos];
				quadpos = 4*quad;

				A.set(subdivtileset[quadpos].pnt);
				nA.set(subdivtileset[quadpos].norm);
				B.set(subdivtileset[quadpos + 1].pnt);
				nB.set(subdivtileset[quadpos + 1].norm);
				C.set(subdivtileset[quadpos + 2].pnt);
				nC.set(subdivtileset[quadpos + 2].norm);
				D.set(subdivtileset[quadpos + 3].pnt);
				nD.set(subdivtileset[quadpos + 3].norm);

				if (possible_phongneartile_B(thispoint, A, B, C, D, nA, nB, nC, nD)) {

					nearset.insert(quad);
					mindistsqr = MIN(mindistsqr, (A - thispoint).length());
					mindistsqr = MIN(mindistsqr, (B - thispoint).length());
					mindistsqr = MIN(mindistsqr, (C - thispoint).length());
					mindistsqr = MIN(mindistsqr, (D - thispoint).length());
				}
			}
		}

		// if none of the above tiles are now possible, check their neighbors
		if (nearset.empty()) {

			for (pos = 0; pos < neighborset2.size(); pos++) {

				int tnum;
				for (int i = 0; i < MAX_TILE_NEIGHBORS; i++) {
					if ((tnum = tneighbors[neighborset2[pos]][i]) != -1) {

						if ((!OldnearsetA.hasmember(tnum)) && (!neighborset.hasmember(tnum))
							                               && (!neighborset2.hasmember(tnum)))
							neighborset3.insert(tnum);
					}		
				}
			}


			// Now check the neighbors
			for (pos = 0; pos < neighborset3.size(); pos++) {

				quad = neighborset3[pos];
				quadpos = 4*quad;

				A.set(subdivtileset[quadpos].pnt);
				nA.set(subdivtileset[quadpos].norm);
				B.set(subdivtileset[quadpos + 1].pnt);
				nB.set(subdivtileset[quadpos + 1].norm);
				C.set(subdivtileset[quadpos + 2].pnt);
				nC.set(subdivtileset[quadpos + 2].norm);
				D.set(subdivtileset[quadpos + 3].pnt);
				nD.set(subdivtileset[quadpos + 3].norm);

				if (possible_phongneartile_B(thispoint, A, B, C, D, nA, nB, nC, nD)) {

					nearset.insert(quad);
					mindistsqr = MIN(mindistsqr, (A - thispoint).length());
					mindistsqr = MIN(mindistsqr, (B - thispoint).length());
					mindistsqr = MIN(mindistsqr, (C - thispoint).length());
					mindistsqr = MIN(mindistsqr, (D - thispoint).length());
				}
			}
		}

	}


/*		DbVector3 testpoint2;
		testpoint2 = DbVector3(1.0, 1.0, 1.0);
		testpoint2.printvals("testpoint2 should be [1.0, 1.0, 1.0] but is:");
 */
	if (nearset.empty()) {

		OldnearsetA.clear();

		for (quad = 0; quad < numsubdivtiles; quad++) {
			for (vert = 0; vert < 4; vert++) {

				testpoint.set(subdivtileset[4*quad + vert].pnt);
				diffvec = testpoint - thispoint;

				if ((distsqr = diffvec.dot(diffvec)) <= mindistsqr) {

					// if strictly less than, dump previous minvals from set
					if (distsqr < mindistsqr) {
						nearset.clear();
						mindistsqr = distsqr;
					}

					nearset.insert(quad);
				}
			}
		}
	}

	// use this to test for bad phong normal intersection caused by being on the
	//   wrong side of the object!
	double mindist = sqrt(mindistsqr);

	// LOAD nearqueue from nearset
	for (pos = 0; pos < nearset.size(); pos++) 
		nearqueue.push(nearset[pos]);

	DbVector3 uvtA, uvtB, uvtC, uvtD;
	DbVector3 pbar, nvec;
	double u, v;
	double rA, rB, rC, rD;
	bool rootfound = false;
	Zerofinder phong_roots;

	Bpoint uvtvals;
	// CHANGED:  In counterclockwise order, vertices run A B C D.  
	//    ---ignore minvert and compute consistently for tile, based on ABC and ACD triangles,
	//    exactly as splitting is done for (u, v, t) interpolation across quad.
	while (!nearqueue.empty()) {

		quad = nearqueue.front();
		quadpos = 4*quad;

		nearqueue.pop();

		A.set(subdivtileset[quadpos].pnt);
		nA.set(subdivtileset[quadpos].norm);
		uvtA.set(subdivtileset[quadpos].u, subdivtileset[quadpos].v, subdivtileset[quadpos].t);
		rA = subdivtileset[quadpos].rad;

		B.set(subdivtileset[quadpos + 1].pnt);
		nB.set(subdivtileset[quadpos + 1].norm);
		uvtB.set(subdivtileset[quadpos + 1].u, subdivtileset[quadpos + 1].v, subdivtileset[quadpos + 1].t);
		rB = subdivtileset[quadpos + 1].rad;

		C.set(subdivtileset[quadpos + 2].pnt);
		nC.set(subdivtileset[quadpos + 2].norm);
		uvtC.set(subdivtileset[quadpos + 2].u, subdivtileset[quadpos + 2].v, subdivtileset[quadpos + 2].t);
		rC = subdivtileset[quadpos + 2].rad;

		D.set(subdivtileset[quadpos + 3].pnt);
		nD.set(subdivtileset[quadpos + 3].norm);
		uvtD.set(subdivtileset[quadpos + 3].u, subdivtileset[quadpos + 3].v, subdivtileset[quadpos + 3].t);
		rD = subdivtileset[quadpos + 3].rad;

		if (possible_phongneartile_B(thispoint, A, B, C, D, nA, nB, nC, nD)) {
			checkedqueue.push(quad);

			// NOTE: phong_roots returns (1 - v, u) for BAC triangle, and (1 - u, v) for DAC triangle
			//   The root-finder uses a different convention than the quad-interpolator.
			//   See the external documentation (Thall dissertation 2001) to make sense of this.
			rootfound = phong_roots.findroot4(&uvtvals, thispoint,
											  C, A, B, nC, nA, nB, uvtC, uvtA, uvtB, rC, rA, rB);

			// compute true uvt for this quad
			if (rootfound) {
				pbar.set(uvtvals.pnt);
				nvec.set(uvtvals.norm);

				// check to make sure are near approximate length, to prevent false nearpoint
				//   at antipodes of model where phong normal is parallel to nearside
				// BETTER TEST --- if length is greater than r for an internal point,
				//    should be looking at other side
				DbVector3 pvec = thispoint - pbar;
				if ((nvec.dot(pvec) < 0) && pvec.length() > uvtvals.rad*1.2)

					rootfound = false;
				else {
					OldnearsetA.clear();
					OldnearsetA.insert(quad);
					break;
				}
			}

			// ELSE try other triangular patch
			rootfound = phong_roots.findroot4(&uvtvals, thispoint,
											  D, A, C, nD, nA, nC, uvtD, uvtA, uvtC, rD, rA, rC);

			// compute true uvt for this quad
			if (rootfound) {
				pbar.set(uvtvals.pnt);
				nvec.set(uvtvals.norm);

				// check to make sure are near approximate length, to prevent false nearpoint
				//   at antipodes of model where phong normal is parallel to nearside
				// BETTER TEST --- if length is greater than r for an internal point,
				//    should be looking at other side
				DbVector3 pvec = thispoint - pbar;
				if ((nvec.dot(pvec) < 0) && pvec.length() > uvtvals.rad*1.2)

					rootfound = false;
				else {
					OldnearsetA.clear();
					OldnearsetA.insert(quad);
					break;
				}
			}
		}

		// If here, no roots found, so add neighbors of this tile to end of queue
		//   if not already in search_set, and add to search_set as well
		int tnum;
		for (int i = 0; i < MAX_TILE_NEIGHBORS; i++) {
			if ((tnum = tneighbors[quad][i]) != -1) {

				if (!nearset.hasmember(tnum)) {
					nearset.insert(tnum);
					nearqueue.push(tnum);
				}
			}
		}
	}

	// if root not found, go through checkedqueue until empty or root-found
	if (!rootfound) {

		while (!checkedqueue.empty()) {

			quad = checkedqueue.front();
			quadpos = 4*quad;
			checkedqueue2.push(quad);

			checkedqueue.pop();

			A.set(subdivtileset[quadpos].pnt);
			nA.set(subdivtileset[quadpos].norm);
			B.set(subdivtileset[quadpos + 1].pnt);
			nB.set(subdivtileset[quadpos + 1].norm);
			C.set(subdivtileset[quadpos + 2].pnt);
			nC.set(subdivtileset[quadpos + 2].norm);
			D.set(subdivtileset[quadpos + 3].pnt);
			nD.set(subdivtileset[quadpos + 3].norm);

			// NOTE: phong_roots returns (1 - v, u) for BAC triangle, and (1 - u, v) for DAC triangle
			//   The root-finder uses a different convention than the quad-interpolator.
			//   See the external documentation (Thall dissertation 2001) to make sense of this.
			rootfound =  phong_roots.findroot(v, u, pbar, nvec, thispoint, B, A, C, nB, nA, nC);
			v = 1.0 - v;
			// Convert (u, v) from triangular patch coordinates to tile coordinates
			//   (see routine interpolate_quad() above for inverse of this, and Thall diss for details)
			// NOTE:  we undo this when we call interpolate_quad2_B()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			u *= v;

			// compute true uvt for this quad
			if (rootfound) {
				interpolate_quad2_B(&uvtvals, u, v, 
								  &subdivtileset[quadpos],
								  &subdivtileset[quadpos + 1],
								  &subdivtileset[quadpos + 2],
								  &subdivtileset[quadpos + 3]);

				// BETTER TEST --- if length is greater than r for an internal point,
				//    should be looking at other side
				DbVector3 pvec = thispoint - pbar;
				if ((nvec.dot(pvec) < 0) && pvec.length() > uvtvals.rad*1.2)

					rootfound = false;
				else {
					OldnearsetA.clear();
					OldnearsetA.insert(quad);
					break;
				}
			}

			// ELSE try other triangular patch

			rootfound = phong_roots.findroot(u, v, pbar, nvec, thispoint, D, A, C, nD, nA, nC);
			u = 1.0 - u;

			// Convert (u, v) from triangular patch coordinates to tile coordinates
			//   (see routine interpolate_quad() above for inverse of this, and Thall diss for details)
			// NOTE:  we undo this when we call interpolate_quad2_B()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			v *= u;

			// compute true uvt for this quad
			if (rootfound) {

				interpolate_quad2_B(&uvtvals, u, v, 
								  &subdivtileset[quadpos],
							  	  &subdivtileset[quadpos + 1],
								  &subdivtileset[quadpos + 2],
								  &subdivtileset[quadpos + 3]);

				// BETTER TEST --- if length is greater than r for an internal point,
				//    should be looking at other side.
				DbVector3 pvec = thispoint - pbar;
				if ((nvec.dot(pvec) < 0) && pvec.length() > uvtvals.rad*1.2)

					rootfound = false;
				else {
					OldnearsetA.clear();
					OldnearsetA.insert(quad);
					break;
				}
			}

		}
	}

	// if root not found, go through checkedqueue until empty or root-found
	if (!rootfound) {

		while (!checkedqueue2.empty()) {

			quad = checkedqueue2.front();
			quadpos = 4*quad;

			checkedqueue2.pop();

			A.set(subdivtileset[quadpos].pnt);
			nA.set(subdivtileset[quadpos].norm);
			B.set(subdivtileset[quadpos + 1].pnt);
			nB.set(subdivtileset[quadpos + 1].norm);
			C.set(subdivtileset[quadpos + 2].pnt);
			nC.set(subdivtileset[quadpos + 2].norm);
			D.set(subdivtileset[quadpos + 3].pnt);
			nD.set(subdivtileset[quadpos + 3].norm);

			// NOTE: phong_roots returns (1 - v, u) for BAC triangle, and (1 - u, v) for DAC triangle
			//   The root-finder uses a different convention than the quad-interpolator.
			//   See the external documentation (Thall dissertation 2001) to make sense of this.
			rootfound =  phong_roots.findroot3(v, u, pbar, nvec, thispoint, B, A, C, nB, nA, nC);
			v = 1.0 - v;
			// Convert (u, v) from triangular patch coordinates to tile coordinates
			//   (see routine interpolate_quad() above for inverse of this, and Thall diss for details)
			// NOTE:  we undo this when we call interpolate_quad2_B()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			u *= v;

			// compute true uvt for this quad
			if (rootfound) {
				interpolate_quad2_B(&uvtvals, u, v, 
								  &subdivtileset[quadpos],
								  &subdivtileset[quadpos + 1],
								  &subdivtileset[quadpos + 2],
								  &subdivtileset[quadpos + 3]);

				// check to make sure are near approximate length, to prevent false nearpoint
				//   at antipodes of model where phong normal is parallel to nearside
				// BETTER TEST --- if length is greater than r for an internal point,
				//    should be looking at other side
				DbVector3 pvec = thispoint - pbar;
				if ((nvec.dot(pvec) < 0) && pvec.length() > uvtvals.rad*1.2)

					rootfound = false;
				else {
					OldnearsetA.clear();
					OldnearsetA.insert(quad);
					break;
				}
			}

			// ELSE try other triangular patch

			rootfound = phong_roots.findroot3(u, v, pbar, nvec, thispoint, D, A, C, nD, nA, nC);
			u = 1.0 - u;

			// Convert (u, v) from triangular patch coordinates to tile coordinates
			//   (see routine interpolate_quad() above for inverse of this, and Thall diss for details)
			// NOTE:  we undo this when we call interpolate_quad2_B()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			v *= u;

			// compute true uvt for this quad
			if (rootfound) {

				interpolate_quad2_B(&uvtvals, u, v, 
								  &subdivtileset[quadpos],
								  &subdivtileset[quadpos + 1],
								  &subdivtileset[quadpos + 2],
								  &subdivtileset[quadpos + 3]);

				// check to make sure are near approximate length, to prevent false nearpoint
				//   at antipodes of model where phong normal is parallel to nearside
				// BETTER TEST --- if length is greater than r for an internal point,
				//    should be looking at other side
				DbVector3 pvec = thispoint - pbar;
				if ((nvec.dot(pvec) < 0) && pvec.length() > uvtvals.rad*1.2)

					rootfound = false;
				else {		
					OldnearsetA.clear();
					OldnearsetA.insert(quad);
					break;
				}
			}

		}
	}

	if (rootfound) {
		pbar.put(nearpoint->pnt);
		nvec.put(nearpoint->norm);
		nearpoint->rad = uvtvals.rad;
		nearpoint->u = uvtvals.u;
		nearpoint->v = uvtvals.v;
		nearpoint->t = uvtvals.t;
		*neardist = (pbar - thispoint).length();

		// if passed in an outsideFLAG pointer, compute sidedness of point
		if (outsideFLAG != NULL)
			*outsideFLAG = (((pbar - thispoint).dot(nvec) < 0) ? true : false);

		return true;
	}
	else {
		fprintf(stderr, "ERR::Pointlist_serverB::subdivNearSurface()::no solution found\n");
		// If NONE of the above tests succeeded, have a failure.  In this case,
		//   return false, and set the Boundary point to the nearest vertex
		OldnearsetA.clear();
		mindistsqr = DBL_MAX;
		int mintilevert = -1;
		DbVector3 mindiffvec;

		for (quad = 0; quad < numsubdivtiles; quad++) {
			for (vert = 0; vert < 4; vert++) {
				testpoint.set(subdivtileset[4*quad + vert].pnt);
				diffvec = testpoint - thispoint;

				if ((distsqr = diffvec.dot(diffvec)) < mindistsqr) {

					// if strictly less than, dump previous minvals from set
					mindistsqr = distsqr;
					mintilevert = 4*quad + vert;
					mindiffvec = diffvec;
				}
			}
		}
		*nearpoint = subdivtileset[mintilevert];
		*neardist = sqrt(mindistsqr);
		// if passed in an outsideFLAG pointer, compute sidedness of point
		if (outsideFLAG != NULL)
			*outsideFLAG = ((mindiffvec.dot(DbVector3(nearpoint->pnt[0], nearpoint->pnt[1], nearpoint->pnt[2])) < 0) 
			               ? true : false);

		return false;
	}
}

inline bool rightsided(bool tpositive, double At, double Bt, double Ct, double Dt)
{
	if (tpositive)
		return (At >= 0.0 && Bt >= 0.0 && Ct >= 0.0 && Dt >= 0.0);
	else
		return (At <= 0.0 && Bt <= 0.0 && Ct <= 0.0 && Dt <= 0.0);
}

/********************************************************************************/
/* subdivNearSurfaceSIDED() -- as above, but only searches tiles on portion		*/
/*    of mesh, as indicated by TOP_p flag.  This will be used in code			*/
/*	  to find involute positions and from them medial atoms corresponding to	*/
/*	  a given boundary position.												*/
/* Returns a bpoint containing XYZ and UVT position on surface, phong normal	*/
/*	  for this position (normalized?  no, I don't think so.) and interpolated	*/
/*	  radius; and return and distance to point in space FROM SURFACE.			*/
/*	  NOTE:  This is not tau, though could be if it	were better so.				*/
/* outsideFLAG:  if NULL, no action												*/
/*				  if non-NULL, set true if point is outside of surface,			*/
/*				  else set false.												*/
/********************************************************************************/
bool Pointlist_serverB::subdivNearSurfaceSIDED(double x0, double y0, double z0, searchSIDE testside,
										  Bpoint *nearpoint, double *neardist, bool *outsideFLAG)
{	
	DbVector3 thispoint(x0, y0, z0), diffvec;
	DbVector3 testpoint;

	IntQueue nearqueue;	   // queue container for int values to be searched
	IntQueue checkedqueue; // second queue to search if no values found during first search of nearqueue
	IntQueue checkedqueue2; // third queue to search if no values found during first two searches

	int pos;
	int quad, quadpos, vert;
	DbVector3 A, B, C, D;
	DbVector3 nA, nB, nC, nD;
	double mindistsqr = DBL_MAX;
	double distsqr;


	bool sol_t_positive;
	switch (testside) {
	case TOPSIDE:
	case TOPCREST:
		sol_t_positive = true;  // solution should have positive t, so can limit tilesearch
		break;
	case BOTTOMSIDE:
	case BOTTOMCREST:
		sol_t_positive = false; // solution should have negative t, so can limit tilesearch
		break;
	default:
		fprintf(stderr, "in subdivNearSurfaceSIDED::ERROR -- not a valid searchSIDE\n");
	}

	if (OldnearsetB.maxsize() != numsubdivtiles) {
		// init all sets for the first time
		OldnearsetB.init(numsubdivtiles);
		nearset.init(numsubdivtiles);
		neighborset.init(numsubdivtiles);
		neighborset2.init(numsubdivtiles);
		neighborset3.init(numsubdivtiles);
	}

	nearset.clear();
	neighborset.clear();
	neighborset2.clear();
	neighborset3.clear();

	// DISABLE CACHING
//	OldnearsetB.clear();

	// if Oldnearset is not empty, use its quads to initialize the routine
	// ---if none is a possible quad, do the brute force search
	if (!OldnearsetB.empty()) {

		for (pos = 0; pos < OldnearsetB.size(); pos++) {

			quad = OldnearsetB[pos];
			quadpos = 4*quad;

			A.set(subdivtileset[quadpos].pnt);
			nA.set(subdivtileset[quadpos].norm);
			B.set(subdivtileset[quadpos + 1].pnt);
			nB.set(subdivtileset[quadpos + 1].norm);
			C.set(subdivtileset[quadpos + 2].pnt);
			nC.set(subdivtileset[quadpos + 2].norm);
			D.set(subdivtileset[quadpos + 3].pnt);
			nD.set(subdivtileset[quadpos + 3].norm);

			if (possible_phongneartile_B(thispoint, A, B, C, D, nA, nB, nC, nD)) {
				if (rightsided(sol_t_positive, subdivtileset[quadpos].t,
											   subdivtileset[quadpos + 1].t, 
											   subdivtileset[quadpos + 2].t,
											   subdivtileset[quadpos + 3].t)) {
					nearset.insert(quad);
					mindistsqr = MIN(mindistsqr, (A - thispoint).length());
					mindistsqr = MIN(mindistsqr, (B - thispoint).length());
					mindistsqr = MIN(mindistsqr, (C - thispoint).length());
					mindistsqr = MIN(mindistsqr, (D - thispoint).length());
				}
			}
		}

		// if none of the above tiles are now possible, check their neighbors
		if (nearset.empty()) {

			for (pos = 0; pos < OldnearsetB.size(); pos++) {

				int tnum;
				for (int i = 0; i < MAX_TILE_NEIGHBORS; i++) {
					if ((tnum = tneighbors[OldnearsetB[pos]][i]) != -1) {

						if (!OldnearsetB.hasmember(tnum))
							neighborset.insert(tnum);
					}		
				}
			}

			// Now check the neighbors
			for (pos = 0; pos < neighborset.size(); pos++) {

				quad = neighborset[pos];
				quadpos = 4*quad;

				A.set(subdivtileset[quadpos].pnt);
				nA.set(subdivtileset[quadpos].norm);
				B.set(subdivtileset[quadpos + 1].pnt);
				nB.set(subdivtileset[quadpos + 1].norm);
				C.set(subdivtileset[quadpos + 2].pnt);
				nC.set(subdivtileset[quadpos + 2].norm);
				D.set(subdivtileset[quadpos + 3].pnt);
				nD.set(subdivtileset[quadpos + 3].norm);

				if (possible_phongneartile_B(thispoint, A, B, C, D, nA, nB, nC, nD)) {

					if (rightsided(sol_t_positive, subdivtileset[quadpos].t,
											   subdivtileset[quadpos + 1].t, 
											   subdivtileset[quadpos + 2].t,
											   subdivtileset[quadpos + 3].t)) {
						nearset.insert(quad);
						mindistsqr = MIN(mindistsqr, (A - thispoint).length());
						mindistsqr = MIN(mindistsqr, (B - thispoint).length());
						mindistsqr = MIN(mindistsqr, (C - thispoint).length());
						mindistsqr = MIN(mindistsqr, (D - thispoint).length());
					}
				}
			}
		}

		// if none of the above tiles are now possible, check their neighbors
		if (nearset.empty()) {

			for (pos = 0; pos < neighborset.size(); pos++) {

				int tnum;
				for (int i = 0; i < MAX_TILE_NEIGHBORS; i++) {
					if ((tnum = tneighbors[neighborset[pos]][i]) != -1) {

						if ((!OldnearsetB.hasmember(tnum)) && (!neighborset.hasmember(tnum)))
							neighborset2.insert(tnum);
					}		
				}
			}

			// Now check the neighbors
			for (pos = 0; pos < neighborset2.size(); pos++) {

				quad = neighborset2[pos];
				quadpos = 4*quad;

				A.set(subdivtileset[quadpos].pnt);
				nA.set(subdivtileset[quadpos].norm);
				B.set(subdivtileset[quadpos + 1].pnt);
				nB.set(subdivtileset[quadpos + 1].norm);
				C.set(subdivtileset[quadpos + 2].pnt);
				nC.set(subdivtileset[quadpos + 2].norm);
				D.set(subdivtileset[quadpos + 3].pnt);
				nD.set(subdivtileset[quadpos + 3].norm);

				if (possible_phongneartile_B(thispoint, A, B, C, D, nA, nB, nC, nD)) {

					if (rightsided(sol_t_positive, subdivtileset[quadpos].t,
											   subdivtileset[quadpos + 1].t, 
											   subdivtileset[quadpos + 2].t,
											   subdivtileset[quadpos + 3].t)) {
						nearset.insert(quad);
						mindistsqr = MIN(mindistsqr, (A - thispoint).length());
						mindistsqr = MIN(mindistsqr, (B - thispoint).length());
						mindistsqr = MIN(mindistsqr, (C - thispoint).length());
						mindistsqr = MIN(mindistsqr, (D - thispoint).length());
					}
				}
			}
		}

		// if none of the above tiles are now possible, check their neighbors
		if (nearset.empty()) {

			for (pos = 0; pos < neighborset2.size(); pos++) {

				int tnum;
				for (int i = 0; i < MAX_TILE_NEIGHBORS; i++) {
					if ((tnum = tneighbors[neighborset2[pos]][i]) != -1) {

						if ((!OldnearsetB.hasmember(tnum)) && (!neighborset.hasmember(tnum))
							                               && (!neighborset2.hasmember(tnum)))
							neighborset3.insert(tnum);
					}		
				}
			}


			// Now check the neighbors
			for (pos = 0; pos < neighborset3.size(); pos++) {

				quad = neighborset3[pos];
				quadpos = 4*quad;

				A.set(subdivtileset[quadpos].pnt);
				nA.set(subdivtileset[quadpos].norm);
				B.set(subdivtileset[quadpos + 1].pnt);
				nB.set(subdivtileset[quadpos + 1].norm);
				C.set(subdivtileset[quadpos + 2].pnt);
				nC.set(subdivtileset[quadpos + 2].norm);
				D.set(subdivtileset[quadpos + 3].pnt);
				nD.set(subdivtileset[quadpos + 3].norm);

				if (possible_phongneartile_B(thispoint, A, B, C, D, nA, nB, nC, nD)) {

					if (rightsided(sol_t_positive, subdivtileset[quadpos].t,
												   subdivtileset[quadpos + 1].t, 
												   subdivtileset[quadpos + 2].t,
												   subdivtileset[quadpos + 3].t)) {
						nearset.insert(quad);
						mindistsqr = MIN(mindistsqr, (A - thispoint).length());
						mindistsqr = MIN(mindistsqr, (B - thispoint).length());
						mindistsqr = MIN(mindistsqr, (C - thispoint).length());
						mindistsqr = MIN(mindistsqr, (D - thispoint).length());
					}
				}
			}
		}

	}

	if (nearset.empty()) {

		OldnearsetB.clear();

		// compute number of tiles and index of first tile in top or bottom grid
		int num_tb_tiles = num_ucoltiles*num_vrowtiles;
		int num_s_tiles_over_2 = (numsubdivtiles - 2*num_tb_tiles)/2;
		int firsttile, tilecount;

		switch (testside) {
		case TOPSIDE:
		case TOPCREST:
			firsttile = 0;
			tilecount = num_tb_tiles;
			testside = TOPSIDE;			// if BOTTOMCREST, just test for TOPSIDE below (2 looping tests)
			break;
		case BOTTOMSIDE:
		case BOTTOMCREST:
			firsttile = num_tb_tiles;
			tilecount = num_tb_tiles;
			testside = BOTTOMSIDE;		// if TOPCREST, just test for TOPSIDE below (2 looping tests)
			break;
//		case TOPCREST:
//			firsttile = 2*num_tb_tiles;
//			tilecount = num_s_tiles_over_2;
//			break;
//		case BOTTOMCREST:
//			firsttile = 2*num_tb_tiles + num_s_tiles_over_2;
//			tilecount = num_s_tiles_over_2;
//			break;

		default:
			fprintf(stderr, "in subdivNearSurfaceSIDED::ERROR -- not a valid searchSIDE\n");
		}

		if (testside == TOPSIDE) {
			for (quad = 0; quad < num_tb_tiles; quad++) {
				for (vert = 0; vert < 4; vert++) {
					testpoint.set(subdivtileset[4*quad + vert].pnt);
					diffvec = testpoint - thispoint;

					if ((distsqr = diffvec.dot(diffvec)) <= mindistsqr) {

						// if strictly less than, dump previous minvals from set
						if (distsqr < mindistsqr) {
							nearset.clear();
							mindistsqr = distsqr;
						}

						nearset.insert(quad);
					}
				}
			}
		}
		if (testside == BOTTOMSIDE) {
			for (quad = num_tb_tiles; quad < 2*num_tb_tiles; quad++) {
				for (vert = 0; vert < 4; vert++) {
					testpoint.set(subdivtileset[4*quad + vert].pnt);
					diffvec = testpoint - thispoint;

					if ((distsqr = diffvec.dot(diffvec)) <= mindistsqr) {

						// if strictly less than, dump previous minvals from set
						if (distsqr < mindistsqr) {
							nearset.clear();
							mindistsqr = distsqr;
						}
						nearset.insert(quad);
					}
				}
			}
		}
		if (testside == TOPSIDE) {
			for (quad = 2*num_tb_tiles; quad < 2*num_tb_tiles + num_s_tiles_over_2; quad++) {
				for (vert = 0; vert < 4; vert++) {
					testpoint.set(subdivtileset[4*quad + vert].pnt);
					diffvec = testpoint - thispoint;

					if ((distsqr = diffvec.dot(diffvec)) <= mindistsqr) {

						// if strictly less than, dump previous minvals from set
						if (distsqr < mindistsqr) {
							nearset.clear();
							mindistsqr = distsqr;
						}	
						nearset.insert(quad);
					}
				}
			}
		}

		firsttile = 2*num_tb_tiles + num_s_tiles_over_2;
		tilecount = num_s_tiles_over_2;
		if (testside == BOTTOMSIDE) {
			for (quad = firsttile; quad < firsttile + tilecount; quad++) {
				for (vert = 0; vert < 4; vert++) {
					testpoint.set(subdivtileset[4*quad + vert].pnt);
					diffvec = testpoint - thispoint;

					if ((distsqr = diffvec.dot(diffvec)) <= mindistsqr) {

						// if strictly less than, dump previous minvals from set
						if (distsqr < mindistsqr) {
							nearset.clear();
							mindistsqr = distsqr;
						}
						nearset.insert(quad);
					}
				}
			}
		}

//		for (quad = firsttile; quad < firsttile + tilecount; quad++) {
//			for (vert = 0; vert < 4; vert++) {
//				testpoint.set(subdivtileset[4*quad + vert].pnt);
//				diffvec = testpoint - thispoint;
//				
//				if ((distsqr = diffvec.dot(diffvec)) <= mindistsqr) {
//						
//						// if strictly less than, dump previous minvals from set
//						if (distsqr < mindistsqr) {
//							nearset.clear();
//							mindistsqr = distsqr;
//						}
//						nearset.insert(quad);
//					}
//			}
//		}

	}

	// use this to test for bad phong normal intersection caused by being on the
	//   wrong side of the object!
	double mindist = sqrt(mindistsqr);

	// LOAD nearqueue from nearset
	for (pos = 0; pos < nearset.size(); pos++) 
		nearqueue.push(nearset[pos]);

	DbVector3 uvtA, uvtB, uvtC, uvtD;
	DbVector3 pbar, nvec;
	double u, v;
	double rA, rB, rC, rD;
	bool rootfound = false;
	Zerofinder phong_roots;

	Bpoint uvtvals;
	// CHANGED:  In counterclockwise order, vertices run A B C D.  
	//    ---ignore minvert and compute consistently for tile, based on ABC and ACD triangles,
	//    exactly as splitting is done for (u, v, t) interpolation across quad.
	while (!nearqueue.empty()) {

		quad = nearqueue.front();
		quadpos = 4*quad;

		nearqueue.pop();

		A.set(subdivtileset[quadpos].pnt);
		nA.set(subdivtileset[quadpos].norm);
		uvtA.set(subdivtileset[quadpos].u, subdivtileset[quadpos].v, subdivtileset[quadpos].t);
		rA = subdivtileset[quadpos].rad;

		B.set(subdivtileset[quadpos + 1].pnt);
		nB.set(subdivtileset[quadpos + 1].norm);
		uvtB.set(subdivtileset[quadpos + 1].u, subdivtileset[quadpos + 1].v, subdivtileset[quadpos + 1].t);
		rB = subdivtileset[quadpos + 1].rad;

		C.set(subdivtileset[quadpos + 2].pnt);
		nC.set(subdivtileset[quadpos + 2].norm);
		uvtC.set(subdivtileset[quadpos + 2].u, subdivtileset[quadpos + 2].v, subdivtileset[quadpos + 2].t);
		rC = subdivtileset[quadpos + 2].rad;

		D.set(subdivtileset[quadpos + 3].pnt);
		nD.set(subdivtileset[quadpos + 3].norm);
		uvtD.set(subdivtileset[quadpos + 3].u, subdivtileset[quadpos + 3].v, subdivtileset[quadpos + 3].t);
		rD = subdivtileset[quadpos + 3].rad;

		if (possible_phongneartile_B(thispoint, A, B, C, D, nA, nB, nC, nD)) {
			checkedqueue.push(quad);

			// NOTE: phong_roots returns (1 - v, u) for BAC triangle, and (1 - u, v) for DAC triangle
			//   The root-finder uses a different convention than the quad-interpolator.
			//   See the external documentation (Thall dissertation 2001) to make sense of this.
			rootfound = phong_roots.findroot4(&uvtvals, thispoint,
											  C, A, B, nC, nA, nB, uvtC, uvtA, uvtB, rC, rA, rB);

			// compute true uvt for this quad (here, already computed in findroot4)
			if (rootfound) {
				pbar.set(uvtvals.pnt);
				nvec.set(uvtvals.norm);

				// BETTER TEST --- if length is greater than r for an internal point,
				//    should be looking at other side---DON'T DO THIS TEST HERE
				//    DON'T WANT OTHER SIDE NO MATTER WHAT...only test t of solution
				if ((sol_t_positive && uvtvals.t < 0) || (!sol_t_positive && uvtvals.t > 0))
					rootfound = false;
				else {
					OldnearsetB.clear();
					OldnearsetB.insert(quad);
					break;
				}
			}

			// ELSE try other triangular patch
			rootfound = phong_roots.findroot4(&uvtvals, thispoint,
											  D, A, C, nD, nA, nC, uvtD, uvtA, uvtC, rD, rA, rC);

			// compute true uvt for this quad (here, already computed in findroot4)
			if (rootfound) {
				pbar.set(uvtvals.pnt);
				nvec.set(uvtvals.norm);

				// BETTER TEST --- if length is greater than r for an internal point,
				//    should be looking at other side---DON'T DO THIS TEST HERE
				//    DON'T WANT OTHER SIDE NO MATTER WHAT...only test t of solution
				if ((sol_t_positive && uvtvals.t < 0) || (!sol_t_positive && uvtvals.t > 0))
					rootfound = false;
				else {
					OldnearsetB.clear();
					OldnearsetB.insert(quad);
					break;
				}
			}
		}

		// If here, no roots found, so add neighbors of this tile to end of queue
		//   if not already in search_set, and add to search_set as well
		int tnum;
		for (int i = 0; i < MAX_TILE_NEIGHBORS; i++) {
			if ((tnum = tneighbors[quad][i]) != -1) {

				if (!nearset.hasmember(tnum)) {
					nearset.insert(tnum);
					nearqueue.push(tnum);
				}
			}

		}
	}

	// if root not found, go through checkedqueue until empty or root-found
	if (!rootfound) {

		while (!checkedqueue.empty()) {

			quad = checkedqueue.front();
			quadpos = 4*quad;
			checkedqueue2.push(quad);

			checkedqueue.pop();

			A.set(subdivtileset[quadpos].pnt);
			nA.set(subdivtileset[quadpos].norm);
			B.set(subdivtileset[quadpos + 1].pnt);
			nB.set(subdivtileset[quadpos + 1].norm);
			C.set(subdivtileset[quadpos + 2].pnt);
			nC.set(subdivtileset[quadpos + 2].norm);
			D.set(subdivtileset[quadpos + 3].pnt);
			nD.set(subdivtileset[quadpos + 3].norm);

			// NOTE: phong_roots returns (1 - v, u) for BAC triangle, and (1 - u, v) for DAC triangle
			//   The root-finder uses a different convention than the quad-interpolator.
			//   See the external documentation (Thall dissertation 2001) to make sense of this.
			rootfound =  phong_roots.findroot(v, u, pbar, nvec, thispoint, B, A, C, nB, nA, nC);
			v = 1.0 - v;
			// Convert (u, v) from triangular patch coordinates to tile coordinates
			//   (see routine interpolate_quad() above for inverse of this, and Thall diss for details)
			// NOTE:  we undo this when we call interpolate_quad2_B()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			u *= v;

			// compute true uvt for this quad (not computed by findroot)
			if (rootfound) {

				interpolate_quad2_B(&uvtvals, u, v, 
								  &subdivtileset[quadpos],
								  &subdivtileset[quadpos + 1],
								  &subdivtileset[quadpos + 2],
								  &subdivtileset[quadpos + 3]);

				// BETTER TEST --- if length is greater than r for an internal point,
				//    should be looking at other side---DON'T DO THIS TEST HERE
				//    DON'T WANT OTHER SIDE NO MATTER WHAT...only test t of solution
				if ((sol_t_positive && uvtvals.t < 0) || (!sol_t_positive && uvtvals.t > 0))
					rootfound = false;
				else {
					OldnearsetB.clear();
					OldnearsetB.insert(quad);
					break;
				}
			}

			// ELSE try other triangular patch

			rootfound = phong_roots.findroot(u, v, pbar, nvec, thispoint, D, A, C, nD, nA, nC);
			u = 1.0 - u;

			// Convert (u, v) from triangular patch coordinates to tile coordinates
			//   (see routine interpolate_quad() above for inverse of this, and Thall diss for details)
			// NOTE:  we undo this when we call interpolate_quad2_B()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			v *= u;

			// compute true uvt for this quad (not computed by findroot)
			if (rootfound) {

				interpolate_quad2_B(&uvtvals, u, v, 
								  &subdivtileset[quadpos],
								  &subdivtileset[quadpos + 1],
								  &subdivtileset[quadpos + 2],
								  &subdivtileset[quadpos + 3]);

				// BETTER TEST --- if length is greater than r for an internal point,
				//    should be looking at other side---DON'T DO THIS TEST HERE
				//    DON'T WANT OTHER SIDE NO MATTER WHAT...only test t of solution
				if ((sol_t_positive && uvtvals.t < 0) || (!sol_t_positive && uvtvals.t > 0))
					rootfound = false;
				else {
					OldnearsetB.clear();
					OldnearsetB.insert(quad);
					break;
				}
			}
		}
	}

	// if root not found, go through checkedqueue until empty or root-found
	if (!rootfound) {

		while (!checkedqueue2.empty()) {

			quad = checkedqueue2.front();
			quadpos = 4*quad;

			checkedqueue2.pop();

			A.set(subdivtileset[quadpos].pnt);
			nA.set(subdivtileset[quadpos].norm);
			B.set(subdivtileset[quadpos + 1].pnt);
			nB.set(subdivtileset[quadpos + 1].norm);
			C.set(subdivtileset[quadpos + 2].pnt);
			nC.set(subdivtileset[quadpos + 2].norm);
			D.set(subdivtileset[quadpos + 3].pnt);
			nD.set(subdivtileset[quadpos + 3].norm);

			// NOTE: phong_roots returns (1 - v, u) for BAC triangle, and (1 - u, v) for DAC triangle
			//   The root-finder uses a different convention than the quad-interpolator.
			//   See the external documentation (Thall dissertation 2001) to make sense of this.
			rootfound =  phong_roots.findroot3(v, u, pbar, nvec, thispoint, B, A, C, nB, nA, nC);
			v = 1.0 - v;
			// Convert (u, v) from triangular patch coordinates to tile coordinates
			//   (see routine interpolate_quad() above for inverse of this, and Thall diss for details)
			// NOTE:  we undo this when we call interpolate_quad2_B()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			u *= v;

			// compute true uvt for this quad (not computed by findroot3)
			if (rootfound) {

				interpolate_quad2_B(&uvtvals, u, v, 
								  &subdivtileset[quadpos],
								  &subdivtileset[quadpos + 1],
								  &subdivtileset[quadpos + 2],
								  &subdivtileset[quadpos + 3]);

				// BETTER TEST --- if length is greater than r for an internal point,
				//    should be looking at other side---DON'T DO THIS TEST HERE
				//    DON'T WANT OTHER SIDE NO MATTER WHAT...only test t of solution
				if ((sol_t_positive && uvtvals.t < 0) || (!sol_t_positive && uvtvals.t > 0))
					rootfound = false;
				else {
					OldnearsetB.clear();
					OldnearsetB.insert(quad);
					break;
				}
			}

			// ELSE try other triangular patch

			rootfound = phong_roots.findroot3(u, v, pbar, nvec, thispoint, D, A, C, nD, nA, nC);
			u = 1.0 - u;

			// Convert (u, v) from triangular patch coordinates to tile coordinates
			//   (see routine interpolate_quad() above for inverse of this, and Thall diss for details)
			// NOTE:  we undo this when we call interpolate_quad2_B()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			v *= u;

			// compute true uvt for this quad (not computed by findroot3)
			if (rootfound) {

				interpolate_quad2_B(&uvtvals, u, v, 
								  &subdivtileset[quadpos],
								  &subdivtileset[quadpos + 1],
								  &subdivtileset[quadpos + 2],
								  &subdivtileset[quadpos + 3]);

				// BETTER TEST --- if length is greater than r for an internal point,
				//    should be looking at other side---DON'T DO THIS TEST HERE
				//    DON'T WANT OTHER SIDE NO MATTER WHAT...only test t of solution
				if ((sol_t_positive && uvtvals.t < 0) || (!sol_t_positive && uvtvals.t > 0))
					rootfound = false;
				else {
					OldnearsetB.clear();
					OldnearsetB.insert(quad);
					break;
				}
			}	
		}
	}

	if (rootfound) {

		pbar.put(nearpoint->pnt);
		nvec.put(nearpoint->norm);
		nearpoint->rad = uvtvals.rad;
		nearpoint->u = uvtvals.u;
		nearpoint->v = uvtvals.v;
		nearpoint->t = uvtvals.t;
		*neardist = (pbar - thispoint).length();

		// if passed in an outsideFLAG pointer, compute sidedness of point
		if (outsideFLAG != NULL)
			*outsideFLAG = (((pbar - thispoint).dot(nvec) < 0) ? true : false);

		return true;
	}
	else {
		// If NONE of the above tests succeeded, have a failure.  
		// In this case, return the NON-sided nearpoint using subdivNearSurface()
		//   and return false.
//		OldnearsetB.clear();

//		subdivNearSurface(x0, y0, z0, nearpoint, neardist, outsideFLAG);

		return false;
	}
}

/********************************************************************************/
/* subdivBpoint_to_Medial() -- using above function in a bisection search,		*/
/*	  find the opposing boundary involute (point on opposite side equidistant	*/
/*    with boundpnt from a medial point along a normal from each) and compute	*/
/*    the medial atom corresponding to these involutes							*/
/********************************************************************************/
/********************************************************************************/
/* NOTE:  this routine fails if involute is along same side, possible if folded	*/
/*    and thick object.  No reason to limit involutes as here, except didn't	*/
/*    want to have to check that weren't getting degeneracies nearby original	*/
/*    point.  Can fix this if need be, but mostly a proof-of-concept right now.	*/
/********************************************************************************/

// search along normal to surface at distance inside in interval
//    [r(1 - EPSmedial), r(1 + EPSmedial)]
const double EPSmedial = 0.1;
const double medTOL = 0.05;
const double paramEPS = 1e-10; // if NEARend and FARend are nearly identical, bad bracketing
                               //   and root not found

// Given alpha in [-1, 1], generate search point in above interval
//   alpha = -1 gives NEAR end, alpha = 1 gives FAR end of interval
inline DbVector3 searchpoint(DbVector3& pvec1, DbVector3& nvec1, double r,
							 double alpha)
{
	return pvec1 - nvec1*(r*(1.0 + alpha*EPSmedial));
}

/*
bool Pointlist_serverB::subdivMedialStuff(Bpoint *boundpnt, Bpoint *boundpnt_opp, Diatom *medialatom)
{
	// Determine which side of object to check---this test may be made arbitrarily
	//    intelligent, searching equi-(u,v) tiles with opposite t={1, -1} accordingly.
	//    Here, just use a brute force search of the opposite t surface.
	// If on a crest, pick the side with opposite t to the Bpoint, to avoid searching
	//    the side nearest the point and thus avoid possible degeneracies by keeping
	//    as large an object angle between top and bottom as possible.
	// NOTE:  if a crest point chosen, the routine calls itself recursively on
	//    the found involute to generate the correct medial frame.  It checks to
	//    see that the initial Bpoint is in the b-norm plane of the generated atom,
	//    else returns FALSE.]

	// Select side for search, opposite to the side of the point
	searchSIDE whichside;
	if (boundpnt->t == 1.0)
		whichside = BOTTOMSIDE;
	else if (boundpnt->t == -1.0)
		whichside = TOPSIDE;
	else if (boundpnt->t > 0.0)
		whichside = BOTTOMCREST;
	else if (boundpnt->t < 0.0)
		whichside = TOPCREST;
	else // boundpnt->t == 0.0
		whichside = BOTTOMCREST; // chose arbitrarily, for now.

	bool crestvalFLAG = (whichside == BOTTOMCREST || whichside == TOPCREST);

	// Set bisection endpoints for search
	DbVector3 pvec, nvec;

	pvec.set(boundpnt->pnt);
	nvec.set(boundpnt->norm);
	nvec.selfnormalize();

	DbVector3 pvec2, nvec2;
	DbVector3 nvec_intercept, nvec2_intercept;

	// THALLMOD
	double NEAR_end = -2.0;
	double FAR_end = 2.0;
	double rval = boundpnt->rad;
	double medTOLscaled = medTOL*rval;
	double dist_side1, dist_newpnt;
	DbVector3 testmedpnt;

	// Clear cache
	OldnearsetB.clear();

	double testval = 0.5;
    bool exitFLAG = false;
	do {
		testmedpnt = searchpoint(pvec, nvec, rval, testval);
		dist_side1 = (testmedpnt - pvec).length();

		if (!subdivNearSurfaceSIDED(testmedpnt.X(), testmedpnt.Y(), testmedpnt.Z(),
									whichside, boundpnt_opp, &dist_newpnt)) {
			fprintf(stderr, "ERROR in Plistserver2::subdivmedialstuff():\n");
			fprintf(stderr, "   no nearpoint found for this point.\n");
			fprintf(stderr, "Handled heuristically.\n");

			// Compute a new medial point by intersecting normals from opposing (u, v, t)
			//   and (u, v, -t) involutes onto the equidistant plane separating the
			//   the two points, then averaging the points of intersection.  Boundary
			//   normals will not be in phong-normal directions, except at Diatom positions
			subdivBposition(boundpnt_opp, boundpnt->u, boundpnt->v, - boundpnt->t, 0);
			pvec2.set(boundpnt_opp->pnt);
			nvec2.set(boundpnt_opp->norm);
			nvec2.selfnormalize();

			medialatom->averageDiatom(pvec, nvec, pvec2, nvec2);

			return true;
		}

		if (fabs(dist_side1 - dist_newpnt) > medTOLscaled) {
			if (dist_side1 > dist_newpnt) {
				FAR_end = testval;
				testval = (testval + NEAR_end)*0.5;
			}
			else {
				NEAR_end = testval;
				testval = (testval + FAR_end)*0.5;
			}
			if (fabs(NEAR_end - FAR_end) < paramEPS) {

				fprintf(stderr, "ERROR in Plistserver2::subdivmedialstuff():\n");
				fprintf(stderr, "Bad bracketing---root not in [-1, 1]*r*medEPS about r\n");
				fprintf(stderr, "Handled heuristically.\n");

				// Compute a new medial point by intersecting normals from opposing (u, v, t)
				//   and (u, v, -t) involutes onto the equidistant plane separating the
				//   the two points, then averaging the points of intersection.  Boundary
				//   normals will not be in phong-normal directions, except at Diatom positions
				subdivBposition(boundpnt_opp, boundpnt->u, boundpnt->v, - boundpnt->t, 0);
				pvec2.set(boundpnt_opp->pnt);
				nvec2.set(boundpnt_opp->norm);
				nvec2.selfnormalize();

				medialatom->averageDiatom(pvec, nvec, pvec2, nvec2);

				return true;
			}
		}
		else
			exitFLAG = true;

	} while (!exitFLAG);
 /*
	do {
		testmedpnt = searchpoint(pvec, nvec, rval, testval);
		dist_side1 = (testmedpnt - pvec).length();

		if (!subdivNearSurfaceSIDED(testmedpnt.X(), testmedpnt.Y(), testmedpnt.Z(),
			whichside, boundpnt_opp, &dist_newpnt)) {

			NEAR_end = testval;
			testval = (testval + FAR_end)*0.5;

			if (fabs(NEAR_end - FAR_end) < paramEPS) {
				fprintf(stderr, "ERROR in Plistserver2::subdivmedialstuff():\n");
				fprintf(stderr, "   no nearpoint found for this point.\n");
				fprintf(stderr, "Oops!  led to Bad bracketing---root not in [-1, 1]*r*medEPS about r\n");
				return false;
			}
		}
		else if (fabs(dist_side1 - dist_newpnt) > medTOLscaled) {
			if (dist_side1 > dist_newpnt) {
				FAR_end = testval;
				testval = (testval + NEAR_end)*0.5;
			}
			else {
				NEAR_end = testval;
				testval = (testval + FAR_end)*0.5;
			}
			if (fabs(NEAR_end - FAR_end) < paramEPS) {
				fprintf(stderr, "Bad bracketing---root not in [-1, 1]*r*medEPS about r\n");
				return false;
			}
		}
		else
			exitFLAG = true;

	} while (!exitFLAG);
	* /
//	if (!crestvalFLAG) {
		DbVector3 newinvolute(boundpnt_opp->pnt[0], boundpnt_opp->pnt[1], boundpnt_opp->pnt[2]);
		*medialatom = Diatom(testmedpnt, (dist_newpnt + dist_side1)*0.5, pvec - testmedpnt, newinvolute - testmedpnt);
		return true;
/*	}
	else {
		Bpoint newbpt = *boundpnt_opp;
		newbpt.print();
		return subdivMedialStuff(&newbpt, boundpnt_opp, medialatom);
	}
 * /
}
*/

/*
 * subdivMedialStuff2() -- computes an approximate medial axis atom based on an averaging
 *		of medial positions implied by P(u, v, t) and its associated P(u, v, -t), by projecting
 *      their normals to the bisecting plane and averaging their points of intersection
 */
/*
bool Pointlist_serverB::subdivMedialStuff2(Bpoint *boundpnt, Bpoint *boundpnt_opp, Diatom *medialatom)
{
	// Set bisection endpoints for search
	DbVector3 pvec, nvec;
	double radius1, radius2;

	pvec.set(boundpnt->pnt);
	nvec.set(boundpnt->norm);
	nvec.selfnormalize();
	radius1 = boundpnt->rad;

	DbVector3 pvec2, nvec2;
	DbVector3 nvec_intercept, nvec2_intercept;

	// Compute a new medial point by intersecting normals from opposing (u, v, t)
	//   and (u, v, -t) involutes onto the equidistant plane separating the
	//   the two points, then averaging the points of intersection.  Boundary
	//   normals will not be in phong-normal directions, except at Diatom positions
	subdivBposition(boundpnt_opp, boundpnt->u, boundpnt->v, -boundpnt->t, 0);
	pvec2.set(boundpnt_opp->pnt);
	nvec2.set(boundpnt_opp->norm);
	nvec2.selfnormalize();
	radius2 = boundpnt_opp->rad;

	medialatom->averageDiatom2(pvec, nvec, radius1, pvec2, nvec2, radius2);

	return true;
}
*/

/********************************************************************************/
/* subdivboundaryinfo() -- pointer to internal list of the (numpoints)			*/
/*	  Bpoint-structs giving vertices of subdivsurf derived boundary.			*/
/* Return NULL if uninitialized													*/
/********************************************************************************/
void Pointlist_serverB::subdivboundaryinfo(int *numpts, Bpoint **bptlistptr)
{
	*numpts = numsubdivpoints;
	*bptlistptr = subdivbpoints;
}

// Added by Tom, does same as above but returns a copy, not the internal
void Pointlist_serverB::copySubdivBoundaryInfo(int * numpts, Bpoint **bptlistptr)
{
    int i;
    Bpoint * listptr;

    *numpts = numsubdivpoints;
    if(subdivbpoints != NULL)
    {
        *bptlistptr = new Bpoint[numsubdivpoints];  
        listptr = *bptlistptr;
        for(i = 0; i < numsubdivpoints; i++)
        {
            listptr->pnt[0] = subdivbpoints[i].pnt[0];
            listptr->pnt[1] = subdivbpoints[i].pnt[1];
            listptr->pnt[2] = subdivbpoints[i].pnt[2];
            listptr->norm[0] = subdivbpoints[i].norm[0];
            listptr->norm[1] = subdivbpoints[i].norm[1];
            listptr->norm[2] = subdivbpoints[i].norm[2];
            listptr->rad = subdivbpoints[i].rad;
            listptr->u = subdivbpoints[i].u;
            listptr->v = subdivbpoints[i].v;
            listptr->t = subdivbpoints[i].t;

            listptr++;
        }
    }
    else
        *bptlistptr = NULL;
}

/********************************************************************************/
/* subdivtileinfo() -- pointer to internal list of the (numtiles) Bpoint quads  */
/*    generated by a Subdivision surface based on the newgrid.					*/
/*	  Returns NULL if uninitialized.											*/
/********************************************************************************/
void Pointlist_serverB::subdivtileinfo(int *ntiles, Bpoint **bptlistptr)
{
	*ntiles = numsubdivtiles;
	*bptlistptr = subdivtileset;
}

// As above, but returns tile vertices as indices into subdivbpoint list
void Pointlist_serverB::subdivtileindexinfo(int *ntiles, int **vertlistptr)
{
	*ntiles = numsubdivtiles;
	*vertlistptr = subdivtileindices;
}

// Added by Tom, does same as above but returns a copy, not the internal
// Also returns number of points not tiles
void Pointlist_serverB::copySubdivTileInfo(int *ntiles, Bpoint **bptlistptr)
{
    int i;
    Bpoint * listptr;

    *ntiles = numsubdivtiles;
    if(subdivtileset != NULL)
    {
        //*bptlistptr = new Bpoint[*ntiles];
        *bptlistptr = new Bpoint[*ntiles*4]; // changed by CL 01/25/02
        listptr = *bptlistptr;
        for(i = 0; i < *ntiles*4; i++)
        {
            listptr->pnt[0] = subdivtileset[i].pnt[0];
            listptr->pnt[1] = subdivtileset[i].pnt[1];
            listptr->pnt[2] = subdivtileset[i].pnt[2];
            listptr->norm[0] = subdivtileset[i].norm[0];
            listptr->norm[1] = subdivtileset[i].norm[1];
            listptr->norm[2] = subdivtileset[i].norm[2];
            listptr->rad = subdivtileset[i].rad;
            listptr->u = subdivtileset[i].u;
            listptr->v = subdivtileset[i].v;
            listptr->t = subdivtileset[i].t;

            listptr++;
        }
    }
    else
        *bptlistptr = NULL;
}

// returns list of quad-tiles as sets of 4 indices
//    into bptlist of vertices returned by subdivboundaryinfo()
void Pointlist_serverB::copySubdivTileIndexInfo(int *ntiles, int **vertlistptr)
{
    *ntiles = numsubdivtiles;
    if(subdivtileindices != NULL) {

        //*vertlistptr = new int[*ntiles];
        *vertlistptr = new int[*ntiles];  // CL 1/25/02
        for(int i = 0; i < *ntiles*4; i++)
			(*vertlistptr)[i] = subdivtileindices[i];
    }
    else
        *vertlistptr = NULL;
}

/********************************************************************************/
/* subdivtileneighbors() -- returns pointer to internal list of Tileneighbors,	*/
/*	 indexed by tilenumbers from above list.  It is an error to call this		*/
/*	 before calling ComputeSubdivBoundaryTiles().								*/
/*   Not of much use without the tilelist, anyway.								*/
/********************************************************************************/
void Pointlist_serverB::subdivtileneighbors(int *ntiles, Tileneighbors **neighbors)
{
	*ntiles = numsubdivtiles;
	*neighbors = tneighbors;
}


/********************************************************************************/
/* subdivvertexneighbors() -- returns pointer to internal list of Tileneighbors,*/
/*	 indexed by vertexnumbers from above list.  It is an error to call this		*/
/*	 before calling ComputeSubdivPointCloud().									*/
/*   Not of much use without the vertexlist, anyway.							*/
/********************************************************************************/
void Pointlist_serverB::subdivvertexneighbors(int *nverts, Tileneighbors **neighbors)
{
	*nverts = numsubdivpoints;
	*neighbors = vneighbors;
}

/*
/******************************************************************************** /
// Distance squared between two centers
/******************************************************************************** /
double d_sqr(double * ctr_0, double * ctr_1) {
	return (ctr_1[0] - ctr_0[0])*(ctr_1[0] - ctr_0[0])
		+ (ctr_1[1] - ctr_0[1])*(ctr_1[1] - ctr_0[1])
		+ (ctr_1[2] - ctr_0[2])*(ctr_1[2] - ctr_0[2]);
}
*/

/********************************************************************************/
/* glRenderSubdivPoints() -- if subdivbpoints != NULL, render surface point cloud*/
/********************************************************************************/
void Pointlist_serverB::glRenderSubdivPoints()
{
#ifdef COLOR_DEBUG
	glColor3d(1.0, 1.0, 1.0);
#endif
	glBegin(GL_POINTS);
	for (int vcnt1 = 0; vcnt1 < numsubdivpoints; vcnt1++) {

		glVertex3d(subdivbpoints[vcnt1].pnt[0], 
				   subdivbpoints[vcnt1].pnt[1],
				   subdivbpoints[vcnt1].pnt[2]);
	}
	glEnd();
}

/********************************************************************************/
/* glRenderSubdivWires() -- if numsubdivtiles != 0, render boundary wire frame	*/
/********************************************************************************/
void Pointlist_serverB::glRenderSubdivWires()
{
#ifdef COLOR_DEBUG
	glColor3d(0.6, 0.6, 0.6);
#endif
	for (int vcnt = 0; vcnt < 4*numsubdivtiles; vcnt += 4) {

		glBegin(GL_LINE_LOOP);
		glVertex3d(subdivtileset[vcnt].pnt[0],
				   subdivtileset[vcnt].pnt[1],
				   subdivtileset[vcnt].pnt[2]);
		glVertex3d(subdivtileset[vcnt + 1].pnt[0],
				   subdivtileset[vcnt + 1].pnt[1],
				   subdivtileset[vcnt + 1].pnt[2]);
		glVertex3d(subdivtileset[vcnt + 2].pnt[0],
				   subdivtileset[vcnt + 2].pnt[1],
				   subdivtileset[vcnt + 2].pnt[2]);
		glVertex3d(subdivtileset[vcnt + 3].pnt[0],
				   subdivtileset[vcnt + 3].pnt[1],
				   subdivtileset[vcnt + 3].pnt[2]);
		glEnd();
	}
}

/********************************************************************************/
/* glRenderSubdivTiles() -- if numsubdivtiles != 0, render boundary tiles		*/
/********************************************************************************/
void Pointlist_serverB::glRenderSubdivTiles(double displacement)
{
	double u, v, t;
	int vcnt;

	// this gives displacement between 1 and zero.
	displacement = (1.0 - cos(displacement * 0.3/M_PI))/2.0;
	displacement *= 0.9;

	glBegin(GL_QUADS);

//	int numtopbottomtiles = 2*num_ucoltiles*num_vrowtiles;
//	for (vcnt = 0; vcnt < 4*numtopbottomtiles;/*numsubdivtiles;*/ vcnt++) {
	for (vcnt = 0; vcnt < 4*numsubdivtiles; vcnt++) {
		u = subdivtileset[vcnt].u;
		v = subdivtileset[vcnt].v;
		t = subdivtileset[vcnt].t;

#ifdef COLOR_DEBUG
#ifndef NORMALIZE_UVT
		glColor3d(u/max_u(), v/max_v(), 1.0 - (t + 1.0)/2.0);
#else
		glColor3d(u, v, t);
#endif
#endif

        glNormal3d(subdivtileset[vcnt].norm[0], subdivtileset[vcnt].norm[1],
            subdivtileset[vcnt].norm[2]);
        glVertex3d(subdivtileset[vcnt].pnt[0] -
            subdivtileset[vcnt].rad * displacement * subdivtileset[vcnt].norm[0],
            subdivtileset[vcnt].pnt[1] -
                subdivtileset[vcnt].rad * displacement * subdivtileset[vcnt].norm[1],
            subdivtileset[vcnt].pnt[2] -
                subdivtileset[vcnt].rad * displacement * subdivtileset[vcnt].norm[2]);
	}
	glEnd();

}


//////////////////////////////////////////////////////////////////////////
// Added by CL 01/21/02.  Add displacements to boundary point list

void Pointlist_serverB::applyDisplacements(int numPts, double * dispVals)
{
	if (numPts != numsubdivpoints)
	{
		fprintf(stderr, "Error in applying boundary displacements in PS2::applyDisplacements.\n");
		return;
	}

	if (dispVals == NULL)
		return;

	// Add displacements to boundary point list
	double offset;
	int i;

	for (i = 0; i < numPts; i++)
	{
		offset = subdivbpoints[i].rad * dispVals[i];
		subdivbpoints[i].pnt[0] += subdivbpoints[i].norm[0] * offset;
		subdivbpoints[i].pnt[1] += subdivbpoints[i].norm[1] * offset;
		subdivbpoints[i].pnt[2] += subdivbpoints[i].norm[2] * offset;
	}

	if(subdivtileindices==NULL)
	{
		fprintf(stderr, "Error in applying boundary displacements in PS2::applyDisplacements.\n");
		return;
	}

	// Add displacement to tile point list
	for (i = 0; i < numsubdivtiles*4; i++)
	{
		offset = subdivtileset[i].rad * dispVals[ subdivtileindices[i] ];
		subdivtileset[i].pnt[0] += subdivtileset[i].norm[0] * offset;
		subdivtileset[i].pnt[1] += subdivtileset[i].norm[1] * offset;
		subdivtileset[i].pnt[2] += subdivtileset[i].norm[2] * offset;
	}

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Pointlist_serverB::ClearBoundaryTilesPointers()
{
	numsubdivtiles=0;
	if (subdivtileset != NULL)
	{
		delete []subdivtileset;
		subdivtileset=NULL;
	}
	if (subdivtileindices != NULL)
	{
		delete []subdivtileindices;
		subdivtileindices=NULL;
	}
	if (tneighbors != NULL)
	{
		delete []tneighbors;
		tneighbors=NULL;
	}
}

void Pointlist_serverB::ClearPointCloudPointers()
{
	numsubdivpoints=0;
	if(subdivbpoints!=NULL)
	{
		delete []subdivbpoints;
		subdivbpoints=NULL;
	}
	if(vneighbors!=NULL)
	{
		delete []vneighbors;
		vneighbors=NULL;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Pointlist_serverB::ComputeSubdivBoundaryTiles(int subLevel)
{
	level=subLevel-SUBDIV_LEVEL_TO_USE;
	if(!subdivSurfInitialized) 
	{
		if(!InitializeSubdivSurf(level))
			return;
		subdivSurfInitialized=true;
	}

	if(currentBoundaryTilesLevel!=level)
	{
		currentBoundaryTilesLevel=level;
//		ClearBoundaryTilesPointers();
	}
	ClearBoundaryTilesPointers();

	// Deal with NULL pointers instead in thisSurf->SubdivTileInfo, GetTileNeighbors and GetTileIndices

	thisSurf->SubdivTileInfo(currentBoundaryTilesLevel, &numsubdivtiles, &subdivtileset);
	thisSurf->GetTileNeighbors(currentBoundaryTilesLevel, &numsubdivtiles, &tneighbors);
	thisSurf->GetTileIndices(currentBoundaryTilesLevel, &subdivtileindices);
}

void Pointlist_serverB::ComputeSubdivPointCloud(int subLevel)
{
	level=subLevel-SUBDIV_LEVEL_TO_USE;
	if(!subdivSurfInitialized) 
	{
		if(!InitializeSubdivSurf(level))
			return;
		subdivSurfInitialized=true;
	}

	if(currentPointCloudLevel!=level)
	{
		currentPointCloudLevel=level;
//		ClearPointCloudPointers();
	}
	ClearPointCloudPointers();

	// deal with the NULL pointers in thisSurf->SubdivBoundaryInfo and GetVertexNeighbors instead

	int numvertneighbors;
	thisSurf->SubdivBoundaryInfo(currentPointCloudLevel, &numsubdivpoints, &subdivbpoints);
	thisSurf->GetVertexNeighbors(currentPointCloudLevel, &numvertneighbors, &vneighbors);

	if (numsubdivpoints!=numvertneighbors)
		fprintf(stderr, "Pointlistserver2::ComputeSubdivPointCloud:  ERROR, bad number of vertices\n");
}

bool NotFoundIn(int sDex, int eDex, int *dex, int key)
{
	int i;
	for(i=sDex; i<=eDex; i++)
		if(key==dex[i])
			return false;
	return true;
}

// To generate the tag list showing which vertices on the boundary are in the blending region
// This is to replace the distanceMap calls in Mask::Mask
void Pointlist_serverB::ComputeFigureMapToRemove(int *numPtsToRemove, short **vertexMapToRemove, \
													   int subLevel, int fId)
{
	// future work to do: update treeId to its real value!
	int treeId=1;
	MyFloat SUB_FIGURE_CUTOFF_U=SUB_FIGURE_CUTOFF_U_ORIGINAL+SUB_FIGURE_CUTOFF_U_RANGE*blendExtents[treeId];
	MyFloat DILATION_STEPS=DILATION_STEPS_ORIGINAL+DILATION_STEPS_RANGE*blendAmounts[treeId];
	MyFloat EPSILON2=EPSILON2_ORIGINAL;
	int i, j, k;
	int interPNum;
	int interPPair[MAX_INTERSECTION_NUMBER*2];
	int interPPair2[MAX_INTERSECTION_NUMBER*2];
	Bpoint2 interPList[MAX_INTERSECTION_NUMBER];
	Bpoint2 interPList2[MAX_INTERSECTION_NUMBER];
	int *polyDexList1, *polyDexList2;

	int nTiles1, nTiles2;
	Bpoint *tile1, *tile2;
	ThallCode::Pointlist_server2 *pListOld[2];

	interPNum=0;						// number of the vertices on the intersection curve
	pListOld[0]=pListOriginal[0]; //tileLists[0];				// Start with the main figure to fill the tile information
	pListOld[1]=pListOriginal[1]; //tileLists[1];				// Need the sub-figure to calculate the intersection
 	pListOld[0]->subdivtileinfo(&nTiles1, &tile1);
	pListOld[1]->subdivtileinfo(&nTiles2, &tile2);

	polyDexList1=new int[nTiles1];
	polyDexList2=new int[nTiles2];
	memset(polyDexList1, 0, nTiles1*sizeof(int));
	memset(polyDexList2, 0, nTiles2*sizeof(int));

	// to check whether we can actually do the checking
	// if the lvl of boundary is even lower than the lvl we did the intersection w/
	// we can not generate the map!
	if(SUBDIV_LEVEL_TO_USE>subLevel)
	{
		// this is the case that is NOT ALLOWED to happen!
		fprintf(stderr, "SUBDIV_LEVEL_TO_USE>subLevel, CASE NOT ALLOWED IN Pointlist_serverB::\
			ComputeFigureVertexMapToRemove!\n");
		return;
	}

	// do the intersection first, find all the intersection segments
	IntersectMeshes(SUB_FIGURE_CUTOFF_U, nTiles1, tile1, nTiles2, tile2, interPPair, interPList, \
		polyDexList1, polyDexList2, &interPNum);	

	if(interPNum>MAX_INTERSECTION_NUMBER || interPNum<=0)
	{
		fprintf(stderr, "TOO MANY or NONE INTERSECTION POINTS found from the TRI-TRI \
			intersections...\n");

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		return;
	}

	// then find the intersection curve
	bool successful=FindIntersectionCurveFromSegments(interPList2, interPPair2, &interPNum, \
		interPList, interPPair, EPSILON2);
	if(!successful)
	{
		fprintf(stderr, "\nThe intersection curve is not closed!!!\nHave to quit right here...\n");
		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;
		return;
	}

	CCMesh *thisMesh, *thisMesh1, *thisMesh2;
	CCVertex *thisVert; //, *thatVert;

	CCPolygon *thisPoly;
	int pSides;
	int curLevel1, curLevel2;

	if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
		curLevel1=SUBDIV_LEVEL_TO_USE+1;//thisMesh1=(pListOld[0]->GetThisSurfPointer())->SubMeshPointer(SUBDIV_LEVEL_TO_USE+1);
	else
		curLevel1=SUBDIV_LEVEL_TO_USE;//thisMesh1=(pListOld[0]->GetThisSurfPointer())->SubMeshPointer(SUBDIV_LEVEL_TO_USE);
	curLevel2=SUBDIV_LEVEL_TO_USE;

	// future work to do: update fId to the exact id of each figure
	// at this moment it's assumed that the entire object complex contains only a 2-figure object
	if(fId==0)
	{
		// main-figure
		// propagete the intersection quads on the boundary of this level
	thisMesh1=(pListOld[0]->GetThisSurfPointer())->SubMeshPointer(curLevel1);

	thisMesh=thisMesh1;
	int vNums=thisMesh->VertList.Length();
	int fNums=thisMesh->FaceList.Length();

	// before everything, need to find the cut curve on the main figure that is connected to the dilated inter-curve
	// then figure out what vertices in the main figure to get rid of!

	int mainVNumEstimate=vNums*1.5;
	int mainFNumEstimate=fNums*1.5;
	int *quadListMain=new int[fNums], quadListMainNum;
	int *quadListMain_2, quadListMainNum_2;
	int *vertListMainQuad=new int[mainVNumEstimate], *vertListMainQuadCnt=new int[mainVNumEstimate], vertListMainQuadNum;
	int *edgeListMainQuad=new int [mainVNumEstimate], edgeListMainQuadNum=0;

	// quadListMainNum counts the total # of the intersected quads
	quadListMainNum=0;
	// the list of all the quads involved with the intersection with the sub-figure
	// index [0, quadListMainNum-1]
	// range [0, fNums-1]
	//memset(quadListMain, 1, mainFNumEstimate*sizeof(int));

	// vertListMainQuadNum counts the # of vertices of all the intersected quads
	vertListMainQuadNum=0;
	// list of vertices of the intersected main-figured
	// index [0, vertListMainQuadNum-1]
	// range [0, vNums-1]
	memset(vertListMainQuad, -1, mainVNumEstimate*sizeof(int));

	// vertListMainQuadCnt counts the # of appearances of each vertex in the intersected 
	//		quads. if it is the valence of the vertex, it should be removed because it's 
	//		already surrouded by quads that should be removed
	// index [0, vNums-1]
	memset(vertListMainQuadCnt, 0, mainVNumEstimate*sizeof(int));

	pListOld[0]->BuildLookUpTableForPDexByUVT();

	// according to the uvt coordinates of the dilated intersection curve, need to recalculate polyDexList1

	// now add the list of quads intersected by the sub-figure into quadListMain
	for(i=0; i<fNums; i++)
		if(polyDexList1[i])
			quadListMain[quadListMainNum++]=i;

	int *edgeListMainQuadDex=new int[mainVNumEstimate];
	memset(edgeListMainQuadDex, -1, sizeof(int)*mainVNumEstimate);

	int curV, curEHd, curEEd;
	for(i=0; i<quadListMainNum; i++)
	{
		thisPoly=thisMesh->FaceList(quadListMain[i]);
		pSides=thisPoly->numsides();
		for(j=0; j<pSides; j++)
		{
			// try to add curV into the vert list
			curV=thisPoly->Pvertices[j];

			// now need to deal with the edges
			//		it's slightly different cause we remove all the shared edges
			curEHd=curV;
			curEEd=thisPoly->Pvertices[(j+1)%pSides];
			k=0;
			while(k<edgeListMainQuadNum)
			{
				// have to reverse the order of curEEd and curEHd because of the topology
				if(curEEd==edgeListMainQuad[2*k] && curEHd==edgeListMainQuad[2*k+1])
					break;
				k++;
			}

			// if it's not, add it to the list
			//		and update the counter by 1 for the first time
			if(k==edgeListMainQuadNum)
			{
				edgeListMainQuadDex[edgeListMainQuadNum]=quadListMain[i];

				edgeListMainQuad[2*edgeListMainQuadNum]=curEHd;
				edgeListMainQuad[2*edgeListMainQuadNum+1]=curEEd;
				edgeListMainQuadNum++;
			}
			// otherwise need to remove the shared edge from the edgeListMainQuad
			else
			{
				if(k<edgeListMainQuadNum-1)
				{
					while(k<edgeListMainQuadNum-1)
					{
						edgeListMainQuadDex[k]=edgeListMainQuadDex[k+1];;

						edgeListMainQuad[2*k]=edgeListMainQuad[2*(k+1)];
						edgeListMainQuad[2*k+1]=edgeListMainQuad[2*(k+1)+1];
						k++;
					}
				}
				edgeListMainQuadNum--;
			}
		}
	}

	// it is time to find the curve on the main figure to connect to the dilated intersection curve!
	// write the codes to find out the cut-curve on the main figure
	int mainInterPNum;

	int *takenMainQuadEdge=new int[mainFNumEstimate];
	memset(takenMainQuadEdge, 0, sizeof(int)*mainFNumEstimate);
	int mainInterCurveDex[MAX_INTERSECTION_NUMBER];
	//mainInterPNum=edgeListMainQuadNum;

	int mainInterCurveDexInEdgeListMainQuad[MAX_INTERSECTION_NUMBER];
	memset(mainInterCurveDexInEdgeListMainQuad, -1, sizeof(int)*MAX_INTERSECTION_NUMBER);

	bool algorithmFailed=false;
	bool timeToCheckCloseness=false;
	int mCurveNum=0;
	MyPoint2D mCurveIndices[MAX_MAIN_CURVE_NUMBER];
	int mCurveStartDex, mCurveEndDex;
	int successfulMCurveID=-1;

	//memset(taken, 0, sizeof(int)*500);
	mainInterCurveDex[0]=edgeListMainQuad[0];
	mainInterCurveDex[1]=edgeListMainQuad[1];
	mainInterCurveDexInEdgeListMainQuad[1-1]=0;
	mCurveStartDex=0;
	j=0;
	i=1;

	takenMainQuadEdge[0]=1;
	//taken[0]=1;
	// generate the whole intersection curve from the separated intersection segments
	//		all the un-ordereed intersection vertices in interPList now will be linked as a curve into interPList2
	successful=false;
	while(i<edgeListMainQuadNum && !successful)
	{
		while(j<=edgeListMainQuadNum && takenMainQuadEdge[j]) 
			j++;
		if(j>edgeListMainQuadNum)
		{
			// another possibility here is that the curve is already closed! (multiple curves could exist)
			// let us figure out whether it is true and whether the current curve is the outer curve

			if(mainInterCurveDex[i]==mainInterCurveDex[mCurveStartDex])
			{
				// it is already a closed curve
				// need to check whether it is indeed the curve we are looking for
				// only the curve that is COUNTER-CLOCKWISE should be kept!
				if(i-mCurveStartDex>edgeListMainQuadNum/2)
				{
					// if the number of vertices in this curve is over the half of number of all the vertices
					// it is a trivial positive!
					mCurveEndDex=i;
					mCurveIndices[mCurveNum].x=mCurveStartDex;
					mCurveIndices[mCurveNum].y=mCurveEndDex;
					successfulMCurveID=mCurveNum;
					mCurveNum++;
					successful=true;
					continue;
				}
				else
				{
					// otherwise, need to continue the searching
					mCurveEndDex=i;
					mCurveIndices[mCurveNum].x=mCurveStartDex;
					mCurveIndices[mCurveNum].y=mCurveEndDex;
					mCurveNum++;

					i++;
					mCurveStartDex=i;
					j=0;
					while(j<=edgeListMainQuadNum && takenMainQuadEdge[j]) 
						j++;
					if(j<=edgeListMainQuadNum)
					{
						mainInterCurveDexInEdgeListMainQuad[i]=j;

						takenMainQuadEdge[j]=1;
						mainInterCurveDex[i]=edgeListMainQuad[2*j];
						mainInterCurveDex[i+1]=edgeListMainQuad[2*j+1];
						j=0;

						//lastTaken=0;
						i++;
					}
					else
					{
						// the algorithm just failed!
						algorithmFailed=true;
					}
					continue;
				}
			}

			if(algorithmFailed)
			{
				// it is time to quit because this is the illegal case!
				bool trap=true;
				fprintf(stderr, "\nCannot figure out the mainInterCurve! time to quit blending codes\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fclose(fDebug);
#endif

			{
				if(polyDexList1!=NULL)
					delete []polyDexList1;
				if(polyDexList2!=NULL)
					delete []polyDexList2;

				if(quadListMain!=NULL)
					delete []quadListMain;
				if(vertListMainQuad!=NULL)
					delete []vertListMainQuad;
				if(vertListMainQuadCnt!=NULL)
					delete []vertListMainQuadCnt;
				if(edgeListMainQuad!=NULL)
					delete []edgeListMainQuad;

				if(edgeListMainQuadDex!=NULL)
					delete []edgeListMainQuadDex;

				if(takenMainQuadEdge!=NULL)
					delete []takenMainQuadEdge;
			}

				return;
			}
		}

		if(mainInterCurveDex[i]==edgeListMainQuad[2*j])
		{
			mainInterCurveDexInEdgeListMainQuad[i]=j;

			takenMainQuadEdge[j]=1;
			i++;
			mainInterCurveDex[i]=edgeListMainQuad[2*j+1];
			j=0;
			if(i==edgeListMainQuadNum)
				timeToCheckCloseness=true;
		}
		else
			if(mainInterCurveDex[i]==edgeListMainQuad[2*j+1])
			{
				// TO REMOVE: these codes are theoretically NOT-REACHABLE!
				mainInterCurveDexInEdgeListMainQuad[i]=j;

				takenMainQuadEdge[j]=1;
				i++;
				mainInterCurveDex[i]=edgeListMainQuad[2*j];
				j=0;
				if(i==edgeListMainQuadNum)
					timeToCheckCloseness=true;
			}
			else j++;

		if(timeToCheckCloseness)
		{
			if(mainInterCurveDex[i]==mainInterCurveDex[mCurveStartDex])
			{
				// it is already a closed curve
				// need to check whether it is indeed the curve we are looking for
				// only the curve that is COUNTER-CLOCKWISE, or the LONGEST should be kept!
				if(i-mCurveStartDex>edgeListMainQuadNum/2)
				{
					mCurveEndDex=i;
					mCurveIndices[mCurveNum].x=mCurveStartDex;
					mCurveIndices[mCurveNum].y=mCurveEndDex;
					successfulMCurveID=mCurveNum;
					mCurveNum++;
					successful=true;
				}
				else
				{
					mCurveEndDex=i;
					mCurveIndices[mCurveNum].x=mCurveStartDex;
					mCurveIndices[mCurveNum].y=mCurveEndDex;
					mCurveNum++;

					successful=true;
				}
			}
			else
			{
				// time to fail the algorithm again, because the curve is not closed!
				fprintf(stderr, "\nCannot figure out the mainInterCurve! time to quit blending codes\n");
			{
				if(polyDexList1!=NULL)
					delete []polyDexList1;
				if(polyDexList2!=NULL)
					delete []polyDexList2;

				if(quadListMain!=NULL)
					delete []quadListMain;
				if(vertListMainQuad!=NULL)
					delete []vertListMainQuad;
				if(vertListMainQuadCnt!=NULL)
					delete []vertListMainQuadCnt;
				if(edgeListMainQuad!=NULL)
					delete []edgeListMainQuad;

				if(edgeListMainQuadDex!=NULL)
					delete []edgeListMainQuadDex;

				if(takenMainQuadEdge!=NULL)
					delete []takenMainQuadEdge;
			}

				return;
			}
		}
	}
	//mainInterPNum=i;

	// by now
	// mCurveNum			holds the total number of the curves from the intersection quads
	// mCurveIndices		holds the start/end indices of the segments of each curve in mainInterCurveDex
	int maxLength=0;
	if(successfulMCurveID<0)
	//if(!successful)
	{
		// then need to figure out which one of the curves is the intersection curve we need
		for(i=0; i<mCurveNum; i++)
		{
			if(maxLength<mCurveIndices[i].y-mCurveIndices[i].x)
			{
				maxLength=mCurveIndices[i].y-mCurveIndices[i].x;
				successfulMCurveID=i;
			}
		}
	}

	// successfulMCurveID	holds the ID of the curve that is the cut curve of the main figure
	mCurveStartDex=mCurveIndices[successfulMCurveID].x;		// right here right now
	mCurveEndDex=mCurveIndices[successfulMCurveID].y;
	mainInterPNum=mCurveEndDex-mCurveStartDex;

	if(mCurveNum<=0 || mainInterPNum<=0)
	{
		fprintf(stderr, "\nCannot find a valid BOUNDARY CURVE on the main-figure!\nHave to quit blending codes...\n");
		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		if(quadListMain!=NULL)
			delete []quadListMain;
		if(vertListMainQuad!=NULL)
			delete []vertListMainQuad;
		if(vertListMainQuadCnt!=NULL)
			delete []vertListMainQuadCnt;
		if(edgeListMainQuad!=NULL)
			delete []edgeListMainQuad;

		if(edgeListMainQuadDex!=NULL)
			delete []edgeListMainQuadDex;

		if(takenMainQuadEdge!=NULL)
			delete []takenMainQuadEdge;
	}

	// now it is time to cut the main figure by this intersection curve
	// need to grow the cut-curve and tag all the traversed quads... 
	// to quadListMain[quadListMainNum...]
	// before that, need two things
	//		1. the list of the indices of quad containing each segment in edgeListMainQuad	done
	//		2. a Queue stucture																done
	int curEdgeDexInEdgeListMainQuad, curEdgeStartDex, curEdgeEndDex;

	// before everything, tag all the edges of the dilated-i-curve related quads to span
	int *quadEdgeToSpan= new int[mainFNumEstimate*4];
	memset(quadEdgeToSpan, 1, mainFNumEstimate*4);

	j=mCurveStartDex;
	for(i=0; i<mainInterPNum; i++)
	{
		curEdgeDexInEdgeListMainQuad=mainInterCurveDexInEdgeListMainQuad[j];
		curEdgeStartDex=edgeListMainQuad[2*curEdgeDexInEdgeListMainQuad];
		curEdgeEndDex=edgeListMainQuad[2*curEdgeDexInEdgeListMainQuad+1];
		thisPoly=thisMesh->FaceList(edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]);
		pSides=thisPoly->numsides();
		for(k=0; k<pSides; k++)
		{
			// need to only span the quad from the edge that was not contributing to the mainInterCurve
			// in other words, only span the quad inward
			if(thisPoly->Pvertices[k]==curEdgeStartDex && thisPoly->Pvertices[(k+1)%pSides]==curEdgeEndDex)
			{
				quadEdgeToSpan[4*edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]+k]=0;
			}
		}
		j++;
	}

	int *visited=new int[mainFNumEstimate];
	memset(visited, 0, sizeof(int)*mainFNumEstimate);
	MyQueue que(mainVNumEstimate);

	// fill in the queue for the first time
	j=mCurveStartDex;
	for(i=0; i<mainInterPNum; i++)
	{
		curEdgeDexInEdgeListMainQuad=mainInterCurveDexInEdgeListMainQuad[j];
		curEdgeStartDex=edgeListMainQuad[2*curEdgeDexInEdgeListMainQuad];
		curEdgeEndDex=edgeListMainQuad[2*curEdgeDexInEdgeListMainQuad+1];
		thisPoly=thisMesh->FaceList(edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]);

	if(!visited[edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]])
	{
		visited[edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]]=1;

		pSides=thisPoly->numsides();
		for(k=0; k<pSides; k++)
		{
			// need to only span the quad from the edge that was not contributing to the mainInterCurve
			// in other words, only span the quad inward
			if(quadEdgeToSpan[4*edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]+k])
			//if(!(thisPoly->Pvertices[k]==curEdgeStartDex && thisPoly->Pvertices[(k+1)%pSides]==curEdgeEndDex))
			{
				if(!visited[thisPoly->Pneighbors[k]])
				{
					if(!polyDexList1[thisPoly->Pneighbors[k]])
					{
						polyDexList1[thisPoly->Pneighbors[k]]=3;

						quadListMain[quadListMainNum++]=thisPoly->Pneighbors[k];
					}

					que.Enque(thisPoly->Pneighbors[k]);
				}
			}
		}
	}
		j++;
	}

	// by traversing the queue, expand quads related to the dilated intersection curve to all the quads
	// that should be removed from the main figure quads
	int curQuadDex;
	while(!que.Empty())
	{
		curQuadDex=que.Deque();
		if(curQuadDex==-1)
		{
			bool trap=true;
		}
		if(!visited[curQuadDex])
		{
			visited[curQuadDex]=1;
			thisPoly=thisMesh->FaceList(curQuadDex);
			pSides=thisPoly->numsides();
			for(k=0; k<pSides; k++)
			{
				if(!visited[thisPoly->Pneighbors[k]])
				{
					if(!polyDexList1[thisPoly->Pneighbors[k]])
					{
						polyDexList1[thisPoly->Pneighbors[k]]=3;

						quadListMain[quadListMainNum++]=thisPoly->Pneighbors[k];
					}

					que.Enque(thisPoly->Pneighbors[k]);
				}
			}
		}
	}
	delete []visited;
	delete []quadEdgeToSpan;

	// till now all the quads that should be removed at level of 'SUBDIV_LEVEL_TO_USE' have been taged
	// now need to follow the tree structure to find all the quads vertices at subLevel to remove!
	if(curLevel1==subLevel)
	{
		// then we are almost done
	}
	else
		if(curLevel1<subLevel)
		{
			// traverse the parent-children trees to find all the child-quads at the successive lvls
			for(i=curLevel1; i<subLevel; i++)
			{
				quadListMain_2=new int[thisMesh1->numfaces()*4];
				quadListMainNum_2=0;
				for(j=0; j<quadListMainNum; j++)
				{
					thisPoly=thisMesh1->FaceList(quadListMain[j]);
					for(k=0; k<thisPoly->numsides(); k++)
						quadListMain_2[quadListMainNum_2++]=thisPoly->childPolygons[k];
				}
				delete []quadListMain;
				quadListMain=quadListMain_2;
				quadListMainNum=quadListMainNum_2;
				// move to the boundary of the next level
				thisMesh1=(pListOld[0]->GetThisSurfPointer())->SubMeshPointer(i+1);
			}
		}
		else
		{
			// this is the case that is NOT ALLOWED to happen!
			// already checked at the beginning
			return;
		}

		// now, time to tag all the vertices in all the tagged quads
		int vert;
		*vertexMapToRemove=new short[thisMesh1->numverts()];
		memset(*vertexMapToRemove, 0, thisMesh1->numverts()*sizeof(short));
		for(i=0; i<quadListMainNum; i++)
		{
			thisPoly=thisMesh1->FaceList(quadListMain[i]);
			for(j=0; j<thisPoly->numsides(); j++)
			{
				vert=thisPoly->Pvertices[j];
				if(!(*vertexMapToRemove)[vert])//if(!(*vertexMapToRemove)[vert] && NotFoundIn(mCurveStartDex, mCurveEndDex, mainInterCurveDex, vert))
					(*vertexMapToRemove)[vert]=1;
			}
		}

		vert=0;
		for(i=0; i<thisMesh1->numverts(); i++)
		{
			if((*vertexMapToRemove)[i])
				vert++;
		}
		*numPtsToRemove=vert;

				if(quadListMain!=NULL)
					delete []quadListMain;
				if(vertListMainQuad!=NULL)
					delete []vertListMainQuad;
				if(vertListMainQuadCnt!=NULL)
					delete []vertListMainQuadCnt;
				if(edgeListMainQuad!=NULL)
					delete []edgeListMainQuad;
				if(edgeListMainQuadDex!=NULL)
					delete []edgeListMainQuadDex;
				if(takenMainQuadEdge!=NULL)
					delete []takenMainQuadEdge;

	}	// main-figure
	else
	{
		// sub-figure
		thisMesh2=(pListOld[1]->GetThisSurfPointer())->SubMeshPointer(subLevel);
		*vertexMapToRemove=new short[thisMesh2->numverts()];
		memset(*vertexMapToRemove, 0, sizeof(short)*thisMesh2->numverts());

		int vert=0;
		MyFloat closestDis;
		MyFloat dis;
		int interPDex;
		for(i=0; i<thisMesh2->numverts(); i++)
		{
			closestDis=999;
			thisVert=thisMesh2->VertList(i);
			for(j=0; j<interPNum; j++)
			{
				dis=(interPList2[j].v_2-thisVert->uvt.v)*(interPList2[j].v_2-thisVert->uvt.v)+\
					(interPList2[j].t_2-thisVert->uvt.t)*(interPList2[j].t_2-thisVert->uvt.t);
				if(dis<closestDis)
				{
					closestDis=dis;
					interPDex=j;
				}
			}
			if(thisVert->uvt.u<interPList2[j].u_2)
			{
				vert++;
				(*vertexMapToRemove)[i]=1;
			}
		}
		*numPtsToRemove=vert;
	}

	if(polyDexList1!=NULL)
		delete []polyDexList1;
	if(polyDexList2!=NULL)
		delete []polyDexList2;

	return;
}

void Pointlist_serverB::ComputeFiguresMapsToRemove(int *numPtsToRemove, short **vertexMapToRemove, \
													   int subLevel)
{
	// future work to do: update treeId to its real value!
	int treeId=1;
	MyFloat SUB_FIGURE_CUTOFF_U=SUB_FIGURE_CUTOFF_U_ORIGINAL+SUB_FIGURE_CUTOFF_U_RANGE*blendExtents[treeId];
	MyFloat DILATION_STEPS=DILATION_STEPS_ORIGINAL+DILATION_STEPS_RANGE*blendAmounts[treeId];
	MyFloat EPSILON2=EPSILON2_ORIGINAL;
	int i, j, k;
	int interPNum;
	int interPPair[MAX_INTERSECTION_NUMBER*2];
	int interPPair2[MAX_INTERSECTION_NUMBER*2];
	Bpoint2 interPList[MAX_INTERSECTION_NUMBER];
	Bpoint2 interPList2[MAX_INTERSECTION_NUMBER];
	int *polyDexList1, *polyDexList2;

	int nTiles1, nTiles2;
	Bpoint *tile1, *tile2;
	ThallCode::Pointlist_server2 *pListOld[2];

	interPNum=0;						// number of the vertices on the intersection curve
	pListOld[0]=pListOriginal[0]; //tileLists[0];				// Start with the main figure to fill the tile information
	pListOld[1]=pListOriginal[1]; //tileLists[1];				// Need the sub-figure to calculate the intersection
 	pListOld[0]->subdivtileinfo(&nTiles1, &tile1);
	pListOld[1]->subdivtileinfo(&nTiles2, &tile2);

	polyDexList1=new int[nTiles1];
	polyDexList2=new int[nTiles2];
	memset(polyDexList1, 0, nTiles1*sizeof(int));
	memset(polyDexList2, 0, nTiles2*sizeof(int));

	// to check whether we can actually do the checking
	// if the lvl of boundary is even lower than the lvl we did the intersection w/
	// we can not generate the map!
	if(SUBDIV_LEVEL_TO_USE>subLevel)
	{
		// this is the case that is NOT ALLOWED to happen!
		fprintf(stderr, "SUBDIV_LEVEL_TO_USE>subLevel, CASE NOT ALLOWED IN Pointlist_serverB::\
			ComputeFigureVertexMapToRemove!\n");
		return;
	}

	// do the intersection first, find all the intersection segments
	IntersectMeshes(SUB_FIGURE_CUTOFF_U, nTiles1, tile1, nTiles2, tile2, interPPair, interPList, \
		polyDexList1, polyDexList2, &interPNum);	

	if(interPNum>MAX_INTERSECTION_NUMBER || interPNum<=0)
	{
		fprintf(stderr, "TOO MANY or NONE INTERSECTION POINTS found from the TRI-TRI \
			intersections...\n");

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		return;
	}

	// then find the intersection curve
	bool successful=FindIntersectionCurveFromSegments(interPList2, interPPair2, &interPNum, \
		interPList, interPPair, EPSILON2);
	if(!successful)
	{
		fprintf(stderr, "\nThe intersection curve is not closed!!!\nHave to quit right here...\n");
		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;
		return;
	}

	CCMesh *thisMesh, *thisMesh1, *thisMesh2;
	CCVertex *thisVert; //, *thatVert;
	CCPolygon *thisPoly;
	int pSides;
	int curLevel1, curLevel2;

	if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
		curLevel1=SUBDIV_LEVEL_TO_USE+1;//thisMesh1=(pListOld[0]->GetThisSurfPointer())->SubMeshPointer(SUBDIV_LEVEL_TO_USE+1);
	else
		curLevel1=SUBDIV_LEVEL_TO_USE;//thisMesh1=(pListOld[0]->GetThisSurfPointer())->SubMeshPointer(SUBDIV_LEVEL_TO_USE);
	curLevel2=SUBDIV_LEVEL_TO_USE;

	// future work to do: update fId to the exact id of each figure
	// at this moment it's assumed that the entire object complex contains only a 2-figure object
	// main-figure
	{
		// main-figure
		// propagete the intersection quads on the boundary of this level
	thisMesh1=(pListOld[0]->GetThisSurfPointer())->SubMeshPointer(curLevel1);

	thisMesh=thisMesh1;
	int vNums=thisMesh->VertList.Length();
	int fNums=thisMesh->FaceList.Length();

	// before everything, need to find the cut curve on the main figure that is connected to the dilated inter-curve
	// then figure out what vertices in the main figure to get rid of!
	int mainVNumEstimate=vNums*1.5;
	int mainFNumEstimate=fNums*1.5;
	int *quadListMain=new int[fNums], quadListMainNum;
	int *quadListMain_2, quadListMainNum_2;
	int *vertListMainQuad=new int[mainVNumEstimate], *vertListMainQuadCnt=new int[mainVNumEstimate], vertListMainQuadNum;
	int *edgeListMainQuad=new int [mainVNumEstimate], edgeListMainQuadNum=0;


	// flag for each vertex of the main-figure quads, 1 means it should be kept otherwise cut
	// index [0, vNums-1]
	// range [0, 1]

	// quadListMainNum counts the total # of the intersected quads
	quadListMainNum=0;
	// the list of all the quads involved with the intersection with the sub-figure
	// index [0, quadListMainNum-1]
	// range [0, fNums-1]
	//memset(quadListMain, 1, mainFNumEstimate*sizeof(int));

	// vertListMainQuadNum counts the # of vertices of all the intersected quads
	vertListMainQuadNum=0;
	// list of vertices of the intersected main-figured
	// index [0, vertListMainQuadNum-1]
	// range [0, vNums-1]
	memset(vertListMainQuad, -1, mainVNumEstimate*sizeof(int));

	// vertListMainQuadCnt counts the # of appearances of each vertex in the intersected 
	//		quads. if it is the valence of the vertex, it should be removed because it's 
	//		already surrouded by quads that should be removed
	// index [0, vNums-1]
	memset(vertListMainQuadCnt, 0, mainVNumEstimate*sizeof(int));

	pListOld[0]->BuildLookUpTableForPDexByUVT();

	// according to the uvt coordinates of the dilated intersection curve, need to recalculate polyDexList1

	// now add the list of quads intersected by the sub-figure into quadListMain
	for(i=0; i<fNums; i++)
		if(polyDexList1[i])
			quadListMain[quadListMainNum++]=i;

	int *edgeListMainQuadDex=new int[mainVNumEstimate];
	memset(edgeListMainQuadDex, -1, sizeof(int)*mainVNumEstimate);

	int curV, curEHd, curEEd;
	for(i=0; i<quadListMainNum; i++)
	{
		thisPoly=thisMesh->FaceList(quadListMain[i]);
		pSides=thisPoly->numsides();
		for(j=0; j<pSides; j++)
		{
			// try to add curV into the vert list
			curV=thisPoly->Pvertices[j];

			// now need to deal with the edges
			//		it's slightly different cause we remove all the shared edges
			curEHd=curV;
			curEEd=thisPoly->Pvertices[(j+1)%pSides];
			k=0;
			while(k<edgeListMainQuadNum)
			{
				// have to reverse the order of curEEd and curEHd because of the topology
				if(curEEd==edgeListMainQuad[2*k] && curEHd==edgeListMainQuad[2*k+1])
					break;
				k++;
			}

			// if it's not, add it to the list
			//		and update the counter by 1 for the first time
			if(k==edgeListMainQuadNum)
			{
				edgeListMainQuadDex[edgeListMainQuadNum]=quadListMain[i];

				edgeListMainQuad[2*edgeListMainQuadNum]=curEHd;
				edgeListMainQuad[2*edgeListMainQuadNum+1]=curEEd;
				edgeListMainQuadNum++;
			}
			// otherwise need to remove the shared edge from the edgeListMainQuad
			else
			{
				if(k<edgeListMainQuadNum-1)
				{
					while(k<edgeListMainQuadNum-1)
					{
						edgeListMainQuadDex[k]=edgeListMainQuadDex[k+1];;

						edgeListMainQuad[2*k]=edgeListMainQuad[2*(k+1)];
						edgeListMainQuad[2*k+1]=edgeListMainQuad[2*(k+1)+1];
						k++;
					}
				}
				edgeListMainQuadNum--;
			}
		}
	}

	// it is time to find the curve on the main figure to connect to the dilated intersection curve!
	// write the codes to find out the cut-curve on the main figure
	int mainInterPNum;

	int *takenMainQuadEdge=new int[mainFNumEstimate];
	memset(takenMainQuadEdge, 0, sizeof(int)*mainFNumEstimate);
	int mainInterCurveDex[MAX_INTERSECTION_NUMBER];
	//mainInterPNum=edgeListMainQuadNum;

	int mainInterCurveDexInEdgeListMainQuad[MAX_INTERSECTION_NUMBER];
	memset(mainInterCurveDexInEdgeListMainQuad, -1, sizeof(int)*MAX_INTERSECTION_NUMBER);

	bool algorithmFailed=false;
	bool timeToCheckCloseness=false;
	int mCurveNum=0;
	MyPoint2D mCurveIndices[MAX_MAIN_CURVE_NUMBER];
	int mCurveStartDex, mCurveEndDex;
	int successfulMCurveID=-1;

	//memset(taken, 0, sizeof(int)*500);
	mainInterCurveDex[0]=edgeListMainQuad[0];
	mainInterCurveDex[1]=edgeListMainQuad[1];
	mainInterCurveDexInEdgeListMainQuad[1-1]=0;
	mCurveStartDex=0;
	j=0;
	//lastTaken=0;
	i=1;

	takenMainQuadEdge[0]=1;
	//taken[0]=1;
	// generate the whole intersection curve from the separated intersection segments
	//		all the un-ordereed intersection vertices in interPList now will be linked as a curve into interPList2
	successful=false;
	while(i<edgeListMainQuadNum && !successful)
	{
		while(j<=edgeListMainQuadNum && takenMainQuadEdge[j]) 
			j++;
		if(j>edgeListMainQuadNum)
		{
			// another possibility here is that the curve is already closed! (multiple curves could exist)
			// let us figure out whether it is true and whether the current curve is the outer curve

			if(mainInterCurveDex[i]==mainInterCurveDex[mCurveStartDex])
			{
				// it is already a closed curve
				// need to check whether it is indeed the curve we are looking for
				// only the curve that is COUNTER-CLOCKWISE should be kept!
				if(i-mCurveStartDex>edgeListMainQuadNum/2)
				{
					// if the number of vertices in this curve is over the half of number of all the vertices
					// it is a trivial positive!
					mCurveEndDex=i;
					mCurveIndices[mCurveNum].x=mCurveStartDex;
					mCurveIndices[mCurveNum].y=mCurveEndDex;
					successfulMCurveID=mCurveNum;
					mCurveNum++;
					successful=true;
					continue;
				}
				else
				{
					// otherwise, need to continue the searching
					mCurveEndDex=i;
					mCurveIndices[mCurveNum].x=mCurveStartDex;
					mCurveIndices[mCurveNum].y=mCurveEndDex;
					mCurveNum++;

					i++;
					mCurveStartDex=i;
					j=0;
					while(j<=edgeListMainQuadNum && takenMainQuadEdge[j]) 
						j++;
					if(j<=edgeListMainQuadNum)
					{
						mainInterCurveDexInEdgeListMainQuad[i]=j;

						takenMainQuadEdge[j]=1;
						mainInterCurveDex[i]=edgeListMainQuad[2*j];
						mainInterCurveDex[i+1]=edgeListMainQuad[2*j+1];
						j=0;

						//lastTaken=0;
						i++;
					}
					else
					{
						// the algorithm just failed!
						algorithmFailed=true;
					}
					continue;
				}
			}

			if(algorithmFailed)
			{
				// it is time to quit because this is the illegal case!
				bool trap=true;
				fprintf(stderr, "\nCannot figure out the mainInterCurve! time to quit blending codes\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fclose(fDebug);
#endif
			{
				if(polyDexList1!=NULL)
					delete []polyDexList1;
				if(polyDexList2!=NULL)
					delete []polyDexList2;


				if(quadListMain!=NULL)
					delete []quadListMain;
				if(vertListMainQuad!=NULL)
					delete []vertListMainQuad;
				if(vertListMainQuadCnt!=NULL)
					delete []vertListMainQuadCnt;
				if(edgeListMainQuad!=NULL)
					delete []edgeListMainQuad;

				if(edgeListMainQuadDex!=NULL)
					delete []edgeListMainQuadDex;

				if(takenMainQuadEdge!=NULL)
					delete []takenMainQuadEdge;
			}

				return;
			}
		}

		if(mainInterCurveDex[i]==edgeListMainQuad[2*j])
		{
			mainInterCurveDexInEdgeListMainQuad[i]=j;

			takenMainQuadEdge[j]=1;
			i++;
			mainInterCurveDex[i]=edgeListMainQuad[2*j+1];
			j=0;
			if(i==edgeListMainQuadNum)
				timeToCheckCloseness=true;
		}
		else
			if(mainInterCurveDex[i]==edgeListMainQuad[2*j+1])
			{
				// TO REMOVE: these codes are theoretically NOT-REACHABLE!
				mainInterCurveDexInEdgeListMainQuad[i]=j;

				takenMainQuadEdge[j]=1;
				i++;
				mainInterCurveDex[i]=edgeListMainQuad[2*j];
				j=0;
				if(i==edgeListMainQuadNum)
					timeToCheckCloseness=true;
			}
			else j++;

		if(timeToCheckCloseness)
		{
			if(mainInterCurveDex[i]==mainInterCurveDex[mCurveStartDex])
			{
				// it is already a closed curve
				// need to check whether it is indeed the curve we are looking for
				// only the curve that is COUNTER-CLOCKWISE, or the LONGEST should be kept!
				if(i-mCurveStartDex>edgeListMainQuadNum/2)
				{
					mCurveEndDex=i;
					mCurveIndices[mCurveNum].x=mCurveStartDex;
					mCurveIndices[mCurveNum].y=mCurveEndDex;
					successfulMCurveID=mCurveNum;
					mCurveNum++;
					successful=true;
				}
				else
				{
					mCurveEndDex=i;
					mCurveIndices[mCurveNum].x=mCurveStartDex;
					mCurveIndices[mCurveNum].y=mCurveEndDex;
					mCurveNum++;

					successful=true;
					//successful=true;
				}
			}
			else
			{
				// time to fail the algorithm again, because the curve is not closed!
				fprintf(stderr, "\nCannot figure out the mainInterCurve! time to quit blending codes\n");
			{
				if(polyDexList1!=NULL)
					delete []polyDexList1;
				if(polyDexList2!=NULL)
					delete []polyDexList2;

				if(quadListMain!=NULL)
					delete []quadListMain;
				if(vertListMainQuad!=NULL)
					delete []vertListMainQuad;
				if(vertListMainQuadCnt!=NULL)
					delete []vertListMainQuadCnt;
				if(edgeListMainQuad!=NULL)
					delete []edgeListMainQuad;

				if(edgeListMainQuadDex!=NULL)
					delete []edgeListMainQuadDex;

				if(takenMainQuadEdge!=NULL)
					delete []takenMainQuadEdge;
			}

				return;
			}
		}
	}
	//mainInterPNum=i;

	// by now
	// mCurveNum			holds the total number of the curves from the intersection quads
	// mCurveIndices		holds the start/end indices of the segments of each curve in mainInterCurveDex
	int maxLength=0;
	if(successfulMCurveID<0)
	//if(!successful)
	{
		// then need to figure out which one of the curves is the intersection curve we need
		for(i=0; i<mCurveNum; i++)
		{
			if(maxLength<mCurveIndices[i].y-mCurveIndices[i].x)
			{
				maxLength=mCurveIndices[i].y-mCurveIndices[i].x;
				successfulMCurveID=i;
			}
		}
	}

	// successfulMCurveID	holds the ID of the curve that is the cut curve of the main figure
	mCurveStartDex=mCurveIndices[successfulMCurveID].x;		// right here right now
	mCurveEndDex=mCurveIndices[successfulMCurveID].y;
	mainInterPNum=mCurveEndDex-mCurveStartDex;

	if(mCurveNum<=0 || mainInterPNum<=0)
	{
		fprintf(stderr, "\nCannot find a valid BOUNDARY CURVE on the main-figure!\nHave to quit blending codes...\n");

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		if(quadListMain!=NULL)
			delete []quadListMain;
		if(vertListMainQuad!=NULL)
			delete []vertListMainQuad;
		if(vertListMainQuadCnt!=NULL)
			delete []vertListMainQuadCnt;
		if(edgeListMainQuad!=NULL)
			delete []edgeListMainQuad;

		if(edgeListMainQuadDex!=NULL)
			delete []edgeListMainQuadDex;

		if(takenMainQuadEdge!=NULL)
			delete []takenMainQuadEdge;
	}

	// now it is time to cut the main figure by this intersection curve
	// need to grow the cut-curve and tag all the traversed quads... 
	// to quadListMain[quadListMainNum...]
	// before that, need two things
	//		1. the list of the indices of quad containing each segment in edgeListMainQuad	done
	//		2. a Queue stucture																done
	int curEdgeDexInEdgeListMainQuad, curEdgeStartDex, curEdgeEndDex;

	// before everything, tag all the edges of the dilated-i-curve related quads to span
	int *quadEdgeToSpan= new int[mainFNumEstimate*4];
	memset(quadEdgeToSpan, 1, mainFNumEstimate*4);

	j=mCurveStartDex;
	for(i=0; i<mainInterPNum; i++)
	{
		curEdgeDexInEdgeListMainQuad=mainInterCurveDexInEdgeListMainQuad[j];
		curEdgeStartDex=edgeListMainQuad[2*curEdgeDexInEdgeListMainQuad];
		curEdgeEndDex=edgeListMainQuad[2*curEdgeDexInEdgeListMainQuad+1];
		thisPoly=thisMesh->FaceList(edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]);
		pSides=thisPoly->numsides();
		for(k=0; k<pSides; k++)
		{
			// need to only span the quad from the edge that was not contributing to the mainInterCurve
			// in other words, only span the quad inward
			if(thisPoly->Pvertices[k]==curEdgeStartDex && thisPoly->Pvertices[(k+1)%pSides]==curEdgeEndDex)
			{
				quadEdgeToSpan[4*edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]+k]=0;
			}
		}
		j++;
	}

	int *visited=new int[mainFNumEstimate];
	memset(visited, 0, sizeof(int)*mainFNumEstimate);
	MyQueue que(mainVNumEstimate);

	// fill in the queue for the first time
	j=mCurveStartDex;
	for(i=0; i<mainInterPNum; i++)
	{
		curEdgeDexInEdgeListMainQuad=mainInterCurveDexInEdgeListMainQuad[j];
		curEdgeStartDex=edgeListMainQuad[2*curEdgeDexInEdgeListMainQuad];
		curEdgeEndDex=edgeListMainQuad[2*curEdgeDexInEdgeListMainQuad+1];
		thisPoly=thisMesh->FaceList(edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]);

	if(!visited[edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]])
	{
		visited[edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]]=1;

		pSides=thisPoly->numsides();
		for(k=0; k<pSides; k++)
		{
			// need to only span the quad from the edge that was not contributing to the mainInterCurve
			// in other words, only span the quad inward
			if(quadEdgeToSpan[4*edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]+k])
			//if(!(thisPoly->Pvertices[k]==curEdgeStartDex && thisPoly->Pvertices[(k+1)%pSides]==curEdgeEndDex))
			{
				if(!visited[thisPoly->Pneighbors[k]])
				{
					if(!polyDexList1[thisPoly->Pneighbors[k]])
					{
						polyDexList1[thisPoly->Pneighbors[k]]=3;

						quadListMain[quadListMainNum++]=thisPoly->Pneighbors[k];
					}

					que.Enque(thisPoly->Pneighbors[k]);
				}
			}
		}
	}
		j++;
	}

	// by traversing the queue, expand quads related to the dilated intersection curve to all the quads
	// that should be removed from the main figure quads
	int curQuadDex;
	while(!que.Empty())
	{
		curQuadDex=que.Deque();
		if(curQuadDex==-1)
		{
			bool trap=true;
		}
		if(!visited[curQuadDex])
		{
			visited[curQuadDex]=1;
			thisPoly=thisMesh->FaceList(curQuadDex);
			pSides=thisPoly->numsides();
			for(k=0; k<pSides; k++)
			{
				if(!visited[thisPoly->Pneighbors[k]])
				{
					if(!polyDexList1[thisPoly->Pneighbors[k]])
					{
						polyDexList1[thisPoly->Pneighbors[k]]=3;

						quadListMain[quadListMainNum++]=thisPoly->Pneighbors[k];
					}

					que.Enque(thisPoly->Pneighbors[k]);
				}
			}
		}
	}
	delete []visited;
	delete []quadEdgeToSpan;

	// till now all the quads that should be removed at level of 'SUBDIV_LEVEL_TO_USE' have been taged
	// now need to follow the tree structure to find all the quads vertices at subLevel to remove!
	if(curLevel1==subLevel)
	{
		// then we are almost done
	}
	else
		if(curLevel1<subLevel)
		{
			// traverse the parent-children trees to find all the child-quads at the successive lvls
			for(i=curLevel1; i<subLevel; i++)
			{
				quadListMain_2=new int[thisMesh1->numfaces()*4];
				quadListMainNum_2=0;
				for(j=0; j<quadListMainNum; j++)
				{
					thisPoly=thisMesh1->FaceList(quadListMain[j]);
					for(k=0; k<thisPoly->numsides(); k++)
						quadListMain_2[quadListMainNum_2++]=thisPoly->childPolygons[k];
				}
				delete []quadListMain;
				quadListMain=quadListMain_2;
				quadListMainNum=quadListMainNum_2;
				// move to the boundary of the next level
				thisMesh1=(pListOld[0]->GetThisSurfPointer())->SubMeshPointer(i+1);
			}
		}
		else
		{
			// this is the case that is NOT ALLOWED to happen!
			// already checked at the beginning
			return;
		}

		// now, time to tag all the vertices in all the tagged quads
		int vert;
		vertexMapToRemove[0]=new short[thisMesh1->numverts()];
		memset(vertexMapToRemove[0], 0, thisMesh1->numverts()*sizeof(short));
		for(i=0; i<quadListMainNum; i++)
		{
			thisPoly=thisMesh1->FaceList(quadListMain[i]);
			for(j=0; j<thisPoly->numsides(); j++)
			{
				vert=thisPoly->Pvertices[j];
				if(!(vertexMapToRemove[0][vert]))//if(!(*vertexMapToRemove)[vert] && NotFoundIn(mCurveStartDex, mCurveEndDex, mainInterCurveDex, vert))
					vertexMapToRemove[0][vert]=1;
			}
		}

		vert=0;
		for(i=0; i<thisMesh1->numverts(); i++)
		{
			if(vertexMapToRemove[0][i])
				vert++;
		}
		numPtsToRemove[0]=vert;

				if(quadListMain!=NULL)
					delete []quadListMain;
				if(vertListMainQuad!=NULL)
					delete []vertListMainQuad;
				if(vertListMainQuadCnt!=NULL)
					delete []vertListMainQuadCnt;
				if(edgeListMainQuad!=NULL)
					delete []edgeListMainQuad;
				if(edgeListMainQuadDex!=NULL)
					delete []edgeListMainQuadDex;
				if(takenMainQuadEdge!=NULL)
					delete []takenMainQuadEdge;
	}	// main-figure is done

	{	// sub-figure
		thisMesh2=(pListOld[1]->GetThisSurfPointer())->SubMeshPointer(subLevel);
		vertexMapToRemove[1]=new short[thisMesh2->numverts()];
		memset(vertexMapToRemove[1], 0, sizeof(short)*thisMesh2->numverts());

		int vert=0;
		MyFloat closestDis;
		MyFloat dis;
		int interPDex;
		for(i=0; i<thisMesh2->numverts(); i++)
		{
			closestDis=999;
			thisVert=thisMesh2->VertList(i);
			for(j=0; j<interPNum; j++)
			{
				dis=(interPList2[j].v_2-thisVert->uvt.v)*(interPList2[j].v_2-thisVert->uvt.v)+\
					(interPList2[j].t_2-thisVert->uvt.t)*(interPList2[j].t_2-thisVert->uvt.t);
				if(dis<closestDis)
				{
					closestDis=dis;
					interPDex=j;
				}
			}
			if(thisVert->uvt.u<interPList2[j].u_2)
			{
				vert++;
				vertexMapToRemove[1][i]=1;
			}
		}
		numPtsToRemove[1]=vert;
	}	// sub-figure is done

	if(polyDexList1!=NULL)
		delete []polyDexList1;
	if(polyDexList2!=NULL)
		delete []polyDexList2;

	return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool Pointlist_serverB::InitializeSubdivSurf(int subLevel)
{
	// object		-> pointer to the entire object complex
	// root			-> root to the current object (single- or multi-figured)
	// subobject	->
	// node			->
	// currentFigureId used to index the current object figure (either in single- or multi-figure
	// object) in the whole complex
	level=subLevel-1;
	subdivSurfInitialized=true;

	return ComputeBlending();
//	return true;
}

bool Pointlist_serverB::InitializeSubdivSurf(M3DObject *obj, int subLevel)
{
	if(subdivSurfInitialized) 
		return true;

	init(obj, subLevel);

	//if(!InitializeSubdivSurf(obj, subLevel)) //if(!InitializeSubdivSurf())
	if(!InitializeSubdivSurf(subLevel)) //if(!InitializeSubdivSurf())
		return false;
	else
		subdivSurfInitialized=true;

	return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Pointlist_serverB::BuildLookUpTableForPDexByUVT()
{
	int i;
	for(i=0; i<figureNum; i++)
		pListOriginal[i]->BuildLookUpTableForPDexByUVT();
}

int Pointlist_serverB::LookUpPDexByUVT(int figureId, Bpoint2 pt)
{
	if(ValidFigureID(figureId))
		return pListOriginal[figureId]->LookUpPDexByUVT(pt);
	else
		return -1;
}

void Pointlist_serverB::glRender(int vNum, int fNum, Bpoint2 *vList, PindexList *fList, RenderStyle renderMode, MyFloat *c)
{
	thisSurf->glRender(vNum, fNum, vList, fList, renderMode, c);
}

void Pointlist_serverB::glRenderPoint(Bpoint2 pnt, MyFloat *c)
{
/*
#ifdef DEBUG_QIONG
	thisSurf->glRender(0, WIREFRAME);
#endif
*/
	thisSurf->glRenderPoint(pnt, c);
}

void Pointlist_serverB::glRenderBlendedObject(int subLvl, RenderStyle renderMode)
{
#ifdef DEBUG_QIONG_SUBDIVISION_02032003
	thisSurf->glRender(-1, WIREFRAME);
#endif
	thisSurf->glRender(subLvl, renderMode);
}

void Pointlist_serverB::glRender(int subLvl, RenderStyle renderMode)
{
#ifdef DEBUG_QIONG_SUBDIVISION_02032003
	thisSurf->glRender(-1, WIREFRAME);
#endif

#ifdef QIONG_BLENDING_NO_DEBUG
	thisSurf->glRenderBlendedObject(subLvl, renderMode);
#endif
}

void Pointlist_serverB::glRenderSegments(int startPDex, int endPDex, Bpoint2 *pntList, MyFloat *c)
{
/*
#ifdef DEBUG_QIONG
	thisSurf->glRender(0, WIREFRAME);
#endif
*/
	thisSurf->glRenderSegments(startPDex, endPDex, pntList, c);
}

void Pointlist_serverB::glRenderCurve(int startPDex, int endPDex, Bpoint2 *pntList, MyFloat *c)
{
/*
#ifdef DEBUG_QIONG
	thisSurf->glRender(0, WIREFRAME);
#endif
*/
	thisSurf->glRenderCurve(startPDex, endPDex, pntList, c);
}

void Pointlist_serverB::glRender(int subLvl, int polyDex, RenderStyle renderMode, MyFloat *c)
{
#ifdef DEBUG_QIONG
	thisSurf->glRender(0, WIREFRAME);
#endif
	thisSurf->glRender(subLvl, polyDex, renderMode, c);
}

void Pointlist_serverB::glRenderBlended(int subLvl, RenderStyle renderMode)
{
	thisSurf->glRender(-1, WIREFRAME);
	thisSurf->glRenderBlended(subLvl, renderMode);
}

//bool Pointlist_serverB::InitializeBlendedSubdivSurf(int vNum, int fNum, Bpoint *vList, PindexList *fList)
bool Pointlist_serverB::InitializeBlendedSubdivSurf(int vNum, int fNum, Bpoint2 *vList, PindexList *fList, int *vTag, int *fTag)
{
#ifndef QIONG_BLENDING_PURE_SMOOTHING
	int i=0;
	CCVertex *thisVertex;
	CCMesh *thisMesh;

	thisSurf=new CCSubdivsurf(NUMLEVELS);

	thisSurf->initBlending(vNum, fNum, vList, fList);
	thisSurf->loadlevelzeroBlending();
	thisSurf->splitandaverageBlending(0, vTag, fTag);

	// after the subdivision in the parameter space, need to find the corresponding vertex-locations 
	//		in the euclidean space!
	//*
	int vertdex;
	thisMesh=thisSurf->submeshes(1);
	Bpoint bp;

	for(vertdex=0; vertdex<thisMesh->numverts(); vertdex++)
	{
		while(thisSurf->vRedex[vertdex]==-1)
			vertdex++;

		thisVertex=thisMesh->VertList(vertdex);

#ifndef DISPLACEMENT_FACTOR
#define DISPLACEMENT_FACTOR 0
#endif
		pListOriginal[thisSurf->vFlag[vertdex]]->subdivBPositionControlMesh(1, &bp, thisVertex->uvt.u, thisVertex->uvt.v, thisVertex->uvt.t);
		thisVertex->loc.set(bp.pnt[0]+DISPLACEMENT_FACTOR*bp.norm[0], bp.pnt[1]+DISPLACEMENT_FACTOR*bp.norm[1], bp.pnt[2]+DISPLACEMENT_FACTOR*bp.norm[2]);
		thisVertex->norm.set(bp.norm[0], bp.norm[1], bp.norm[2]);
		thisVertex->rad=bp.rad;
	}


	//*/
	// Lastly, compute the limit mesh
	thisSurf->computelimitmeshBlending(1);

	thisSurf->splitandaverageBlending_1(1);
	thisSurf->splitandaverageBlending_2(2);
	//thisSurf->splitandaverageBlending_1(2);
	return true;

#else

	thisSurf=new CCSubdivsurf(NUMLEVELS);

	thisSurf->initBlending(vNum, fNum, vList, fList);
	thisSurf->loadlevelzeroPureSmoothing();

	for(int slevel = 0; slevel < MIN_SUBDIV_LEVEL_TO_INTERPOLATE_BP; slevel++) 
		thisSurf->splitandaverage(slevel);
	return true;
#endif
}

bool Pointlist_serverB::InitializeBlendedSubdivSurf_1(int subLevel, int lvlToCalculateIntersection, \
	int vNum, int fNum, Bpoint2 *vList, PindexList *fList, int *vTag, int *fTag)
{
#ifndef QIONG_BLENDING_PURE_SMOOTHING
	int i=0;
	CCVertex *thisVertex;
	CCMesh *thisMesh;

	thisSurf=new CCSubdivsurf(NUMLEVELS);

	thisSurf->initBlending(vNum, fNum, vList, fList);
	thisSurf->loadlevelzeroBlending();

	//*
	int vertdex;
	thisMesh=thisSurf->submeshes(0);
	Bpoint bp;
	for(vertdex=0; vertdex<thisMesh->numverts(); vertdex++)
	{
		thisVertex=thisMesh->VertList(vertdex);

#ifndef DISPLACEMENT_FACTOR
#define DISPLACEMENT_FACTOR 0
#endif
		if(vTag[vertdex]>1)
			pListOriginal[0]->subdivBPositionControlMesh(1, &bp, thisVertex->uvt.u, thisVertex->uvt.v, thisVertex->uvt.t);
		else
			pListOriginal[1]->subdivBPositionControlMesh(1, &bp, thisVertex->uvt.u, thisVertex->uvt.v, thisVertex->uvt.t);

		thisVertex->loc.set(bp.pnt[0]+DISPLACEMENT_FACTOR*bp.norm[0], bp.pnt[1]+DISPLACEMENT_FACTOR*bp.norm[1], \
			bp.pnt[2]+DISPLACEMENT_FACTOR*bp.norm[2]);

		if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
		{
			if(vTag[vertdex]<=1)
			{
				thisVertex->norm.set(-bp.norm[0], -bp.norm[1], -bp.norm[2]);
				thisVertex->limitnorm.set(-bp.norm[0], -bp.norm[1], -bp.norm[2]);
			}
			else
			{
				thisVertex->norm.set(bp.norm[0], bp.norm[1], bp.norm[2]);
				thisVertex->limitnorm.set(bp.norm[0], bp.norm[1], bp.norm[2]);
			}
		}
		else
		{
			thisVertex->norm.set(bp.norm[0], bp.norm[1], bp.norm[2]);
			thisVertex->limitnorm.set(bp.norm[0], bp.norm[1], bp.norm[2]);
		}

		thisVertex->rad=bp.rad;
	}
	//*/

	// Lastly, compute the limit mesh
	// instead of calculating all the subdivision levels until NUMLEVELS-1, only subdivide it 
	// once and display the surface
	level=MIN_SUBDIV_LEVEL_TO_INTERPOLATE_BP>subLevel?MIN_SUBDIV_LEVEL_TO_INTERPOLATE_BP:subLevel;

//	printf("level=%d\n", level);

	for(int slevel = 0; slevel < level; slevel++) //for(int slevel = 0; slevel < 1; slevel++) //for(int slevel = 0; slevel < NUMLEVELS - 1; slevel++)
	{
		thisSurf->splitandaverage(slevel);

		for(int vDex=0; vDex<thisSurf->SubMeshPointer(slevel+1)->numverts(); vDex++)
			thisSurf->SubMeshPointer(slevel+1)->VertList(vDex)->norm.selfnormalize();
	}
	subdivSurfInitialized=true;

	return true;

#else

	thisSurf=new CCSubdivsurf(NUMLEVELS);

	//*
	thisSurf->initBlending(vNum, fNum, vList, fList);
	thisSurf->loadlevelzeroPureSmoothing();
	//*/

	for(int slevel=0; slevel < MIN_SUBDIV_LEVEL_TO_INTERPOLATE_BP; slevel++) // for(int slevel = 0; slevel < NUMLEVELS - 1; slevel++)
		thisSurf->splitandaverage(slevel);
	subdivSurfInitialized=true;
	return true;
#endif
}

bool Pointlist_serverB::InitializeBlendedSubdivSurf_New(int subLevel, int lvlToCalculateIntersection, int vNum, int fNum, \
	Bpoint2 *vList, PindexList *fList, int *vTag, int *fTag)
{
	switch(lvlToCalculateIntersection)
	{
	case 0:
	case 1:
		return InitializeBlendedSubdivSurf_1(subLevel, lvlToCalculateIntersection, vNum, fNum, vList, fList, vTag, fTag);
		break;
	default:
		return 0;
		break;
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/********************************************************************************/
/* subdivBposition() -- given (u, v, t), return Bpoint on the tiled	subdivision */
/*    generated boundary corresponding to it.  u, v, and t are double in ranges	*/
/*	  [0, max_u()],[0, max_v()], [-1, 1].										*/
/*	If successfully computed point (if legal address), return true, else false	*/
/********************************************************************************/
//bool Pointlist_serverB::subdivBPositionNormal(Bpoint *BPosition, Bpoint *BNormal, double u, double v, double t)
#define MAX 99999
bool Pointlist_serverB::WithinPolygonUVT(CCMesh *thisMesh, CCPolygon *thisPoly, double u, double v, double t)
{
	int dex;
	CCVertex *thisVert;
	double uMin=MAX, uMax=-MAX, vMin=MAX, vMax=-MAX, tMin=MAX, tMax=-MAX;

	for(dex=0; dex<thisPoly->numsides(); dex++)
	{
		thisVert=thisMesh->VertList(thisPoly->Pvertices[dex]);
		if(uMin>thisVert->uvt.u) uMin=thisVert->uvt.u;
		if(uMax<thisVert->uvt.u) uMax=thisVert->uvt.u;
		if(vMin>thisVert->uvt.v) vMin=thisVert->uvt.v;
		if(vMax<thisVert->uvt.v) vMax=thisVert->uvt.v;
		if(tMin>thisVert->uvt.t) tMin=thisVert->uvt.t;
		if(tMax<thisVert->uvt.t) tMax=thisVert->uvt.t;
	}
	if(u>=uMin && u<=uMax && v>=vMin && v<=vMax && t>=tMin && t<=tMax)
		return true;
	else
		return false;
}

bool Pointlist_serverB::WithinPolygonUVT(int subLvl, CCPolygon *thisPoly, double u, double v, double t)
{
	int dex;
	CCMesh *thisMesh=thisSurf->SubMeshPointer(subLvl);
	CCVertex *thisVert;
	double uMin=MAX, uMax=-MAX, vMin=MAX, vMax=-MAX, tMin=MAX, tMax=-MAX;

	for(dex=0; dex<thisPoly->numsides(); dex++)
	{
		thisVert=thisMesh->VertList(thisPoly->Pvertices[dex]);
		if(uMin>thisVert->uvt.u) uMin=thisVert->uvt.u;
		if(uMax<thisVert->uvt.u) uMax=thisVert->uvt.u;
		if(vMin>thisVert->uvt.v) vMin=thisVert->uvt.v;
		if(vMax<thisVert->uvt.v) vMax=thisVert->uvt.v;
		if(tMin>thisVert->uvt.t) tMin=thisVert->uvt.t;
		if(tMax<thisVert->uvt.t) tMax=thisVert->uvt.t;
	}
	if(u>=uMin && u<=uMax && v>=vMin && v<=vMax && t>=tMin && t<=tMax)
		return true;
	else
		return false;
}

int Pointlist_serverB::WithinPolyDexAtSublvlUVT(int pDex, int subLvl, double u, double v, double t)
{
	if(subLvl<1)
		return pDex;
	else
	{
		CCMesh *lvlMesh0;
		//CCVertex *quadVert[4];
		lvlMesh0=thisSurf->SubMeshPointer(0);
		CCPolygon *thisPoly=lvlMesh0->FaceList(pDex);
		return WithinPolygonAtSublvlUVT(thisPoly, subLvl, u, v, t);
	}
}

int Pointlist_serverB::WithinPolygonAtSublvlUVT(CCPolygon *thisPoly, int subLvl, double u, double v, double t)
{
	bool found;
	int lvl, dex, polyDex=-1;
	CCPolygon *childPoly;

	CCMesh *lvlMesh[NUMLEVELS];
	//CCVertex *quadVert[4];
	for(dex=1; dex<4; dex++)
	//for(dex=1; dex<NUMLEVELS; dex++)
		lvlMesh[dex]=thisSurf->SubMeshPointer(dex);

	for(lvl=1; lvl<=subLvl; lvl++)
	{
		found=false;
		for(dex=0; dex<thisPoly->numsides(); dex++)
		{
			childPoly=lvlMesh[lvl]->FaceList(thisPoly->childPolygons[dex]);
			if(WithinPolygonUVT(lvl, childPoly, u, v, t))
			{
				found=true;
				break;
			}
		}
		if(found)
		{
			polyDex=thisPoly->childPolygons[dex];
			thisPoly=childPoly;
		}
	}
	if(found)
	{
		return polyDex;
	}
	else
		return -1;
}

bool Pointlist_serverB::subdivBposition(Bpoint *Bpnt, double u, double v, double t, double w)
{
	if(subdivSurfInitialized) 
		subdivBPosition(Bpnt, u, v, t, w);
	else
		return false;

	return true;
	//printf("%s\n", "\nPointlist_serverB::subdivBposition shoule all be replaced by subdivBPosition!");
	//return false;
}

// need a new parameter for the lvl needed
bool Pointlist_serverB::subdivBPosition(Bpoint *Bpnt, double u, double v, double t, double w)
{
	// you can't call this function if the subdivision surfaces have not been calculated yet!
	if(!subdivSurfInitialized) 
		return false;
	// if (u, v, t) coords out of bounds, return false
//	if (u < 0.0 || u > max_u() || v < 0.0 || v > max_v() || t < -1.0 || t > 1.0)
//		return false;

	CCMesh *lvlMesh[NUMLEVELS];
	CCPolygon *thisPoly, *childPoly;
	//CCVertex *quadVert[4];

	int dex;
	int subLvl;
	bool found;

	// if (u, v, t) coords out of bounds, return false
//	if (u < 0.0 || u > max_u() || v < 0.0 || v > max_v() || t < -1.0 || t > 1.0)
//		return false;

	for(dex=0; dex<=MIN_SUBDIV_LEVEL_TO_INTERPOLATE_BP; dex++) //for(dex=0; dex<NUMLEVELS; dex++)
		lvlMesh[dex]=thisSurf->SubMeshPointer(dex);

	// from level 0, find the root polygon for the given (u, v, t)
	dex=0;
	found=false;
	while(dex<lvlMesh[0]->numfaces() && !found)
	{
		thisPoly=lvlMesh[0]->FaceList(dex);
		if(WithinPolygonUVT(lvlMesh[0], thisPoly, u, v, t))
		{
			found=true;
		}
		dex++;
	}

	// from now on, we check the child polygons of the current polygon
	//		after the loop, childPoly is the pointer to the found polygon
	//		containing the current given paramters of (u, v, t)
	for(subLvl=1; subLvl<=MIN_SUBDIV_LEVEL_TO_INTERPOLATE_BP; subLvl++) //for(subLvl=1; subLvl<NUMLEVELS; subLvl++)
	{
		//CCPolygon *thispoly = oldlevel->FaceList(polydex);
		int nsides=thisPoly->numsides();
		dex=0;
		found=false;
		while(dex<nsides && !found)
		{
			childPoly=lvlMesh[subLvl]->FaceList(thisPoly->childPolygons[dex]);
			if(WithinPolygonUVT(lvlMesh[subLvl], childPoly, u, v, t))
			{
				found=true;
			}
			dex++;
		}
		if(found)
			thisPoly=childPoly;
	}


	// now it's time to actually do the interpolation
	subLvl=MIN_SUBDIV_LEVEL_TO_INTERPOLATE_BP; //subLvl=NUMLEVELS-1;


	CCVertex BPnt;
	if(!InterpolateQuad(&BPnt, u, v, t, 0, thisPoly, subLvl))
		return false;
	BPnt.limitnorm.selfnormalize();
	Bpnt->pnt[0]=BPnt.limitloc.x();
	Bpnt->pnt[1]=BPnt.limitloc.y();
	Bpnt->pnt[2]=BPnt.limitloc.z();
	Bpnt->norm[0]=BPnt.limitnorm.x();
	Bpnt->norm[1]=BPnt.limitnorm.y();
	Bpnt->norm[2]=BPnt.limitnorm.z();
	Bpnt->rad=BPnt.rad;

	return true;
}

bool Pointlist_serverB::InterpolateQuad(CCVertex *vert, double u, double v, double t, double w, CCPolygon *poly, int subLvl)
{
	// you can't call this function if the subdivision surfaces have not been calculated yet!
	if(!subdivSurfInitialized) 
		return false;

	// if (u, v, t) coords out of bounds, return false
//	if (u < 0.0 || u > max_u() || v < 0.0 || v > max_v() || t < -1.0 || t > 1.0)
//		return false;

	int startPos, dex;
	CCMesh *mesh=thisSurf->SubMeshPointer(subLvl);
	CCVertex *quadVert[4];

	CCVertex *thisVert, *thatVert;
	thisVert=mesh->VertList(poly->Pvertices[0]);
	thatVert=mesh->VertList(poly->Pvertices[1]);

	// here, need to figure out how to interpolate different part of the patches
	// because 
	//		t is 1 or -1
	//			interpolate by (u, v)
	//		else
	//			interpolate by (v, t) or (u, t)
	int caseNum;
	startPos=0;
	if(t==1.0 || t==-1.0)
	{
		if(thisVert->uvt.u!=thatVert->uvt.u) 
		{
			startPos=1;
		}

		caseNum=TOP_BOTTOM;
	}
	else
	{
		if(mesh->VertList(poly->Pvertices[0])->uvt.u==mesh->VertList(poly->Pvertices[1])->uvt.u &&
		   mesh->VertList(poly->Pvertices[1])->uvt.u==mesh->VertList(poly->Pvertices[2])->uvt.u &&
		   mesh->VertList(poly->Pvertices[2])->uvt.u==mesh->VertList(poly->Pvertices[3])->uvt.u   )
		{
			if(thisVert->uvt.v!=thatVert->uvt.v)
			{
				startPos=1;
			}
			caseNum=CREST_V;
		}
		else
		{
			if(thisVert->uvt.u!=thatVert->uvt.u)
			{
				startPos=1;
			}
			caseNum=CREST_U;
		}
	}

	// I want to schedule the quadVert as follows:
	for(dex=0; dex<4; dex++)
	{
		quadVert[dex]=mesh->VertList(poly->Pvertices[startPos++%4]);
	}

	/*
			2 - - - 1
			|		|
			|		|  vv
			|		|  |
			3 - - - 0  |
				uu---
	*/

	// re-calculate the parameters uu/vv to make the interpolation much more uniform
	double uu, vv;
	switch(caseNum)
	{
		case TOP_BOTTOM:
			uu=(u-quadVert[0]->uvt.u)/(quadVert[3]->uvt.u-quadVert[0]->uvt.u);
			vv=(v-quadVert[0]->uvt.v)/(quadVert[1]->uvt.v-quadVert[0]->uvt.v);
			break;
		case CREST_U:
			uu=(u-quadVert[0]->uvt.u)/(quadVert[3]->uvt.u-quadVert[0]->uvt.u);
			vv=(t-quadVert[0]->uvt.t)/(quadVert[1]->uvt.t-quadVert[0]->uvt.t);
			break;
		case CREST_V:
			uu=(v-quadVert[0]->uvt.v)/(quadVert[3]->uvt.v-quadVert[0]->uvt.v);
			vv=(t-quadVert[0]->uvt.t)/(quadVert[1]->uvt.t-quadVert[0]->uvt.t);
			break;
	}
	uu=uu<0?-uu:uu;
	vv=vv<0?-vv:vv;

	// the bi-linear interpolation looks like
	/*
			  thatVert

			2 - - - - - 1
			|			|
			|			|
			|	vert	|  vv
			|			|  |
			|			|  |
			3 - - - - - 0  |
					uu---

			  thisVert
	*/


	//*
	CCVertex thisVert2;
	CCVertex thatVert2;
	InterpolateVertex(&thisVert2, uu, quadVert[0], quadVert[3]);
	InterpolateVertex(&thatVert2, uu, quadVert[1], quadVert[2]);
	InterpolateVertex(vert, vv, &thisVert2, &thatVert2);
	//*/

	return true;
}

bool Pointlist_serverB::InterpolateVertex(CCVertex *vert, double p, CCVertex *vert1, CCVertex *vert2)
{
	vert->limitloc=(1.0-p)*vert1->limitloc+p*vert2->limitloc;
	vert->limitnorm=(1.0-p)*vert1->limitnorm+p*vert2->limitnorm;
	vert->rad=(1.0-p)*vert1->rad+p*vert2->rad;
	return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// blending codes are here

// Now the correspondence is defiend as
//	dilation along the sub-figure	blendExtents
//	dilation along the main-figure	blendAmounts

using namespace std;

// ----------------------------  Start of non-class functions  ----------------------------



//		pointlist server is also called here to generate the figure surface first
//		qiong caution: these parameters should somehow be determined dynamically
//					   if i have more time i will play with these

bool Pointlist_serverB::SamePointByDistance(MyFloat epsilon, Bpoint2 p1, Bpoint2 p2)
{
	MyFloat mag=0;
	bool res;
	mag+=(p1.pnt[0]-p2.pnt[0])*(p1.pnt[0]-p2.pnt[0]);
	mag+=(p1.pnt[1]-p2.pnt[1])*(p1.pnt[1]-p2.pnt[1]);
	mag+=(p1.pnt[2]-p2.pnt[2])*(p1.pnt[2]-p2.pnt[2]);
	res=sqrt(mag)<epsilon;
	return res;
}
bool Pointlist_serverB::SamePointByDistanceSquare(MyFloat epsilon, Bpoint2 p1, Bpoint2 p2)
{
	MyFloat mag=0;
	bool res;
	mag+=(p1.pnt[0]-p2.pnt[0])*(p1.pnt[0]-p2.pnt[0]);
	mag+=(p1.pnt[1]-p2.pnt[1])*(p1.pnt[1]-p2.pnt[1]);
	mag+=(p1.pnt[2]-p2.pnt[2])*(p1.pnt[2]-p2.pnt[2]);
	res=mag<epsilon;
	return res;
}
bool Pointlist_serverB::SamePoint(MyFloat epsilon, Bpoint2 p1, Bpoint2 p2)
{
	bool res;
	res=p1.pnt[0]-p2.pnt[0]>=-epsilon&&p1.pnt[0]-p2.pnt[0]<=epsilon&&
		p1.pnt[1]-p2.pnt[1]>=-epsilon&&p1.pnt[1]-p2.pnt[1]<=epsilon&&
		p1.pnt[2]-p2.pnt[2]>=-epsilon&&p1.pnt[2]-p2.pnt[2]<=epsilon;
	return res;
}
bool Pointlist_serverB::SamePoint(MyFloat epsilon, Bpoint p1, Bpoint p2)
{
	bool res;
	res=p1.pnt[0]-p2.pnt[0]>=-epsilon&&p1.pnt[0]-p2.pnt[0]<=epsilon&&
		p1.pnt[1]-p2.pnt[1]>=-epsilon&&p1.pnt[1]-p2.pnt[1]<=epsilon&&
		p1.pnt[2]-p2.pnt[2]>=-epsilon&&p1.pnt[2]-p2.pnt[2]<=epsilon;
	return res;
}

#define NUMERICAL_ZERO 	1.0e-010
bool NumericallyEqual(MyFloat v1, MyFloat v2)
{
	if(fabs(v1-v2)<NUMERICAL_ZERO) 
		return true;
	else
		return false;
}

// the UVTR/UVTR_2 should be carried as well as the positions
void Pointlist_serverB::RegularizeDilatedIntersectionCurve(Bpoint2 *newC, Bpoint2 *oldC, int *pNum)
{
	int i, j, k;
	int interPNum=*pNum;
	double dist;

	*pNum=0;
	for(k=0; k<3; k++)
	{
		//newC[0].pnt[k]=oldC[0].pnt[k];
		newC[0]=oldC[0];
	}

	for(i=0; i<interPNum; i++)
	{
		// (i, j) will be an active segment
		dist=0;
		j=i+1;

		// dist SHOULD be RE-CALCULATED between newC[*pNum] and oldC[j]
		// NOT between oldC[i] and oldC[j]!
		for(k=0; k<3; k++)
//			dist+=(oldC[i].pnt[k]-oldC[j].pnt[k])*(oldC[i].pnt[k]-oldC[j].pnt[k]);
			dist+=(newC[*pNum].pnt[k]-oldC[j].pnt[k])*(newC[*pNum].pnt[k]-oldC[j].pnt[k]);

		// only keep the segment longer than the sqrt(MINIMAL_SEGMENT_LENGTH_SQUARE)
		if(dist>=MINIMAL_SEGMENT_LENGTH_SQUARE)
		{
			/*
				might need to be done:	later on I might also use the angle between adjacent segments 
										to eliminate colinear sub-segments
				if(...)
				{
				}

				the real implementation could be found below, outside this loop!
			*/

			(*pNum)++;
			//for(k=0; k<3; k++)
			{
				//newC[*pNum].pnt[k]=oldC[j].pnt[k];
				newC[*pNum]=oldC[j];
			}
		}
		else
			if(i==interPNum-1)
			{
				//for(k=0; k<3; k++)
				{
					//newC[*pNum].pnt[k]=oldC[j].pnt[k];
					newC[*pNum]=oldC[j];
				}
			}
			else
				for(k=0; k<2; k++)
				{
//					newPairs[2*(*pNum)+k]=oldPairs[2*j+k];
					//newPairs[2*(*pNum)+1]=oldPairs[2*i+1];
				}
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_REGULARIZED
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fprintf(fDebug, "distance between dilated inter-point %d and %d is %f\n", i, j, dist);
	#else
				fprintf(stderr, "distance between dilated inter-point %d and %d is %f\n", i, j, dist);
	#endif
#endif
	}
}
// the UVTR/UVTR_2 should be carried as well as the positions
void Pointlist_serverB::RegularizeDilatedIntersectionCurve(int numPoints, Bpoint2 *newC, Bpoint2 *oldC, int *pNum)
{
	int i, j, k;
	int interPNum=*pNum;
	double dist;

	// Calculate the total length of the intersection 'curve'
	double totalDist=0;
	for(i=0; i<interPNum; i++)
	{
		dist=0;
		for(k=0; k<3; k++)
			dist+=(oldC[i+1].pnt[k]-oldC[i].pnt[k])*(oldC[i+1].pnt[k]-oldC[i].pnt[k]);
		totalDist+=sqrt(dist);
	}

	// Determine the minimal segment legnth square by the (totalDist/numPoints)^2
	double minLenSquare=totalDist/numPoints;
	minLenSquare*=minLenSquare;

	*pNum=0;
	for(k=0; k<3; k++)
	{
		//newC[0].pnt[k]=oldC[0].pnt[k];
		newC[0]=oldC[0];
	}

	for(i=0; i<interPNum; i++)
	{
		// (i, j) will be an active segment
		dist=0;
		j=i+1;

		// dist SHOULD be RE-CALCULATED between newC[*pNum] and oldC[j]
		// NOT between oldC[i] and oldC[j]!
		for(k=0; k<3; k++)
//			dist+=(oldC[i].pnt[k]-oldC[j].pnt[k])*(oldC[i].pnt[k]-oldC[j].pnt[k]);
			dist+=(newC[*pNum].pnt[k]-oldC[j].pnt[k])*(newC[*pNum].pnt[k]-oldC[j].pnt[k]);

		// only keep the segment longer than the sqrt(MINIMAL_SEGMENT_LENGTH_SQUARE)
		if(dist>=minLenSquare)	//if(dist>=MINIMAL_SEGMENT_LENGTH_SQUARE)
		{
			/*
				might need to be done:	later on I might also use the angle between adjacent segments 
										to eliminate colinear sub-segments
				if(...)
				{
				}

				the real implementation could be found below, outside this loop!
			*/

			(*pNum)++;
			//for(k=0; k<3; k++)
			{
				//newC[*pNum].pnt[k]=oldC[j].pnt[k];
				newC[*pNum]=oldC[j];
			}
		}
		else
			if(i==interPNum-1)
			{
				//for(k=0; k<3; k++)
				{
					//newC[*pNum].pnt[k]=oldC[j].pnt[k];
					newC[*pNum]=oldC[j];
				}
			}
			else
				for(k=0; k<2; k++)
				{
//					newPairs[2*(*pNum)+k]=oldPairs[2*j+k];
					//newPairs[2*(*pNum)+1]=oldPairs[2*i+1];
				}
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_REGULARIZED
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fprintf(fDebug, "distance between dilated inter-point %d and %d is %f\n", i, j, dist);
	#else
				fprintf(stderr, "distance between dilated inter-point %d and %d is %f\n", i, j, dist);
	#endif
#endif
	}
}


// to deal with the invalid vt value of the intersection vertex, which should be removed
bool ValidVTCoordinatesInSubFigure(MyFloat vMaxS, Bpoint2 bp)
{
	if(NumericallyEqual(bp.t_2, 1.0) || NumericallyEqual(bp.t_2, -1.0))
		return true;
	else
	{
		if(NumericallyEqual(bp.v_2, 0.0) || NumericallyEqual(bp.v_2, vMaxS))
			return true;
		else
			return false;
	}
}

void Pointlist_serverB::RegularizeIntersectionCurve(Bpoint2 *newC, Bpoint2 *oldC, int *newPairs, int *oldPairs, int *pNum)
{
	int i, j, k;
	int interPNum=*pNum;
	double dist;

	*pNum=0;
	for(k=0; k<3; k++)
	{
		//newC[0].pnt[k]=oldC[0].pnt[k];
		newC[0]=oldC[0];
	}
	for(k=0; k<2; k++)
		newPairs[k]=oldPairs[k];

	for(i=0; i<interPNum; i++)
	{
		// (i, j) will be an active segment
		dist=0;

		j=i+1;
		// dist SHOULD be RE-CALCULATED between newC[*pNum] and oldC[j]
		// NOT between oldC[i] and oldC[j]!
		for(k=0; k<3; k++)
//			dist+=(oldC[i].pnt[k]-oldC[j].pnt[k])*(oldC[i].pnt[k]-oldC[j].pnt[k]);
			dist+=(newC[*pNum].pnt[k]-oldC[j].pnt[k])*(newC[*pNum].pnt[k]-oldC[j].pnt[k]);

		// only keep the segment longer than the sqrt(MINIMAL_SEGMENT_LENGTH_SQUARE)
		if(dist>=MINIMAL_SEGMENT_LENGTH_SQUARE)
		{
			(*pNum)++;
			//for(k=0; k<3; k++)
			{
				//newC[*pNum].pnt[k]=oldC[j].pnt[k];
				newC[*pNum]=oldC[j];
			}
			for(k=0; k<2; k++)
				newPairs[2*(*pNum)+k]=oldPairs[2*j+k];

		}
		else
			if(i==interPNum-1)
			{
				//for(k=0; k<3; k++)
				{
					//newC[*pNum].pnt[k]=oldC[j].pnt[k];
					newC[*pNum]=oldC[j];
				}
				for(k=0; k<2; k++)
					newPairs[2*(*pNum)+k]=oldPairs[2*j+k];

			}
			else
				for(k=0; k<2; k++)
				{
					newPairs[2*(*pNum)+k]=oldPairs[2*j+k];
					//newPairs[2*(*pNum)+1]=oldPairs[2*i+1];
				}
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_REGULARIZED
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
				//fprintf(fDebug, "distance between inter-point %d and %d is %f\n", i, j, dist);
	#else
				//fprintf(stderr, "distance between inter-point %d and %d is %f\n", i, j, dist);
	#endif
#endif
	}

}
// to deal with the invalid vt value of the intersection vertex, which should be removed
void Pointlist_serverB::RegularizeIntersectionCurve(MyFloat vMaxS, Bpoint2 *newC, Bpoint2 *oldC, int *newPairs, int *oldPairs, int *pNum)
{
	int i, j, k;
	int interPNum=*pNum;
	double dist;

	*pNum=0;
	for(k=0; k<3; k++)
	{
		//newC[0].pnt[k]=oldC[0].pnt[k];
		newC[0]=oldC[0];
	}
	for(k=0; k<2; k++)
		newPairs[k]=oldPairs[k];

	for(i=0; i<interPNum; i++)
	{
		// (i, j) will be an active segment
		dist=0;
		j=i+1;
		for(k=0; k<3; k++)
			dist+=(oldC[i].pnt[k]-oldC[j].pnt[k])*(oldC[i].pnt[k]-oldC[j].pnt[k]);

		// only keep the segment longer than the sqrt(MINIMAL_SEGMENT_LENGTH_SQUARE)
		if(dist>=MINIMAL_SEGMENT_LENGTH_SQUARE && ValidVTCoordinatesInSubFigure(vMaxS, oldC[j]))
		{
			(*pNum)++;
			//for(k=0; k<3; k++)
			{
				//newC[*pNum].pnt[k]=oldC[j].pnt[k];
				newC[*pNum]=oldC[j];
			}
			for(k=0; k<2; k++)
				newPairs[2*(*pNum)+k]=oldPairs[2*j+k];

		}
		else
			if(i==interPNum-1)
			{
				//for(k=0; k<3; k++)
				{
					//newC[*pNum].pnt[k]=oldC[j].pnt[k];
					newC[*pNum]=oldC[j];
				}
				for(k=0; k<2; k++)
					newPairs[2*(*pNum)+k]=oldPairs[2*j+k];

			}
			else
				for(k=0; k<2; k++)
				{
					newPairs[2*(*pNum)+k]=oldPairs[2*j+k];
					//newPairs[2*(*pNum)+1]=oldPairs[2*i+1];
				}
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_REGULARIZED
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		//fprintf(fDebug, "distance between inter-point %d and %d is %f\n", i, j, dist);
	#else
		//fprintf(stderr, "distance between inter-point %d and %d is %f\n", i, j, dist);
	#endif
#endif
	}
}
// the UVTR/UVTR_2 should be carried as well as the positions
void Pointlist_serverB::RegularizeIntersectionCurve(int numPoints, Bpoint2 *newC, Bpoint2 *oldC, int *newPairs, int *oldPairs, int *pNum)
{
	int i, j, k;
	int interPNum=*pNum;
	double dist;

	// Calculate the total length of the intersection 'curve'
	double totalDist=0;
	for(i=0; i<interPNum; i++)
	{
		dist=0;
		for(k=0; k<3; k++)
			dist+=(oldC[i+1].pnt[k]-oldC[i].pnt[k])*(oldC[i+1].pnt[k]-oldC[i].pnt[k]);
		totalDist+=sqrt(dist);
	}

	// Determine the minimal segment legnth square by the (totalDist/numPoints)^2
	double minLenSquare=totalDist/numPoints;
	minLenSquare*=minLenSquare;

	*pNum=0;
	for(k=0; k<3; k++)
	{
		//newC[0].pnt[k]=oldC[0].pnt[k];
		newC[0]=oldC[0];
	}
	for(k=0; k<2; k++)
		newPairs[k]=oldPairs[k];

	for(i=0; i<interPNum; i++)
	{
		// (i, j) will be an active segment
		dist=0;

		j=i+1;
		// dist SHOULD be RE-CALCULATED between newC[*pNum] and oldC[j]
		// NOT between oldC[i] and oldC[j]!
		for(k=0; k<3; k++)
//			dist+=(oldC[i].pnt[k]-oldC[j].pnt[k])*(oldC[i].pnt[k]-oldC[j].pnt[k]);
			dist+=(newC[*pNum].pnt[k]-oldC[j].pnt[k])*(newC[*pNum].pnt[k]-oldC[j].pnt[k]);

		// only keep the segment longer than the sqrt(MINIMAL_SEGMENT_LENGTH_SQUARE)
		if(dist>=minLenSquare)	//if(dist>=MINIMAL_SEGMENT_LENGTH_SQUARE)
		{
			(*pNum)++;
			//for(k=0; k<3; k++)
			{
				//newC[*pNum].pnt[k]=oldC[j].pnt[k];
				newC[*pNum]=oldC[j];
			}
			for(k=0; k<2; k++)
				newPairs[2*(*pNum)+k]=oldPairs[2*j+k];

		}
		else
			if(i==interPNum-1)
			{
				//for(k=0; k<3; k++)
				{
					//newC[*pNum].pnt[k]=oldC[j].pnt[k];
					newC[*pNum]=oldC[j];
				}
				for(k=0; k<2; k++)
					newPairs[2*(*pNum)+k]=oldPairs[2*j+k];

			}
			else
				for(k=0; k<2; k++)
				{
					newPairs[2*(*pNum)+k]=oldPairs[2*j+k];
					//newPairs[2*(*pNum)+1]=oldPairs[2*i+1];
				}
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_REGULARIZED
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
				//fprintf(fDebug, "distance between inter-point %d and %d is %f\n", i, j, dist);
	#else
				//fprintf(stderr, "distance between inter-point %d and %d is %f\n", i, j, dist);
	#endif
#endif
	}
}

void Pointlist_serverB::RegularizeIntersectionCurveUVT(double uMaxM, double vMaxM, double uMaxS, double vMaxS, Bpoint2 *curv, int pNum)
{
	int i;
	for(i=0; i<pNum+1; i++)
	{
		if(NumericallyEqual(curv[i].t, 1.0))
		{
			// top/bottom vertex
			curv[i].t=1.0;
		}
		else 
			if(NumericallyEqual(curv[i].t, -1.0))
			{
				// top/bottom vertex
				curv[i].t=-1.0;
			}
			else
			{
				// crest vertex
				if(NumericallyEqual(curv[i].u, 0))
					curv[i].u=0;
				if(NumericallyEqual(curv[i].u, uMaxM))
					curv[i].u=uMaxM;
				if(NumericallyEqual(curv[i].v, 0))
					curv[i].v=0;
				if(NumericallyEqual(curv[i].v, vMaxM))
					curv[i].v=vMaxM;
			}

		if(NumericallyEqual(curv[i].t_2, 1.0))
		{
			// top/bottom vertex
			curv[i].t_2=1.0;
		}
		else 
			if(NumericallyEqual(curv[i].t_2, -1.0))
			{
				// top/bottom vertex
				curv[i].t_2=-1.0;
			}
			else
			{
				// crest vertex
				if(NumericallyEqual(curv[i].u_2, 0))
					curv[i].u_2=0;
				if(NumericallyEqual(curv[i].u_2, uMaxS))
					curv[i].u_2=uMaxS;
				if(NumericallyEqual(curv[i].v_2, 0))
					curv[i].v_2=0;
				if(NumericallyEqual(curv[i].v_2, vMaxS))
					curv[i].v_2=vMaxS;
			}
	}
}

void Pointlist_serverB::ConvertCoordinates(MyPoint2D *c2D, Bpoint2 *curv, MyFloat maxU, MyFloat maxV, int pNum)
{
	int i=0, j=0;
	//int pos;
	int startPos=-1;

	// -2 for not-defined
	// 1  for top
	// 0  for crest
	// -1 for bottom
	int type=0;

	//
	//	   | 1 |
	//	-----------
	//	   |   |
	//	 4 | 0 | 3
	//	-----------
	//	   | 2 |
	//
	int area=0;
	int subCase=0;
	int oppositeSide=0;

	while(i<pNum)
	{

		//if(NumericallyEqual(curv[i].t, 1.0) || NumericallyEqual(curv[i].t, -1.0))
		if(curv[i].t==1||curv[i].t==-1)
		{
			//if(NumericallyEqual(curv[i].t, 1.0))
			if(curv[i].t==1.0)
			{
				type=1;
				oppositeSide=-1;
				startPos=i;
			}
			else
			{
				type=-1;
				oppositeSide=1;
				startPos=i;
			}
			break;
		}
		i++;
	}

	// if none of the vertices is on the top/bottom
	// the whole curve is on one of the crest!
	//if(i==pNum)
	{
		//type=0;
		// then need to figure which crest the curve is on

	}

	// now we only deal with the case of starting from TOP/BOTTOM
	// except the case that the whole intersection cruve falls in the crest region!
	if(startPos==-1) startPos=0;

	// calculate all the points from 0 to pNum-1
	j=startPos%pNum;

if(type==1)
{
	for(i=0; i<pNum; i++)
	{
		// this isn't right cause '+2' does not move the point accordingly in the uvt 2d space
		// instead, I need to calculate the 'alpha' value 
		//if(NumericallyEqual(curv[j].t, -oppositeSide))
		if(curv[j].t==-oppositeSide || curv[j].t==oppositeSide)
		{
			// top
			if(curv[j].t==-oppositeSide)
			{
				c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
				c2D[j].y=curv[j].v+2*FACTOR_T_2_V;
				area=0;
			}
			// bottom
			else
				// when to use subCase?
				switch(area)
				{
				case 1:
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
					c2D[j].y=(maxV-curv[j].v)+maxV+4*FACTOR_T_2_V;
					break;
				case 2:
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
					c2D[j].y=-curv[j].v;
					break;
				case 3:
					c2D[j].y=curv[j].v+2*FACTOR_T_2_V;
					c2D[j].x=(maxU-curv[j].u)+maxU+4*FACTOR_T_2_U;
					break;
				case 4:
					c2D[j].y=curv[j].v+2*FACTOR_T_2_V;
					c2D[j].x=-curv[j].u;
					break;
				default:
					fprintf(stderr, "ConvertCorordinates failed!\n");
					bool trap=true;
					break;
				}
		}
		else
		{
	//
	//	   | 1 |
	//	-----------
	//	   |   |
	//	 4 | 0 | 3
	//	-----------
	//	   | 2 |
	//
			//if(area==0)
			//{
				//if(!NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t!=oppositeSide)
				//{
				// now to determine the whole scheme of mapping
				//		in other words, need to figure out which area to map the crest/bottom points
				if(subCase!=area) subCase=area;
				if(curv[j].v==maxV) area=1;
				else
					if(curv[j].v==0) area=2;
					else
						if(curv[j].u==maxU) area=3;
						else
							if(curv[j].u==0) area=4;
							else
							{
								fprintf(stderr, "ConvertCorordinates failed!\n");
								bool trap=true;
							}
				//}
				/*
				else
				{
					// the current vertex falls into the other side of the model, top(bottom) -> bottom(top)
					switch(subCase)
					{
					case 1:
						area=subCase;
						break;
					case 2:
						area=subCase;
						break;
					case 3:
						area=subCase;
						break;
					case 4:
						area=subCase;
						break;
					default:
						fprintf(stderr, "ConvertCorordinates failed!\n");
						bool trap=true;
						break;
					}
				}
				*/
			//}
			switch(area)
			{
			case 1:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t==oppositeSide)
				{
				}
				//else
				{
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
					c2D[j].y=(3-curv[j].t)*FACTOR_T_2_V+maxV;
				}
				break;
			case 2:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t==oppositeSide)
				{
				}
				//else
				{
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
					c2D[j].y=(1+curv[j].t)*FACTOR_T_2_V;
				}
				break;
			case 3:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t==oppositeSide)
				{
				}
				//else
				{
					c2D[j].y=curv[j].v+2*FACTOR_T_2_V;
					c2D[j].x=(3-curv[j].t)*FACTOR_T_2_U+maxU;
				}
				break;
			case 4:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t==oppositeSide)
				{
				}
				//else
				{
					c2D[j].y=curv[j].v+2*FACTOR_T_2_V;
					c2D[j].x=(1+curv[j].t)*FACTOR_T_2_U;
				}
				break;
			default:
				fprintf(stderr, "ConvertCorordinates failed!\n");
				break;
				//return;
			}
		}

		// time for the next vertex on the curve
		if(j==pNum-1)
			j=0;
		else
			j++;
	}
}
else
if(type==-1)
{
	for(i=0; i<pNum; i++)
	{
		// this isn't right cause '+2' does not move the point accordingly in the uvt 2d space
		//if(NumericallyEqual(curv[j].t, -oppositeSide))
		if(curv[j].t==-oppositeSide || curv[j].t==oppositeSide)
		{
			// top
			if(curv[j].t==-oppositeSide)
			{
				c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
				c2D[j].y=(maxV-curv[j].v)+2*FACTOR_T_2_V;
				area=0;
			}
			// bottom
			else
				// when to use subCase?
				switch(area)
				{
				case 1:
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;				// same x
					c2D[j].y=curv[j].v+maxV+4*FACTOR_T_2_V;
					break;
				case 2:
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;				// same x
					c2D[j].y=curv[j].v-maxV;
					break;
				case 3:	
					c2D[j].y=(maxV-curv[j].v)+2*FACTOR_T_2_V;
					c2D[j].x=(maxU-curv[j].u)+maxU+4*FACTOR_T_2_U;	// same x
					break;
				case 4:
					c2D[j].y=(maxV-curv[j].v)+2*FACTOR_T_2_V;
					c2D[j].x=-curv[j].u;							// same x
					break;
				default:
					fprintf(stderr, "ConvertCorordinates failed!\n");
					bool trap=true;
					break;
				}
		}
		else
		{
	//
	//	   | 1 |
	//	-----------
	//	   |   |
	//	 4 | 0 | 3
	//	-----------
	//	   | 2 |
	//
			//if(area==0)
			//{
				//if(!NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t!=oppositeSide)
				//{
				// now to determine the whole scheme of mapping
				//		in other words, need to figure out which area to map the crest/bottom points
				if(subCase!=area) subCase=area;
				if(curv[j].v==0) area=1;
				else
					if(curv[j].v==maxV) area=2;
					else
						if(curv[j].u==maxU) area=3;
						else
							if(curv[j].u==0) area=4;
							else
							{
								fprintf(stderr, "ConvertCorordinates failed!\n");
								bool trap=true;
							}
				//}
				/*
				else
				{
					// the current vertex falls into the other side of the model, top(bottom) -> bottom(top)
					switch(subCase)
					{
					case 1:
						area=subCase;
						break;
					case 2:
						area=subCase;
						break;
					case 3:
						area=subCase;
						break;
					case 4:
						area=subCase;
						break;
					default:
						fprintf(stderr, "ConvertCorordinates failed!\n");
						bool trap=true;
						break;
					}
				}
				*/
			//}
			switch(area)
			{
			case 1:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t==oppositeSide)
				{
				}
				//else
				{
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
					c2D[j].y=(3+curv[j].t)*FACTOR_T_2_V+maxV;
				}
				break;
			case 2:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t==oppositeSide)
				{
				}
				//else
				{
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
					c2D[j].y=(1-curv[j].t)*FACTOR_T_2_V;
				}
				break;
			case 3:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t==oppositeSide)
				{
				}
				//else
				{
					c2D[j].y=(maxV-curv[j].v)+2*FACTOR_T_2_V;
					c2D[j].x=(3+curv[j].t)*FACTOR_T_2_U+maxU;
				}
				break;
			case 4:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t==oppositeSide)
				{
				}
				//else
				{
					c2D[j].y=(maxV-curv[j].v)+2*FACTOR_T_2_V;
					c2D[j].x=(1-curv[j].t)*FACTOR_T_2_U;
				}
				break;
			default:
				fprintf(stderr, "ConvertCorordinates failed!\n");
				break;
				//return;
			}
		}

		// time for the next vertex on the curve
		if(j==pNum-1)
			j=0;
		else
			j++;
	}//---------------------------------------------------------------------- TEMPORARILY COMMENTED
/*
	for(i=0; i<pNum; i++)
	{
		// this ain't right cause '+2' does not move the point accordingly in the uvt 2d space
		//*
		//if(NumericallyEqual(curv[j].t, -oppositeSide))
		if(curv[j].t==-oppositeSide)
		{
			// top
			c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
			c2D[j].y=curv[j].v+2*FACTOR_T_2_V;
			area=0;
		}
		else
		{
			if(area==0)
			{
				//if(!NumericallyEqual(curv[j].t, oppositeSide))
				if(curv[j].t!=oppositeSide)
				{
				// now to determine the whole scheme of mapping
				//		in other words, need to figure out which area to map the crest/bottom points
				if(curv[j].v==maxV) area=1;
				else
					if(curv[j].v==0) area=2;
					else
						if(curv[j].u==maxU) area=3;
						else
							if(curv[j].u==0) area=4;
							else
								bool trap=true;
				if(subCase<0) subCase=area;
				}
				else
				{
					// the current vertex falls into the other side of the model, top(bottom) -> bottom(top)
					// the current vertex falls into the other side of the model, top(bottom) -> bottom(top)
					switch(subCase)
					{
					case 1:
						area=subCase;
						break;
					case 2:
						area=subCase;
						break;
					case 3:
						area=subCase;
						break;
					case 4:
						area=subCase;
						break;
					default:
						fprintf(stderr, "ConvertCorordinates failed!\n");
						bool trap=true;
					}
				}
			}
			switch(area)
			{
			case 1:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				if(curv[j].t==oppositeSide)
				{
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
					c2D[j].y=(maxV-curv[j].v)+maxV+4*FACTOR_T_2_V;
				}
				else
				{
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
					c2D[j].y=(1+curv[j].t)+maxV+2*FACTOR_T_2_V;
				}
				break;
			case 2:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				if(curv[j].t==oppositeSide)
				{
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
					c2D[j].y=-curv[j].v;
				}
				else
				{
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
					c2D[j].y=1-curv[j].t;
				}
				break;
			case 3:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				if(curv[j].t==oppositeSide)
				{
					c2D[j].y=curv[j].v+2*FACTOR_T_2_V;
					c2D[j].x=(maxU-curv[j].u)+maxU+4*FACTOR_T_2_U;
				}
				else
				{
					c2D[j].y=curv[j].v+2*FACTOR_T_2_V;
					c2D[j].x=(1+curv[j].t)+maxU+2*FACTOR_T_2_U;
				}
				break;
			case 4:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				if(curv[j].t==oppositeSide)
				{
					c2D[j].y=curv[j].v+2*FACTOR_T_2_V;
					c2D[j].x=-curv[j].u;
				}
				else
				{
					c2D[j].y=curv[j].v+2*FACTOR_T_2_V;
					c2D[j].x=1-curv[j].t;
				}
				break;
			default:
				fprintf(stderr, "ConvertCorordinates failed!\n");
				break;
				//return;
			}
		}
		if(j==pNum-1)
			j=0;
		else
			j++;
	}
*/
}
	else
	{
		// then the whole intersection curve stays on the crest of the shape
		for(i=0; i<pNum; i++)
		{
	//
	//	   | 1 |
	//	-----------
	//	   |   |
	//	 4 | 0 | 3
	//	-----------
	//	   | 2 |
	//
			//if(area==0)
			//{
				//if(!NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t!=oppositeSide)
				//{
				// now to determine the whole scheme of mapping
				//		in other words, need to figure out which area to map the crest/bottom points
				if(subCase!=area) subCase=area;
				if(curv[j].v==maxV) area=1;
				else
					if(curv[j].v==0) area=2;
					else
						if(curv[j].u==maxU) area=3;
						else
							if(curv[j].u==0) area=4;
							else
							{
								fprintf(stderr, "ConvertCorordinates failed!\n");
								bool trap=true;
							}
				//}
				/*
				else
				{
					// the current vertex falls into the other side of the model, top(bottom) -> bottom(top)
					switch(subCase)
					{
					case 1:
						area=subCase;
						break;
					case 2:
						area=subCase;
						break;
					case 3:
						area=subCase;
						break;
					case 4:
						area=subCase;
						break;
					default:
						fprintf(stderr, "ConvertCorordinates failed!\n");
						bool trap=true;
						break;
					}
				}
				*/
			//}
			switch(area)
			{
			case 1:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t==oppositeSide)
				{
				}
				//else
				{
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
					c2D[j].y=(3-curv[j].t)*FACTOR_T_2_V+maxV;
				}
				break;
			case 2:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t==oppositeSide)
				{
				}
				//else
				{
					c2D[j].x=curv[j].u+2*FACTOR_T_2_U;
					c2D[j].y=(1+curv[j].t)*FACTOR_T_2_V;
				}
				break;
			case 3:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t==oppositeSide)
				{
				}
				//else
				{
					c2D[j].y=curv[j].v+2*FACTOR_T_2_V;
					c2D[j].x=(3-curv[j].t)*FACTOR_T_2_U+maxU;
				}
				break;
			case 4:
				//if(NumericallyEqual(curv[j].t, oppositeSide))
				//if(curv[j].t==oppositeSide)
				{
				}
				//else
				{
					c2D[j].y=curv[j].v+2*FACTOR_T_2_V;
					c2D[j].x=(1+curv[j].t)*FACTOR_T_2_U;
				}
				break;
			default:
				fprintf(stderr, "ConvertCorordinates failed!\n");
				break;
				//return;
			}

			// time for the next vertex on the curve
			if(j==pNum-1)
				j=0;
			else
				j++;
		}
	}

	// then assign back 0 to pNum
	c2D[pNum].x=c2D[0].x;
	c2D[pNum].y=c2D[0].y;

	c2D[pNum+1].x=subCase;
	c2D[pNum+1].y=type;
	//c2D[pNum+1].x=area;
	//c2D[pNum+1].y=-oppositeSide;
}


void Pointlist_serverB::DilatePoint(Bpoint2 *pnt, MyFloat maxU, MyFloat maxV, MyFloat deltaU, MyFloat deltaV)
{
	int theFlag;
	MyFloat deltaU1, deltaU2, deltaV1, deltaV2, deltaT;
	//MyFloat u, v, t;
	//MyFloat alpha=1;

	//link = figureTreeNode->getLink(i);

	/*
	if((pnt->u == 0 || pnt->u == maxU) && (centerV == 0 || centerV == maxV))
		theFlag = 1;
	if ((pnt->v == 0 || pnt->v == maxU) && (centerU == 0 || centerU == maxV))
		theFlag = 1;
	if(link == NULL)
		continue;
	*/

	// theFlag will be the flag for all the 'corner' points
	if((pnt->u == 0 || pnt->u == maxU) && (pnt->v == 0 || pnt->v == maxV))
		theFlag = 1;

	/*
	int flag=0;
	//int flag;	is this right? no initialization at all?
	if(pnt->t>-1)
		if(pnt->u==0 || pnt->u==maxU)
			flag=1;
		else 
			if(pnt->v==0||pnt->v == maxV)
				//else if (link->v == 0 || link->v == maxU)
				flag = 2;
			else
				flag = 0;	
	if(pnt->t<=-1 )
	{// on the bottom face
		if (flag == 1) // u crest;
		{
			deltaU=-deltaU;
			//deltaV=deltaV;
		}
		else 
			if (flag == 2) // v crest;
			{
				//deltaU=deltaU;
				deltaV=-deltaV;
			}
			else
			{
				deltaU=-deltaU;
				deltaV=-deltaV;
			}
	}
			/*
			if (link->t >-1 ){
				deltaU=theDeltaU;
				deltaV=theDeltaV;
			}
			else
			{// on the bottom face
				if (flag == 1) // u crest;
				{
					deltaU=-theDeltaU;
					deltaV=theDeltaV;
				}
				else if (flag == 2) // v crest;
				{
					deltaU=theDeltaU;
					deltaV=-theDeltaV;
				}
				else
				{
					deltaU=-theDeltaU;
					deltaV=-theDeltaV;
				}
			}
			//*/

	//u=pnt->u; v=pnt->v; t=pnt->t;
	//printf("before: u= %f, v= %f t= %f deltaU= %f deltaV= %f\n", pnt->u, pnt->v, pnt->t, deltaU, deltaV);

	if (pnt->t == -1 || pnt->t == 1)
	{
		// Neither u nor v is in the crest region.
		if (pnt->u + deltaU <= maxU && pnt->v + deltaV <= maxV && 			
			pnt->u + deltaU >= 0 && pnt->v + deltaV >= 0)
		{
			//Neither u+deltaU nor v+deltaV is in the crest region.
			pnt->u+=deltaU; 
			pnt->v+=deltaV;
		}
		else
		{
			if(pnt->u + deltaU > maxU)
			{
				// u+deltaU is in the crest region.
				deltaU1=maxU-pnt->u;
				deltaU2=deltaU-deltaU1;
				pnt->u+=deltaU1; 
				deltaT=deltaU2;// /alpha;
				//vDelta=pnt->v/maxV*(maxV+2*alpha)-pnt->v-alpha;
				//deltaT=deltaU2/sqrt(alpha*alpha+vDelta*vDelta);
				if (pnt->t == 1)
					pnt->t-=fabs(deltaT);
				else
					pnt->t+=fabs(deltaT);
			}
			else 
			{
				if (pnt->u + deltaU < 0)
				{
					// u+deltaU is in the crest region.
					deltaU1=0-pnt->u;
					deltaU2=deltaU-deltaU1;
					pnt->u+=deltaU1; 
					deltaT=deltaU2;// /alpha;
					//vDelta=pnt->v/maxV*(maxV+2*alpha)-pnt->v-alpha;
					//deltaT=deltaU2/sqrt(alpha*alpha+vDelta*vDelta);
					if (pnt->t == 1)
						pnt->t-=fabs(deltaT);
					else
						pnt->t+=fabs(deltaT);
				}
				else
					pnt->u+=deltaU;
			}

			if ( pnt->v + deltaV > maxV)			
			{ 		
				//pnt->v+deltaV is in the crest region.
				deltaV1=maxV-pnt->v;
				deltaV2=deltaV-deltaV1;
				pnt->v+=deltaV1;
				deltaT=deltaV2;// /alpha;
				//uDelta=pnt->u/maxU*(maxU+2*alpha)-pnt->u-alpha;
				//deltaT=deltaV2/sqrt(alpha*alpha+uDelta*uDelta);
				if (pnt->t == 1)
					pnt->t-=fabs(deltaT);
				else
					pnt->t+=fabs(deltaT); 
			}
			else 
			{
				if (pnt->v + deltaV < 0)
				{						
					//pnt->v + deltaV < 0) 
					//pnt->v+deltaV is in the crest region.
					deltaV1=0-pnt->v;
					deltaV2=deltaV-deltaV1;
					pnt->v+=deltaV1;
					deltaT=deltaV2;// /alpha;
					//uDelta=pnt->u/maxU*(maxU+2*alpha)-pnt->u-alpha;
					//deltaT=deltaV2/sqrt(alpha*alpha+uDelta*uDelta);
					if (pnt->t == 1)
						pnt->t-=fabs(deltaT);
					else
						pnt->t+=fabs(deltaT); 
				}
				else
					pnt->v+=deltaV;
			}
		}
	} // Neither u nor v is but either u+deltaU or v+deltaV is in the crest region.
	else 
		// u/v is in the crest
		// need to figure out which crest it's in
		if ((pnt->u == 0) || (pnt->u == maxU))
		{
			// u is in the crest region.				
//if (deltaU != 0)
			{ 
				// moving in u direction first!
				if (theFlag == 1)
				{
					// actually moving in v direction
					if ( pnt->v == 0)
						pnt->v=fabs(deltaU);
					else
						pnt->v-=fabs(deltaU);
				}
				else
				{
					//vDelta=pnt->v/maxV*(maxV+2*alpha)-pnt->v-alpha;
					//if (pnt->t >= 0)
					if (pnt->u == maxU)
						deltaT=-deltaU;// /alpha;
					//deltaT=-deltaU/sqrt(alpha*alpha+vDelta*vDelta);
					else 
						deltaT=deltaU;// /alpha;
					//deltaT=deltaU/sqrt(alpha*alpha+vDelta*vDelta);
					if (pnt->t+deltaT >= -1 && pnt->t+deltaT <= 1)
					{
						// still ends in the crest region.
						pnt->t+=deltaT; 
					}
					else
					{
						// ends in the non-crest region.
						if (pnt->t+deltaT < -1)
						{
							deltaT=-1-pnt->t;
							pnt->t=-1;

//deltaU1=deltaT;// *alpha;

							//deltaU1=deltaT*sqrt(alpha*alpha+vDelta*vDelta);
							// this segment of codes are pretty weird
//deltaU2=deltaU-deltaU1;
							if (pnt->u == 0)
//pnt->u+=fabs(deltaU-deltaU1);
								pnt->u=fabs(deltaT);
							else
//pnt->u-=fabs(deltaU+deltaU1);
								pnt->u-=fabs(deltaT);
						}
						else
						{
							deltaT=1-pnt->t;
							pnt->t=1;
//deltaU1=deltaT;// *alpha;
							//deltaU1=deltaT*sqrt(alpha*alpha+vDelta*vDelta);
							//deltaU2=deltaU-deltaU1;
							if (pnt->u == 0)
//pnt->u+=fabs(deltaU-deltaU1);
								pnt->u=fabs(deltaT);
							else
								pnt->u-=fabs(deltaT);						
						}
					}
				}
			}
//else 
			{ 
				// then moving in v direction.
				// to add codes here to actaully deal with moving in v direction after in u direction
				if (pnt->v+deltaV >= 0 && pnt->v+deltaV <= maxV)
				{
					// still in the range of V, normal movement
					pnt->v=pnt->v+deltaV;
				}
				else
				{
					if (pnt->v+deltaV < 0)
					{
						deltaV1=0-pnt->v;
						pnt->v=0;
						deltaU2=deltaV-deltaV1;
						if ( pnt->u == 0)
							pnt->u=fabs(deltaU2);
						else
							pnt->u=pnt->u-fabs(deltaU2);
					}
					else
					{
						deltaV1=maxV-pnt->v;
						pnt->v=maxV;
						deltaU2=deltaV-deltaV1;
						if ( pnt->u == 0)
							pnt->u=fabs(deltaU2);
						else
							pnt->u=pnt->u-fabs(deltaU2);
					}
				}
			}
		} // for the u cresit

		else
		{
			// v is in the crest region.
//if (deltaV != 0)
			{ 
				// moving in v direction first!
				if (theFlag == 1)
				{
					// actually moving in u direction
					if ( pnt->u == 0)
						pnt->u=fabs(deltaV);
					else
						pnt->u=pnt->u-fabs(deltaV);
				}				
				else
				{
					deltaT=deltaV;// /alpha;
					//uDelta=pnt->u/maxU*(maxU+2*alpha)-pnt->u-alpha;
					//deltaT=deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					//if (pnt->t >= 0)
					if (pnt->v == maxV)
						deltaT=-deltaV;// /alpha;
						//deltaT=-deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					else 
						deltaT=deltaV;// /alpha;
						//deltaT=deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					if (pnt->t+deltaT >= -1 && pnt->t+deltaT <= 1)
					{
						// still ends in the crest region.
						pnt->t=pnt->t+deltaT;
					}
					else
					{
						// ends in the non-crest region.
						if (pnt->t+deltaT < -1)
						{
							deltaT=-1-pnt->t;
							pnt->t=-1;
							deltaV1=deltaT;// *alpha;
							//deltaV1=deltaT*sqrt(alpha*alpha+uDelta*uDelta);
							//deltaV2=deltaV-deltaV1;
							if (pnt->v == 0)
								pnt->v=pnt->v+fabs(deltaV-deltaV1);
							else
								pnt->v=pnt->v-fabs(deltaV+deltaV1);
						}
						else
						{
							deltaT=1-pnt->t;
							pnt->t=1;
							deltaV1=deltaT;// *alpha;
							//deltaV1=deltaT*sqrt(alpha*alpha+uDelta*uDelta);
							//deltaV2=deltaV-deltaV1;
							if (pnt->v == 0)
								pnt->v=pnt->v+fabs(deltaV-deltaV1);
							else
								pnt->v=pnt->v-fabs(deltaV+deltaV1);						
						}
					}
				}
			}
//else 
			{ 
				// 110602
				// to add codes here to actaully deal with moving in u direction after in v direction

				// then moving in u direction.
				if (pnt->u+deltaU <= maxU && pnt->u+deltaU >= 0)
				{
					// still in the range of U, normal movement
					pnt->u=pnt->u+deltaU;
				}
				else
				{
					if (pnt->u+deltaU < 0)
					{
						deltaU1=0-pnt->u;
						pnt->u=0;
						deltaV2=deltaU-deltaU1;
						if ( pnt->v == 0)
							pnt->v=fabs(deltaV2);
						else
							pnt->v=pnt->v-fabs(deltaV2);
					}
					else
					{ 
						// pnt->u+deltaU <= maxU
						deltaU1=maxU-pnt->u;
						pnt->u=maxU;
						deltaV2=deltaU-deltaU1;
						if ( pnt->v == 0)
							pnt->v=fabs(deltaV2);
						else
							pnt->v=pnt->v-fabs(deltaV2);
					}
				}
			}
		} // for the v crest
//	printf("after: u= %f, v= %f t= %f\n", pnt->u, pnt->v, pnt->t);
}


#define MAGNITUDE_EPSILON_2D_VECTOR 0.000001
#define MAGNITUDE_EPSILON_2D_DILATION_AMOUNT 0.001
bool Normalize2D(MyPoint2D *p)
{
	MyFloat m=p->MagSquare();
	if(m<MAGNITUDE_EPSILON_2D_VECTOR)
	{
#ifndef QIONG_BLENDING_NO_DEBUG
		fprintf(stderr, "ERROR M3DBlendedRend::Normalize2D: vector is too close to 0 to normalize!\n");
#endif
		return false;
	}
	m=sqrt(m);
	*p=*p/m;
	//p->x=p->x/m;
	//p->y=p->y/m;
	return true;
}
bool Valid2D(MyPoint2D p)
{
	MyFloat m=p.MagSquare();
#ifndef QIONG_BLENDING_NO_DEBUG
	fprintf(stderr, "The dilation amount is u/v: %f/%f\n", p.x, p.y);
#endif
	if(m<MAGNITUDE_EPSILON_2D_DILATION_AMOUNT)
	{
#ifndef QIONG_BLENDING_NO_DEBUG
		fprintf(stderr, "ERROR M3DBlendedRend::Valid2D: vector is too close to 0!\n");
#endif
		return false;
	}
	return true;
}
MyFloat DotProduct2D(MyPoint2D p1, MyPoint2D p2)
{
	return p1*p2;
}


// returned value:
//					false: DO NOT need to reverse the intersection curve
//					true:  Do need to reverse the intersection curve
bool Pointlist_serverB::CheckIntersectionCurveOrientation(Bpoint2 *subC, int subPNum, Bpoint2 *intC, int intPNum)
{
	int orientation;
	MyFloat tagTotal=0;
	MyFloat theSign;

	int i;//, j;
	for(i=0; i<subPNum; i++)
	{
		/*
		if(subC[i].t_2<0)
			theSign=-1;
		else
			if(subC[i+1].t_2)<=0
				theSign=-1;
			else
				theSign=1;
		*/
		theSign=(subC[i].t+subC[i+1].t)/2;
		if(theSign<0) theSign=-1;		
		else theSign=1;

		tagTotal+=(subC[i+1].v-subC[i].v)*theSign;
	}
	if(tagTotal>0)
		orientation=1;
	else
		if(tagTotal<0)
			orientation=-1;
		else
		{
			fprintf(stderr, "\nThen this is the case that should never happen in CheckIntersectionCurveOrientation!\n");
			return false;
		}

	tagTotal=0;
	for(i=0; i<intPNum; i++)
	{
		theSign=(intC[i].t_2+intC[i+1].t_2)/2;
		if(theSign<0) theSign=-1;		
		else theSign=1;
		tagTotal+=(intC[i+1].v_2-intC[i].v_2)*theSign;
	}
	tagTotal*=orientation;
	if(tagTotal>0)
	{
		// then the intersection curve has the same orientation as the sub-cut-curve
		// which means we DO NOT need to reverse the intersection curve!
		return false;
	}
	else
		if(tagTotal<0)
		{
			// then the intersection curve has the different orientation as the sub-cut-curve
			// which means we DO need to reverse the intersection curve!
			return true;
		}
		else
		{
			fprintf(stderr, "\nThen this is the case that should never happen in CheckIntersectionCurveOrientation!\n");
			return false;
		}
}

bool IsCounterClockwise(MyPoint2D *c2D, int pNum)
{
	int i, j;
	MyFloat theSign, dot, angle;
	MyPoint2D center, vec_0, vec_1;
	center.Set(0, 0);
	angle=0;
	for(i=0; i<pNum; i++)
	{
		center=center+c2D[i];
	}
	center=center/pNum;

	for(i=0; i<pNum; i++)
	{
		if(i==pNum-1)
			j=0;
		else
			j=i+1;
		vec_0=c2D[i]-center;
		vec_1=c2D[j]-center;
		theSign=vec_0.x*vec_1.y-vec_0.y*vec_1.x; //c3=a1b2-a2b1
		if(theSign>0)
			theSign=1;
		else
			if(theSign==0)
				theSign=0;
			else
				theSign=-1;
		Normalize2D(&vec_0);
		Normalize2D(&vec_1);
		//vec_0.Normalize();
		//vec_1.Normalize();
		dot=vec_0*vec_1;
		angle+=theSign*acos(dot);
	}
	return (angle>0);
}

int CrestRegion(int type, Bpoint2 bp, MyFloat maxU, MyFloat maxV)
{
	switch(type)
	{
	case 1:
	if(bp.v==maxV || bp.v==0)
	{
		// area 1/2
		if(bp.v==maxV)
		{
			// area 1
			return 1;
		}
		else
		{
			// area 2
			return 2;
		}
	}
	else
		if(bp.u==maxU || bp.u==0)
		{
			// area 3/4
			if(bp.u==maxU)
			{
				// area 3
				return 3;
			}
			else
			{
				// area 4
				return 4;
			}
		}
		else
		{
			// this is the case that is not allowed!
			// area ???
			return 0;
		}
	return true;

	case -1:
	if(bp.v==maxV || bp.v==0)
	{
		// area 1/2
		if(bp.v==0)
		{
			// area 1
			return 1;		
		}
		else
		{
			// area 2
			return 2;
		}
	}
	else
		if(bp.u==maxU || bp.u==0)
		{
			// area 3/4
			if(bp.u==maxU)
			{
				// area 3
				return 3;			
			}
			else
			{
				// area 4
				return 4;
			}
		}
		else
		{
			// this is the case that is not allowed!
			// area ???
			return 0;
		}

	default:
	// this is the case that is not allowed!
	return 0;
	}
}

bool Pointlist_serverB::ConvertUVT2XY_Bottom(int type, int area, Bpoint2 curv, MyPoint2D *curv2D, MyFloat maxU, MyFloat maxV)
{
	switch(type)
	{
	case 1:
		switch(area)
		{
		case 1:
			curv2D->x=curv.u+2*FACTOR_T_2_U;
			curv2D->y=(maxV-curv.v)+maxV+4*FACTOR_T_2_V;
			break;
		case 2:
			curv2D->x=curv.u+2*FACTOR_T_2_U;
			curv2D->y=-curv.v;
			break;

		case 3:
			curv2D->y=curv.v+2*FACTOR_T_2_V;
			curv2D->x=(maxU-curv.u)+maxU+4*FACTOR_T_2_U;
			break;
		case 4:
			curv2D->y=curv.v+2*FACTOR_T_2_V;
			curv2D->x=-curv.u;
			break;
		default:
			return false;
		}
		return true;

	case -1:
		switch(area)
		{
		case 1:
			curv2D->x=curv.u+2*FACTOR_T_2_U;				// same x
			curv2D->y=curv.v+maxV+4*FACTOR_T_2_V;
			break;
		case 2:
			curv2D->x=curv.u+2*FACTOR_T_2_U;				// same x
			curv2D->y=curv.v-maxV;
			break;
		case 3:	
			curv2D->y=(maxV-curv.v)+2*FACTOR_T_2_V;
			curv2D->x=(maxU-curv.u)+maxU+4*FACTOR_T_2_U;	// same x
			break;
		case 4:
			curv2D->y=(maxV-curv.v)+2*FACTOR_T_2_V;
			curv2D->x=-curv.u;								// same x
			break;
		default:
			return false;
		}
		return true;

	default:
	// this is the case that is not allowed!
		return false;
	}
}

bool Pointlist_serverB::ConvertUVT2XY_Crest(int type, Bpoint2 curv, MyPoint2D *curv2D, MyFloat maxU, MyFloat maxV)
{
	switch(type)
	{
	case 1:
	if(curv.v==maxV || curv.v==0)
	{
		// area 1/2
		if(curv.v==maxV)
		{
			// area 1
			curv2D->x=curv.u+2*FACTOR_T_2_U;
			curv2D->y=(3-curv.t)*FACTOR_T_2_V+maxV;
		}
		else
		{
			// area 2
			curv2D->x=curv.u+2*FACTOR_T_2_U;
			curv2D->y=(1+curv.t)*FACTOR_T_2_V;
		}
	}
	else
		if(curv.u==maxU || curv.u==0)
		{
			// area 3/4
			if(curv.u==maxU)
			{
				// area 3
				curv2D->x=(3-curv.t)*FACTOR_T_2_U+maxU;
				curv2D->y=curv.v+2*FACTOR_T_2_V;
			}
			else
			{
				// area 4
				curv2D->x=(1+curv.t)*FACTOR_T_2_U;
				curv2D->y=curv.v+2*FACTOR_T_2_V;
			}
		}
		else
		{
			// this is the case that is not allowed!
			// area ???
			return false;
		}
	return true;

	case -1:
	if(curv.v==maxV || curv.v==0)
	{
		// area 1/2
		if(curv.v==0)
		{
			// area 1
			curv2D->x=curv.u+2*FACTOR_T_2_U;
			curv2D->y=(1+curv.t)*FACTOR_T_2_V;
		}
		else
		{
			// area 2
			curv2D->x=curv.u+2*FACTOR_T_2_U;
			curv2D->y=(-1-curv.t)*FACTOR_T_2_V-maxV;
		}
	}
	else
		if(curv.u==maxU || curv.u==0)
		{
			// area 3/4
			if(curv.u==maxU)
			{
				// area 3
				curv2D->x=(3-curv.t)*FACTOR_T_2_U+maxU;
				curv2D->y=-curv.v;
			}
			else
			{
				// area 4
				curv2D->x=(1+curv.t)*FACTOR_T_2_U;
				curv2D->y=-curv.v;
			}
		}
		else
		{
			// this is the case that is not allowed!
			// area ???
			return false;
		}
	return true;

	default:
	// this is the case that is not allowed!
	return false;
	}
}

bool Pointlist_serverB::DilateCurve(ThallCode::Pointlist_server2 *pList, Bpoint2 *curv, int pNum, MyFloat dt, int steps)
{
	//bool need2ReverseIntersectionCurve;
	MyFloat maxU=pList->max_u(), maxV=pList->max_v();
	MyPoint2D curv2D[200];
	MyPoint2D delta[200];
	int i, s, j;
	int area, type;//oppositeSide;
	MyPoint2D tmpP, tmpP2, centerP;
	//MyFloat tmpX, tmpY;
	MyFloat theSign;
	MyPoint2D tmpPSum, tmpPSumNormalized;

	int pDex[3];
	//int pDex0, pDex1, pDex2; 

	bool algorithmFailed=false;

	// first, need to figure out which area the curve lies in
	i=0;
	while(i<pNum)
	{
		if(NumericallyEqual(curv[i].t, 1.0))
		{
				// the intersection curve lies at least partially on the top
			type=1;
				//oppositeSide=-1;
				//startPos=i;
			break;
		}
		i++;
	}
	if(i==pNum)
	{
		// either the intersection curve is on the bottom
		// or on the crest
		type=0;
	}

	for(s=0; s<steps; s++)
	{
		// instead of using the function 'ConvertCoordinates'
		// I process the intersection curve chunk by chunk
		// so far it's the only way to deal with the corner problem
		// therefore we have to know the orientation of the curve first

#ifdef QIONG_BLENDING_DILATION
		// convert the UVT coordinates into the Euclidean 2D space
		ConvertCoordinates(curv2D, curv, maxU, maxV, pNum);
		MyFloat color[3];
		color[0]=0;
		color[1]=1;
		color[2]=0;
		Bpoint2 dCurv[200];
		if(s==0)
		{
			for(i=0; i<pNum+1; i++)
			{
				dCurv[i].pnt[0]=curv2D[i].x*2;
				dCurv[i].pnt[1]=curv2D[i].y*2;
				dCurv[i].pnt[2]=.65;
			}
			if(s==0 || s==steps-1)
				//pList->glRenderCurve(0, pNum, dCurv, color);
				glRenderCurve(0, pNum, dCurv, color);
		}
#endif

// used one extra integer to carry the area information, this should be updated later!
//		area=curv2D[pNum+1].x;
//		type=curv2D[pNum+1].y; //oppositeSide=curv2D[pNum+1].y;

		// calculate the diffusion direction for each single point on the curve
		//		using the coordinates in the Euclidean 2D
//		centerP.x=curv2D[0].x;
//		centerP.y=curv2D[0].y;
//		for(i=1; i<pNum; i++)
//		{
//			centerP=centerP+curv2D[i];
			//centerP.x+=curv2D[i].x;
			//centerP.y+=curv2D[i].y;
//		}
//		centerP=centerP/pNum;
		//centerP.x/=pNum;
		//centerP.y/=pNum;

		for(i=0; i<pNum; i++)
		{
			if(i==0)
				pDex[0]=pNum-1;
			else
				pDex[0]=i-1;
			pDex[1]=i;
			if(i==pNum-1)
				pDex[2]=0;
			else
				pDex[2]=i+1;

			// first, need to figure out which area the current vertex lies in
			// then convert the uvt to the 2D xy space 
			if(curv[pDex[1]].t==1.0)
			{
				// on the top
				curv2D[pDex[1]].x=curv[pDex[1]].u+2*FACTOR_T_2_U;
				curv2D[pDex[1]].y=curv[pDex[1]].v+2*FACTOR_T_2_V;
				for(j=0; j<3; j+=2)
				{
					if(curv[pDex[j]].t==-1)				
					{
						// this is the case that is not allowed!
						algorithmFailed=true;
					}
					else
						if(curv[pDex[j]].t==1.0)
						{
							curv2D[pDex[j]].x=curv[pDex[j]].u+2*FACTOR_T_2_U;
							curv2D[pDex[j]].y=curv[pDex[j]].v+2*FACTOR_T_2_V;
						}
						else
						{
							if(!ConvertUVT2XY_Crest(1, curv[pDex[j]], &curv2D[pDex[j]], maxU, maxV))
							{
								// the algorithm failed!
								algorithmFailed=true;
							}
						}
				}
			}
			else
				if(curv[pDex[1]].t==-1.0)
				{
					// on the bottom
					curv2D[pDex[1]].x=curv[pDex[1]].u+2*FACTOR_T_2_U;
					curv2D[pDex[1]].y=-curv[pDex[1]].v;
					for(j=0; j<3; j+=2)
					{
						if(curv[pDex[j]].t==1.0)
						{
							// this is the case that is not allowed!
							algorithmFailed=true;
						}
						else
							if(curv[pDex[j]].t==-1.0)
							{
								curv2D[pDex[j]].x=curv[pDex[j]].u+2*FACTOR_T_2_U;
								curv2D[pDex[j]].y=-curv[pDex[j]].v;
							}
							else
							{
								if(!ConvertUVT2XY_Crest(-1, curv[pDex[j]], &curv2D[pDex[j]], maxU, maxV))
								{
									// the algorithm failed!
									algorithmFailed=true;
								}
							}
					}
				}
				else
				{
					// on the crest
					if(!(area=CrestRegion(1, curv[pDex[1]], maxU, maxV)))
					{
						// this is the case that is not allowed!
						algorithmFailed=true;
					}
					if(!ConvertUVT2XY_Crest(1, curv[pDex[1]], &curv2D[pDex[1]], maxU, maxV))
					{
						// this is the case that is not allowed!
						algorithmFailed=true;
					}
					for(j=0; j<3; j+=2)
					{
						if(curv[pDex[j]].t==1.0)
						{
							curv2D[pDex[j]].x=curv[pDex[j]].u+2*FACTOR_T_2_U;
							curv2D[pDex[j]].y=curv[pDex[j]].v+2*FACTOR_T_2_V;
						}
						else
							if(curv[pDex[j]].t==-1.0)
							{
								if(!ConvertUVT2XY_Bottom(1, area, curv[pDex[j]], &curv2D[pDex[j]], maxU, maxV))
								{
									// the algorithm failed!
									algorithmFailed=true;
								}
							}
							else
							{
								// to do: update the XY for pDex[0] and pDex[2] according to 'area'
								// this is vital to solve the corner problem!!!!!!!!!!!!!!!!!!!!!!!
// right here right now
								if(!ConvertUVT2XY_Crest(1, curv[pDex[j]], &curv2D[pDex[j]], maxU, maxV))
								{
									// the algorithm failed!
									algorithmFailed=true;
								}
							}
					}
				}

			if(algorithmFailed)
				return false;

#define MAGNITUDE_EPSILON_2D_COLINEAR .01
#define EPSILON_2D_PERPENDICULAR 0.1
			// time to determine what
			// curv2D[pDex0, pDex1 & pDex2] are

			// need to check the flatness of the neighbors of pDex1 here by 'pDex0---pDex1---pDex2'
			// if it's too flat, need to assign the dilation direction perpendicular to the segments
			tmpP=curv2D[pDex[1]]-curv2D[pDex[0]];
			//tmpP.x=curv2D[pDex1].x-curv2D[pDex0].x;
			//tmpP.y=curv2D[pDex1].y-curv2D[pDex0].y;
			tmpP2=curv2D[pDex[1]]-curv2D[pDex[2]];
			//tmpP2.x=curv2D[pDex1].x-curv2D[pDex2].x;
			//tmpP2.y=curv2D[pDex1].y-curv2D[pDex2].y;
			Normalize2D(&tmpP);
			Normalize2D(&tmpP2);

			// this is purely weird
			theSign=tmpP.x*tmpP2.y-tmpP.y*tmpP2.x; //c3=a1b2-a2b1

			// no need to reverse the intersection curve later
			if(theSign>0)
			{
				// convex
				tmpPSum=tmpP+tmpP2;
				tmpPSumNormalized=tmpPSum;
				Normalize2D(&tmpPSumNormalized);
				tmpPSumNormalized=tmpPSumNormalized*2;
				delta[pDex[1]]=tmpPSumNormalized-tmpPSum;
			}
			else
				if(theSign==0)
				{
					// flat
					// delta[pDex1]=tmpP*i;
					delta[pDex[1]].x=-tmpP.y;
					delta[pDex[1]].y=tmpP.x;
					delta[pDex[1]]=delta[pDex[1]]*2;
				}
				else
				{
					// concave
					tmpPSum=tmpP+tmpP2;
					tmpPSumNormalized=tmpPSum;
					Normalize2D(&tmpPSumNormalized);
					tmpPSumNormalized=tmpPSumNormalized*(-2);
					delta[pDex[1]]=tmpPSumNormalized-tmpPSum;
				}

			// this part is NOT useful any longer!
			// check whether the dilation vector is pointing away from the 'center' of the curve
//			tmpP=curv2D[pDex1]-centerP;
			// if it is, reverse it!
			// 110702 this algorithm is quite buggy! center doesn't mean everything
			//			realy need to come up with a much more robust way to test the dilation direction!
			// already fixed!
//			if(DotProduct2D(tmpP, delta[pDex1])<-EPSILON_2D_PERPENDICULAR)
//			{
//				delta[pDex1].x=-delta[pDex1].x;
//				delta[pDex1].y=-delta[pDex1].y;
//			}

			if(curv[pDex[1]].t==-1.0)
			{
				//delta[pDex[1]].x=-delta[pDex[1]].x;
				delta[pDex[1]].y=-delta[pDex[1]].y;
			}

		}

		// dialate the curve in the UVT orordinates 2D space
		//		by the direction calculated in Euclidean 2D
		for(i=0; i<pNum; i++)
		{
			//DilatePoint(&curv[i], maxU, maxV, 0, 0.1);//dt*delta[i].x, dt*delta[i].y);
			DilatePoint(&curv[i], maxU, maxV, dt*delta[i].x, dt*delta[i].y);
#ifdef QIONG_BLENDING_DILATION_DEBUG
			fprintf(stderr, "The dilation amount is u/v: %f/%f\n", delta[i].x, delta[i].y);
#endif
		}

		/*
		// I don't actually need this part now, because my dialation still happens in UVT coordinates
		// convert the Euclidean 2D space back to the UVT coordinates 
		// ConvertCoordinates(curv, curv2D, maxU, maxV, interPNum);
		*/
#ifdef QIONG_BLENDING_DILATION
		color[0]=0;
		color[1]=1;
		color[2]=0;
		//if(s%2==0)
		{
		for(i=0; i<pNum+1; i++)
		{
			dCurv[i].pnt[0]=curv2D[i].x*2;
			dCurv[i].pnt[1]=curv2D[i].y*2;
			dCurv[i].pnt[2]=.65;
		}
		if(s==0 || s==steps-1)
			//pList->glRenderCurve(0, pNum, dCurv, color);
			glRenderCurve(0, pNum, dCurv, color);
		}
#endif
	}

	// after all the mess, need to map the points in UVT coordinates back to 3D
	Bpoint bp;
	for(i=0; i<pNum; i++)
	{
		// map each point back by calling the UVT->Boundary point routin in the Pointlist_server2
		//pList->subdivBPosition(&curv[i], curv[i].u, curv[i].v, curv[i].t);
		pList->subdivBPosition(&bp, curv[i].u, curv[i].v, curv[i].t);
		//		to assign a Bpoint object to a Bpoint2 object, only the values of
		//		pnt/norm/rad are tranferred so far!
		curv[i]=bp;
	}

	curv[pNum]=curv[0];
	//return need2ReverseIntersectionCurve;
	return true;
}


// right here right now
void MeshCurves(int *triNum, int *triList, int *c2StartPos, int pNum1, Bpoint2 *c1, int pNum2, Bpoint2 *c2)
{
	// connect the 2 curves in mainInterCurve[mainInterPNum] and interCurveD[interPNum]
	int i, k;

	bool goWithC1;
	bool c1Done=false, c2Done=false;
	bool done=false;
	int c1Dex, c2Dex;
	int c1Next, c2Next;
	int c1Picked, c2Picked;

	//int c2StartPos=0;
	*c2StartPos=0;

	//MyFloat diff1;
	MyFloat dis1, dis2, dis;
	MyFloat angleC1, angleC2;
	MyFloat angle, angle1, angle2;
	//MyFloat dis3;

	dis1=9999;
	for(i=0; i<pNum2; i++)
	{
		dis2=(c1[0].pnt[0]-c2[i].pnt[0])*(c1[0].pnt[0]-c2[i].pnt[0])+
			 (c1[0].pnt[1]-c2[i].pnt[1])*(c1[0].pnt[1]-c2[i].pnt[1])+
			 (c1[0].pnt[2]-c2[i].pnt[2])*(c1[0].pnt[2]-c2[i].pnt[2]);
//		dis2=(curv1[mCurveStartDex].pnt[0]-curv2[i].pnt[0])*(curv1[mCurveStartDex].pnt[0]-curv2[i].pnt[0])+
//			 (curv1[mCurveStartDex].pnt[1]-curv2[i].pnt[1])*(curv1[mCurveStartDex].pnt[1]-curv2[i].pnt[1])+
//			 (curv1[mCurveStartDex].pnt[2]-curv2[i].pnt[2])*(curv1[mCurveStartDex].pnt[2]-curv2[i].pnt[2]);
		if(dis1>dis2)
		{
			dis1=dis2;
			*c2StartPos=i;
		}
	}


	Bpoint interDir, interDir1, interDir2, interDir3;
	MyFloat dot;//, dot1, dot2;
	c1Dex=0;
	c2Dex=*c2StartPos;
	c1Picked=0;
	c2Picked=0;

	*triNum=0;

	bool secondThoughtWithC1;
	// start the intersection curve from the closest point to the first poin on the main-figure curve
	while(!(c1Done&&c2Done))
	{
		//int mainPrev=INTMOD(mainDex-1, mainInterPNum);
		//int interPrev=INTMOD(interDex-1, interPNum);
		//int interNext2=INTMOD(interDex+2, interPNum);

		c1Next=(c1Dex+1)%pNum1;
		c2Next=(c2Dex+1)%pNum2;
		if(!c1Done && !c2Done)
		{
			#define TRIANGULATION_CRETERION .9
			for(k=0; k<3; k++)
				interDir1.pnt[k]=(c2[c2Dex].pnt[k]+c1[c1Next].pnt[k])/2;
			for(k=0; k<3; k++)
				interDir.pnt[k]=c2[c2Dex].pnt[k]-interDir1.pnt[k];
			dis=sqrt(interDir.pnt[0]*interDir.pnt[0]+interDir.pnt[1]*interDir.pnt[1]+interDir.pnt[2]*interDir.pnt[2]);
			for(k=0; k<3; k++)
				interDir.pnt[k]=c2[c2Next].pnt[k]-interDir1.pnt[k];
			dis1=sqrt(interDir.pnt[0]*interDir.pnt[0]+interDir.pnt[1]*interDir.pnt[1]+interDir.pnt[2]*interDir.pnt[2]);
			for(k=0; k<3; k++)
				interDir.pnt[k]=c1[c1Dex].pnt[k]-interDir1.pnt[k];
			dis2=sqrt(interDir.pnt[0]*interDir.pnt[0]+interDir.pnt[1]*interDir.pnt[1]+interDir.pnt[2]*interDir.pnt[2]);
			if(dis1<dis2)
				angleC1=dis1/dis;
			else 
				angleC1=dis2/dis;

			for(k=0; k<3; k++)
				interDir2.pnt[k]=(c2[c2Next].pnt[k]+c1[c1Dex].pnt[k])/2;
			for(k=0; k<3; k++)
				interDir.pnt[k]=c2[c2Next].pnt[k]-interDir2.pnt[k];
			dis=sqrt(interDir.pnt[0]*interDir.pnt[0]+interDir.pnt[1]*interDir.pnt[1]+interDir.pnt[2]*interDir.pnt[2]);
			for(k=0; k<3; k++)
				interDir.pnt[k]=c2[c2Dex].pnt[k]-interDir2.pnt[k];
			dis1=sqrt(interDir.pnt[0]*interDir.pnt[0]+interDir.pnt[1]*interDir.pnt[1]+interDir.pnt[2]*interDir.pnt[2]);
			for(k=0; k<3; k++)
				interDir.pnt[k]=c1[c1Next].pnt[k]-interDir2.pnt[k];
			dis2=sqrt(interDir.pnt[0]*interDir.pnt[0]+interDir.pnt[1]*interDir.pnt[1]+interDir.pnt[2]*interDir.pnt[2]);
			if(dis1<dis2)
				angleC2=dis1/dis;
			else 
				angleC2=dis2/dis;

			if(angleC1<TRIANGULATION_CRETERION && angleC2<TRIANGULATION_CRETERION)
			{
			for(k=0; k<3; k++)
				interDir1.pnt[k]=(c1[c1Next].pnt[k]-c1[c1Dex].pnt[k]);
			dot=sqrt(interDir1.pnt[0]*interDir1.pnt[0]+interDir1.pnt[1]*interDir1.pnt[1]+interDir1.pnt[2]*interDir1.pnt[2]);
			for(k=0; k<3; k++)
				interDir1.pnt[k]=interDir1.pnt[k]/dot;
			for(k=0; k<3; k++)
				interDir2.pnt[k]=(c2[c2Dex].pnt[k]-c1[c1Next].pnt[k]);
			dot=sqrt(interDir2.pnt[0]*interDir2.pnt[0]+interDir2.pnt[1]*interDir2.pnt[1]+interDir2.pnt[2]*interDir2.pnt[2]);
			for(k=0; k<3; k++)
				interDir2.pnt[k]=interDir2.pnt[k]/dot;
			for(k=0; k<3; k++)
				interDir3.pnt[k]=(c1[c1Dex].pnt[k]-c2[c2Dex].pnt[k]);
			dot=sqrt(interDir3.pnt[0]*interDir3.pnt[0]+interDir3.pnt[1]*interDir3.pnt[1]+interDir3.pnt[2]*interDir3.pnt[2]);
			for(k=0; k<3; k++)
				interDir3.pnt[k]=interDir3.pnt[k]/dot;
			angle1=fabs(interDir1.pnt[0]*interDir2.pnt[0]+interDir1.pnt[1]*interDir2.pnt[1]+interDir1.pnt[2]*interDir2.pnt[2]);
			//angle1=fabs(interDir2.pnt[0]*interDir3.pnt[0]+interDir2.pnt[1]*interDir3.pnt[1]+interDir2.pnt[2]*interDir3.pnt[2]);
			angle=fabs(interDir2.pnt[0]*interDir3.pnt[0]+interDir2.pnt[1]*interDir3.pnt[1]+interDir2.pnt[2]*interDir3.pnt[2]);
			if(angle1<angle) angle1=angle;
			//angle1=fabs(interDir3.pnt[0]*interDir1.pnt[0]+interDir3.pnt[1]*interDir1.pnt[1]+interDir3.pnt[2]*interDir1.pnt[2]);
			angle=fabs(interDir3.pnt[0]*interDir1.pnt[0]+interDir3.pnt[1]*interDir1.pnt[1]+interDir3.pnt[2]*interDir1.pnt[2]);
			if(angle1<angle) angle1=angle;
			////////////////////////////
			for(k=0; k<3; k++)
				interDir1.pnt[k]=(c1[c1Next].pnt[k]-c2[c2Next].pnt[k]);
			dot=sqrt(interDir1.pnt[0]*interDir1.pnt[0]+interDir1.pnt[1]*interDir1.pnt[1]+interDir1.pnt[2]*interDir1.pnt[2]);
			for(k=0; k<3; k++)
				interDir1.pnt[k]=interDir1.pnt[k]/dot;
			for(k=0; k<3; k++)
				interDir3.pnt[k]=(c2[c2Next].pnt[k]-c2[c2Dex].pnt[k]);
			dot=sqrt(interDir3.pnt[0]*interDir3.pnt[0]+interDir3.pnt[1]*interDir3.pnt[1]+interDir3.pnt[2]*interDir3.pnt[2]);
			for(k=0; k<3; k++)
				interDir3.pnt[k]=interDir3.pnt[k]/dot;
			angle=fabs(interDir1.pnt[0]*interDir2.pnt[0]+interDir1.pnt[1]*interDir2.pnt[1]+interDir1.pnt[2]*interDir2.pnt[2]);
			if(angle1<angle) angle1=angle;
			angle=fabs(interDir2.pnt[0]*interDir3.pnt[0]+interDir2.pnt[1]*interDir3.pnt[1]+interDir2.pnt[2]*interDir3.pnt[2]);
			if(angle1<angle) angle1=angle;
			angle=fabs(interDir3.pnt[0]*interDir1.pnt[0]+interDir3.pnt[1]*interDir1.pnt[1]+interDir3.pnt[2]*interDir1.pnt[2]);
			if(angle1<angle) angle1=angle;

			for(k=0; k<3; k++)
				interDir1.pnt[k]=(c2[c2Next].pnt[k]-c1[c1Dex].pnt[k]);
			dot=sqrt(interDir1.pnt[0]*interDir1.pnt[0]+interDir1.pnt[1]*interDir1.pnt[1]+interDir1.pnt[2]*interDir1.pnt[2]);
			for(k=0; k<3; k++)
				interDir1.pnt[k]=interDir1.pnt[k]/dot;
			for(k=0; k<3; k++)
				interDir2.pnt[k]=(c2[c2Dex].pnt[k]-c2[c2Next].pnt[k]);
			dot=sqrt(interDir2.pnt[0]*interDir2.pnt[0]+interDir2.pnt[1]*interDir2.pnt[1]+interDir2.pnt[2]*interDir2.pnt[2]);
			for(k=0; k<3; k++)
				interDir2.pnt[k]=interDir2.pnt[k]/dot;
			for(k=0; k<3; k++)
				interDir3.pnt[k]=(c1[c1Dex].pnt[k]-c2[c2Dex].pnt[k]);
			dot=sqrt(interDir3.pnt[0]*interDir3.pnt[0]+interDir3.pnt[1]*interDir3.pnt[1]+interDir3.pnt[2]*interDir3.pnt[2]);
			for(k=0; k<3; k++)
				interDir3.pnt[k]=interDir3.pnt[k]/dot;
			angle2=fabs(interDir1.pnt[0]*interDir2.pnt[0]+interDir1.pnt[1]*interDir2.pnt[1]+interDir1.pnt[2]*interDir2.pnt[2]);
			angle=fabs(interDir2.pnt[0]*interDir3.pnt[0]+interDir2.pnt[1]*interDir3.pnt[1]+interDir2.pnt[2]*interDir3.pnt[2]);
			// 2003/06/01
			//if(angle<angle) angle2=angle;
			if(angle2<angle) angle2=angle;
			angle=fabs(interDir3.pnt[0]*interDir1.pnt[0]+interDir3.pnt[1]*interDir1.pnt[1]+interDir3.pnt[2]*interDir1.pnt[2]);
			if(angle2<angle) angle2=angle;
			////////////////////////////
			for(k=0; k<3; k++)
				interDir2.pnt[k]=(c1[c1Next].pnt[k]-c2[c2Next].pnt[k]);
			dot=sqrt(interDir2.pnt[0]*interDir2.pnt[0]+interDir2.pnt[1]*interDir2.pnt[1]+interDir2.pnt[2]*interDir2.pnt[2]);
			for(k=0; k<3; k++)
				interDir2.pnt[k]=interDir2.pnt[k]/dot;
			for(k=0; k<3; k++)
				interDir3.pnt[k]=(c1[c1Dex].pnt[k]-c1[c1Next].pnt[k]);
			dot=sqrt(interDir3.pnt[0]*interDir3.pnt[0]+interDir3.pnt[1]*interDir3.pnt[1]+interDir3.pnt[2]*interDir3.pnt[2]);
			for(k=0; k<3; k++)
				interDir3.pnt[k]=interDir3.pnt[k]/dot;
			angle=fabs(interDir1.pnt[0]*interDir2.pnt[0]+interDir1.pnt[1]*interDir2.pnt[1]+interDir1.pnt[2]*interDir2.pnt[2]);
			if(angle2<angle) angle2=angle;
			angle=fabs(interDir2.pnt[0]*interDir3.pnt[0]+interDir2.pnt[1]*interDir3.pnt[1]+interDir2.pnt[2]*interDir3.pnt[2]);
			if(angle2<angle) angle2=angle;
			angle=fabs(interDir3.pnt[0]*interDir1.pnt[0]+interDir3.pnt[1]*interDir1.pnt[1]+interDir3.pnt[2]*interDir1.pnt[2]);
			if(angle2<angle) angle2=angle;

			if(angle1<angle2)
				secondThoughtWithC1=true;
			else
				secondThoughtWithC1=false;
			}

			if(angleC1>angleC2)
			{
				if(angleC1<TRIANGULATION_CRETERION)
					goWithC1=secondThoughtWithC1;
				// this means both triangulation is invalid, need to pick the relatively better one
				// with the bigger minimal angle
				else
					goWithC1=true;
			}
			else
			{
				if(angleC2<TRIANGULATION_CRETERION)
					goWithC1=secondThoughtWithC1;
				// this means both triangulation is invalid, need to pick the relatively better one
				// with the bigger minimal angle
				else
					goWithC1=false;
			}
			if(goWithC1)
			{
				// this validity check might still be necessary for a real robust implementation
				// check the validity of the segment of Inter(dex)->Main(next) by Inter(dex)->Inter(next)
				// with Inter(dex)->Inter(prev)
				// pick m[mNext] and i[iDex]

				//if(mainPicked<4)
				{
					goWithC1=true;
					c1Picked++;
				}
			}
			else
			{
				// pick m[mDex] and i[iNext]
				//if(interPicked<4)
				{
					goWithC1=false;
					c2Picked++;
				}
			}
		}
		else
			if(c1Done)
			{
				// pick m[mNext] and i[iDex]
				goWithC1=false;
			}
			else
				if(c2Done)
				{
					// pick m[mNext] and i[iDex]
					goWithC1=true;
				}

		if(goWithC1)
		{
			// pick c1[mNext] and c2[iDex]
//			fList[fOffsetT+fOffset5].Create(3);
//			fList[fOffsetT+fOffset5][0]=reDex[mainInterCurveDex[mainDex]];
//			fList[fOffsetT+fOffset5][1]=vOffset+interDex%interPNum;//reDex[subVertOldDex[0]];
//			fList[fOffsetT+fOffset5][2]=reDex[mainInterCurveDex[mainNext]];
//			fTag[fOffsetT+fOffset5]=3;
//			fOffset5++;

			triList[*triNum]=1;
			(*triNum)++;

			c1Dex=c1Next;
			if(c1Dex==0)
				c1Done=true;
		}
		else
		{
			// pick c1[mDex] and c2[iNext]
//			fList[fOffsetT+fOffset5].Create(3);
//			fList[fOffsetT+fOffset5][0]=reDex[mainInterCurveDex[mainDex]];
//			fList[fOffsetT+fOffset5][1]=vOffset+interDex%interPNum;//reDex[subVertOldDex[0]];
//			fList[fOffsetT+fOffset5][2]=vOffset+interNext%interPNum;//reDex[subVertOldDex[0]];
//			fTag[fOffsetT+fOffset5]=3;
//			fOffset5++;

			triList[*triNum]=2;
			(*triNum)++;

			c2Dex=c2Next;
			if(c2Dex==(*c2StartPos))
				c2Done=true;
		}
	}
}

void Pointlist_serverB::IntersectMeshes(MyFloat uCutOff, int nTiles1, Bpoint *tile1, int nTiles2, \
	Bpoint *tile2, int *interPPair, Bpoint2 *interPList, int *polyDexList1, int *polyDexList2, int *intPNum)
{
	int i, j, coPlanar;
	MyFloat minU, minU2;				// 
	MyFloat maxU, maxU2;
	Bpoint2 interP1, interP2;			// 
	int firstLoop=1, interPNum=0;

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	MyFloat minV, minV2;
#endif

if(uCutOff>=0)
{
	for(i=0; i<4*nTiles1; i+=4)
	{
		for(j=0; j<4*nTiles2; j+=4)
		{
			if(tile2[j].u<tile2[j+1].u)
				minU=tile2[j].u;
			else
				minU=tile2[j+1].u;
			if(tile2[j+2].u<tile2[j+3].u)
				minU2=tile2[j+2].u;
			else
				minU2=tile2[j+3].u;
			if(minU>minU2)
				minU=minU2;

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
			if(tile2[j].v<tile2[j+1].v)
				minV=tile2[j].v;
			else
				minV=tile2[j+1].v;
			if(tile2[j+2].v<tile2[j+3].v)
				minV2=tile2[j+2].v;
			else
				minV2=tile2[j+3].v;
			if(minV>minV2)
				minV=minV2;
			if(firstLoop) 
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fprintf(fDebug, "Poly %d: u/v, %f/%f\n", j/4, minU, minV);
#else
				fprintf(stderr, "Poly %d: u/v, %f/%f\n", j/4, minU, minV);
#endif

#endif

			// so far, polyReDex2 has the index of the sub-figure quad corresponding to the intersected 
			// segments below and this quad will be cut, therefore all the quads connected to this with 
			// smaller u should be redirected to this quad

			if((minU<uCutOff) &&
				//tri_tri_intersect_with_isectline_(tile1[i].pnt, tile1[i+1].pnt, tile1[i+2].pnt, 
				//								 tile2[j].pnt, tile2[j+1].pnt, tile2[j+2].pnt, 
				//								 &coPlanar, interP1, interP2))
				tri_tri_intersect_with_isectline_w_UVTRad(tile1[i], tile1[i+1], tile1[i+2], \
												tile2[j], tile2[j+1], tile2[j+2], \
												&coPlanar, &interP1, &interP2))
			{
				if(!polyDexList1[i/4])
					polyDexList1[i/4]=1;
				if(!polyDexList2[j/4])
					polyDexList2[j/4]=1;
				interPPair[2*interPNum]=i/4;
				//interPPair[2*interPNum+1]=index1/4;
				interPPair[2*interPNum+1]=j/4;

				interPList[2*interPNum]=interP1;
				interPList[2*interPNum+1]=interP2;

				interPNum++;
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fprintf(fDebug, "Intersection between tile %d and %d\n", i/4, j/4);
	#else
				fprintf(stderr, "Intersection between tile %d and %d\n", i/4, j/4);
	#endif
#endif
			}
			if((minU<uCutOff) && 
				//tri_tri_intersect_with_isectline(tile1[i].pnt, tile1[i+2].pnt, tile1[i+3].pnt, 
				//								 tile2[j].pnt, tile2[j+2].pnt, tile2[j+3].pnt, 
				//								 &coPlanar, interP1, interP2))
				tri_tri_intersect_with_isectline_w_UVTRad(tile1[i], tile1[i+2], tile1[i+3], \
												tile2[j], tile2[j+2], tile2[j+3], \
												&coPlanar, &interP1, &interP2))
			{
				if(!polyDexList1[i/4])
					polyDexList1[i/4]=1;
				if(!polyDexList2[j/4])
					polyDexList2[j/4]=1;
				interPPair[2*interPNum]=i/4;
				//interPPair[2*interPNum+1]=index1/4;
				interPPair[2*interPNum+1]=j/4;

				interPList[2*interPNum]=interP1;
				interPList[2*interPNum+1]=interP2;

				interPNum++;
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fprintf(fDebug, "Intersection between tile %d and %d\n", i/4, j/4);
	#else
				fprintf(stderr, "Intersection between tile %d and %d\n", i/4, j/4);
	#endif
#endif
			}
			if((minU<uCutOff) && 
				//tri_tri_intersect_with_isectline(tile1[i].pnt, tile1[i+1].pnt, tile1[i+2].pnt, 
				//								 tile2[j].pnt, tile2[j+2].pnt, tile2[j+3].pnt, 
				//								 &coPlanar, interP1, interP2))
				tri_tri_intersect_with_isectline_w_UVTRad(tile1[i], tile1[i+1], tile1[i+2], \
												tile2[j], tile2[j+2], tile2[j+3], \
												&coPlanar, &interP1, &interP2))
			{
				if(!polyDexList1[i/4])
					polyDexList1[i/4]=1;
				if(!polyDexList2[j/4])
					polyDexList2[j/4]=1;
				interPPair[2*interPNum]=i/4;
				//interPPair[2*interPNum+1]=index1/4;
				interPPair[2*interPNum+1]=j/4;

				interPList[2*interPNum]=interP1;
				interPList[2*interPNum+1]=interP2;

				interPNum++;
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fprintf(fDebug, "Intersection between tile %d and %d\n", i/4, j/4);
	#else
				fprintf(stderr, "Intersection between tile %d and %d\n", i/4, j/4);
	#endif
#endif
			}
			if((minU<uCutOff) && 
				//tri_tri_intersect_with_isectline(tile1[i].pnt, tile1[i+2].pnt, tile1[i+3].pnt, 
				//						 tile2[j].pnt, tile2[j+1].pnt, tile2[j+2].pnt, 
				//						 &coPlanar, interP1, interP2))
				tri_tri_intersect_with_isectline_w_UVTRad(tile1[i], tile1[i+2], tile1[i+3], \
												tile2[j], tile2[j+1], tile2[j+2], \
												&coPlanar, &interP1, &interP2))
			{
				if(!polyDexList1[i/4])
					polyDexList1[i/4]=1;
				if(!polyDexList2[j/4])
					polyDexList2[j/4]=1;
				interPPair[2*interPNum]=i/4;
				//interPPair[2*interPNum+1]=index1/4;
				interPPair[2*interPNum+1]=j/4;

				interPList[2*interPNum]=interP1;
				interPList[2*interPNum+1]=interP2;

				interPNum++;
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fprintf(fDebug, "Intersection between tile %d and %d\n", i/4, j/4);
	#else
				fprintf(stderr, "Intersection between tile %d and %d\n", i/4, j/4);
	#endif
#endif
			}
		}
		if(firstLoop) firstLoop=0;
	}
	*intPNum=interPNum;
}
else
{
	uCutOff=-uCutOff;
	for(i=0; i<4*nTiles1; i+=4)
	{
		for(j=0; j<4*nTiles2; j+=4)
		{
			if(tile2[j].u<tile2[j+1].u)
				maxU=tile2[j].u;
			else
				maxU=tile2[j+1].u;
			if(tile2[j+2].u<tile2[j+3].u)
				maxU2=tile2[j+2].u;
			else
				maxU2=tile2[j+3].u;
			if(maxU>maxU2)
				maxU=maxU2;
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
			int minV, minV2;
			if(tile2[j].v<tile2[j+1].v)
				minV=tile2[j].v;
			else
				minV=tile2[j+1].v;
			if(tile2[j+2].v<tile2[j+3].v)
				minV2=tile2[j+2].v;
			else
				minV2=tile2[j+3].v;
			if(minV>minV2)
				minV=minV2;
			if(firstLoop) 
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fprintf(fDebug, "Poly %d: u/v, %f/%f\n", j/4, maxU, minV);
#else
				fprintf(stderr, "Poly %d: u/v, %f/%f\n", j/4, maxU, minV);
#endif

#endif

			// so far, polyReDex2 has the index of the sub-figure quad corresponding to the intersected 
			// segments below and this quad will be cut, therefore all the quads connected to this with 
			// smaller u should be redirected to this quad

			if((maxU>uCutOff) &&
				//tri_tri_intersect_with_isectline_(tile1[i].pnt, tile1[i+1].pnt, tile1[i+2].pnt, 
				//								 tile2[j].pnt, tile2[j+1].pnt, tile2[j+2].pnt, 
				//								 &coPlanar, interP1, interP2))
				tri_tri_intersect_with_isectline_w_UVTRad(tile1[i], tile1[i+1], tile1[i+2], \
												tile2[j], tile2[j+1], tile2[j+2], \
												&coPlanar, &interP1, &interP2))
			{
				if(!polyDexList1[i/4])
					polyDexList1[i/4]=1;
				if(!polyDexList2[j/4])
					polyDexList2[j/4]=1;
				interPPair[2*interPNum]=i/4;
				//interPPair[2*interPNum+1]=index1/4;
				interPPair[2*interPNum+1]=j/4;

				interPList[2*interPNum]=interP1;
				interPList[2*interPNum+1]=interP2;
				interPNum++;
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fprintf(fDebug, "Intersection between tile %d and %d\n", i/4, j/4);
	#else
				fprintf(stderr, "Intersection between tile %d and %d\n", i/4, j/4);
	#endif
#endif
			}
			if((maxU>uCutOff) &&
				//tri_tri_intersect_with_isectline(tile1[i].pnt, tile1[i+2].pnt, tile1[i+3].pnt, 
				//								 tile2[j].pnt, tile2[j+2].pnt, tile2[j+3].pnt, 
				//								 &coPlanar, interP1, interP2))
				tri_tri_intersect_with_isectline_w_UVTRad(tile1[i], tile1[i+2], tile1[i+3], \
												tile2[j], tile2[j+2], tile2[j+3], \
												&coPlanar, &interP1, &interP2))
			{
				if(!polyDexList1[i/4])
					polyDexList1[i/4]=1;
				if(!polyDexList2[j/4])
					polyDexList2[j/4]=1;
				interPPair[2*interPNum]=i/4;
				//interPPair[2*interPNum+1]=index1/4;
				interPPair[2*interPNum+1]=j/4;

				interPList[2*interPNum]=interP1;
				interPList[2*interPNum+1]=interP2;
				interPNum++;
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fprintf(fDebug, "Intersection between tile %d and %d\n", i/4, j/4);
	#else
				fprintf(stderr, "Intersection between tile %d and %d\n", i/4, j/4);
	#endif
#endif
			}
			if((maxU>uCutOff) &&
				//tri_tri_intersect_with_isectline(tile1[i].pnt, tile1[i+1].pnt, tile1[i+2].pnt, 
				//								 tile2[j].pnt, tile2[j+2].pnt, tile2[j+3].pnt, 
				//								 &coPlanar, interP1, interP2))
				tri_tri_intersect_with_isectline_w_UVTRad(tile1[i], tile1[i+1], tile1[i+2], \
												tile2[j], tile2[j+2], tile2[j+3], \
												&coPlanar, &interP1, &interP2))
			{
				if(!polyDexList1[i/4])
					polyDexList1[i/4]=1;
				if(!polyDexList2[j/4])
					polyDexList2[j/4]=1;
				interPPair[2*interPNum]=i/4;
				//interPPair[2*interPNum+1]=index1/4;
				interPPair[2*interPNum+1]=j/4;

				interPList[2*interPNum]=interP1;
				interPList[2*interPNum+1]=interP2;
				interPNum++;
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fprintf(fDebug, "Intersection between tile %d and %d\n", i/4, j/4);
	#else
				fprintf(stderr, "Intersection between tile %d and %d\n", i/4, j/4);
	#endif
#endif
			}
			if((maxU>uCutOff) &&
				//tri_tri_intersect_with_isectline(tile1[i].pnt, tile1[i+2].pnt, tile1[i+3].pnt, 
				//						 tile2[j].pnt, tile2[j+1].pnt, tile2[j+2].pnt, 
				//						 &coPlanar, interP1, interP2))
				tri_tri_intersect_with_isectline_w_UVTRad(tile1[i], tile1[i+2], tile1[i+3], \
												tile2[j], tile2[j+1], tile2[j+2], \
												&coPlanar, &interP1, &interP2))
			{
				if(!polyDexList1[i/4])
					polyDexList1[i/4]=1;
				if(!polyDexList2[j/4])
					polyDexList2[j/4]=1;
				interPPair[2*interPNum]=i/4;
				interPPair[2*interPNum+1]=j/4;

				interPList[2*interPNum]=interP1;
				interPList[2*interPNum+1]=interP2;
				interPNum++;
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fprintf(fDebug, "Intersection between tile %d and %d\n", i/4, j/4);
	#else
				fprintf(stderr, "Intersection between tile %d and %d\n", i/4, j/4);
	#endif
#endif
			}
		}
		if(firstLoop) firstLoop=0;
	}
	*intPNum=interPNum;
}
}

bool Pointlist_serverB::FindIntersectionCurveFromSegments(Bpoint2 *interPList2, int *interPPair2, int *intPNum, \
											   Bpoint2 *interPList, int *interPPair, MyFloat EPSILON2)
{
	int i, j;
	int zeroLenSegmentCounter;
	int badSegmentCounter;
	int taken[MAX_INTERSECTION_NUMBER];
	int interPNum=*intPNum;

	bool successful;
	MyFloat EPSILON2_TEMP=EPSILON2;

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\n\n\ntrying to find the closed intersection curve: 1st trial from %d points,\n", interPNum);
	#else
	fprintf(stderr, "\n\n\ntrying to find the closed intersection curve: 1st trial from %d points,\n", interPNum);
	#endif
#endif

	successful=true;
	badSegmentCounter=0;
	zeroLenSegmentCounter=0;
	memset(taken, 0, sizeof(int)*MAX_INTERSECTION_NUMBER);

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "trying to find the segment which is actually a vertex\n");
	#endif
#endif
	for(i=0; i<interPNum; i++)
	{
		if(SamePoint(EPSILON2, interPList[2*i], interPList[2*i+1]))
		{
			taken[i]=1;
			zeroLenSegmentCounter++;
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			fprintf(fDebug, "seg %d: %f %f %f, %f %f %f\n", i, \
				interPList[2*i].pnt[0], interPList[2*i].pnt[1], interPList[2*i].pnt[2], \
				interPList[2*i+1].pnt[0], interPList[2*i+1].pnt[1], interPList[2*i+1].pnt[2]);
	#else
			fprintf(stderr, "seg %d is vertex-segment\n", i);
	#endif
#endif
		}
	}
	i=0;
	while(taken[i]) 
		i++;

	interPList2[0]=interPList[2*i];
	interPList2[1]=interPList[2*i+1];
	interPPair2[0]=interPPair[2*i];
	interPPair2[1]=interPPair[2*i+1];
	taken[i]=1;

	j=i;
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\n0th point found at seg %d: %f %f %f, %f %f %f\n", j, \
		interPList[2*j].pnt[0], interPList[2*j].pnt[1], interPList[2*j].pnt[2], \
	interPList[2*j+1].pnt[0], interPList[2*j+1].pnt[1], interPList[2*j+1].pnt[2]);
	#else
	fprintf(stderr, "\n0th point found at %dth segment\n", j);
	#endif
#endif

	// should start from exactly 1 instead of the index of the next segment
	// i++;
	i=1;

	j=0;
	// generate the whole intersection curve from the separated intersection segments
	//		all the un-ordereed intersection vertices in interPList now will be linked as a curve into interPList2
	while(i<interPNum-zeroLenSegmentCounter)
	{
		//if(lastTaken==j) lastTaken=-1;
		while(taken[j] && j<=interPNum) j++;

		if(j>interPNum)
		{
			//		this does not necessarily mean the curve is not closed!
			//		it could be the case that the segment actually has zero-length
			//		therefore all such segments will be tagged and counted before entering this cycle
			bool trap=true;

#ifdef QIONG_DEBUG_FILE_TO_WRITE
			fprintf(fDebug, "\n1st trial failed!\nneed to print out all the 'bad' segments:\n");
			//fprintf(fDebug, "\n\n\n%d\n", interPNum);
			for(int ii=0; ii<interPNum; ii++)
			{
				if(!taken[ii])
					fprintf(fDebug, "seg %d: %f %f %f, %f %f %f\n", ii, interPList[2*ii].pnt[0], interPList[2*ii].pnt[1], interPList[2*ii].pnt[2], \
					interPList[2*ii+1].pnt[0], interPList[2*ii+1].pnt[1], interPList[2*ii+1].pnt[2]);
			}
#endif

			// the curve could be CLOSED before TRAVERSING all the segments!
			if(SamePoint(EPSILON2_TEMP, interPList2[i], interPList2[0]))
			{
				for(i=0; i<interPNum; i++)
				{
					if(!taken[i])
						badSegmentCounter++;
				}

				// which means the curve IS ALREADY closed!
				break;
			}
			else
			{
				successful=false;
				break;
			}
			// if the algorithm fails, i need to check whether the curve itself has already been closed or not
			// if it is, need to update the actual length of the curve
			//return;
		}

		//if(lastTaken==-1) lastTaken=j;
		//if(interPList2[i].pnt[0]-interPList[2*j].pnt[0]>=-EPSILON2&&interPList2[i].pnt[0]-interPList[2*j].pnt[0]<=EPSILON2)
		//if(SamePoint(EPSILON2*10, interPList2[i], interPList[2*j]))
		if(SamePoint(EPSILON2_TEMP, interPList2[i], interPList[2*j]))
		{
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			fprintf(fDebug, "%dth point found at seg %d: %f %f %f, %f %f %f\n", i, j, interPList[2*j].pnt[0], interPList[2*j].pnt[1], interPList[2*j].pnt[2], \
			interPList[2*j+1].pnt[0], interPList[2*j+1].pnt[1], interPList[2*j+1].pnt[2]);
	#else
			fprintf(stderr, "%dth point found at %dth segment\n", i, j);
	#endif
#endif

			// need to check whether the segment j and the last taken segment are ACTUALLY THE SAME!
			if(i>0 && SamePoint(EPSILON2_TEMP, interPList2[i-1], interPList[2*j+1]))
			{
				// which means the segment j and the last taken segment are ACTUALLY THE SAME!
				badSegmentCounter++;
				taken[j]=1;
			}
			else
			{
			taken[j]=1;
			interPPair2[2*i]=interPPair[2*j];
			interPPair2[2*i+1]=interPPair[2*j+1];
			i++;
			interPList2[i]=interPList[2*j+1];
			//for(int k=0; k<3; k++)
			//{
			//	interPList2[i].pnt[k]=interPList[2*j+1].pnt[k];
			//}
			j=0;
			}
		}
		else
			//if(interPList2[i].pnt[0]-interPList[2*j+1].pnt[0]>=-EPSILON2&&interPList2[i].pnt[0]-interPList[2*j+1].pnt[0]<=EPSILON2)
			//if(SamePoint(EPSILON2*10, interPList2[i], interPList[2*j+1]))
			if(SamePoint(EPSILON2_TEMP, interPList2[i], interPList[2*j+1]))
			{
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			fprintf(fDebug, "%dth point found at seg %d: %f %f %f, %f %f %f\n", i, j, interPList[2*j].pnt[0], interPList[2*j].pnt[1], interPList[2*j].pnt[2], \
			interPList[2*j+1].pnt[0], interPList[2*j+1].pnt[1], interPList[2*j+1].pnt[2]);
	#else
			fprintf(stderr, "%dth point found at %dth segment\n", i, j);
	#endif
#endif

				// need to check whether the segment j and the last taken segment are ACTUALLY THE SAME!
				if(i>0 && SamePoint(EPSILON2_TEMP, interPList2[i-1], interPList[2*j]))
				{
					// which means the segment j and the last taken segment are ACTUALLY THE SAME!
					badSegmentCounter++;
					taken[j]=1;
				}
				else
				{
				taken[j]=1;
				interPPair2[2*i]=interPPair[2*j];
				interPPair2[2*i+1]=interPPair[2*j+1];
				i++;

				interPList2[i]=interPList[2*j];
				//for(int k=0; k<3; k++)
				//{
				//	interPList2[i].pnt[k]=interPList[2*j].pnt[k];
				//}
				j=0;
				}
			}
			else j++;
	}

	EPSILON2_TEMP=4*EPSILON2_TEMP*EPSILON2_TEMP;

if(!successful)
{
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\n\n\ntrying to find the closed intersection curve: 2nd trial from %d points,\n", interPNum);
	//fprintf(fDebug, "0th point found at seg %d: %f %f %f, %f %f %f\n", j, interPList[2*j].pnt[0], interPList[2*j].pnt[1], interPList[2*j].pnt[2], \
	//interPList[2*j+1].pnt[0], interPList[2*j+1].pnt[1], interPList[2*j+1].pnt[2]);
	#else
	fprintf(stderr, "\n\n\ntrying to find the closed intersection curve: 2nd trial from %d points,\n", interPNum);
	//fprintf(stderr, "%dth point found at %dth segment\n", i, j);
	#endif
#endif

	successful=true;
	zeroLenSegmentCounter=0;
	badSegmentCounter=0;
	memset(taken, 0, sizeof(int)*MAX_INTERSECTION_NUMBER);

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "trying to find the segment which is actually a vertex\n");
	#endif
#endif

	for(i=0; i<interPNum; i++)
	{
		if(SamePoint(EPSILON2, interPList[2*i], interPList[2*i+1]))
		{
			taken[i]=1;
			zeroLenSegmentCounter++;

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			fprintf(fDebug, "seg %d: %f %f %f, %f %f %f\n", i, interPList[2*i].pnt[0], interPList[2*i].pnt[1], interPList[2*i].pnt[2], \
			interPList[2*i+1].pnt[0], interPList[2*i+1].pnt[1], interPList[2*i+1].pnt[2]);
	#else
			fprintf(stderr, "seg %d is vertex-segment\n", i);
	#endif
#endif
		}
	}
	i=0;
	while(taken[i]) 
		i++;

	interPList2[0]=interPList[2*i];
	interPList2[1]=interPList[2*i+1];
	interPPair2[0]=interPPair[2*i];
	interPPair2[1]=interPPair[2*i+1];
	taken[i]=1;

	j=i;

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	//fprintf(fDebug, "\n\n\ntrying to find the closed intersection curve: 2nd trial from %d points,\n", interPNum);
	fprintf(fDebug, "\n0th point found at seg %d: %f %f %f, %f %f %f\n", j, interPList[2*j].pnt[0], interPList[2*j].pnt[1], interPList[2*j].pnt[2], \
	interPList[2*j+1].pnt[0], interPList[2*j+1].pnt[1], interPList[2*j+1].pnt[2]);
	#else
	//fprintf(stderr, "\n\n\ntrying to find the closed intersection curve: 2nd trial from %d points,\n", interPNum);
	fprintf(stderr, "\n0th point found at %dth segment\n", j);
	#endif
#endif

	// should start from exactly 1 instead of the index of the next segment
	// i++;
	i=1;

	j=0;
	// generate the whole intersection curve from the separated intersection segments
	//		all the un-ordereed intersection vertices in interPList now will be linked as a curve 
	//		into interPList2
	while(i<interPNum-zeroLenSegmentCounter)
	{
		//if(lastTaken==j) lastTaken=-1;
		while(taken[j] && j<=interPNum) j++;

		if(j>interPNum)
		{
			//		this does not necessarily mean the curve is not closed!
			//		it could be the case that the segment actually has zero-length
			//		therefore all such segments will be tagged and counted before entering this cycle
			bool trap=true;

#ifdef QIONG_DEBUG_FILE_TO_WRITE
			fprintf(fDebug, "\n2nd trial failed!\nneed to print out all the 'bad' segments:\n");
			//fprintf(fDebug, "\n\n\n%d\n", interPNum);
			for(int ii=0; ii<interPNum; ii++)
			{
				if(!taken[ii])
					fprintf(fDebug, "seg %d: %f %f %f, %f %f %f\n", ii, interPList[2*ii].pnt[0], interPList[2*ii].pnt[1], interPList[2*ii].pnt[2], \
					interPList[2*ii+1].pnt[0], interPList[2*ii+1].pnt[1], interPList[2*ii+1].pnt[2]);
			}
#endif

			// the curve could be CLOSED before TRAVERSING all the segments!
			if(SamePoint(EPSILON2, interPList2[i], interPList2[0]))
			{
				for(i=0; i<interPNum; i++)
				{
					if(!taken[i])
						badSegmentCounter++;
				}

				// which means the curve IS ALREADY closed!
				break;
			}
			else
			{
				successful=false;
				break;
			}
			// if the algorithm fails, i need to check whether the curve itself has already been closed or not
			// if it is, need to update the actual length of the curve
			//return;
		}

		//if(lastTaken==-1) lastTaken=j;
		//if(interPList2[i].pnt[0]-interPList[2*j].pnt[0]>=-EPSILON2&&interPList2[i].pnt[0]-interPList[2*j].pnt[0]<=EPSILON2)
		//if(SamePoint(EPSILON2*10, interPList2[i], interPList[2*j]))
		if(SamePointByDistanceSquare(EPSILON2_TEMP, interPList2[i], interPList[2*j]))
		{
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			fprintf(fDebug, "%dth point found at seg %d: %f %f %f, %f %f %f\n", i, j, interPList[2*j].pnt[0], interPList[2*j].pnt[1], interPList[2*j].pnt[2], \
			interPList[2*j+1].pnt[0], interPList[2*j+1].pnt[1], interPList[2*j+1].pnt[2]);
	#else
			fprintf(stderr, "%dth point found at %dth segment\n", i, j);
	#endif
#endif

				// very tricky bug here (maybe not jus a BUG)
				// need to check whether the segment j and the last taken segment are ACTUALLY THE SAME!
				if(i>0 && SamePointByDistanceSquare(EPSILON2_TEMP, interPList2[i-1], interPList[2*j+1]))
				{
					// which means the segment j and the last taken segment are ACTUALLY THE SAME!
					badSegmentCounter++;
					taken[j]=1;
				}
				else
				{
					taken[j]=1;
					interPPair2[2*i]=interPPair[2*j];
					interPPair2[2*i+1]=interPPair[2*j+1];
					i++;

					interPList2[i]=interPList[2*j+1];
					//for(int k=0; k<3; k++)
					//{
					//	interPList2[i].pnt[k]=interPList[2*j+1].pnt[k];
					//}
					j=0;
				}
		}
		else
			//if(interPList2[i].pnt[0]-interPList[2*j+1].pnt[0]>=-EPSILON2&&interPList2[i].pnt[0]-interPList[2*j+1].pnt[0]<=EPSILON2)
			//if(SamePoint(EPSILON2*10, interPList2[i], interPList[2*j+1]))
			if(SamePointByDistanceSquare(EPSILON2_TEMP, interPList2[i], interPList[2*j+1]))
			{
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			fprintf(fDebug, "%dth point found at seg %d: %f %f %f, %f %f %f\n", i, j, interPList[2*j].pnt[0], interPList[2*j].pnt[1], interPList[2*j].pnt[2], \
			interPList[2*j+1].pnt[0], interPList[2*j+1].pnt[1], interPList[2*j+1].pnt[2]);
	#else
			fprintf(stderr, "%dth point found at %dth segment\n", i, j);
	#endif
#endif

				// need to check whether the segment j and the last taken segment are ACTUALLY THE SAME!
				if(i>0 && SamePointByDistanceSquare(EPSILON2_TEMP, interPList2[i-1], interPList[2*j]))
				{
					// which means the segment j and the last taken segment are ACTUALLY THE SAME!
					badSegmentCounter++;
					taken[j]=1;
				}
				else
				{
				taken[j]=1;
				interPPair2[2*i]=interPPair[2*j];
				interPPair2[2*i+1]=interPPair[2*j+1];
				i++;

				interPList2[i]=interPList[2*j];
				//for(int k=0; k<3; k++)
				//{
				//	interPList2[i].pnt[k]=interPList[2*j].pnt[k];
				//}
				j=0;
				}
			}
			else j++;
	}
}

	//interPNum-=zeroLenSegmentCounter;
	interPNum-=zeroLenSegmentCounter+badSegmentCounter;
	*intPNum=interPNum;
	return successful;
}

bool Pointlist_serverB::FindCurveFromSegments()
{
	return true;
}

inline double CurveLength(Bpoint2 *curv, int pNum)
{
	int localDex1,
		localDex2;
	double	len=0,
			dist;
	for(localDex1=0; localDex1<pNum; localDex1++)
	{
		dist=0;
		for(localDex2=0; localDex2<3; localDex2++)
			dist+=(curv[localDex1+1].pnt[localDex2]-curv[localDex1].pnt[localDex2])*
				  (curv[localDex1+1].pnt[localDex2]-curv[localDex1].pnt[localDex2]);
		len+=sqrt(dist);
	}
	return len;
}

bool Pointlist_serverB::ComputeBlending(int treeId)
{
#ifndef PARAMETERS_FOR_TORUS_OBJECT
//	if(rMode==WIREFRAME)
//		glLineWidth(1.5);

	MyFloat SUB_FIGURE_CUTOFF_U=SUB_FIGURE_CUTOFF_U_ORIGINAL+SUB_FIGURE_CUTOFF_U_RANGE*blendExtents[treeId];
	MyFloat DILATION_STEPS=DILATION_STEPS_ORIGINAL+DILATION_STEPS_RANGE*blendAmounts[treeId];

	// blending amount testing parameters
	//MyFloat SUB_FIGURE_CUTOFF_U=SUB_FIGURE_CUTOFF_U_ORIGINAL+SUB_FIGURE_CUTOFF_U_RANGE*subBlending;
	//MyFloat DILATION_STEPS=DILATION_STEPS_ORIGINAL+DILATION_STEPS_RANGE*mainBlending;

	// Since the parameters controling the blending amount has been implemented,
	// I will use them instead. Their names should be changed accordingly later.
	//MyFloat SUB_FIGURE_CUTOFF_U=SUB_FIGURE_CUTOFF_U_ORIGINAL+SUB_FIGURE_CUTOFF_U_RANGE*(blendExtent==0?.01:blendExtent);
	//MyFloat DILATION_STEPS=DILATION_STEPS_ORIGINAL+DILATION_STEPS_RANGE*(blendAmount==0?.01:blendAmount);
	MyFloat EPSILON2=EPSILON2_ORIGINAL;

	int i, j, k;
#ifndef QIONG_BLENDING_NO_DEBUG
	MyFloat color[3];
#endif

	int interPNum;
    //int ret;
	//int taken[MAX_INTERSECTION_NUMBER];
	int interPPair[MAX_INTERSECTION_NUMBER*2];
	int interPPair2[MAX_INTERSECTION_NUMBER*2];
	Bpoint2 interPList[MAX_INTERSECTION_NUMBER];
	Bpoint2 interPList2[MAX_INTERSECTION_NUMBER];
	Bpoint2 interCurve[MAX_INTERSECTION_NUMBER];
	Bpoint2 interCurveD[MAX_INTERSECTION_NUMBER];

	int *polyDexList1, *polyDexList2;
	//int *polyDexList1=new int[nTiles1], *polyDexList2=new int[nTiles2];
	//int polyDexList1[MAX_INTERSECTION_NUMBER];
	//int polyDexList2[MAX_INTERSECTION_NUMBER];

	int nTiles1, nTiles2;
	Bpoint *tile1, *tile2;
    list<M3DFigureTreeNode *> children;
	ThallCode::Pointlist_server2 *pListOld[2];

#ifdef DEBUG
	cout << "Pointlist_serverB::renderFigure() called with treeId = " << treeId << '\n';
#endif

#ifdef QIONG_DEBUG_FILE_TO_WRITE
	//FILE *fDebug;
	char fName[200];
	sprintf(fName, "%s", "debug.dat");
	fDebug=fopen(fName, "w+");
#endif

	interPNum=0;						// number of the vertices on the intersection curve
	//int coPlanar;						// ??

	//MyFloat interP1[3], interP2[3];	// ??
	//Bpoint2 interP1, interP2;
	//MyFloat interUVT[3];
	MyFloat minU; //, minU2;				// ??
	//MyFloat minV, minV2;

//	currentFigure = treeId;				// ??
    //M3DFigureTreeNode *root=object->getFigureTreeRoot(treeId);

	/*
	pList1=tileLists[0];				// Start with the main figure to fill the tile information
	pList2=tileLists[1];				// Need the sub-figure to calculate the intersection
 	pList1->subdivtileinfo(&nTiles1, &tile1);
	pList2->subdivtileinfo(&nTiles2, &tile2);
	*/
	pListOld[0]=pListOriginal[0]; //tileLists[0];				// Start with the main figure to fill the tile information
	pListOld[1]=pListOriginal[1]; //tileLists[1];				// Need the sub-figure to calculate the intersection
 	pListOld[0]->subdivtileinfo(&nTiles1, &tile1);
	pListOld[1]->subdivtileinfo(&nTiles2, &tile2);

	polyDexList1=new int[nTiles1];
	polyDexList2=new int[nTiles2];
	memset(polyDexList1, 0, nTiles1*sizeof(int));
	memset(polyDexList2, 0, nTiles2*sizeof(int));
	//memset(polyDexList1, 0, MAX_INTERSECTION_NUMBER*sizeof(int));
	//memset(polyDexList2, 0, MAX_INTERSECTION_NUMBER*sizeof(int));

#ifndef QIONG_BLENDING_NO_DEBUG
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nThe re-indexing information of level %d sub-figure \
		quads are removed.\n", SUBDIV_LEVEL_TO_USE);
	#else
	fprintf(stderr, "\nThe re-indexing information of level %d sub-figure \
		quads are removed.\n", SUBDIV_LEVEL_TO_USE);
	#endif
#endif

// .1 To intersect the two meshes here!
	IntersectMeshes(SUB_FIGURE_CUTOFF_U, nTiles1, tile1, nTiles2, tile2, interPPair, interPList, \
		polyDexList1, polyDexList2, &interPNum);

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	for(i=0; i<interPNum; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "Intersection point %d: u/v/t in main figure %f/%f/%f\n", \
			i, interPList[i].u, interPList[i].v, interPList[i].t);
	#else
		fprintf(stderr, "Intersection point %d: u/v/t in main figure %f/%f/%f\n", \
			i, interPList[i].u, interPList[i].v, interPList[i].t);
	#endif
	}
#endif

//#ifndef QIONG_BLENDING_NO_DEBUG
	if(interPNum>MAX_INTERSECTION_NUMBER || interPNum<=0)
	{
		fprintf(stderr, "TOO MANY or NONE INTERSECTION POINTS found from the TRI-TRI intersections...\n");

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		return false;
	}
//#endif

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	color[0]=0;
	color[1]=1;
	color[2]=0;
	glRenderSegments(0, interPNum-1, interPList, color);
#endif

#ifdef QIONG_BLENDING_DILATION_DEBUG
	for(i=0; i<interPNum; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "seg %d: %f %f %f, %f %f %f\n", i, \
			interPList[2*i].pnt[0], interPList[2*i].pnt[1], interPList[2*i].pnt[2], \
			interPList[2*i+1].pnt[0], interPList[2*i+1].pnt[1], interPList[2*i+1].pnt[2]);
	#else
		fprintf(stderr, "seg %d: %f %f %f, %f %f %f\n", i, \
			interPList[2*i].pnt[0], interPList[2*i].pnt[1], interPList[2*i].pnt[2], \
			interPList[2*i+1].pnt[0], interPList[2*i+1].pnt[1], interPList[2*i+1].pnt[2]);
	#endif
	}
#endif

	// at this stage, all the intersection vertices are found!

// .2 Need to find the closed intersection curve!
	bool successful;
	successful=FindIntersectionCurveFromSegments(interPList2, interPPair2, &interPNum, \
		interPList, interPPair, EPSILON2);

	if(!successful)
	{
		fprintf(stderr, "\nThe intersection curve is not closed!!!\nHave to quit right here...\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fclose(fDebug);
#endif

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		return false;
	}


#ifdef QIONG_BLENDING_INTERSECTION_QUAD
	color[0]=1;
	color[1]=0;
	color[2]=0;
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "The list of main-figure quads that was intersected by sub-figure (1)\n \
		or intersected by the dilated i-curve (2)...\n");
	#endif
if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
{
	for(i=0; i<nTiles1; i++)
	{
		if(polyDexList1[i])
		{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			if(i%8==0)
				fprintf(fDebug, "\n");
			fprintf(fDebug, "case %d:%d\t", polyDexList1[i], i);
	#endif
			glRender(SUBDIV_LEVEL_TO_USE+1, i, rMode, color);
		}
	}
}
else
{
	for(i=0; i<nTiles1; i++)
	{
		if(polyDexList1[i])
		{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			if(i%8==0)
				fprintf(fDebug, "\n");
			fprintf(fDebug, "case %d:%d\t", polyDexList1[i], i);
	#endif
			glRender(SUBDIV_LEVEL_TO_USE, i, rMode, color);
		}
	}
}
	for(i=0; i<nTiles2; i++)
	{
		if(polyDexList2[i])
			glRender(SUBDIV_LEVEL_TO_USE, i, rMode, color);
	}
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\n");
	#endif
#endif	

	//	make sure the intersection curve is closed somehow
	//if(!SamePoint(EPSILON2, interCurve[0], interCurve[interPNum]))
	if(!SamePointByDistanceSquare(EPSILON2, interPList2[0], interPList2[interPNum]))
	//if(!SamePointByDistanceSquare(EPSILON2, interCurve[0], interCurve[interPNum]))
	{
		fprintf(stderr, "\nThe intersection curve is not closed!!!\nHave to quit right here...\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fclose(fDebug);
#endif
		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		return false;
	}

#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nThe main-figure u/v/t of all the %d vertices on the original i-curve...\n", interPNum);
	for(i=0; i<=interPNum; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f/%f\t", i, interPList2[i].u, interPList2[i].v, interPList2[i].t);
	}
	fprintf(fDebug, "\nThe sub-figure u/v/t of all the %d vertices on the original i-curve...\n", interPNum);
	for(i=0; i<=interPNum; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f/%f\t", i, interPList2[i].u_2, interPList2[i].v_2, interPList2[i].t_2);
	}
	fprintf(fDebug, "\n");
#endif

	//		now the task is to trim all the unnecessary segments including: too short and colinear segments
	//		the results are stored in 'interCurve and interPPair'
	//		interPPair[2*i]   stores the information of which quad the current vertex i is on main-figure
	//		interPPair[2*i+1] stores the information of which quad the current segment [i, i+1] is on sub-figure
	// the minimal-length of the segment to keep should be dynamically determined by the sub-figure itself
	// the # of points left should be compatible to the # of points on the sub-figure CUT curve
	//int numPointsOnCutCurve=(int)(2*(pListOld[1]->max_v()+2)*pow(2, SUBDIV_LEVEL_TO_USE)*1.1);	// hard-coded 1.1 here
	RegularizeIntersectionCurve(numPointsOnCutCurve, interCurve, interPList2, interPPair, interPPair2, &interPNum);
	//RegularizeIntersectionCurve(interCurve, interPList2, interPPair, interPPair2, &interPNum);
	if(interPNum>MAX_INTERSECTION_NUMBER || interPNum<=0)
	{
		fprintf(stderr, "NONE INTERSECTION POINTS left after the RegularizeIntersectionCurve!\n");

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		return false;
	}

#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nThe v/t of all the %d vertices on the regularized i-curve...\n", interPNum);
	for(i=0; i<=interPNum; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f\t", i, interCurve[i].v_2, interCurve[i].t_2);
	}
	fprintf(fDebug, "\n");
#endif

	//		the numerical problem has to be dealt with somewhere
	//		after the interpolation of the barycentric coordinates of the intersection curve, those UVT.Ts which
	//		were supposed to be 1/-1 or in the crest the UVT.UVs to be 0, uMax, 0, vMax are not exactly these
	//		values, then we have to PUT THEM BACK to the values they are assumed to be!
	RegularizeIntersectionCurveUVT(pListOld[0]->max_u(), pListOld[0]->max_v(), \
		pListOld[1]->max_u(), pListOld[1]->max_v(),	interCurve, interPNum);

#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nThe v/t of all the %d vertices on the regularized UVT i-curve...\n", interPNum);
	for(i=0; i<=interPNum; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f\t", i, interCurve[i].v_2, interCurve[i].t_2);
	}
	fprintf(fDebug, "\n");
#endif

	//		make sure the intersection curve is closed somehow

	// all the UVTRad information should have been done when calculating the intersections!
	// Need to figure out the UVT coordinates for the intersection points on the main-figure
	//		the implementation in the Pointlist_server2 is an approximation
	//MyFloat interUVT[3];
	//for(i=0; i<interPNum; i++)
	//{
	//	pListOld[0]->SubdivBUVT(interUVT, 0, interPPair[2*i], interCurve[i].pnt);
		//CalculateUVT(tile1[i].pnt, tile1[i+1].pnt, tile1[i+2].pnt, tile1[i+3].pnt, interP1, interP1);
	//	interCurve[i].u=interUVT[0];
	//	interCurve[i].v=interUVT[1];
	//	interCurve[i].t=interUVT[2];
	//}
	//interCurve[interPNum].u=interCurve[0].u;
	//interCurve[interPNum].v=interCurve[0].v;
	//interCurve[interPNum].t=interCurve[0].t;
	//interCurve[interPNum]=interCurve[0];

	// this is not necessary any longer since interCurve[interPNum] is not used anyway!
	//interCurve[interPNum]=interCurve[0];

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_REGULARIZED
	color[0]=1;
	color[1]=0;
	color[2]=0;
	//glRenderCurve(0, interPNum, interCurve, color);
	glRenderCurve(0, interPNum, interCurve, color);
#endif


	// use the 'interCurve' as temporary storage, the final dilated intersection curve should still be
	// stored in interCurveD
	//memcpy(interCurveD, interCurve, (interPNum+5)*sizeof(Bpoint2));

#ifdef QIONG_DEBUG_20030510
	// try to check all the t coordinates of all the intersection vertices
	// they are supposed to be 1/-1 if uv belongs to (0, uMax) and (0, vMax) respectively
	for(i=0; i<=interPNum; i++)
		//fprintf(stderr, "%8.7f\t%d\t", interCurveD[i].t, interCurveD[i].t==-1.0);
		fprintf(stderr, "%d\t", interCurveD[i].t==-1.0);
	fprintf(stderr, "\n");
#endif

// time to re-mesh the entire surface boundary!
	CCMesh *thisMesh, *thisMesh1, *thisMesh2;
	CCVertex *thisVert;//, *thatVert;
	Bpoint tempVert;

	CCPolygon *thisPoly;
	int pSides;

	//thisMesh1=pListOld[0]->thisSurf->SubMeshPointer(SUBDIV_LEVEL_TO_USE);
	//thisMesh2=pListOld[1]->thisSurf->SubMeshPointer(SUBDIV_LEVEL_TO_USE);
	if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
		thisMesh1=(pListOld[0]->GetThisSurfPointer())->SubMeshPointer(SUBDIV_LEVEL_TO_USE+1);
	else
		thisMesh1=(pListOld[0]->GetThisSurfPointer())->SubMeshPointer(SUBDIV_LEVEL_TO_USE);

	thisMesh2=(pListOld[1]->GetThisSurfPointer())->SubMeshPointer(SUBDIV_LEVEL_TO_USE);

// 1.0 cut the subfigure!
	// for sub-figure, re-relabling vertices

	int vNumEstimate=thisMesh1->VertList.Length();
	int fNumEstimate=thisMesh1->FaceList.Length();
	vNumEstimate+=thisMesh2->VertList.Length();
	fNumEstimate+=thisMesh2->FaceList.Length();
	vNumEstimate*=1.5;
	fNumEstimate*=1.5;

	thisMesh=thisMesh2;
	int vNums=thisMesh->VertList.Length();
	int fNums=thisMesh->FaceList.Length();

	int vOffsetT=0, fOffsetT=0;
	int vOffset=0, fOffset=0;			// for remained and cut vertices and patches on sub-figure

	//int  fOffset2=0;					// for sewing part of sub-figure
	//int  fOffset3=0;					// for remained patches on main-figure

	int *reDex=new int[vNumEstimate];
	int reduction=0;
	Bpoint2 *vList=new Bpoint2[vNumEstimate];
	PindexList *fList=new PindexList[fNumEstimate];
	int *vTag=new int[vNumEstimate];
	int *fTag=new int[fNumEstimate];
	memset(vTag, -1, vNumEstimate*sizeof(int));
	memset(fTag, -1, fNumEstimate*sizeof(int));
	for(i=0; i<vNumEstimate; i++) 
		reDex[i]=i;

	// reindex the vertices for the sub-figure first
	int vDex;
	for(vDex=0; vDex<vNums; vDex++)
	{
		// reDex is the important array to re-index all the orginal sub-figure vertices to the new array in vList
		reDex[vDex]-=reduction;
		thisVert=thisMesh->VertList(vDex);
		//		the codes picking vertices to remove are not CORRECT!
		if(thisVert->uvt.u<SUB_FIGURE_CUTOFF_U)
		{
			//if(thisVert->uvt.u-SUB_FIGURE_CUTOFF_U<-1)// || (thisVert->uvt.u<SUB_FIGURE_CUTOFF_U&&thisVert->uvt.t==0))
			if( (thisVert->uvt.u-SUB_FIGURE_CUTOFF_U<=-UVT_UNIT) || 
				(thisVert->uvt.u==0 && thisVert->uvt.t>-1 && thisVert->uvt.t<1 && thisVert->uvt.v>0 && thisVert->uvt.v<pListOld[1]->max_v()) )
			    //(thisVert->uvt.u-SUB_FIGURE_CUTOFF_U>-1.0/(SUBDIV_LEVEL_TO_USE+1) && (thisVert->uvt.t>-1 && thisVert->uvt.t<1) && thisVert->uvt.v>0 && thisVert->uvt.v<pListOld[1]->max_v()) )
			{
				reDex[vDex]=-1;
				reduction++;
				//vOffset--;
			}
			else
			{
				// need to cut the sub-figure, generate new vertices here with new U
				/*
				for(int ii=0; ii<vNums; ii++)
				{
					thatVert=thisMesh->VertList(ii);
					if(thatVert->uvt.u>SUB_FIGURE_CUTOFF_U&&thatVert->uvt.u-SUB_FIGURE_CUTOFF_U<1&&
					   thatVert->uvt.v==thisVert->uvt.v   &&thatVert->uvt.t==thisVert->uvt.t)
							break;
				}
				*/
				//double t=SUB_FIGURE_CUTOFF_U-thisVert->uvt.u;
				tempVert.u=SUB_FIGURE_CUTOFF_U+SUB_FIGURE_CUTOFF_TRANSFORM;
				tempVert.v=thisVert->uvt.v;
				tempVert.t=thisVert->uvt.t;
				pListOld[1]->subdivBPosition(&tempVert, tempVert.u, tempVert.v, tempVert.t);
				//tempVert.limitloc=(1-t)*thisVert->limitloc+t*thatVert->limitloc;
				//tempVert.limitnorm=(1-t)*thisVert->limitnorm+t*thatVert->limitnorm;
				//tempVert.rad=(1-t)*thisVert->rad+t*thatVert->rad;
				//tempVert.uvt.u=SUB_FIGURE_CUTOFF_U;
				//tempVert.uvt.v=thisVert->uvt.v;
				//tempVert.uvt.t=thisVert->uvt.t;

				// remove garbage later! 110802
				vList[vOffset].setvals(tempVert.pnt[0], tempVert.pnt[1], tempVert.pnt[2],
									   tempVert.norm[0], tempVert.norm[1], tempVert.norm[2], 
									   tempVert.rad, tempVert.u, tempVert.v, tempVert.t);
				vTag[vOffset]=1;
				vOffset++;
			}
		}
		else
		{
			// need to lift up this curve to guarantee the cut sub-figure is more uniform
			if(thisVert->uvt.u-SUB_FIGURE_CUTOFF_U<UVT_UNIT)
			//if(thisVert->uvt.u==int(SUB_FIGURE_CUTOFF_U+.5))
			{
				Bpoint bp;
				pListOld[1]->subdivBPosition(&bp, thisVert->uvt.u+SUB_FIGURE_CUTOFF_TRANSFORM, thisVert->uvt.v, thisVert->uvt.t);
				vList[vOffset].setvals(bp.pnt[0], bp.pnt[1], bp.pnt[2],
									   bp.norm[0], bp.norm[1], bp.norm[2], 
									   bp.rad, thisVert->uvt.u+SUB_FIGURE_CUTOFF_TRANSFORM, thisVert->uvt.v, thisVert->uvt.t);
				//vList[vOffset].setvals(thisVert->limitloc.x(), thisVert->limitloc.y(), thisVert->limitloc.z()+0.035,
				//					   thisVert->limitnorm.x(), thisVert->limitnorm.y(), thisVert->limitnorm.x(), 
				//					   thisVert->rad, thisVert->uvt.u, thisVert->uvt.v, thisVert->uvt.t);
			}
			else
				vList[vOffset].setvals(thisVert->limitloc.x(), thisVert->limitloc.y(), thisVert->limitloc.z(),
									   thisVert->limitnorm.x(), thisVert->limitnorm.y(), thisVert->limitnorm.x(), 
									   thisVert->rad, thisVert->uvt.u, thisVert->uvt.v, thisVert->uvt.t);
			vTag[vOffset]=0;
			vOffset++;
		}
	}

	// working on a more rubust algorithm to find the cut-curve on the sub-figure
	// keep all the remaining patches of the sub-figure first

	//int maxU
	double maxU;

	// all the variables to calculate the cut-curve on the sub-figure
	int quadListSub[MAX_INTERSECTION_NUMBER*2], quadListSubNum=0;
	//int *quadListSub=new int[mainFNumEstimate], quadListSubNum;
	//int *vertListMainQuad=new int[mainVNumEstimate], *vertListMainQuadCnt=new int[mainVNumEstimate], vertListMainQuadNum;
	int edgeListSubQuad[MAX_INTERSECTION_NUMBER*2], edgeListSubQuadNum=0;
	//int *edgeListSubQuad=new int [mainVNumEstimate], edgeListSubQuadNum=0;
	//int *vertListMain=new int[mainVNumEstimate];
	int takenSubQuadEdge[MAX_INTERSECTION_NUMBER];
	memset(takenSubQuadEdge, 0, sizeof(int)*MAX_INTERSECTION_NUMBER);
	//int *takenSubQuadEdge=new int[mainFNumEstimate];

	int fDex;
	for(fDex=0; fDex<fNums; fDex++)
	{
		int fvdex;
		thisPoly=thisMesh->FaceList(fDex);
		minU=999;
		maxU=-1;
		for(fvdex=0; fvdex<thisPoly->numsides(); fvdex++)
		{
			thisVert=thisMesh->VertList(thisPoly->Pvertices[fvdex]);
			//thisVert=thisMesh->VertList(reDex[thisPoly->Pvertices[fvdex]]);
			if(minU>thisVert->uvt.u) minU=thisVert->uvt.u;
			if(maxU<thisVert->uvt.u) maxU=thisVert->uvt.u;
		}
		if(minU>SUB_FIGURE_CUTOFF_U)
		{
			fList[fOffset].Create(thisPoly->numsides());
			for(fvdex=0; fvdex<thisPoly->numsides(); fvdex++)
			{
				fList[fOffset][fvdex]=reDex[thisPoly->Pvertices[fvdex]];
			}
			fTag[fOffset]=0;
			fOffset++;
		}
		else
			if(maxU>=SUB_FIGURE_CUTOFF_U)
			{
				// 1. re-write the codes to connect the cut-curve on the subfigure and the dilated-i-curve
				quadListSub[quadListSubNum++]=fDex;

				fList[fOffset].Create(thisPoly->numsides());
				for(fvdex=0; fvdex<thisPoly->numsides(); fvdex++)
				{
					fList[fOffset][fvdex]=reDex[thisPoly->Pvertices[fvdex]];
				}
				fTag[fOffset]=1;
				fOffset++;
			}
	}
	vOffsetT+=vOffset;
	fOffsetT+=fOffset;

#ifdef QIONG_BLENDING_SUB_FIGURE
	color[0]=0;
	color[1]=0;
	color[2]=1;
	glRender(vOffsetT, fOffsetT, vList, fList, rMode, color);
#endif
	// done: so far vOffset vertices and fOffset faces in sub-figure have been inserted into vList and fList

//*
// 2.0 insert all the dilated intersection point into the vList!
// 2.1 before that, need to calculate the cut-curve of the sub-figure
	// now need to find the cut-curve of the subfigure, and connect the dilated intersection curve to it
	// then the sub-figure
	// all the variables to calculate the cut-curve on the sub-figure
	//int quadListSub[MAX_INTERSECTION_NUMBER*2], quadListSubNum=0;
	//int edgeListSubQuad[MAX_INTERSECTION_NUMBER*2], edgeListSubQuadNum=0;
	//int takenSubQuadEdge[MAX_INTERSECTION_NUMBER];

	// first find all the edges contributing to the cut-curve on the sub-figure
	for(i=0; i<quadListSubNum; i++)
	{
		thisPoly=thisMesh->FaceList(quadListSub[i]);
		pSides=thisPoly->numsides();
		for(k=0; k<pSides; k++)
		{
			thisVert=thisMesh->VertList(thisPoly->Pvertices[k]);
			if(thisVert->uvt.u<SUB_FIGURE_CUTOFF_U)
			{
				thisVert=thisMesh->VertList(thisPoly->Pvertices[(k+1)%pSides]);
				if(thisVert->uvt.u<SUB_FIGURE_CUTOFF_U)
				{
					edgeListSubQuad[2*edgeListSubQuadNum]=thisPoly->Pvertices[k];
					edgeListSubQuad[2*edgeListSubQuadNum+1]=thisPoly->Pvertices[(k+1)%pSides];
					edgeListSubQuadNum++;
				}
			}
		}
	}

	// time to connect the edges to form the cut-curve
	int subInterPNum=edgeListSubQuadNum;
	int subInterCurveDex[MAX_INTERSECTION_NUMBER];
	Bpoint2 subInterCurve[MAX_INTERSECTION_NUMBER];
	subInterCurveDex[0]=edgeListSubQuad[0];
	subInterCurveDex[1]=edgeListSubQuad[1];
	takenSubQuadEdge[0]=1;
	i=1;
	j=1;
	while(i<subInterPNum)
	{
		while(takenSubQuadEdge[j]) 
			j++;
		if(j>=edgeListSubQuadNum)
		{
			// the cut-curve can not be connected at all
			// time to fail the algorithm
			// i.e., this should NOT possibly happen
			fprintf(stderr, "\nThe intersection curve is not closed!!!\nHave to quit right here...\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
			fclose(fDebug);
#endif

			if(polyDexList1!=NULL)
				delete []polyDexList1;
			if(polyDexList2!=NULL)
				delete []polyDexList2;

			if(reDex!=NULL)
				delete []reDex;
			if(vList!=NULL)
				delete []vList;
			if(fList!=NULL)
				delete []fList;
			if(vTag!=NULL)
				delete []vTag;
			if(fTag!=NULL)
				delete []fTag;

			return false;
		}
		if(subInterCurveDex[i]==edgeListSubQuad[2*j])
		{
			i++;
			subInterCurveDex[i]=edgeListSubQuad[2*j+1];
			takenSubQuadEdge[j]=1;
			j=0;
		}
		else
			j++;
	}

	if(subInterPNum<=0)
	{
		// the cut-curve can not be connected at all
		// time to fail the algorithm
		// i.e., this should NOT possibly happen if the topology is valid, must be something wrong w/ the u/v/t
		fprintf(stderr, "\nNo CUT CURVE found on the sub-figure!\nHave to quit right here...\n");

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		if(reDex!=NULL)
			delete []reDex;
		if(vList!=NULL)
			delete []vList;
		if(fList!=NULL)
			delete []fList;
		if(vTag!=NULL)
			delete []vTag;
		if(fTag!=NULL)
			delete []fTag;

		return false;
	}

	// for indentation model, need to reverse this curve!
	if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
	{
		for(i=0; i<subInterPNum/2; i++)
		{
			j=subInterCurveDex[i];
			subInterCurveDex[i]=subInterCurveDex[subInterPNum-1-i];
			subInterCurveDex[subInterPNum-1-i]=j;
		}
		subInterCurveDex[subInterPNum]=subInterCurveDex[0];
	}

	// assign the values for the cut-curve
	for(i=0; i<=subInterPNum; i++)
	{
		subInterCurve[i]=vList[reDex[subInterCurveDex[i]]];
	}

#ifdef QIONG_DEBUG_20030510
	color[0]=1;
	color[1]=1;
	color[2]=1;
	for(i=0; i<vOffsetT; i++)
		if(vTag[i]==1)
		{
			glRenderPoint(vList[i], color);
		}
#endif

#ifdef QIONG_BLENDING_SUB_FIGURE_INTERSECTION_CURVE
	for(i=0; i<=subInterPNum; i++) //mainInterPNum; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "%d\t", subInterCurveDex[i]);
	#else
		fprintf(stderr, "%d\t", subInterCurveDex[i]);
	#endif
	}
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "\n");
	#else
		fprintf(stderr, "\n");
	#endif

	color[0]=0;
	color[1]=1;
	color[2]=0;
	glRenderCurve(0, subInterPNum, subInterCurve, color);
	//glRenderCurve(0, mainInterPNum, mainInterCurve, color);

	//pList1->glRenderCurve(1, 2, interCurveD, color);
#endif

	// still, need to check the closeness of the cut-curve to make it robust!
	if(subInterCurveDex[subInterPNum]!=subInterCurveDex[0])
	{
		// the cut-curve is not closed
		// time to fail the algorithm
		// i.e., this should NOT possibly happen
		fprintf(stderr, "\nThe cut-curve of the sub-figure is not closed!!!\nHave to quit right here...\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fclose(fDebug);
#endif

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		if(reDex!=NULL)
			delete []reDex;
		if(vList!=NULL)
			delete []vList;
		if(fList!=NULL)
			delete []fList;
		if(vTag!=NULL)
			delete []vTag;
		if(fTag!=NULL)
			delete []fTag;

		return false;
	}

// 2.2 after calculation of the cut-curve on the sub-figure
//     I will dilate the intersection curve in the main-figure's uvt space

	// NEW FUNCTION to check the orientation of the intersection curve
	// instead of using the DilateCurve to check the orientaion (later I won't map the whole curve 
	// into the 2D space then do the dilation any longer. instead, I will deal with the intersection 
	// curve chunk by chunk
	bool need2ReverseIntersectionCurve=CheckIntersectionCurveOrientation(subInterCurve, subInterPNum, interCurve, interPNum);
//	fprintf(stderr, "Need to reverse the intersection curve? %s\n", need2ReverseIntersectionCurve?"Yes.":"No");

#ifdef QIONG_DEBUG_20030510
	// try to check all the t coordinates of all the intersection vertices
	// they are supposed to be 1/-1 if uv belongs to (0, uMax) and (0, vMax) respectively
	for(i=0; i<interPNum; i++)
		//fprintf(stderr, "%8.7f\t%d\t", interCurveD[i].t, interCurveD[i].t==-1.0);
		fprintf(stderr, "%d\t", interCurveD[i].t==-1.0);
	fprintf(stderr, "\n");
#endif

	// now it might be time to reverse the whole intersection curve including all the 
	// other related arrays
	if(need2ReverseIntersectionCurve)
	{
		Bpoint2 tmpPnt;
		int tmpPair;
		for(i=0; i<interPNum/2; i++)
		{
			j=interPNum-1-i;

			tmpPnt=interCurve[i];
			interCurve[i]=interCurve[j];
			interCurve[j]=tmpPnt;

			tmpPair=interPPair[2*i];
			interPPair[2*i]=interPPair[2*j];
			interPPair[2*j]=tmpPair;

			tmpPair=interPPair[2*i+1];
			interPPair[2*i+1]=interPPair[2*j+1];
			interPPair[2*j+1]=tmpPair;
		}
		interCurve[interPNum]=interCurve[0];
		interPPair[2*interPNum]=interPPair[0];
		interPPair[2*interPNum+1]=interPPair[1];
	}

	double iCurveLengthBeforeDilation=CurveLength(interCurve, interPNum);
	if(!DilateCurve(pListOld[0], interCurve, interPNum, DILATION_DT, DILATION_STEPS))
	{
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "\nThe dilation of the intersection curve just FAILED, QUITTING...\n");
#endif
		fprintf(stderr, "\nThe dilation of the intersection curve just FAILED, QUITTING...\n");

#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fclose(fDebug);
#endif

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		if(reDex!=NULL)
			delete []reDex;
		if(vList!=NULL)
			delete []vList;
		if(fList!=NULL)
			delete []fList;
		if(vTag!=NULL)
			delete []vTag;
		if(fTag!=NULL)
			delete []fTag;

		return false;
	}
	double iCurveLengthAfterDilation=CurveLength(interCurve, interPNum);
	if(iCurveLengthBeforeDilation > iCurveLengthAfterDilation)
	{
		fprintf(stderr, "\nINCORRECT ORIENTATION of the INTERSECTION CURVE in DILATION, QUITTING...\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fclose(fDebug);
#endif

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		if(reDex!=NULL)
			delete []reDex;
		if(vList!=NULL)
			delete []vList;
		if(fList!=NULL)
			delete []fList;
		if(vTag!=NULL)
			delete []vTag;
		if(fTag!=NULL)
			delete []fTag;

		return false;
	}


	// Re-regularize the intersection curve after the dilation
	// after the dilation step, still need to re-regularize the dilated curve because at high-curvature
	// part, vertices tend to be crunched together
	//memcpy(interCurveD, interCurve, (interPNum+5)*sizeof(Bpoint2));
	RegularizeDilatedIntersectionCurve(numPointsOnCutCurve, interCurveD, interCurve, &interPNum);
	//RegularizeDilatedIntersectionCurve(interCurveD, interCurve, &interPNum);
	if(interPNum>MAX_INTERSECTION_NUMBER || interPNum<=0)
	{
		fprintf(stderr, "NONE DILATED INTERSECTION POINTS left after the RegularizeDilatedIntersectionCurve!\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fclose(fDebug);
#endif

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		if(reDex!=NULL)
			delete []reDex;
		if(vList!=NULL)
			delete []vList;
		if(fList!=NULL)
			delete []fList;
		if(vTag!=NULL)
			delete []vTag;
		if(fTag!=NULL)
			delete []fTag;

		return false;
	}

#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nThe v/t of all the %d vertices on the regularized dilated i-curve...\n", interPNum);
	for(i=0; i<=interPNum; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f\t", i, interCurveD[i].v_2, interCurveD[i].t_2);
	}
	fprintf(fDebug, "\n");
#endif

	//bool wrappedAround;
#ifdef QIONG_BLENDING_INTERSECTION_CURVE_DEBUG	
	bool wrappedAround;
	for(i=0; i<=interPNum; i++)
		fprintf(stderr, "%d, %d\n", interPPair[2*i], interPPair[2*i+1]);
	// by this way we can find the segment on same curve patch one by one
	// we start with the main figure
	fprintf(stderr, "%d segments in total!\n", interPNum);
	fprintf(stderr, "main-figure now!\n", interPNum);

	i=0;
	j=0;
	while(interPPair[2*i]==interPPair[2*(i+1)] && i<interPNum) i++;
	if(i<interPNum-1)
	{
	wrappedAround=false;
	while(j<interPNum)
	{
		j=i+1;
		while(interPPair[2*(j%interPNum)]==interPPair[2*((j+1)%interPNum)] && !wrappedAround) 
		{
			j++;
			if(j%interPNum==i+1) 
				wrappedAround=true;
		}
		// then everytime, i+1->j+1 is the new set of segments on the same patch of main figure
		fprintf(stderr, "segment from %d to %d...\n", i+1, j+1);
		i=j;
	}
	}

	wrappedAround=false;

	// then the sub-figure
#ifndef QIONG_BLENDING_NO_DEBUG
	fprintf(stderr, "sub-figure now!\n", interPNum);
#endif

	i=0;
	j=0;
	while(interPPair[2*i+1]==interPPair[2*(i+1)+1] && i<interPNum) i++;
	if(i<interPNum-1)
	{
	while(j<interPNum)
	{
		j=i+1;
		while(interPPair[2*(j%interPNum)+1]==interPPair[2*((j+1)%interPNum)+1] && !wrappedAround) 
		{
			j++;
			if(j%interPNum==i+1) 
				wrappedAround=true;
		}
		// then everytime, i+1->j+1 is the new set of segments on the same patch of sub figure
#ifndef QIONG_BLENDING_NO_DEBUG
		fprintf(stderr, "segment from %d to %d...\n", i+1, j+1);
#endif
		i=j;
	}
	}
#endif

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_DILATED
	color[0]=0;
	color[1]=1;
	color[2]=0;
	//glRenderCurve(0, interPNum, interCurveD, color);
	glRenderCurve(0, interPNum, interCurveD, color);
#endif

#ifdef QIONG_DEBUG_20030625
	color[0]=.5;
	color[1]=1;
	color[2]=.5;
	//glRenderPoint(interCurveD[11], color);
	glRenderPoint(interCurveD[11], color);
#endif

// 2.5 finally, insert all the vertices on the dilated intersection curve to the vList
	int vOffset2=0, fOffset2=0;

	for(i=0; i<interPNum; i++)
	{
		vList[vOffsetT+vOffset2].setvals(interCurveD[i].pnt[0], interCurveD[i].pnt[1], interCurveD[i].pnt[2], 
										 interCurveD[i].norm[0], interCurveD[i].norm[1], interCurveD[i].norm[2], 
										 interCurveD[i].rad, interCurveD[i].u, interCurveD[i].v, interCurveD[i].t,
										 interCurveD[i].rad_2, interCurveD[i].u_2, interCurveD[i].v, interCurveD[i].t_2);
		vTag[vOffsetT+vOffset2]=2;
		vOffset2++;
	}
	vOffsetT+=vOffset2;
	fOffsetT+=fOffset2;

// 3.0 connect the cut sub-figure with the dilated intersection curve
#ifndef QIONG_BLENDING_NO_DEBUG
	fprintf(stderr, "sub-figure now!\n");
#endif

	int vOffset3=0, fOffset3=0;
	int fvdex;

	//int subVertOldDex[2];

//	int polyDex2;
	// need to re-write all the codes here to connect the dilated-intersection-curve and 
	// the subfigure cut-curve
	//MyFloat tClosest=999, vClosest=999;
	int triNum=0, triList[MAX_INTERSECTION_NUMBER];
	memset(triList, 0, MAX_INTERSECTION_NUMBER*sizeof(int));
	int subDex, subNext;
	int interDex, interNext;
	int interStartPos;

	interStartPos=0;
	MeshCurves(&triNum, triList, &interStartPos, subInterPNum, subInterCurve, interPNum, interCurveD);

	subDex=0;
	interDex=interStartPos%interPNum;
	for(i=0; i<triNum; i++)
	{
		switch(triList[i])
		{
		case 1:
			subNext=(subDex+1)%subInterPNum;
				fList[fOffsetT+fOffset3].Create(3);
				fList[fOffsetT+fOffset3][0]=vOffset+interDex;
				fList[fOffsetT+fOffset3][1]=reDex[subInterCurveDex[subNext]];
				fList[fOffsetT+fOffset3][2]=reDex[subInterCurveDex[subDex]];
				fOffset3++;
			subDex=subNext;
			break;
		case 2:
			interNext=(interDex+1)%interPNum;
				fList[fOffsetT+fOffset3].Create(3);
				fList[fOffsetT+fOffset3][0]=vOffset+interDex;
				fList[fOffsetT+fOffset3][1]=vOffset+interNext;
				fList[fOffsetT+fOffset3][2]=reDex[subInterCurveDex[subDex]];
				fOffset3++;
			interDex=interNext;
			break;
		default:
			fprintf(stderr, "\nUNEXPECTED EXCEPTION FROM MeshCurves(...)!\n");
			return false;
		}
	}

#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\n");
#endif

#ifdef QIONG_BLENDING_SUB_FIGURE_INTERSECTION
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nThe v/t of all the %d vertices on the cut-curve of the sub-figure...\n", subInterPNum);
	for(i=0; i<=subInterPNum; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f\t", i, subInterCurve[i].v, subInterCurve[i].t);
	}
	fprintf(fDebug, "\nThe v/t of all the %d vertices on the dilated-i-curve...\n", interPNum);
	for(i=0; i<=interPNum; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f\t", i, interCurveD[i].v_2, interCurveD[i].t_2);
	}
	fprintf(fDebug, "\n");
	#endif
#endif

	vOffsetT+=vOffset3;
	fOffsetT+=fOffset3;

#ifdef QIONG_BLENDING_SUB_FIGURE_INTERSECTION
	color[0]=1;
	color[1]=1;
	color[2]=0;
	glRender(vOffsetT, fOffset3, vList, &fList[fOffset+fOffset2], rMode, color);
#endif

#ifdef QIONG_BLENDING_SUB_FIGURE_QUAD_LIST
	color[0]=1;
	color[1]=0;
	color[2]=0;
	glRender(SUBDIV_LEVEL_TO_USE, 39, rMode, color);
//New InterCurvD-subPolyDex: 3    4       4       5       5       38      38
//38      39      39      35      20      19      18      59      59      58
//58      58      3
#endif

//*
// 4.0 find the cut version of the main-figure
	// qiong: re-write this part to re-mesh the blending area
	// for main-figure, re-relabling vertices first
	thisMesh=thisMesh1;
	vNums=thisMesh->VertList.Length();
	fNums=thisMesh->FaceList.Length();
	int vOffset4=0, fOffset4=0;

	// before everything, need to find the cut curve on the main figure that is connected to the dilated inter-curve
	// then figure out what vertices in the main figure to get rid of!

	int mainVNumEstimate=vNums*1.5;
	int mainFNumEstimate=fNums*1.5;
	int *quadListMain=new int[mainFNumEstimate], quadListMainNum;
	int *vertListMainQuad=new int[mainVNumEstimate], *vertListMainQuadCnt=new int[mainVNumEstimate], vertListMainQuadNum;
	int *edgeListMainQuad=new int [mainVNumEstimate], edgeListMainQuadNum=0;
	int *vertListMain=new int[mainVNumEstimate];

	//int quadListMain[200], quadListNum;
	//int vertListMainQuad[500], vertListMainQuadCnt[500], vertListMainQuadNum;
	//int edgeListMainQuad[500], edgeListMainQuadNum=0;

	// flag for each vertex of the main-figure quads, 1 means it should be kept otherwise cut
	// index [0, vNums-1]
	// range [0, 1]
	memset(vertListMain, 1, mainVNumEstimate*sizeof(int));

	// quadListMainNum counts the total # of the intersected quads
	quadListMainNum=0;
	// the list of all the quads involved with the intersection with the sub-figure
	// index [0, quadListMainNum-1]
	// range [0, fNums-1]
	memset(quadListMain, 1, mainFNumEstimate*sizeof(int));

	// vertListMainQuadNum counts the # of vertices of all the intersected quads
	vertListMainQuadNum=0;
	// list of vertices of the intersected main-figured
	// index [0, vertListMainQuadNum-1]
	// range [0, vNums-1]
	memset(vertListMainQuad, -1, mainVNumEstimate*sizeof(int));

	// vertListMainQuadCnt counts the # of appearances of each vertex in the intersected 
	//		quads. if it is the valence of the vertex, it should be removed because it's 
	//		already surrouded by quads that should be removed
	// index [0, vNums-1]
	memset(vertListMainQuadCnt, 0, mainVNumEstimate*sizeof(int));

	pListOld[0]->BuildLookUpTableForPDexByUVT();

	// according to the uvt coordinates of the dilated intersection curve, need to recalculate polyDexList1
	int pDexTmp; //, subLvl;
	//CCPolygon *thisPoly;
	for(i=0; i<interPNum; i++)
	{
		//pDexTmp=pListOld[0]->lookUpTable[LookUpPDexByUVT(pListOld[0], interCurveD[i])];
		pDexTmp=pListOld[0]->LookUpPDexByUVT(interCurveD[i]);//lookUpTable[LookUpPDexByUVT(pListOld[0], interCurveD[i])];
if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
		pDexTmp=pListOld[0]->WithinPolyDexAtSublvlUVT(pDexTmp, SUBDIV_LEVEL_TO_USE+1, interCurveD[i].u, interCurveD[i].v, interCurveD[i].t);
else
		pDexTmp=pListOld[0]->WithinPolyDexAtSublvlUVT(pDexTmp, SUBDIV_LEVEL_TO_USE, interCurveD[i].u, interCurveD[i].v, interCurveD[i].t);

		//for(subLvl=1; subLvl<=SUBDIV_LEVEL_TO_USE; subLvl++)
		//{
		//	thisPoly->
		//	//thisPoly=pListOld[0]->
		//	//pDexTmp=...;
		//}
		if(!polyDexList1[pDexTmp]) 
			polyDexList1[pDexTmp]=2;
	}
	//Pointlist_serverB::LookUpPolygonIndexByUVT(Pointlist_server2 *pList, interCurveD[i])

#ifdef QIONG_DEBUG_20030510
	for(i=0; i<fNums; i++)
		if(polyDexList1[i])
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			fprintf(fDebug, "%d\t", i);
	#else
			fprintf(stderr, "%d\t", i);
	#endif
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\n");
	#else
	fprintf(stderr, "\n");
	#endif
#endif

#ifdef QIONG_BLENDING_MAIN_FIGURE_QUAD_LIST
	color[0]=1;
	color[1]=0;
	color[2]=0;

if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
{
	glRender(SUBDIV_LEVEL_TO_USE+1, 8, rMode, color);

}
else
{
	glRender(SUBDIV_LEVEL_TO_USE, 15, rMode, color);
}
//	7	10	12	17	18	20	25	26	
#endif

	// now add the list of quads intersected by the dilated intersection curve into quadListMain
	for(i=0; i<fNums; i++)
		if(polyDexList1[i])
			quadListMain[quadListMainNum++]=i;

#ifndef QIONG_BLENDING_NO_DEBUG
#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nRecording all the 'curV's for vertListMainQuadCnt\n");
#endif
#endif

	int *edgeListMainQuadDex=new int[mainVNumEstimate];
	memset(edgeListMainQuadDex, -1, sizeof(int)*mainVNumEstimate);

	int curV, curEHd, curEEd;
	for(i=0; i<quadListMainNum; i++)
	{
		thisPoly=thisMesh->FaceList(quadListMain[i]);
		pSides=thisPoly->numsides();
		for(j=0; j<pSides; j++)
		{
			// try to add curV into the vert list
			curV=thisPoly->Pvertices[j];
/*
			//objectFound=false;
			k=0;
			// check whether curV is already in the list
			while(k<vertListMainQuadNum)
			{
				if(curV==vertListMainQuad[k])
				{
					//objectFound=true;
					break;
				}
				k++;
			}
			// if it's not, add it to the list
			//		and update the counter by 1 for the first time
			if(k==vertListMainQuadNum)
			{
				vertListMainQuad[vertListMainQuadNum++]=curV;
				vertListMainQuadCnt[curV]++;
#ifndef QIONG_BLENDING_NO_DEBUG
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				if(curV>=mainVNumEstimate)
					fprintf(fDebug, "%d\t", curV);
#endif
#endif
			}
			// otherwise need to update the counter and find the shared vertices
			else
			{
				vertListMainQuadCnt[curV]++;
#ifndef QIONG_BLENDING_NO_DEBUG
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				if(curV>=mainVNumEstimate)
					fprintf(fDebug, "%d\t", curV);
#endif
#endif
			}
*/

			// now need to deal with the edges
			//		it's slightly different cause we remove all the shared edges
			curEHd=curV;
			curEEd=thisPoly->Pvertices[(j+1)%pSides];
			k=0;
			while(k<edgeListMainQuadNum)
			{
				// have to reverse the order of curEEd and curEHd because of the topology
				if(curEEd==edgeListMainQuad[2*k] && curEHd==edgeListMainQuad[2*k+1])
					break;
				k++;
			}

			// if it's not, add it to the list
			//		and update the counter by 1 for the first time
			if(k==edgeListMainQuadNum)
			{
				edgeListMainQuadDex[edgeListMainQuadNum]=quadListMain[i];

				edgeListMainQuad[2*edgeListMainQuadNum]=curEHd;
				edgeListMainQuad[2*edgeListMainQuadNum+1]=curEEd;
				edgeListMainQuadNum++;
			}
			// otherwise need to remove the shared edge from the edgeListMainQuad
			else
			{
				if(k<edgeListMainQuadNum-1)
				{
					while(k<edgeListMainQuadNum-1)
					{
						edgeListMainQuadDex[k]=edgeListMainQuadDex[k+1];;

						edgeListMainQuad[2*k]=edgeListMainQuad[2*(k+1)];
						edgeListMainQuad[2*k+1]=edgeListMainQuad[2*(k+1)+1];
						k++;
					}
				}
				edgeListMainQuadNum--;
			}
		}
	}

#ifndef QIONG_BLENDING_NO_DEBUG
#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nRecording all the 'vertListMainQuad[i]'s for vertListMainQuadCnt\n");
#endif
#endif

	// it is time to find the curve on the main figure to connect to the dilated intersection curve!
	// write the codes to find out the cut-curve on the main figure
	int mainInterPNum;

	int *takenMainQuadEdge=new int[mainFNumEstimate];
	memset(takenMainQuadEdge, 0, sizeof(int)*mainFNumEstimate);
	Bpoint2 mainInterCurve[MAX_INTERSECTION_NUMBER];
	int mainInterCurveDex[MAX_INTERSECTION_NUMBER];
	//mainInterPNum=edgeListMainQuadNum;

	int mainInterCurveDexInEdgeListMainQuad[MAX_INTERSECTION_NUMBER];
	memset(mainInterCurveDexInEdgeListMainQuad, -1, sizeof(int)*MAX_INTERSECTION_NUMBER);

	bool algorithmFailed=false;
	bool timeToCheckCloseness=false;
	int mCurveNum=0;
	MyPoint2D mCurveIndices[MAX_MAIN_CURVE_NUMBER];
	int mCurveStartDex, mCurveEndDex;
	int successfulMCurveID=-1;

	//memset(taken, 0, sizeof(int)*500);
	mainInterCurveDex[0]=edgeListMainQuad[0];
	mainInterCurveDex[1]=edgeListMainQuad[1];
	mainInterCurveDexInEdgeListMainQuad[1-1]=0;
	mCurveStartDex=0;
	j=0;
	//lastTaken=0;
	i=1;

	takenMainQuadEdge[0]=1;
	//taken[0]=1;
	// generate the whole intersection curve from the separated intersection segments
	//		all the un-ordereed intersection vertices in interPList now will be linked as a curve into interPList2
	successful=false;
	while(i<edgeListMainQuadNum && !successful)
	{
		while(j<=edgeListMainQuadNum && takenMainQuadEdge[j]) 
			j++;
		if(j>edgeListMainQuadNum)
		{
			// another possibility here is that the curve is already closed! (multiple curves could exist)
			// let us figure out whether it is true and whether the current curve is the outer curve

			if(mainInterCurveDex[i]==mainInterCurveDex[mCurveStartDex])
			{
				// it is already a closed curve
				// need to check whether it is indeed the curve we are looking for
				// only the curve that is COUNTER-CLOCKWISE should be kept!
				if(i-mCurveStartDex>edgeListMainQuadNum/2)
				{
					// if the number of vertices in this curve is over the half of number of all the vertices
					// it is a trivial positive!
					mCurveEndDex=i;
					mCurveIndices[mCurveNum].x=mCurveStartDex;
					mCurveIndices[mCurveNum].y=mCurveEndDex;
					successfulMCurveID=mCurveNum;
					mCurveNum++;
					successful=true;
					continue;
				}
				else
				{
					// otherwise, need to continue the searching
					mCurveEndDex=i;
					mCurveIndices[mCurveNum].x=mCurveStartDex;
					mCurveIndices[mCurveNum].y=mCurveEndDex;
					mCurveNum++;

					i++;
					mCurveStartDex=i;
					j=0;
					while(j<=edgeListMainQuadNum && takenMainQuadEdge[j]) 
						j++;
					if(j<=edgeListMainQuadNum)
					{
						mainInterCurveDexInEdgeListMainQuad[i]=j;

						takenMainQuadEdge[j]=1;
						mainInterCurveDex[i]=edgeListMainQuad[2*j];
						mainInterCurveDex[i+1]=edgeListMainQuad[2*j+1];
						j=0;

						//lastTaken=0;
						i++;
					}
					else
					{
						// the algorithm just failed!
						algorithmFailed=true;
					}
					continue;
				}
			}

			if(algorithmFailed)
			{
				// it is time to quit because this is the illegal case!
				bool trap=true;
				fprintf(stderr, "\nCannot figure out the mainInterCurve! time to quit blending codes\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fclose(fDebug);
#endif
			{
				if(polyDexList1!=NULL)
					delete []polyDexList1;
				if(polyDexList2!=NULL)
					delete []polyDexList2;

				if(reDex!=NULL)
					delete []reDex;
				if(vList!=NULL)
					delete []vList;
				if(fList!=NULL)
					delete []fList;
				if(vTag!=NULL)
					delete []vTag;
				if(fTag!=NULL)
					delete []fTag;

				if(quadListMain!=NULL)
					delete []quadListMain;
				if(vertListMainQuad!=NULL)
					delete []vertListMainQuad;
				if(vertListMainQuadCnt!=NULL)
					delete []vertListMainQuadCnt;
				if(edgeListMainQuad!=NULL)
					delete []edgeListMainQuad;
				if(vertListMain!=NULL)
					delete []vertListMain;

				if(edgeListMainQuadDex!=NULL)
					delete []edgeListMainQuadDex;

				if(takenMainQuadEdge!=NULL)
					delete []takenMainQuadEdge;
			}

				return false;
			}
		}

		if(mainInterCurveDex[i]==edgeListMainQuad[2*j])
		{
			mainInterCurveDexInEdgeListMainQuad[i]=j;

			takenMainQuadEdge[j]=1;
			i++;
			mainInterCurveDex[i]=edgeListMainQuad[2*j+1];
			j=0;
			if(i==edgeListMainQuadNum)
				timeToCheckCloseness=true;
		}
		else
			if(mainInterCurveDex[i]==edgeListMainQuad[2*j+1])
			{
				// TO REMOVE: these codes are theoretically NOT-REACHABLE!
				mainInterCurveDexInEdgeListMainQuad[i]=j;

				takenMainQuadEdge[j]=1;
				i++;
				mainInterCurveDex[i]=edgeListMainQuad[2*j];
				j=0;
				if(i==edgeListMainQuadNum)
					timeToCheckCloseness=true;
			}
			else j++;

		if(timeToCheckCloseness)
		{
			if(mainInterCurveDex[i]==mainInterCurveDex[mCurveStartDex])
			{
				// it is already a closed curve
				// need to check whether it is indeed the curve we are looking for
				// only the curve that is COUNTER-CLOCKWISE, or the LONGEST should be kept!
				if(i-mCurveStartDex>edgeListMainQuadNum/2)
				{
					mCurveEndDex=i;
					mCurveIndices[mCurveNum].x=mCurveStartDex;
					mCurveIndices[mCurveNum].y=mCurveEndDex;
					successfulMCurveID=mCurveNum;
					mCurveNum++;
					successful=true;
				}
				else
				{
					mCurveEndDex=i;
					mCurveIndices[mCurveNum].x=mCurveStartDex;
					mCurveIndices[mCurveNum].y=mCurveEndDex;
					mCurveNum++;

					successful=true;
					//successful=true;
				}
			}
			else
			{
				// time to fail the algorithm again, because the curve is not closed!
				fprintf(stderr, "\nCannot figure out the mainInterCurve! time to quit blending codes\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fclose(fDebug);
#endif
			{
				if(polyDexList1!=NULL)
					delete []polyDexList1;
				if(polyDexList2!=NULL)
					delete []polyDexList2;

				if(reDex!=NULL)
					delete []reDex;
				if(vList!=NULL)
					delete []vList;
				if(fList!=NULL)
					delete []fList;
				if(vTag!=NULL)
					delete []vTag;
				if(fTag!=NULL)
					delete []fTag;

				if(quadListMain!=NULL)
					delete []quadListMain;
				if(vertListMainQuad!=NULL)
					delete []vertListMainQuad;
				if(vertListMainQuadCnt!=NULL)
					delete []vertListMainQuadCnt;
				if(edgeListMainQuad!=NULL)
					delete []edgeListMainQuad;
				if(vertListMain!=NULL)
					delete []vertListMain;

				if(edgeListMainQuadDex!=NULL)
					delete []edgeListMainQuadDex;

				if(takenMainQuadEdge!=NULL)
					delete []takenMainQuadEdge;
			}

				return false;
			}
		}
	}
	//mainInterPNum=i;

	// by now
	// mCurveNum			holds the total number of the curves from the intersection quads
	// mCurveIndices		holds the start/end indices of the segments of each curve in mainInterCurveDex
	int maxLength=0;
	if(successfulMCurveID<0)
	//if(!successful)
	{
		// then need to figure out which one of the curves is the intersection curve we need
		for(i=0; i<mCurveNum; i++)
		{
			if(maxLength<mCurveIndices[i].y-mCurveIndices[i].x)
			{
				maxLength=mCurveIndices[i].y-mCurveIndices[i].x;
				successfulMCurveID=i;
			}
		}
	}

	// successfulMCurveID	holds the ID of the curve that is the cut curve of the main figure
	mCurveStartDex=mCurveIndices[successfulMCurveID].x;		// right here right now
	mCurveEndDex=mCurveIndices[successfulMCurveID].y;
	mainInterPNum=mCurveEndDex-mCurveStartDex;

	if(mCurveNum<=0 || mainInterPNum<=0)
	{
		fprintf(stderr, "\nCannot find a valid BOUNDARY CURVE on the main-figure!\nHave to quit now...\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fclose(fDebug);
#endif
		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		if(reDex!=NULL)
			delete []reDex;
		if(vList!=NULL)
			delete []vList;
		if(fList!=NULL)
			delete []fList;
		if(vTag!=NULL)
			delete []vTag;
		if(fTag!=NULL)
			delete []fTag;

		if(quadListMain!=NULL)
			delete []quadListMain;
		if(vertListMainQuad!=NULL)
			delete []vertListMainQuad;
		if(vertListMainQuadCnt!=NULL)
			delete []vertListMainQuadCnt;
		if(edgeListMainQuad!=NULL)
			delete []edgeListMainQuad;
		if(vertListMain!=NULL)
			delete []vertListMain;

		if(edgeListMainQuadDex!=NULL)
			delete []edgeListMainQuadDex;

		if(takenMainQuadEdge!=NULL)
			delete []takenMainQuadEdge;
	}

#ifdef QIONG_BLENDING_DILATED_ICURVE_RELATED_QUADS
	color[0]=0;
	color[1]=0;
	color[2]=1;
	//*
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "The list of main-figure quads related to the dilated-i-curve...\n");
	#endif
	j=mCurveStartDex;
if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
{
	for(i=0; i<mainInterPNum; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		if(i%8==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d\t", edgeListMainQuadDex[mainInterCurveDexInEdgeListMainQuad[j]]);
	#endif
		glRender(SUBDIV_LEVEL_TO_USE+1, edgeListMainQuadDex[mainInterCurveDexInEdgeListMainQuad[j]], rMode, color);
		j++;
	}
}
else
{
	for(i=0; i<mainInterPNum; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		if(i%8==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d\t", edgeListMainQuadDex[mainInterCurveDexInEdgeListMainQuad[j]]);
	#endif
		glRender(SUBDIV_LEVEL_TO_USE, edgeListMainQuadDex[mainInterCurveDexInEdgeListMainQuad[j]], rMode, color);
		j++;
	}
}
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\n");
	#endif
#endif

	// now it is time to cut the main figure by this intersection curve
	// need to grow the cut-curve and tag all the traversed quads... 
	// to quadListMain[quadListMainNum...]
	// before that, need two things
	//		1. the list of the indices of quad containing each segment in edgeListMainQuad	done
	//		2. a Queue stucture																done
	int curEdgeDexInEdgeListMainQuad, curEdgeStartDex, curEdgeEndDex;

	// before everything, tag all the edges of the dilated-i-curve related quads to span
	int *quadEdgeToSpan= new int[mainFNumEstimate*4];
	memset(quadEdgeToSpan, 1, mainFNumEstimate*4);

	j=mCurveStartDex;
	for(i=0; i<mainInterPNum; i++)
	{
		curEdgeDexInEdgeListMainQuad=mainInterCurveDexInEdgeListMainQuad[j];
		curEdgeStartDex=edgeListMainQuad[2*curEdgeDexInEdgeListMainQuad];
		curEdgeEndDex=edgeListMainQuad[2*curEdgeDexInEdgeListMainQuad+1];
		thisPoly=thisMesh->FaceList(edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]);
		pSides=thisPoly->numsides();
		for(k=0; k<pSides; k++)
		{
			// need to only span the quad from the edge that was not contributing to the mainInterCurve
			// in other words, only span the quad inward
			if(thisPoly->Pvertices[k]==curEdgeStartDex && thisPoly->Pvertices[(k+1)%pSides]==curEdgeEndDex)
			{
				quadEdgeToSpan[4*edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]+k]=0;
			}
		}
		j++;
	}

	int *visited=new int[mainFNumEstimate];
	memset(visited, 0, sizeof(int)*mainFNumEstimate);
	MyQueue que(mainVNumEstimate);

	// fill in the queue for the first time
	j=mCurveStartDex;
	for(i=0; i<mainInterPNum; i++)
	{
		curEdgeDexInEdgeListMainQuad=mainInterCurveDexInEdgeListMainQuad[j];
		curEdgeStartDex=edgeListMainQuad[2*curEdgeDexInEdgeListMainQuad];
		curEdgeEndDex=edgeListMainQuad[2*curEdgeDexInEdgeListMainQuad+1];
		thisPoly=thisMesh->FaceList(edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]);

		if(!visited[edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]])
		{
			visited[edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]]=1;

			pSides=thisPoly->numsides();
			for(k=0; k<pSides; k++)
			{
				// need to only span the quad from the edge that was not contributing to the mainInterCurve
				// in other words, only span the quad inward
				if(quadEdgeToSpan[4*edgeListMainQuadDex[curEdgeDexInEdgeListMainQuad]+k])
				//if(!(thisPoly->Pvertices[k]==curEdgeStartDex && thisPoly->Pvertices[(k+1)%pSides]==curEdgeEndDex))
				{
#ifdef QIONG_BLENDING_EXPANDED_MAIN_FIGURE_QUADS_TO_REMOVE
					if(thisPoly->Pneighbors[k]==15)
					{
						bool ttt=true;
					}
#endif
					if(!visited[thisPoly->Pneighbors[k]])
					{
						if(!polyDexList1[thisPoly->Pneighbors[k]])
						{
							polyDexList1[thisPoly->Pneighbors[k]]=3;

							quadListMain[quadListMainNum++]=thisPoly->Pneighbors[k];
						}

						que.Enque(thisPoly->Pneighbors[k]);
					}
				}
			}
		}
		j++;
	}

	// by traversing the queue, expand quads related to the dilated intersection curve to all the quads
	// that should be removed from the main figure quads
	int curQuadDex;
	while(!que.Empty())
	{
		curQuadDex=que.Deque();
		if(curQuadDex==-1)
		{
			bool trap=true;
		}
		if(!visited[curQuadDex])
		{
			visited[curQuadDex]=1;
			thisPoly=thisMesh->FaceList(curQuadDex);
			pSides=thisPoly->numsides();
			for(k=0; k<pSides; k++)
			{
#ifdef QIONG_BLENDING_EXPANDED_MAIN_FIGURE_QUADS_TO_REMOVE
				if(thisPoly->Pneighbors[k]==15)
				{
					bool ttt=true;
				}
#endif
				if(!visited[thisPoly->Pneighbors[k]])
				{
					if(!polyDexList1[thisPoly->Pneighbors[k]])
					{
						polyDexList1[thisPoly->Pneighbors[k]]=3;

						quadListMain[quadListMainNum++]=thisPoly->Pneighbors[k];
					}

					que.Enque(thisPoly->Pneighbors[k]);
				}
			}
		}
	}
	delete []visited;
	delete []quadEdgeToSpan;

#ifdef QIONG_BLENDING_EXPANDED_MAIN_FIGURE_QUADS_TO_REMOVE
	color[0]=0;
	color[1]=1;
	color[2]=0;
	//*
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "The list of main-figure quads that was found by spanning the dilated-i-curve quads...\n");
	#endif
if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
{
	for(i=0; i<nTiles1; i++)
	{
		if(polyDexList1[i]==3)
		{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			if(i%8==0)
				fprintf(fDebug, "\n");
			fprintf(fDebug, "%d\t", i);
	#endif
			glRender(SUBDIV_LEVEL_TO_USE+1, i, rMode, color);
		}
	}
}
else
{
	for(i=0; i<nTiles1; i++)
	{
		if(polyDexList1[i]==3)
		{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			if(i%8==0)
				fprintf(fDebug, "\n");
			fprintf(fDebug, "%d\t", i);
	#endif
			glRender(SUBDIV_LEVEL_TO_USE, i, rMode, color);
		}
	}
}

	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\n");
	#endif
#endif

	// update the quad-tags indicating that they are in the cut region
	// also tag all the vertices that are entirely surrounded by these cut quads
	for(i=0; i<quadListMainNum; i++)
	{
		thisPoly=thisMesh->FaceList(quadListMain[i]);
		pSides=thisPoly->numsides();
		for(j=0; j<pSides; j++)
		{
			// try to add curV into the vert list
			curV=thisPoly->Pvertices[j];
			k=0;
			// check whether curV is already in the list
			while(k<vertListMainQuadNum)
			{
				if(curV==vertListMainQuad[k])
					break;
				k++;
			}
			// if it's not, add it to the list
			//		and update the counter by 1 for the first time
			if(k==vertListMainQuadNum)
			{
				vertListMainQuad[vertListMainQuadNum++]=curV;
				vertListMainQuadCnt[curV]++;
#ifndef QIONG_BLENDING_NO_DEBUG
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				if(curV>=mainVNumEstimate)
					fprintf(fDebug, "%d\t", curV);
#endif
#endif
			}
			// otherwise need to update the counter and find the shared vertices
			else
			{
				vertListMainQuadCnt[curV]++;
#ifndef QIONG_BLENDING_NO_DEBUG
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				if(curV>=mainVNumEstimate)
					fprintf(fDebug, "%d\t", curV);
#endif
#endif
			}	
		}
	}

	// tag all the vertices on the main figure to remove
	for(i=0; i<vertListMainQuadNum; i++)
	{
		thisVert=thisMesh->VertList(vertListMainQuad[i]);
		if(thisVert->v_valence==vertListMainQuadCnt[vertListMainQuad[i]])
		//if(thisVert->v_valence==vertListMainQuadCnt[i])
			vertListMain[vertListMainQuad[i]]=0;
			//vertListMain[i]=0;
#ifndef QIONG_BLENDING_NO_DEBUG
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		if(vertListMainQuad[i]>=mainVNumEstimate)
			fprintf(fDebug, "%d\t", vertListMainQuad[i]);
#endif
#endif
	}

	//int subPNum;
	//Bpoint subFCurve[100];
	//int reDex[200];
	//int reduction;
	reduction=0;
	// reDex is the important array to re-index all the orginal main-figure vertices to the new array in vList
//	for(i=0; i<MAX_INTERSECTION_NUMBER; i++)
	for(i=0; i<vNumEstimate; i++)
		reDex[i]=i;
//*
	// reindex the vertices for the main-figure first
	//		need to remove all the vertices according to vertListMain[]
	//		0 means to be removed
	for(vDex=0; vDex<vNums; vDex++)
	{
		if(!vertListMain[vDex])
		{
			reDex[vDex]=-1;
			reduction++;
			continue;
		}
		//reDex[vDex]-=reduction;
		//thisVert=thisMesh->VertList(vDex);
		else
		// 110802 need to get rid of all the points inside the dilated intersection curve
		//if(polyDexList1[])
		//{
		//	reDex[vDex]=-1;
		//	reduction++;
		//}
		//else
		{
			thisVert=thisMesh->VertList(vDex);
			reDex[vDex]=vOffsetT+vOffset4;
			//if(vDex==26||vDex==17||vDex==21||vDex==25||vDex==22||vDex==18)
			//vList[vOffsetT+vOffset4].setvals(thisVert->limitloc.x(), thisVert->limitloc.y(), thisVert->limitloc.z()+.2,
			//						thisVert->limitnorm.x(), thisVert->limitnorm.y(), thisVert->limitnorm.x(), 
			//						thisVert->rad, thisVert->uvt.u, thisVert->uvt.v, thisVert->uvt.t);
			//else
			vList[vOffsetT+vOffset4].setvals(thisVert->limitloc.x(), thisVert->limitloc.y(), thisVert->limitloc.z(),
									thisVert->limitnorm.x(), thisVert->limitnorm.y(), thisVert->limitnorm.x(), 
									thisVert->rad, thisVert->uvt.u, thisVert->uvt.v, thisVert->uvt.t);
			vTag[vOffsetT+vOffset4]=3;
			vOffset4++;
		}
	}

	// keep all the remaining patches of the main-figure first
	for(fDex=0; fDex<fNums; fDex++)
	{
		CCPolygon *thisPoly=thisMesh->FaceList(fDex);
		minU=999;
		double maxU=-1;
		if(!polyDexList1[fDex])
		{
			fList[fOffsetT+fOffset4].Create(thisPoly->numsides());
			for(fvdex=0; fvdex<thisPoly->numsides(); fvdex++)
			{
				fList[fOffsetT+fOffset4][fvdex]=reDex[thisPoly->Pvertices[fvdex]];
			}
			fTag[fOffsetT+fOffset4]=4;
			fOffset4++;
		}
	}

	vOffsetT+=vOffset4;
	fOffsetT+=fOffset4;

#ifdef QIONG_BLENDING_MAIN_FIGURE
	color[0]=1;
	color[1]=0;
	color[2]=1;
	glRender(vOffsetT, fOffset4, vList, &fList[fOffset+fOffset2+fOffset3], rMode, color);
#endif

	// generate the cut curve on the main figure, stored in mainInterCurve[0...mainInterPNum]
	for(i=mCurveStartDex; i<=mCurveEndDex; i++)
	//for(i=0; i<=mainInterPNum; i++)
	{
		thisVert=thisMesh->VertList(mainInterCurveDex[i]);

		// this is in fact NOT necessary
		//pListOld[0]->subdivBPosition(&tempPnt, thisVert->uvt.u, thisVert->uvt.v, thisVert->uvt.t);

		//		to assign a Bpoint object to a Bpoint2 object, only the values of
		//		pnt/norm/rad are tranferred so far!

		mainInterCurve[i-mCurveStartDex].pnt[0]=thisVert->limitloc.x();
		mainInterCurve[i-mCurveStartDex].pnt[1]=thisVert->limitloc.y();
		mainInterCurve[i-mCurveStartDex].pnt[2]=thisVert->limitloc.z();
		//mainInterCurve[i].norm[0]=thisVert->limitnorm.x();
		//mainInterCurve[i].norm[1]=thisVert->limitnorm.y();
		//mainInterCurve[i].norm[2]=thisVert->limitnorm.z();
		//mainInterCurve[i].rad=thisVert->rad;
		//mainInterCurve[i].=tempPnt;
		//mainInterCurve[i].u=thisVert->uvt.u;
		//mainInterCurve[i].v=thisVert->uvt.v;
		//mainInterCurve[i].t=thisVert->uvt.t;
	}
	/*
	Bpoint tempPnt;
	for(i=0; i<=mainInterPNum; i++)
	{
		thisVert=thisMesh->VertList(mainInterCurveDex[i]);
		//pListOld[0]->subdivBPosition(&mainInterCurve[i], thisVert->uvt.u, thisVert->uvt.v, thisVert->uvt.t);
		pListOld[0]->subdivBPosition(&tempPnt, thisVert->uvt.u, thisVert->uvt.v, thisVert->uvt.t);

		//		to assign a Bpoint object to a Bpoint2 object, only the values of
		//		pnt/norm/rad are tranferred so far!
		mainInterCurve[i]=tempPnt;
		//mainInterCurve[i].u=thisVert->uvt.u;
		//mainInterCurve[i].v=thisVert->uvt.v;
		//mainInterCurve[i].t=thisVert->uvt.t;
	}
	//*/

#ifdef QIONG_BLENDING_MAIN_FIGURE_INTERSECTION_CURVE
	for(i=mCurveStartDex; i<=mCurveEndDex; i++) //mainInterPNum; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "%d\t", mainInterCurveDex[i]);
	#else
		fprintf(stderr, "%d\t", mainInterCurveDex[i]);
	#endif
	}
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "\n");
	#else
		fprintf(stderr, "\n");
	#endif

	color[0]=0;
	color[1]=1;
	color[2]=0;
	glRenderCurve(0, mainInterPNum, mainInterCurve, color);
#endif

//*
// 5.0 connect the rest of main-figure to the intersection curve
	int vOffset5=0, fOffset5=0;

	//fprintf(stderr, "main-figure now!\n", interPNum);

	// connect the 2 curves in mainInterCurve[mainInterPNum] and interCurveD[interPNum]
//	bool goWithMain;
//	bool mDone=false, iDone=false;
//	bool done=false;
//	int triNum=0, triList[MAX_INTERSECTION_NUMBER];
	memset(triList, 0, MAX_INTERSECTION_NUMBER*sizeof(int));
	int mainDex, mainNext;
//	int mainNext, interNext;
//	int mainPicked, interPicked;
	//int interStartPos=0;
	interStartPos=0;

	MeshCurves(&triNum, triList, &interStartPos, mainInterPNum, mainInterCurve, interPNum, interCurveD);
	mainDex=0;
	interDex=interStartPos%interPNum;
	for(i=0; i<triNum; i++)
	{
		switch(triList[i])
		{
		case 1:
			mainNext=(mainDex+1)%mainInterPNum;
			fList[fOffsetT+fOffset5].Create(3);
			fList[fOffsetT+fOffset5][0]=reDex[mainInterCurveDex[mainNext]];
			fList[fOffsetT+fOffset5][1]=vOffset+interDex%interPNum;			//reDex[subVertOldDex[0]];
			fList[fOffsetT+fOffset5][2]=reDex[mainInterCurveDex[mainDex]];
			fTag[fOffsetT+fOffset5]=3;
			fOffset5++;
			mainDex=mainNext;
			break;
		case 2:
			interNext=(interDex+1)%interPNum;
			fList[fOffsetT+fOffset5].Create(3);
			fList[fOffsetT+fOffset5][0]=vOffset+interNext;					//reDex[subVertOldDex[0]];
			fList[fOffsetT+fOffset5][1]=vOffset+interDex;					//reDex[subVertOldDex[0]];
			fList[fOffsetT+fOffset5][2]=reDex[mainInterCurveDex[mainDex]];
			fTag[fOffsetT+fOffset5]=3;
			fOffset5++;
			interDex=interNext;
			break;
		default:
			fprintf(stderr, "\nUNEXPECTED EXCEPTION FROM MeshCurves(...)!\n");
			return false;
		}
	}

	vOffsetT+=vOffset5;
	fOffsetT+=fOffset5;

#ifdef QIONG_BLENDING_INTERSECTION_MAIN_FIGURE
	color[0]=0;
	color[1]=1;
	color[2]=1;
	glRender(vOffsetT, fOffset5, vList, &fList[fOffset+fOffset2+fOffset3+fOffset4], rMode, color);
	//glRender(vOffsetT, 1, vList, &fList[fOffset+fOffset2+fOffset3+fOffset4+fOffset5/2+9], FILLEDPOLY, color);
#endif

	//for(i=0; i<vOffsetT-vOffset5-vOffset4; i++)

#ifdef QIONG_DEBUG_20030510
	color[0]=1;
	color[1]=1;
	color[2]=1;
	//for(i=55; i<65; i++)
	//	glRenderPoint(vList[59], color);
	bool trap1=true;
#endif

//*
#ifdef QIONG_BLENDING_NO_DEBUG
//	fprintf(stderr, "Time to calculate the blended OBJECT......");
	//ThallCode::Pointlist_server2 *pList=new Pointlist_server2();
//	if(blendedPList!=NULL)
//		delete blendedPList;
//	blendedPList=new Pointlist_server2();
//	blendedPList->InitializeBlendedSubdivSurf_New(SUBDIV_LEVEL_TO_USE, vOffsetT, fOffsetT, vList, fList, vTag, fTag, pListOld, IS_IT_AN_INDENTATION_MODEL);
	InitializeBlendedSubdivSurf_New(level, SUBDIV_LEVEL_TO_USE, vOffsetT, fOffsetT, vList, fList, vTag, fTag);
	cout << '.' << flush;
#endif
//*/

#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fclose(fDebug);
#endif

#ifdef DEBUG
	cout << "Tree " << treeId << " yielded " << vertexLists[treeId].size() << " verticies" << endl;
#endif

	if(polyDexList1!=NULL)
		delete []polyDexList1;
	if(polyDexList2!=NULL)
		delete []polyDexList2;

	if(reDex!=NULL)
		delete []reDex;
	if(vList!=NULL)
		delete []vList;
	if(fList!=NULL)
		delete []fList;
	if(vTag!=NULL)
		delete []vTag;
	if(fTag!=NULL)
		delete []fTag;

	if(quadListMain!=NULL)
		delete []quadListMain;
	if(vertListMainQuad!=NULL)
		delete []vertListMainQuad;
	if(vertListMainQuadCnt!=NULL)
		delete []vertListMainQuadCnt;
	if(edgeListMainQuad!=NULL)
		delete []edgeListMainQuad;
	if(vertListMain!=NULL)
		delete []vertListMain;

	if(edgeListMainQuadDex!=NULL)
		delete []edgeListMainQuadDex;

	if(takenMainQuadEdge!=NULL)
		delete []takenMainQuadEdge;
	return true;
#else
	return ComputeBlending_temp(treeId);
#endif
}

bool Pointlist_serverB::ComputeBlending_temp(int treeId)
{
#ifdef PARAMETERS_FOR_TORUS_OBJECT
	if(rMode==WIREFRAME)
		glLineWidth(1.5);

	MyFloat SUB_FIGURE_CUTOFF_U=SUB_FIGURE_CUTOFF_U_ORIGINAL+SUB_FIGURE_CUTOFF_U_RANGE*blendExtents[treeId];
	MyFloat DILATION_STEPS=DILATION_STEPS_ORIGINAL+DILATION_STEPS_RANGE*blendAmounts[treeId];

	MyFloat SUB_FIGURE_CUTOFF_U_2=SUB_FIGURE_CUTOFF_U_ORIGINAL_2+SUB_FIGURE_CUTOFF_U_RANGE_2*blendExtents[treeId];
	MyFloat DILATION_STEPS_2=DILATION_STEPS_ORIGINAL_2+DILATION_STEPS_RANGE_2*blendAmounts[treeId];

	// blending amount testing parameters
	//MyFloat SUB_FIGURE_CUTOFF_U=SUB_FIGURE_CUTOFF_U_ORIGINAL+SUB_FIGURE_CUTOFF_U_RANGE*subBlending;
	//MyFloat DILATION_STEPS=DILATION_STEPS_ORIGINAL+DILATION_STEPS_RANGE*mainBlending;

	// Since the parameters controling the blending amount has been implemented,
	// I will use them instead. Their names should be changed accordingly later.
	//MyFloat SUB_FIGURE_CUTOFF_U=SUB_FIGURE_CUTOFF_U_ORIGINAL+SUB_FIGURE_CUTOFF_U_RANGE*(blendExtent==0?.01:blendExtent);
	//MyFloat DILATION_STEPS=DILATION_STEPS_ORIGINAL+DILATION_STEPS_RANGE*(blendAmount==0?.01:blendAmount);
	MyFloat EPSILON2=EPSILON2_ORIGINAL;

	int i, j, k;

#ifndef QIONG_BLENDING_NO_DEBUG
	MyFloat color[3];
#endif

	int interPNum[2];
    //int ret;
	//int taken[2][MAX_INTERSECTION_NUMBER];
	int interPPair[2][MAX_INTERSECTION_NUMBER*2];
	int interPPair2[2][MAX_INTERSECTION_NUMBER*2];
	Bpoint2 interPList[2][MAX_INTERSECTION_NUMBER];
	Bpoint2 interPList2[2][MAX_INTERSECTION_NUMBER];
	Bpoint2 interCurve[2][MAX_INTERSECTION_NUMBER];
	Bpoint2 interCurveD[2][MAX_INTERSECTION_NUMBER];

	int **polyDexList1, **polyDexList2;
	polyDexList1=new int*[2];
	polyDexList2=new int*[2];
	//int *polyDexList1=new int[nTiles1], *polyDexList2=new int[nTiles2];
	//int polyDexList1[MAX_INTERSECTION_NUMBER];
	//int polyDexList2[MAX_INTERSECTION_NUMBER];

	int nTiles1, nTiles2;
	Bpoint *tile1, *tile2;
    list<M3DFigureTreeNode *> children;
	ThallCode::Pointlist_server2 *pListOld[2];

#ifdef DEBUG
	cout << "Pointlist_serverB::renderFigure() called with treeId = " << treeId << '\n';
#endif

#ifdef QIONG_DEBUG_FILE_TO_WRITE
	//FILE *fDebug;
	char fName[200];
	sprintf(fName, "%s", "debug.dat");
	fDebug=fopen(fName, "w+");
#endif

	interPNum[0]=0;						// number of the vertices on the intersection curve
	interPNum[1]=0;

	MyFloat minU; //, minU2;				// ??
	MyFloat maxU; //, maxU2;

	pListOld[0]=pListOriginal[0]; //tileLists[0];				// Start with the main figure to fill the tile information
	pListOld[1]=pListOriginal[1]; //tileLists[1];				// Need the sub-figure to calculate the intersection
 	pListOld[0]->subdivtileinfo(&nTiles1, &tile1);
	pListOld[1]->subdivtileinfo(&nTiles2, &tile2);

	for(i=0; i<2; i++)
	{
		polyDexList1[i]=new int[nTiles1];
		polyDexList2[i]=new int[nTiles2];
		memset(polyDexList1[i], 0, nTiles1*sizeof(int));
		memset(polyDexList2[i], 0, nTiles2*sizeof(int));
	}

#ifndef QIONG_BLENDING_NO_DEBUG
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nThe re-indexing information of level %d sub-figure quads are removed.\n", SUBDIV_LEVEL_TO_USE);
	#else
	fprintf(stderr, "\nThe re-indexing information of level %d sub-figure quads are removed.\n", SUBDIV_LEVEL_TO_USE);
	#endif
#endif

/// .1.
	// To intersect the two meshes here for the first time
	IntersectMeshes(SUB_FIGURE_CUTOFF_U, nTiles1, tile1, nTiles2, tile2, interPPair[0], interPList[0], \
		polyDexList1[0], polyDexList2[0], &interPNum[0]);

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	for(i=0; i<interPNum[0]; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "Intersection point %d: u/v/t in main figure %f/%f/%f\n", i, interPList[0][i].u, interPList[0][i].v, interPList[0][i].t);
	#else
		fprintf(stderr, "Intersection point %d: u/v/t in main figure %f/%f/%f\n", i, interPList[0][i].u, interPList[0][i].v, interPList[0][i].t);
	#endif
	}
#endif

/// .2.
	// To intersect the two meshes here for the second time
	//int firstLoop=1;
	IntersectMeshes(-SUB_FIGURE_CUTOFF_U_2, nTiles1, tile1, nTiles2, tile2, interPPair[1], interPList[1], \
		polyDexList1[1], polyDexList2[1], &interPNum[1]);

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	for(i=0; i<interPNum[1]; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "Intersection point %d: u/v/t in main figure %f/%f/%f\n", i, interPList[1][i].u, interPList[1][i].v, interPList[1][i].t);
	#else
		fprintf(stderr, "Intersection point %d: u/v/t in main figure %f/%f/%f\n", i, interPList[1][i].u, interPList[1][i].v, interPList[1][i].t);
	#endif
	}
#endif

//#ifndef QIONG_BLENDING_NO_DEBUG
	if(interPNum[0]>MAX_INTERSECTION_NUMBER || interPNum[1]>MAX_INTERSECTION_NUMBER || \
		interPNum[0]<=0 || interPNum[1]<=0)
	{
		fprintf(stderr, "TOO MANY or NONE INTERSECTION POINTS found from the TRI-TRI intersections...\n");

		for(i=0; i<2; i++)
		{
			if(polyDexList1[i]!=NULL)
				delete [](polyDexList1[i]);
			if(polyDexList2[i]!=NULL)
				delete [](polyDexList2[i]);
		}
		delete []polyDexList1;
		delete []polyDexList2;

		return false;
	}
//#endif

#ifndef QIONG_BLENDING_NO_DEBUG
	if(interPNum[0]>MAX_INTERSECTION_NUMBER || interPNum[1]>MAX_INTERSECTION_NUMBER)
	{
		fprintf(stderr, "TOO MANY INTERSECTION POINTS TO HANDLE FROM the TRI-TRI intersections...\n");

		for(i=0; i<2; i++)
		{
			if(polyDexList1[i]!=NULL)
				delete [](polyDexList1[i]);
			if(polyDexList2[i]!=NULL)
				delete [](polyDexList2[i]);
		}
		delete []polyDexList1;
		delete []polyDexList2;

		return false;
	}
#endif

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_ORIGINAL
	color[0]=0;
	color[1]=1;
	color[2]=0;
	glRenderSegments(0, interPNum[0]-1, interPList[0], color);
	glRenderSegments(0, interPNum[1]-1, interPList[1], color);
#endif

#ifdef QIONG_BLENDING_DILATION_DEBUG
	for(i=0; i<interPNum[0]; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "seg %d: %f %f %f, %f %f %f\n", i, \
			interPList[0][2*i].pnt[0], interPList[0][2*i].pnt[1], interPList[0][2*i].pnt[2], \
			interPList[0][2*i+1].pnt[0], interPList[0][2*i+1].pnt[1], interPList[0][2*i+1].pnt[2]);
	#else
		fprintf(stderr, "seg %d: %f %f %f, %f %f %f\n", i, \
			interPList[0][2*i].pnt[0], interPList[0][2*i].pnt[1], interPList[0][2*i].pnt[2], \
			interPList[0][2*i+1].pnt[0], interPList[0][2*i+1].pnt[1], interPList[0][2*i+1].pnt[2]);
	#endif
	}
	for(i=0; i<interPNum[1]; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "seg %d: %f %f %f, %f %f %f\n", i, \
			interPList[1][2*i].pnt[0], interPList[1][2*i].pnt[1], interPList[1][2*i].pnt[2], \
			interPList[1][2*i+1].pnt[0], interPList[1][2*i+1].pnt[1], interPList[1][2*i+1].pnt[2]);
	#else
		fprintf(stderr, "seg %d: %f %f %f, %f %f %f\n", i, \
			interPList[1][2*i].pnt[0], interPList[1][2*i].pnt[1], interPList[1][2*i].pnt[2], \
			interPList[1][2*i+1].pnt[0], interPList[1][2*i+1].pnt[1], interPList[1][2*i+1].pnt[2]);
	#endif
	}
#endif

	// at this stage, all the intersection vertices are found!

/// .3. try to connect the first intersection curve
/// .4. try to connect the second intersection curve
	int cDex;

	bool successful;

for(cDex=0; cDex<2; cDex++)
{
	successful=FindIntersectionCurveFromSegments(interPList2[cDex], interPPair2[cDex], &interPNum[cDex], \
		interPList[cDex], interPPair[cDex], EPSILON2);

	if(!successful)
	{
		fprintf(stderr, "\nThe intersection curve is not closed!!!\nHave to quit right here...\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fclose(fDebug);
#endif
		//if(polyReDex2!=NULL)
		//	delete []polyReDex2;

		if(polyDexList1[cDex]!=NULL)
			delete [](polyDexList1[0]);
		if(polyDexList2[cDex]!=NULL)
			delete [](polyDexList2[0]);

		if(polyDexList1[1]!=NULL)
			delete [](polyDexList1[1]);
		if(polyDexList2[1]!=NULL)
			delete [](polyDexList2[1]);

		delete []polyDexList1;
		delete []polyDexList2;

		return false;
	}
}


#ifdef QIONG_BLENDING_INTERSECTION_QUAD
	color[0]=1;
	color[1]=0;
	color[2]=0;
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "The list of main-figure quads that was intersected by sub-figure (1)\n  or intersected by the dilated i-curve (2)...\n");
	#endif

if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
{
	for(i=0; i<nTiles1; i++)
	{
		if(polyDexList1[0][i])
		{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			if(i%8==0)
				fprintf(fDebug, "\n");
			fprintf(fDebug, "case %d:%d\t", polyDexList1[0][i], i);
	#endif
		}
	}
	for(i=0; i<nTiles1; i++)
	{
		if(polyDexList1[1][i])
		{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			if(i%8==0)
				fprintf(fDebug, "\n");
			fprintf(fDebug, "case %d:%d\t", polyDexList1[1][i], i);
	#endif
			switch(surfaceStyle)
			{
			case M3D_NONE:			
				break;
			case M3D_POINTCLOUD:
				glRender(SUBDIV_LEVEL_TO_USE+1, i, VPOINTS, color);
				break;
			case M3D_WIREFRAME:
				glRender(SUBDIV_LEVEL_TO_USE+1, i, WIREFRAME, color);
				break;
			case M3D_SOLID:
				glRender(SUBDIV_LEVEL_TO_USE+1, i, FILLEDPOLY, color);
				break;
			}
		}
	}
}
else
{
	for(i=0; i<nTiles1; i++)
	{
		if(polyDexList1[0][i])
		{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			if(i%8==0)
				fprintf(fDebug, "\n");
			fprintf(fDebug, "case %d:%d\t", polyDexList1[0][i], i);
	#endif
			switch(surfaceStyle)
			{
			case M3D_NONE:			
				break;
			case M3D_POINTCLOUD:
				glRender(SUBDIV_LEVEL_TO_USE, i, VPOINTS, color);
				break;
			case M3D_WIREFRAME:
				glRender(SUBDIV_LEVEL_TO_USE, i, WIREFRAME, color);
				break;
			case M3D_SOLID:
				glRender(SUBDIV_LEVEL_TO_USE, i, FILLEDPOLY, color);
				break;
			}
		}
	}
	for(i=0; i<nTiles1; i++)
	{
		if(polyDexList1[1][i])
		{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			if(i%8==0)
				fprintf(fDebug, "\n");
			fprintf(fDebug, "case %d:%d\t", polyDexList1[1][i], i);
	#endif
			switch(surfaceStyle)
			{
			case M3D_NONE:			
				break;
			case M3D_POINTCLOUD:
				glRender(SUBDIV_LEVEL_TO_USE, i, VPOINTS, color);
				break;
			case M3D_WIREFRAME:
				glRender(SUBDIV_LEVEL_TO_USE, i, WIREFRAME, color);
				break;
			case M3D_SOLID:
				glRender(SUBDIV_LEVEL_TO_USE, i, FILLEDPOLY, color);
				break;
			}
		}
	}
}

	color[0]=0;
	color[1]=1;
	color[2]=0;
	for(i=0; i<nTiles2; i++)
	{
		if(polyDexList2[0][i])
			switch(surfaceStyle)
			{
			case M3D_NONE:			
				break;
			case M3D_POINTCLOUD:
				pListOld[1]->glRender(SUBDIV_LEVEL_TO_USE, i, VPOINTS, color);
				break;
			case M3D_WIREFRAME:
				pListOld[1]->glRender(SUBDIV_LEVEL_TO_USE, i, WIREFRAME, color);
				break;
			case M3D_SOLID:
				pListOld[1]->glRender(SUBDIV_LEVEL_TO_USE, i, FILLEDPOLY, color);
				break;
			}
		if(polyDexList2[1][i])
			switch(surfaceStyle)
			{
			case M3D_NONE:			
				break;
			case M3D_POINTCLOUD:
				pListOld[1]->glRender(SUBDIV_LEVEL_TO_USE, i, VPOINTS, color);
				break;
			case M3D_WIREFRAME:
				pListOld[1]->glRender(SUBDIV_LEVEL_TO_USE, i, WIREFRAME, color);
				break;
			case M3D_SOLID:
				pListOld[1]->glRender(SUBDIV_LEVEL_TO_USE, i, FILLEDPOLY, color);
				break;
			}
	}
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\n");
	#endif
#endif	

/// .5.
	//	make sure the intersection curve is closed somehow
	//if(interCurve[0].pnt[0]-interCurve[interPNum].pnt[0]<-EPSILON2||interCurve[0].pnt[0]-interCurve[interPNum].pnt[0]>EPSILON2)
	//if(!SamePoint(EPSILON2, interCurve[0], interCurve[interPNum]))
	if(!SamePointByDistanceSquare(EPSILON2, interPList2[0][0], interPList2[0][interPNum[0]]) ||
	   !SamePointByDistanceSquare(EPSILON2, interPList2[1][0], interPList2[1][interPNum[1]])   )
	//if(!SamePointByDistanceSquare(EPSILON2, interCurve[0], interCurve[interPNum]))
	{
		fprintf(stderr, "\nThe intersection curve is not closed!!!\nHave to quit right here...\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fclose(fDebug);
#endif
		for(i=0; i<2; i++)
		{
			if(polyDexList1[i]!=NULL)
				delete polyDexList1[i];
			if(polyDexList2[i]!=NULL)
				delete polyDexList2[i];
		}

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		return false;
	}

#ifdef QIONG_DEBUG_FILE_TO_WRITE
for(cDex=0; cDex<2; cDex++)
{
	fprintf(fDebug, "\nThe main-figure u/v/t of all the %d vertices on the original i-curve %d...\n", interPNum[cDex], cDex);
	for(i=0; i<=interPNum[cDex]; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f/%f\t", i, interPList2[cDex][i].u, interPList2[cDex][i].v, interPList2[cDex][i].t);
	}
	fprintf(fDebug, "\nThe sub-figure u/v/t of all the %d vertices on the original i-curve %d...\n", interPNum[cDex], cDex);
	for(i=0; i<=interPNum[cDex]; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f/%f\t", i, interPList2[cDex][i].u_2, interPList2[cDex][i].v_2, interPList2[cDex][i].t_2);
	}
	fprintf(fDebug, "\n");
}
#endif

	//		now the task is to trim all the unnecessary segments including: too short and colinear segments
	//		the results are stored in 'interCurve and interPPair'
	//		interPPair[2*i]   stores the information of which quad the current vertex i is on main-figure
	//		interPPair[2*i+1] stores the information of which quad the current segment [i, i+1] is on sub-figure
	//RegularizeIntersectionCurve(pListOld[1]->max_v(), interCurve, interPList2, interPPair, interPPair2, &interPNum);

for(cDex=0; cDex<2; cDex++)
{
	RegularizeIntersectionCurve(interCurve[cDex], interPList2[cDex], interPPair[cDex], interPPair2[cDex], &interPNum[cDex]);

	if(interPNum[cDex]>MAX_INTERSECTION_NUMBER || interPNum[cDex]<=0)
	{
		fprintf(stderr, "NONE INTERSECTION POINTS left after the REGULARIZATION!\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fclose(fDebug);
#endif
		for(i=0; i<2; i++)
		{
			if(polyDexList1[i]!=NULL)
				delete polyDexList1[i];
			if(polyDexList2[i]!=NULL)
				delete polyDexList2[i];
		}

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		return false;
	}
#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nThe v/t of all the %d vertices on the regularized i-curve %d...\n", interPNum, cDex);
	for(i=0; i<=interPNum[cDex]; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f\t", i, interCurve[cDex][i].v_2, interCurve[cDex][i].t_2);
	}
	fprintf(fDebug, "\n");
#endif

	//		the numerical problem has to be dealt with somewhere
	//		after the interpolation of the barycentric coordinates of the intersection curve, those UVT.Ts which
	//		were supposed to be 1/-1 or in the crest the UVT.UVs to be 0, uMax, 0, vMax are not exactly these
	//		values, then we have to PUT THEM BACK to the values they are assumed to be!

	RegularizeIntersectionCurveUVT(pListOld[0]->max_u(), pListOld[0]->max_v(), \
		pListOld[1]->max_u(), pListOld[1]->max_v(),	interCurve[cDex], interPNum[cDex]);
#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nThe v/t of all the %d vertices on the regularized UVT i-curve...\n", interPNum);
	for(i=0; i<=interPNum[cDex]; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f\t", i, interCurve[cDex][i].v_2, interCurve[cDex][i].t_2);
	}
	fprintf(fDebug, "\n");
#endif

	//		make sure the intersection curve is closed somehow

	// all the UVTRad information should have been done when calculating the intersections!
	// Need to figure out the UVT coordinates for the intersection points on the main-figure
	//		the implementation in the Pointlist_server2 is an approximation

	// this is not necessary any longer since interCurve[interPNum] is not used anyway!

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_REGULARIZED
	color[0]=1;
	color[1]=0;
	color[2]=0;
	glRenderCurve(0, interPNum[cDex], interCurve[cDex], color);
#endif

	// use the 'interCurve' as temporary storage, the final dilated intersection curve should still be
	// stored in interCurveD

#ifdef QIONG_DEBUG_20030510
	// try to check all the t coordinates of all the intersection vertices
	// they are supposed to be 1/-1 if uv belongs to (0, uMax) and (0, vMax) respectively
	for(i=0; i<=interPNum[cDex]; i++)
		//fprintf(stderr, "%8.7f\t%d\t", interCurveD[i].t, interCurveD[i].t==-1.0);
		fprintf(stderr, "%d\t", interCurveD[cDex][i].t==-1.0);
	fprintf(stderr, "\n");
#endif
}


/// 0.0. time to re-mesh the entire surface boundary!
	CCMesh *thisMesh, *thisMesh1, *thisMesh2;
	CCVertex *thisVert;//, *thatVert;
	Bpoint tempVert;

	CCPolygon *thisPoly;
	int pSides;

	if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
		thisMesh1=(pListOld[0]->GetThisSurfPointer())->SubMeshPointer(SUBDIV_LEVEL_TO_USE+1);
	else
		thisMesh1=(pListOld[0]->GetThisSurfPointer())->SubMeshPointer(SUBDIV_LEVEL_TO_USE);

	thisMesh2=(pListOld[1]->GetThisSurfPointer())->SubMeshPointer(SUBDIV_LEVEL_TO_USE);

/// 1.0. cut the subfigure!

	// for sub-figure, re-relabling vertices
	int vNumEstimate=thisMesh1->VertList.Length();
	int fNumEstimate=thisMesh1->FaceList.Length();
	vNumEstimate+=thisMesh2->VertList.Length();
	fNumEstimate+=thisMesh2->FaceList.Length();
	vNumEstimate*=1.5;
	fNumEstimate*=1.5;

	thisMesh=thisMesh2;
	int vNums=thisMesh->VertList.Length();
	int fNums=thisMesh->FaceList.Length();

	int vOffsetT=0, fOffsetT=0;
	int vOffset=0, fOffset=0;			// for remained and cut vertices and patches on sub-figure

	int **reDex=new int*[2];
	reDex[0]=new int[vNumEstimate];
	reDex[1]=new int[vNumEstimate];

	int reduction=0;
	Bpoint2 *vList=new Bpoint2[vNumEstimate];
	PindexList *fList=new PindexList[fNumEstimate];
	int *vTag=new int[vNumEstimate];
	int *fTag=new int[fNumEstimate];
	memset(vTag, -1, vNumEstimate*sizeof(int));
	memset(fTag, -1, fNumEstimate*sizeof(int));
	for(i=0; i<vNumEstimate; i++) reDex[0][i]=i;

	// reindex the vertices for the sub-figure first
	for(int vDex=0; vDex<vNums; vDex++)
	{
		// reDex is the important array to re-index all the orginal sub-figure vertices to the new array in vList
		reDex[0][vDex]-=reduction;
		thisVert=thisMesh->VertList(vDex);
		//		the codes picking vertices to remove are not CORRECT!
		if(thisVert->uvt.u<SUB_FIGURE_CUTOFF_U)
		{
			//if(thisVert->uvt.u-SUB_FIGURE_CUTOFF_U<-1)// || (thisVert->uvt.u<SUB_FIGURE_CUTOFF_U&&thisVert->uvt.t==0))
			if( (thisVert->uvt.u-SUB_FIGURE_CUTOFF_U<=-UVT_UNIT) || 
				(thisVert->uvt.u==0 && thisVert->uvt.t>-1 && thisVert->uvt.t<1 && thisVert->uvt.v>0 && thisVert->uvt.v<pListOld[1]->max_v()) )
			    //(thisVert->uvt.u-SUB_FIGURE_CUTOFF_U>-1.0/(SUBDIV_LEVEL_TO_USE+1) && (thisVert->uvt.t>-1 && thisVert->uvt.t<1) && thisVert->uvt.v>0 && thisVert->uvt.v<pListOld[1]->max_v()) )
			{
				reDex[0][vDex]=-1;
				reduction++;
				//vOffset--;
			}
			else
			{
				// need to cut the sub-figure, generate new vertices here with new U
				//double t=SUB_FIGURE_CUTOFF_U-thisVert->uvt.u;
				tempVert.u=SUB_FIGURE_CUTOFF_U+SUB_FIGURE_CUTOFF_TRANSFORM;
				tempVert.v=thisVert->uvt.v;
				tempVert.t=thisVert->uvt.t;
				pListOld[1]->subdivBPosition(&tempVert, tempVert.u, tempVert.v, tempVert.t);

				// remove garbage later! 110802
				vList[vOffset].setvals(tempVert.pnt[0], tempVert.pnt[1], tempVert.pnt[2],
									   tempVert.norm[0], tempVert.norm[1], tempVert.norm[2], 
									   tempVert.rad, tempVert.u, tempVert.v, tempVert.t);
				vTag[vOffset]=1;
				vOffset++;
			}
		}
		else 
			if(thisVert->uvt.u>SUB_FIGURE_CUTOFF_U_2)
			{
				if( (thisVert->uvt.u-SUB_FIGURE_CUTOFF_U_2>=UVT_UNIT) || 
					(thisVert->uvt.u==pListOld[1]->max_u() && thisVert->uvt.t>-1 && thisVert->uvt.t<1 && thisVert->uvt.v>0 && thisVert->uvt.v<pListOld[1]->max_v()) )
				{
					reDex[0][vDex]=-1;
					reduction++;
				}
				else
				{
					// need to cut the sub-figure, generate new vertices here with new U
					//double t=SUB_FIGURE_CUTOFF_U-thisVert->uvt.u;
					tempVert.u=SUB_FIGURE_CUTOFF_U_2+SUB_FIGURE_CUTOFF_TRANSFORM_2;
					tempVert.v=thisVert->uvt.v;
					tempVert.t=thisVert->uvt.t;
					pListOld[1]->subdivBPosition(&tempVert, tempVert.u, tempVert.v, tempVert.t);

					// remove garbage later! 110802
					vList[vOffset].setvals(tempVert.pnt[0], tempVert.pnt[1], tempVert.pnt[2],
						tempVert.norm[0], tempVert.norm[1], tempVert.norm[2], 
						tempVert.rad, tempVert.u, tempVert.v, tempVert.t);
					vTag[vOffset]=1;
					vOffset++;
				}
			}
			else
			{
				// need to lift up this curve to guarantee the cut sub-figure is more uniform
				if(thisVert->uvt.u-SUB_FIGURE_CUTOFF_U<UVT_UNIT)
				{
					Bpoint bp;
					pListOld[1]->subdivBPosition(&bp, thisVert->uvt.u+SUB_FIGURE_CUTOFF_TRANSFORM, thisVert->uvt.v, thisVert->uvt.t);
					vList[vOffset].setvals(bp.pnt[0], bp.pnt[1], bp.pnt[2],
						bp.norm[0], bp.norm[1], bp.norm[2], 
						bp.rad, thisVert->uvt.u+SUB_FIGURE_CUTOFF_TRANSFORM, thisVert->uvt.v, thisVert->uvt.t);
				}
				else
					if(SUB_FIGURE_CUTOFF_U_2-thisVert->uvt.u<UVT_UNIT)
					{
						Bpoint bp;
						pListOld[1]->subdivBPosition(&bp, thisVert->uvt.u+SUB_FIGURE_CUTOFF_TRANSFORM_2, thisVert->uvt.v, thisVert->uvt.t);
						vList[vOffset].setvals(bp.pnt[0], bp.pnt[1], bp.pnt[2],
							bp.norm[0], bp.norm[1], bp.norm[2], 
							bp.rad, thisVert->uvt.u+SUB_FIGURE_CUTOFF_TRANSFORM_2, thisVert->uvt.v, thisVert->uvt.t);
					}
					else
					{
						vList[vOffset].setvals(thisVert->limitloc.x(), thisVert->limitloc.y(), thisVert->limitloc.z(),
							thisVert->limitnorm.x(), thisVert->limitnorm.y(), thisVert->limitnorm.x(), 
							thisVert->rad, thisVert->uvt.u, thisVert->uvt.v, thisVert->uvt.t);
					}

				vTag[vOffset]=0;
				vOffset++;
		}
	}

/// 1.5. prepare to find the cut-curve(s)
	// working on a more rubust algorithm to find the cut-curve on the sub-figure
	// keep all the remaining patches of the sub-figure first

	// all the variables to calculate the cut-curve on the sub-figure
	int quadListSub[2][MAX_INTERSECTION_NUMBER*2], quadListSubNum[2];
	int edgeListSubQuad[2][MAX_INTERSECTION_NUMBER*2], edgeListSubQuadNum[2];
	int takenSubQuadEdge[2][MAX_INTERSECTION_NUMBER];

for(cDex=0; cDex<2; cDex++)
{
	memset(takenSubQuadEdge[cDex], 0, sizeof(int)*MAX_INTERSECTION_NUMBER);
	quadListSubNum[cDex]=0;
	edgeListSubQuadNum[cDex]=0;
}

	for(int fDex=0; fDex<fNums; fDex++)
	{
		thisPoly=thisMesh->FaceList(fDex);
		minU=999;
		maxU=-1;
		for(int fvdex=0; fvdex<thisPoly->numsides(); fvdex++)
		{
			thisVert=thisMesh->VertList(thisPoly->Pvertices[fvdex]);
			if(minU>thisVert->uvt.u) minU=thisVert->uvt.u;
			if(maxU<thisVert->uvt.u) maxU=thisVert->uvt.u;
		}
		if(minU>SUB_FIGURE_CUTOFF_U && maxU<SUB_FIGURE_CUTOFF_U_2)
		{
			fList[fOffset].Create(thisPoly->numsides());
			for(fvdex=0; fvdex<thisPoly->numsides(); fvdex++)
			{
				fList[fOffset][fvdex]=reDex[0][thisPoly->Pvertices[fvdex]];
			}
			fTag[fOffset]=0;
			fOffset++;
		}
		else
		{
			if(minU<=SUB_FIGURE_CUTOFF_U && maxU>=SUB_FIGURE_CUTOFF_U)
			{
				// 1. re-write the codes to connect the cut-curve on the subfigure and the dilated-i-curve
				quadListSub[0][quadListSubNum[0]++]=fDex;

				fList[fOffset].Create(thisPoly->numsides());
				for(fvdex=0; fvdex<thisPoly->numsides(); fvdex++)
				{
					fList[fOffset][fvdex]=reDex[0][thisPoly->Pvertices[fvdex]];
				}
				fTag[fOffset]=1;
				fOffset++;
			}
			else
				if(maxU>=SUB_FIGURE_CUTOFF_U_2 && minU<=SUB_FIGURE_CUTOFF_U_2)
				{
					// 1. re-write the codes to connect the cut-curve on the subfigure and the dilated-i-curve
					quadListSub[1][quadListSubNum[1]++]=fDex;

					fList[fOffset].Create(thisPoly->numsides());
					for(fvdex=0; fvdex<thisPoly->numsides(); fvdex++)
					{
						fList[fOffset][fvdex]=reDex[0][thisPoly->Pvertices[fvdex]];
					}
					fTag[fOffset]=1;
					fOffset++;
				}
		}
	}
	vOffsetT+=vOffset;
	fOffsetT+=fOffset;

#ifdef QIONG_BLENDING_SUB_FIGURE
	color[0]=0;
	color[1]=0;
	color[2]=1;
	glRender(vOffsetT, fOffsetT, vList, fList, rMode, color);
#endif

	// done: so far vOffset vertices and fOffset faces in sub-figure have been inserted into vList and fList

//*
/// 2.0. insert all the dilated intersection point into the vList!
/// 2.1. before that, need to calculate the cut-curve of the sub-figure
	// now need to find the cut-curve of the subfigure, and connect the dilated intersection curve to it
	// then the sub-figure
	// all the variables to calculate the cut-curve on the sub-figure
	int subInterPNum[2];
	int subInterCurveDex[2][MAX_INTERSECTION_NUMBER];
	Bpoint2 subInterCurve[2][MAX_INTERSECTION_NUMBER];

for(cDex=0; cDex<2; cDex++)
{
	// first find all the edges contributing to the cut-curve on the sub-figure
	for(i=0; i<quadListSubNum[cDex]; i++)
	{
		thisPoly=thisMesh->FaceList(quadListSub[cDex][i]);
		pSides=thisPoly->numsides();
		for(k=0; k<pSides; k++)
		{
			thisVert=thisMesh->VertList(thisPoly->Pvertices[k]);
			switch(cDex)
			{
			case 0:
				if(thisVert->uvt.u<SUB_FIGURE_CUTOFF_U)
				{
					thisVert=thisMesh->VertList(thisPoly->Pvertices[(k+1)%pSides]);
					if(thisVert->uvt.u<SUB_FIGURE_CUTOFF_U)
					{
						edgeListSubQuad[cDex][2*edgeListSubQuadNum[cDex]]=thisPoly->Pvertices[k];
						edgeListSubQuad[cDex][2*edgeListSubQuadNum[cDex]+1]=thisPoly->Pvertices[(k+1)%pSides];
						edgeListSubQuadNum[cDex]++;
					}
				}
				break;
			case 1:
				if(thisVert->uvt.u>SUB_FIGURE_CUTOFF_U_2)
				{
					thisVert=thisMesh->VertList(thisPoly->Pvertices[(k+1)%pSides]);
					if(thisVert->uvt.u>SUB_FIGURE_CUTOFF_U_2)
					{
						edgeListSubQuad[cDex][2*edgeListSubQuadNum[cDex]]=thisPoly->Pvertices[k];
						edgeListSubQuad[cDex][2*edgeListSubQuadNum[cDex]+1]=thisPoly->Pvertices[(k+1)%pSides];
						edgeListSubQuadNum[cDex]++;
					}
				}
				break;
			default:
				break;
			}
		}
	}

	// time to connect the edges to form the cut-curve
	subInterPNum[cDex]=edgeListSubQuadNum[cDex];
	subInterCurveDex[cDex][0]=edgeListSubQuad[cDex][0];
	subInterCurveDex[cDex][1]=edgeListSubQuad[cDex][1];
	takenSubQuadEdge[cDex][0]=1;
	i=1;
	j=1;
	while(i<subInterPNum[cDex])
	{
		while(takenSubQuadEdge[cDex][j]) 
			j++;
		if(j>=edgeListSubQuadNum[cDex])
		{
			// the cut-curve can not be connected at all
			// time to fail the algorithm
			// i.e., this should NOT possibly happen
			fprintf(stderr, "\nThe intersection curve is not closed!!!\nHave to quit right here...\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
			fclose(fDebug);
#endif
			for(i=0; i<2; i++)
			{
				if(polyDexList1[i]!=NULL)
					delete polyDexList1[i];
				if(polyDexList2[i]!=NULL)
					delete polyDexList2[i];
				if(reDex[i]!=NULL)
					delete reDex[i];
			}

			if(polyDexList1!=NULL)
				delete []polyDexList1;
			if(polyDexList2!=NULL)
				delete []polyDexList2;

			if(reDex!=NULL)
				delete []reDex;
			if(vList!=NULL)
				delete []vList;
			if(fList!=NULL)
				delete []fList;
			if(vTag!=NULL)
				delete []vTag;
			if(fTag!=NULL)
				delete []fTag;

			return false;
		}
		if(subInterCurveDex[cDex][i]==edgeListSubQuad[cDex][2*j])
		{
			i++;
			subInterCurveDex[cDex][i]=edgeListSubQuad[cDex][2*j+1];
			takenSubQuadEdge[cDex][j]=1;
			j=0;
		}
		else
			j++;
	}

	if(subInerPNum[cDex]<=0)
	{
		// the cut-curve can not be connected at all
		// time to fail the algorithm
		// i.e., this should NOT possibly happen if the topology is valid, must be something wrong w/ the u/v/t
		fprintf(stderr, "\nNo CUT CURVE found on the sub-figure!\nHave to quit right here...\n");

#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fclose(fDebug);
#endif
		for(i=0; i<2; i++)
		{
			if(polyDexList1[i]!=NULL)
				delete polyDexList1[i];
			if(polyDexList2[i]!=NULL)
				delete polyDexList2[i];
			if(reDex[i]!=NULL)
				delete reDex[i];
		}

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		if(reDex!=NULL)
			delete []reDex;
		if(vList!=NULL)
			delete []vList;
		if(fList!=NULL)
			delete []fList;
		if(vTag!=NULL)
			delete []vTag;
		if(fTag!=NULL)
			delete []fTag;

		return false;
	}

	// for 'cup'-shaped object, need to treat it differently
	// for indentation model, need to reverse this curve!
	//if(IS_IT_AN_INDENTATION_MODEL)
	{
		for(i=0; i<subInterPNum[cDex]/2; i++)
		{
			j=subInterCurveDex[cDex][i];
			subInterCurveDex[cDex][i]=subInterCurveDex[cDex][subInterPNum[cDex]-1-i];
			subInterCurveDex[cDex][subInterPNum[cDex]-1-i]=j;
		}
		subInterCurveDex[cDex][subInterPNum[cDex]]=subInterCurveDex[cDex][0];
	}

	// assign the values for the cut-curve
	for(i=0; i<=subInterPNum[cDex]; i++)
	{
		subInterCurve[cDex][i]=vList[reDex[0][subInterCurveDex[cDex][i]]];
	}

#ifdef QIONG_BLENDING_SUB_FIGURE_INTERSECTION_CURVE
	for(i=0; i<=subInterPNum[cDex]; i++) //mainInterPNum; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "%d\t", subInterCurveDex[cDex][i]);
	#else
		fprintf(stderr, "%d\t", subInterCurveDex[cDex][i]);
	#endif
	}
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "\n");
	#else
		fprintf(stderr, "\n");
	#endif

	color[0]=0;
	color[1]=1;
	color[2]=0;
	glRenderCurve(0, subInterPNum[cDex], subInterCurve[cDex], color);
#endif

	// still, need to check the closeness of the cut-curve to make it robust!
	if(subInterCurveDex[cDex][subInterPNum[cDex]]!=subInterCurveDex[cDex][0])
	{
		// the cut-curve is not closed
		// time to fail the algorithm
		// i.e., this should NOT possibly happen
		fprintf(stderr, "\nThe cut-curve of the sub-figure is not closed!!!\nHave to quit right here...\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fclose(fDebug);
#endif
		for(i=0; i<2; i++)
		{
			if(polyDexList1[i]!=NULL)
				delete polyDexList1[i];
			if(polyDexList2[i]!=NULL)
				delete polyDexList2[i];
			if(reDex[i]!=NULL)
				delete reDex[i];
		}

		if(polyDexList1!=NULL)
			delete []polyDexList1;
		if(polyDexList2!=NULL)
			delete []polyDexList2;

		if(reDex!=NULL)
			delete []reDex;
		if(vList!=NULL)
			delete []vList;
		if(fList!=NULL)
			delete []fList;
		if(vTag!=NULL)
			delete []vTag;
		if(fTag!=NULL)
			delete []fTag;

		return false;
	}
}

#ifdef QIONG_DEBUG_20030510
	color[0]=1;
	color[1]=1;
	color[2]=1;
	for(i=0; i<vOffsetT; i++)
		if(vTag[i]==1)
		{
			glRenderPoint(vList[i], color);
		}
#endif

//*
/// 2.2. after calculation of the cut-curve on the sub-figure
//		I will dilate the intersection curve in the main-figure's uvt space

	// NEW FUNCTION to check the orientation of the intersection curve
	// instead of using the DilateCurve to check the orientaion (later I won't map the whole curve 
	// into the 2D space then do the dilation any longer. instead, I will deal with the intersection 
	// curve chunk by chunk
	bool need2ReverseIntersectionCurve;
for(cDex=0; cDex<2; cDex++)
{
	need2ReverseIntersectionCurve=CheckIntersectionCurveOrientation(subInterCurve[cDex], subInterPNum[cDex], interCurve[cDex], interPNum[cDex]);
	fprintf(stderr, "Need to reverse the intersection curve? %s\n", need2ReverseIntersectionCurve?"Yes.":"No");

#ifdef QIONG_DEBUG_20030510
	// try to check all the t coordinates of all the intersection vertices
	// they are supposed to be 1/-1 if uv belongs to (0, uMax) and (0, vMax) respectively
	for(i=0; i<interPNum[cDex]; i++)
		//fprintf(stderr, "%8.7f\t%d\t", interCurveD[i].t, interCurveD[i].t==-1.0);
		fprintf(stderr, "%d\t", interCurveD[cDex][i].t==-1.0);
	fprintf(stderr, "\n");
#endif

	// now it might be time to reverse the whole intersection curve including all the 
	// other related arrays
	if(need2ReverseIntersectionCurve)
	{
		Bpoint2 tmpPnt;
		int tmpPair;
		for(i=0; i<interPNum[cDex]/2; i++)
		{
			j=interPNum[cDex]-1-i;

			tmpPnt=interCurve[cDex][i];
			interCurve[cDex][i]=interCurve[cDex][j];
			interCurve[cDex][j]=tmpPnt;

			tmpPair=interPPair[cDex][2*i];
			interPPair[cDex][2*i]=interPPair[cDex][2*j];
			interPPair[cDex][2*j]=tmpPair;

			tmpPair=interPPair[cDex][2*i+1];
			interPPair[cDex][2*i+1]=interPPair[cDex][2*j+1];
			interPPair[cDex][2*j+1]=tmpPair;
		}
		interCurve[cDex][interPNum[cDex]]=interCurve[cDex][0];
		interPPair[cDex][2*interPNum[cDex]]=interPPair[cDex][0];
		interPPair[cDex][2*interPNum[cDex]+1]=interPPair[cDex][1];
	}

	// implement the curve dialation and remeshing codes!
	MyFloat ACTUAL_DILATION_STEPS;
	switch(cDex)
	{
	case 0:
		ACTUAL_DILATION_STEPS=DILATION_STEPS;
		break;
	case 1:
		ACTUAL_DILATION_STEPS=DILATION_STEPS_2;
		break;
	default:
		break;
	}

	if(!DilateCurve(pListOld[0], interCurve[cDex], interPNum[cDex], DILATION_DT, ACTUAL_DILATION_STEPS))
	{
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "\nThe dilation of the intersection curve just FAILED, QUITTING...\n");
#endif
		fprintf(stderr, "\nThe dilation of the intersection curve just FAILED, QUITTING...\n");
		return false;
	}

	// Re-regularize the intersection curve after the dilation
	// after the dilation step, still need to re-regularize the dilated curve because at high-curvature
	// part, vertices tend to be crunched together
	//memcpy(interCurveD, interCurve, (interPNum+5)*sizeof(Bpoint2));
	RegularizeDilatedIntersectionCurve(interCurveD[cDex], interCurve[cDex], &interPNum[cDex]);
#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nThe v/t of all the %d vertices on the regularized dilated i-curve...\n", interPNum[cDex]);
	for(i=0; i<=interPNum[cDex]; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f\t", i, interCurveD[cDex][i].v_2, interCurveD[cDex][i].t_2);
	}
	fprintf(fDebug, "\n");
#endif

#ifdef QIONG_BLENDING_INTERSECTION_CURVE_DILATED
	color[0]=0;
	color[1]=1;
	color[2]=0;
	//glRenderCurve(0, interPNum[cDex], interCurveD[cDex], color);
	glRenderCurve(0, interPNum[cDex], interCurveD[cDex], color);
#endif

#ifdef QIONG_DEBUG_20030625
	color[0]=.5;
	color[1]=1;
	color[2]=.5;
	//glRenderPoint(interCurveD[cDex][11], color);
	glRenderPoint(interCurveD[cDex][11], color);
#endif
}

//*
/// 2.5. finally, insert all the vertices on the dilated intersection curve to the vList
	int vOffset2=0, fOffset2=0;
for(cDex=0; cDex<2; cDex++)
{
	for(i=0; i<interPNum[cDex]; i++)
	{
		vList[vOffsetT+vOffset2].setvals(interCurveD[cDex][i].pnt[0], interCurveD[cDex][i].pnt[1], interCurveD[cDex][i].pnt[2], 
										 interCurveD[cDex][i].norm[0], interCurveD[cDex][i].norm[1], interCurveD[cDex][i].norm[2], 
										 interCurveD[cDex][i].rad, interCurveD[cDex][i].u, interCurveD[cDex][i].v, interCurveD[cDex][i].t,
										 interCurveD[cDex][i].rad_2, interCurveD[cDex][i].u_2, interCurveD[cDex][i].v, interCurveD[cDex][i].t_2);
		vTag[vOffsetT+vOffset2]=2;
		vOffset2++;
	}
}

	vOffsetT+=vOffset2;
	fOffsetT+=fOffset2;

/// 3.0. connect the cut sub-figure with the dilated intersection curve
#ifndef QIONG_BLENDING_NO_DEBUG
	fprintf(stderr, "sub-figure now!\n");
#endif

	int vOffset3=0, fOffset3=0;

	int fvdex;
	//CCVertex *subVert[2];
	int subVertOldDex[2];

	// need to re-write all the codes here to connect the dilated-intersection-curve and 
	// the subfigure cut-curve

	//MyFloat tClosest=999, vClosest=999;
	MyPoint2D vtTemp, curTemp, nxtTemp, preTemp;
	int startI;
	int actualVOffset, fOffset3Temp;

for(cDex=0; cDex<2; cDex++)
{
	if(cDex==0)
		actualVOffset=vOffset;
	else
		actualVOffset=vOffset+interPNum[0];

	// find the first matching chunk of segments on the intersection curve to the first 
	// segment on the cut-curve of the sub-figure
	vtTemp.Set(9999, 9999);
	curTemp.Set(subInterCurve[cDex][0].v, subInterCurve[cDex][0].t);
	for(k=0; k<interPNum[cDex]; k++)
	{
		nxtTemp.Set(interCurveD[cDex][k].v_2, interCurveD[cDex][k].t_2);
		if(vtTemp.CircularDistance(curTemp)>nxtTemp.CircularDistance(curTemp))
		{
			startI=k;
			vtTemp.Set(interCurveD[cDex][k].v_2, interCurveD[cDex][k].t_2);
		}

	}

	//MyFloat tRange, vRange;
	i=startI; 
	j=startI; //j=startI+1;
#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "%d points on the cut-curve of the sub-figure\n", subInterPNum);
#endif
	for(fvdex=0; fvdex<subInterPNum[cDex]; fvdex++)
	{
		vtTemp.Set (subInterCurve[cDex][fvdex+1].v, subInterCurve[cDex][fvdex+1].t);
		preTemp.Set(interCurveD[cDex][INTMOD(j-1, interPNum[cDex])].v_2, interCurveD[cDex][INTMOD(j-1, interPNum[cDex])].t_2);
		curTemp.Set(interCurveD[cDex][INTMOD(j,   interPNum[cDex])].v_2, interCurveD[cDex][INTMOD(j,   interPNum[cDex])].t_2);
		nxtTemp.Set(interCurveD[cDex][INTMOD(j+1, interPNum[cDex])].v_2, interCurveD[cDex][INTMOD(j+1, interPNum[cDex])].t_2);
		while(!(curTemp.CircularDistance(vtTemp)<=preTemp.CircularDistance(vtTemp) && 
			    curTemp.CircularDistance(vtTemp)<=nxtTemp.CircularDistance(vtTemp)   ))
		{
			j++;
			preTemp=curTemp;
			curTemp=nxtTemp;
			nxtTemp.Set(interCurveD[cDex][INTMOD(j+1, interPNum[cDex])].v_2, interCurveD[cDex][INTMOD(j+1, interPNum[cDex])].t_2);
		}

		subVertOldDex[0]=subInterCurveDex[cDex][fvdex];
		subVertOldDex[1]=subInterCurveDex[cDex][fvdex+1];

//
//	subVertOldDex[0]			subVertOldDex[1]
//
//	i   i+1   i+2 ...                      ... j
//

		if((j-i+1)%2==1)
		{
			for(k=i; k<(i)+(j-i)/2; k++)
			{
				fList[fOffsetT+fOffset3].Create(3);
				fList[fOffsetT+fOffset3][2]=actualVOffset+(k%interPNum[cDex]);
				fList[fOffsetT+fOffset3][1]=reDex[0][subVertOldDex[0]];
				fList[fOffsetT+fOffset3][0]=actualVOffset+((k+1)%interPNum[cDex]);
				fTag[fOffsetT+fOffset3]=2;
				fOffset3++;
			}
			for(k=(i)+(j-i)/2; k<j; k++)
			{
				fList[fOffsetT+fOffset3].Create(3);
				fList[fOffsetT+fOffset3][2]=actualVOffset+(k%interPNum[cDex]);
				fList[fOffsetT+fOffset3][1]=reDex[0][subVertOldDex[1]];
				fList[fOffsetT+fOffset3][0]=actualVOffset+((k+1)%interPNum[cDex]);
				fTag[fOffsetT+fOffset3]=2;
				fOffset3++;
			}
			fList[fOffsetT+fOffset3].Create(3);
			fList[fOffsetT+fOffset3][2]=reDex[0][subVertOldDex[0]];
			fList[fOffsetT+fOffset3][1]=reDex[0][subVertOldDex[1]];
			fList[fOffsetT+fOffset3][0]=actualVOffset+((i)+(j-i)/2)%interPNum[cDex];
			fTag[fOffsetT+fOffset3]=2;
			fOffset3++;
		}
		else
		{
			for(k=i; k<(i)+(j-i+1)/2; k++)
			{
				fList[fOffsetT+fOffset3].Create(3);
				fList[fOffsetT+fOffset3][2]=actualVOffset+(k%interPNum[cDex]);
				fList[fOffsetT+fOffset3][1]=reDex[0][subVertOldDex[0]];
				fList[fOffsetT+fOffset3][0]=actualVOffset+((k+1)%interPNum[cDex]);
				fTag[fOffsetT+fOffset3]=2;
				fOffset3++;
			}
			for(k=(i)+(j-i+1)/2; k<j; k++)
			{
				fList[fOffsetT+fOffset3].Create(3);
				fList[fOffsetT+fOffset3][2]=actualVOffset+(k%interPNum[cDex]);
				fList[fOffsetT+fOffset3][1]=reDex[0][subVertOldDex[1]];
				fList[fOffsetT+fOffset3][0]=actualVOffset+((k+1)%interPNum[cDex]);
				fTag[fOffsetT+fOffset3]=2;
				fOffset3++;
			}
			fList[fOffsetT+fOffset3].Create(3);
			fList[fOffsetT+fOffset3][2]=reDex[0][subVertOldDex[0]];
			fList[fOffsetT+fOffset3][1]=reDex[0][subVertOldDex[1]];
			fList[fOffsetT+fOffset3][0]=actualVOffset+((i)+(j-i+1)/2)%interPNum[cDex];
			fTag[fOffsetT+fOffset3]=2;
			fOffset3++;
		}
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		if(fvdex%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %d->%d\t", fvdex, i%interPNum[cDex], j%interPNum[cDex]);
#endif
		i=j;
		//j++;
	}
#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\n");
#endif

#ifdef QIONG_BLENDING_SUB_FIGURE_INTERSECTION
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nThe v/t of all the %d vertices on the cut-curve %d of the sub-figure...\n", subInterPNum[cDex], cDex);
	for(i=0; i<=subInterPNum[cDex]; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f\t", i, subInterCurve[cDex][i].v, subInterCurve[cDex][i].t);
	}
	fprintf(fDebug, "\nThe v/t of all the %d vertices on the dilated-i-curve %d...\n", interPNum[cDex], cDex);
	for(i=0; i<=interPNum[cDex]; i++)
	{
		if(i%4==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d: %f/%f\t", i, interCurveD[cDex][i].v_2, interCurveD[cDex][i].t_2);
	}
	fprintf(fDebug, "\n");
	#endif
#endif
	if(cDex==0)
		fOffset3Temp=fOffset3;
}

	vOffsetT+=vOffset3;
	fOffsetT+=fOffset3;

#ifdef QIONG_BLENDING_SUB_FIGURE_INTERSECTION
	color[0]=1;
	color[1]=1;
	color[2]=0;
	glRender(vOffsetT, fOffset3, vList, &fList[fOffset+fOffset2], rMode, color);
#endif

#ifdef QIONG_BLENDING_SUB_FIGURE_QUAD_LIST
	color[0]=1;
	color[1]=0;
	color[2]=0;
	pListOld[1]->glRender(SUBDIV_LEVEL_TO_USE, 39, rMode, color);
#endif

//*
/// 4.0. find the cut version of the main-figure
	// re-write this part to re-mesh the blending area
	// for main-figure, re-relabling vertices first
	thisMesh=thisMesh1;
	vNums=thisMesh->VertList.Length();
	fNums=thisMesh->FaceList.Length();
	int vOffset4=0, fOffset4=0;

	// before everything, need to find the cut curve on the main figure that is connected to the dilated inter-curve
	// then figure out what vertices in the main figure to get rid of!

	int mainVNumEstimate=vNums*1.5;
	int mainFNumEstimate=fNums*1.5;

	int **quadListMain=new int*[2];
	int quadListMainNum[2];
	int **vertListMainQuad=new int*[2];
	int **vertListMainQuadCnt=new int*[2];
	int vertListMainQuadNum[2];
	int **edgeListMainQuad=new int*[2];
	int edgeListMainQuadNum[2];
	int **vertListMain=new int*[2];

	pListOld[0]->BuildLookUpTableForPDexByUVT();

	int **edgeListMainQuadDex=new int*[2];
	int mainInterPNum[2];

	int **takenMainQuadEdge=new int*[2];

	Bpoint2 mainInterCurve[2][MAX_INTERSECTION_NUMBER];
	int mainInterCurveDex[2][MAX_INTERSECTION_NUMBER];

	int mainInterCurveDexInEdgeListMainQuad[2][MAX_INTERSECTION_NUMBER];
	MyPoint2D mCurveIndices[2][MAX_MAIN_CURVE_NUMBER];

	int mCurveNum[2];
	int mCurveStartDex[2], mCurveEndDex[2];
	int successfulMCurveID[2];

for(cDex=0; cDex<2; cDex++)
{
	quadListMain[cDex]=new int[mainFNumEstimate];
	//quadListMain[1]=new int[mainFNumEstimate]
	vertListMainQuad[cDex]=new int[mainVNumEstimate];
	//vertListMainQuad[1]=new int[mainVNumEstimate];
	vertListMainQuadCnt[cDex]=new int[mainVNumEstimate];
	//vertListMainQuad[1=new int[mainVNumEstimate];
	edgeListMainQuad[cDex]=new int[mainVNumEstimate];
	//edgeListMainQuad[1]=new int[mainVNumEstimate];
	vertListMain[cDex]=new int[mainVNumEstimate];

	// flag for each vertex of the main-figure quads, 1 means it should be kept otherwise cut
	// index [0, vNums-1]
	// range [0, 1]
	memset(vertListMain[cDex], 1, mainVNumEstimate*sizeof(int));

	// quadListMainNum counts the total # of the intersected quads
	quadListMainNum[cDex]=0;
	// the list of all the quads involved with the intersection with the sub-figure
	// index [0, quadListMainNum-1]
	// range [0, fNums-1]
	memset(quadListMain[cDex], 1, mainFNumEstimate*sizeof(int));

	// vertListMainQuadNum counts the # of vertices of all the intersected quads
	vertListMainQuadNum[cDex]=0;
	// list of vertices of the intersected main-figured
	// index [0, vertListMainQuadNum-1]
	// range [0, vNums-1]
	memset(vertListMainQuad[cDex], -1, mainVNumEstimate*sizeof(int));

	// vertListMainQuadCnt counts the # of appearances of each vertex in the intersected 
	//		quads. if it is the valence of the vertex, it should be removed because it's 
	//		already surrouded by quads that should be removed
	// index [0, vNums-1]
	memset(vertListMainQuadCnt[cDex], 0, mainVNumEstimate*sizeof(int));

	edgeListMainQuadNum[cDex]=0;

	// according to the uvt coordinates of the dilated intersection curve, need to recalculate polyDexList1
	int pDexTmp; //, subLvl;
	//CCPolygon *thisPoly;
	for(i=0; i<interPNum[cDex]; i++)
	{
		//pDexTmp=pListOld[0]->lookUpTable[LookUpPDexByUVT(pListOld[0], interCurveD[i])];
		pDexTmp=pListOld[0]->LookUpPDexByUVT(interCurveD[cDex][i]);//lookUpTable[LookUpPDexByUVT(pListOld[0], interCurveD[i])];
if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
		pDexTmp=pListOld[0]->WithinPolyDexAtSublvlUVT(pDexTmp, SUBDIV_LEVEL_TO_USE+1, \
		interCurveD[cDex][i].u, interCurveD[cDex][i].v, interCurveD[cDex][i].t);
else
		pDexTmp=pListOld[0]->WithinPolyDexAtSublvlUVT(pDexTmp, SUBDIV_LEVEL_TO_USE, \
		interCurveD[cDex][i].u, interCurveD[cDex][i].v, interCurveD[cDex][i].t);

		if(!polyDexList1[cDex][pDexTmp]) 
			polyDexList1[cDex][pDexTmp]=2;
	}
	//Pointlist_serverB::LookUpPolygonIndexByUVT(Pointlist_server2 *pList, interCurveD[i])

#ifdef QIONG_DEBUG_20030510
	for(i=0; i<fNums; i++)
		if(polyDexList1[cDex][i])
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			fprintf(fDebug, "%d\t", i);
	#else
			fprintf(stderr, "%d\t", i);
	#endif
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\n");
	#else
	fprintf(stderr, "\n");
	#endif
#endif

	// now add the list of quads intersected by the dilated intersection curve into quadListMain
	for(i=0; i<fNums; i++)
		if(polyDexList1[cDex][i])
			quadListMain[cDex][quadListMainNum[cDex]++]=i;

#ifndef QIONG_BLENDING_NO_DEBUG
#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nRecording all the 'curV's for vertListMainQuadCnt\n");
#endif
#endif

	edgeListMainQuadDex[cDex]=new int[mainVNumEstimate];
	memset(edgeListMainQuadDex[cDex], -1, sizeof(int)*mainVNumEstimate);

	int curV, curEHd, curEEd;
	for(i=0; i<quadListMainNum[cDex]; i++)
	{
		thisPoly=thisMesh->FaceList(quadListMain[cDex][i]);
		pSides=thisPoly->numsides();
		for(j=0; j<pSides; j++)
		{
			// try to add curV into the vert list
			curV=thisPoly->Pvertices[j];

			// now need to deal with the edges
			//		it's slightly different cause we remove all the shared edges
			curEHd=curV;
			curEEd=thisPoly->Pvertices[(j+1)%pSides];
			k=0;
			while(k<edgeListMainQuadNum[cDex])
			{
				// have to reverse the order of curEEd and curEHd because of the topology
				if(curEEd==edgeListMainQuad[cDex][2*k] && curEHd==edgeListMainQuad[cDex][2*k+1])
					break;
				k++;
			}

			// if it's not, add it to the list
			//		and update the counter by 1 for the first time
			if(k==edgeListMainQuadNum[cDex])
			{
				edgeListMainQuadDex[cDex][edgeListMainQuadNum[cDex]]=quadListMain[cDex][i];

				edgeListMainQuad[cDex][2*edgeListMainQuadNum[cDex]]=curEHd;
				edgeListMainQuad[cDex][2*edgeListMainQuadNum[cDex]+1]=curEEd;
				edgeListMainQuadNum[cDex]++;
			}
			// otherwise need to remove the shared edge from the edgeListMainQuad
			else
			{
				if(k<edgeListMainQuadNum[cDex]-1)
				{
					while(k<edgeListMainQuadNum[cDex]-1)
					{
						edgeListMainQuadDex[cDex][k]=edgeListMainQuadDex[cDex][k+1];;

						edgeListMainQuad[cDex][2*k]=edgeListMainQuad[cDex][2*(k+1)];
						edgeListMainQuad[cDex][2*k+1]=edgeListMainQuad[cDex][2*(k+1)+1];
						k++;
					}
				}
				edgeListMainQuadNum[cDex]--;
			}
		}
	}

#ifndef QIONG_BLENDING_NO_DEBUG
#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\nRecording all the 'vertListMainQuad[cDex][i]'s for vertListMainQuadCnt\n");
#endif
#endif

	// it is time to find the curve on the main figure to connect to the dilated intersection curve!
	// write the codes to find out the cut-curve on the main figure
	//mainInterPNum=edgeListMainQuadNum;
	takenMainQuadEdge[cDex]=new int[mainFNumEstimate];
	memset(takenMainQuadEdge[cDex], 0, sizeof(int)*mainFNumEstimate);

	memset(mainInterCurveDexInEdgeListMainQuad[cDex], -1, sizeof(int)*MAX_INTERSECTION_NUMBER);

	bool algorithmFailed=false;
	bool timeToCheckCloseness=false;

	mCurveNum[cDex]=0;
	successfulMCurveID[cDex]=-1;

	mainInterCurveDex[cDex][0]=edgeListMainQuad[cDex][0];
	mainInterCurveDex[cDex][1]=edgeListMainQuad[cDex][1];
	mainInterCurveDexInEdgeListMainQuad[cDex][1-1]=0;
	mCurveStartDex[cDex]=0;
	j=0;
	i=1;

	takenMainQuadEdge[cDex][0]=1;

	// generate the whole intersection curve from the separated intersection segments
	//		all the un-ordereed intersection vertices in interPList now will be linked as a curve into interPList2
	successful=false;
	while(i<edgeListMainQuadNum[cDex] && !successful)
	{
		while(j<=edgeListMainQuadNum[cDex] && takenMainQuadEdge[cDex][j]) 
			j++;
		if(j>edgeListMainQuadNum[cDex])
		{
			// another possibility here is that the curve is already closed! (multiple curves could exist)
			// let us figure out whether it is true and whether the current curve is the outer curve

			if(mainInterCurveDex[cDex][i]==mainInterCurveDex[cDex][mCurveStartDex[cDex]])
			{
				// it is already a closed curve
				// need to check whether it is indeed the curve we are looking for
				// only the curve that is COUNTER-CLOCKWISE should be kept!
				if(i-mCurveStartDex[cDex]>edgeListMainQuadNum[cDex]/2)
				{
					// if the number of vertices in this curve is over the half of number of all the vertices
					// it is a trivial positive!
					mCurveEndDex[cDex]=i;
					mCurveIndices[cDex][mCurveNum[cDex]].x=mCurveStartDex[cDex];
					mCurveIndices[cDex][mCurveNum[cDex]].y=mCurveEndDex[cDex];
					successfulMCurveID[cDex]=mCurveNum[cDex];
					mCurveNum[cDex]++;
					successful=true;
					continue;
				}
				else
				{
					// otherwise, need to continue the searching
					mCurveEndDex[cDex]=i;
					mCurveIndices[cDex][mCurveNum[cDex]].x=mCurveStartDex[cDex];
					mCurveIndices[cDex][mCurveNum[cDex]].y=mCurveEndDex[cDex];
					mCurveNum[cDex]++;

					i++;
					mCurveStartDex[cDex]=i;
					j=0;
					while(j<=edgeListMainQuadNum[cDex] && takenMainQuadEdge[cDex][j]) 
						j++;
					if(j<=edgeListMainQuadNum[cDex])
					{
						mainInterCurveDexInEdgeListMainQuad[cDex][i]=j;

						takenMainQuadEdge[cDex][j]=1;
						mainInterCurveDex[cDex][i]=edgeListMainQuad[cDex][2*j];
						mainInterCurveDex[cDex][i+1]=edgeListMainQuad[cDex][2*j+1];
						j=0;

						//lastTaken=0;
						i++;
					}
					else
					{
						// the algorithm just failed!
						algorithmFailed=true;
					}
					continue;
				}
			}

			if(algorithmFailed)
			{
				// it is time to quit because this is the illegal case!
				bool trap=true;
				fprintf(stderr, "\nCannot figure out the mainInterCurve! time to quit blending codes\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fclose(fDebug);
#endif
			{
				for(i=0; i<2; i++)
				{
					if(polyDexList1[i]!=NULL)
						delete polyDexList1[i];
					if(polyDexList2[i]!=NULL)
						delete polyDexList2[i];
					if(reDex[i]!=NULL)
						delete reDex[i];

					if(quadListMain[i]!=NULL)
						delete quadListMain[i];
					if(vertListMainQuad[i]!=NULL)
						delete vertListMainQuad[i];
					if(vertListMainQuadCnt[i]!=NULL)
						delete vertListMainQuadCnt[i];
					if(edgeListMainQuad[i]!=NULL)			
						delete edgeListMainQuad[i];
					if(vertListMain[i]!=NULL)
						delete vertListMain[i];

					if(edgeListMainQuadDex[i]!=NULL)
						delete edgeListMainQuadDex[i];

					if(takenMainQuadEdge[i]!=NULL)
						delete takenMainQuadEdge[i];
				}

				if(polyDexList1!=NULL)
					delete []polyDexList1;
				if(polyDexList2!=NULL)
					delete []polyDexList2;

				if(reDex!=NULL)
					delete []reDex;
				if(vList!=NULL)
					delete []vList;
				if(fList!=NULL)
					delete []fList;
				if(vTag!=NULL)
					delete []vTag;
				if(fTag!=NULL)
					delete []fTag;

				if(quadListMain!=NULL)
					delete []quadListMain;
				if(vertListMainQuad!=NULL)
					delete []vertListMainQuad;
				if(vertListMainQuadCnt!=NULL)
					delete []vertListMainQuadCnt;
				if(edgeListMainQuad!=NULL)
					delete []edgeListMainQuad;
				if(vertListMain!=NULL)
					delete []vertListMain;

				if(edgeListMainQuadDex!=NULL)
					delete []edgeListMainQuadDex;

				if(takenMainQuadEdge!=NULL)
					delete []takenMainQuadEdge;			}

				return false;
			}
		}

		if(mainInterCurveDex[cDex][i]==edgeListMainQuad[cDex][2*j])
		{
			mainInterCurveDexInEdgeListMainQuad[cDex][i]=j;

			takenMainQuadEdge[cDex][j]=1;
			i++;
			mainInterCurveDex[cDex][i]=edgeListMainQuad[cDex][2*j+1];
			j=0;
			if(i==edgeListMainQuadNum[cDex])
				timeToCheckCloseness=true;
		}
		else
			if(mainInterCurveDex[cDex][i]==edgeListMainQuad[cDex][2*j+1])
			{
				// TO REMOVE: these codes are theoretically NOT-REACHABLE!
				mainInterCurveDexInEdgeListMainQuad[cDex][i]=j;

				takenMainQuadEdge[cDex][j]=1;
				i++;
				mainInterCurveDex[cDex][i]=edgeListMainQuad[cDex][2*j];
				j=0;
				if(i==edgeListMainQuadNum[cDex])
					timeToCheckCloseness=true;
			}
			else j++;

		if(timeToCheckCloseness)
		{
			if(mainInterCurveDex[cDex][i]==mainInterCurveDex[cDex][mCurveStartDex[cDex]])
			{
				// it is already a closed curve
				// need to check whether it is indeed the curve we are looking for
				// only the curve that is COUNTER-CLOCKWISE, or the LONGEST should be kept!
				if(i-mCurveStartDex[cDex]>edgeListMainQuadNum[cDex]/2)
				{
					mCurveEndDex[cDex]=i;
					mCurveIndices[cDex][mCurveNum[cDex]].x=mCurveStartDex[cDex];
					mCurveIndices[cDex][mCurveNum[cDex]].y=mCurveEndDex[cDex];
					successfulMCurveID[cDex]=mCurveNum[cDex];
					mCurveNum[cDex]++;
					successful=true;
				}
				else
				{
					mCurveEndDex[cDex]=i;
					mCurveIndices[cDex][mCurveNum[cDex]].x=mCurveStartDex[cDex];
					mCurveIndices[cDex][mCurveNum[cDex]].y=mCurveEndDex[cDex];
					mCurveNum[cDex]++;

					successful=true;
					//successful=true;
				}
			}
			else
			{
				// time to fail the algorithm again, because the curve is not closed!
				fprintf(stderr, "\nCannot figure out the mainInterCurve! time to quit blending codes\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				fclose(fDebug);
#endif
			{
				for(i=0; i<2; i++)
				{
					if(polyDexList1[i]!=NULL)
						delete polyDexList1[i];
					if(polyDexList2[i]!=NULL)
						delete polyDexList2[i];
					if(reDex[i]!=NULL)
						delete reDex[i];

					if(quadListMain[i]!=NULL)
						delete quadListMain[i];
					if(vertListMainQuad[i]!=NULL)
						delete vertListMainQuad[i];
					if(vertListMainQuadCnt[i]!=NULL)
						delete vertListMainQuadCnt[i];
					if(edgeListMainQuad[i]!=NULL)			
						delete edgeListMainQuad[i];
					if(vertListMain[i]!=NULL)
						delete vertListMain[i];

					if(edgeListMainQuadDex[i]!=NULL)
						delete edgeListMainQuadDex[i];

					if(takenMainQuadEdge[i]!=NULL)
						delete takenMainQuadEdge[i];
				}

				if(polyDexList1!=NULL)
					delete []polyDexList1;
				if(polyDexList2!=NULL)
					delete []polyDexList2;

				if(reDex!=NULL)
					delete []reDex;
				if(vList!=NULL)
					delete []vList;
				if(fList!=NULL)
					delete []fList;
				if(vTag!=NULL)
					delete []vTag;
				if(fTag!=NULL)
					delete []fTag;

				if(quadListMain!=NULL)
					delete []quadListMain;
				if(vertListMainQuad!=NULL)
					delete []vertListMainQuad;
				if(vertListMainQuadCnt!=NULL)
					delete []vertListMainQuadCnt;
				if(edgeListMainQuad!=NULL)
					delete []edgeListMainQuad;
				if(vertListMain!=NULL)
					delete []vertListMain;

				if(edgeListMainQuadDex!=NULL)
					delete []edgeListMainQuadDex;

				if(takenMainQuadEdge!=NULL)
					delete []takenMainQuadEdge;
			}

				return false;
			}
		}
	}
	//mainInterPNum=i;

	// by now
	// mCurveNum			holds the total number of the curves from the intersection quads
	// mCurveIndices		holds the start/end indices of the segments of each curve in mainInterCurveDex
	int maxLength=0;
	if(successfulMCurveID[cDex]<0)
	//if(!successful)
	{
		// then need to figure out which one of the curves is the intersection curve we need
		for(i=0; i<mCurveNum[cDex]; i++)
		{
			if(maxLength<mCurveIndices[cDex][i].y-mCurveIndices[cDex][i].x)
			{
				maxLength=mCurveIndices[cDex][i].y-mCurveIndices[cDex][i].x;
				successfulMCurveID[cDex]=i;
			}
		}
	}

	// by now
	// successfulMCurveID	holds the ID of the curve that is the cut curve of the main figure
	mCurveStartDex[cDex]=mCurveIndices[cDex][successfulMCurveID[cDex]].x;
	mCurveEndDex[cDex]=mCurveIndices[cDex][successfulMCurveID[cDex]].y;
	mainInterPNum[cDex]=mCurveEndDex[cDex]-mCurveStartDex[cDex];

	if(mCurveNum[cDex]<=0 || mainInterPNum[cDex]<=0)
	{
		fprintf(stderr, "\nCannot find a valid BOUNDARY CURVE on the main-figure!\nHave to quit blending codes...\n");
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fclose(fDebug);
#endif
			{
				for(i=0; i<2; i++)
				{
					if(polyDexList1[i]!=NULL)
						delete polyDexList1[i];
					if(polyDexList2[i]!=NULL)
						delete polyDexList2[i];
					if(reDex[i]!=NULL)
						delete reDex[i];

					if(quadListMain[i]!=NULL)
						delete quadListMain[i];
					if(vertListMainQuad[i]!=NULL)
						delete vertListMainQuad[i];
					if(vertListMainQuadCnt[i]!=NULL)
						delete vertListMainQuadCnt[i];
					if(edgeListMainQuad[i]!=NULL)			
						delete edgeListMainQuad[i];
					if(vertListMain[i]!=NULL)
						delete vertListMain[i];

					if(edgeListMainQuadDex[i]!=NULL)
						delete edgeListMainQuadDex[i];

					if(takenMainQuadEdge[i]!=NULL)
						delete takenMainQuadEdge[i];
				}

				if(polyDexList1!=NULL)
					delete []polyDexList1;
				if(polyDexList2!=NULL)
					delete []polyDexList2;

				if(reDex!=NULL)
					delete []reDex;
				if(vList!=NULL)
					delete []vList;
				if(fList!=NULL)
					delete []fList;
				if(vTag!=NULL)
					delete []vTag;
				if(fTag!=NULL)
					delete []fTag;

				if(quadListMain!=NULL)
					delete []quadListMain;
				if(vertListMainQuad!=NULL)
					delete []vertListMainQuad;
				if(vertListMainQuadCnt!=NULL)
					delete []vertListMainQuadCnt;
				if(edgeListMainQuad!=NULL)
					delete []edgeListMainQuad;
				if(vertListMain!=NULL)
					delete []vertListMain;

				if(edgeListMainQuadDex!=NULL)
					delete []edgeListMainQuadDex;

				if(takenMainQuadEdge!=NULL)
					delete []takenMainQuadEdge;
			}
		return false;
	}

#ifdef QIONG_BLENDING_DILATED_ICURVE_RELATED_QUADS
	color[0]=0;
	color[1]=0;
	color[2]=1;
	//*
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "The list of main-figure quads related to the dilated-i-curve...\n");
	#endif
	j=mCurveStartDex[cDex];
if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
{
	for(i=0; i<mainInterPNum[cDex]; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		if(i%8==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d\t", edgeListMainQuadDex[cDex][mainInterCurveDexInEdgeListMainQuad[cDex][j]]);
	#endif
		glRender(SUBDIV_LEVEL_TO_USE+1, edgeListMainQuadDex[cDex][mainInterCurveDexInEdgeListMainQuad[cDex][j]], rMode, color);
		j++;
	}
}
else
{
	for(i=0; i<mainInterPNum[cDex]; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		if(i%8==0)
			fprintf(fDebug, "\n");
		fprintf(fDebug, "%d\t", edgeListMainQuadDex[cDex][mainInterCurveDexInEdgeListMainQuad[cDex][j]]);
	#endif
		glRender(SUBDIV_LEVEL_TO_USE, edgeListMainQuadDex[cDex][mainInterCurveDexInEdgeListMainQuad[cDex][j]], rMode, color);
		j++;
	}
}
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\n");
	#endif
#endif

	// now it is time to cut the main figure by this intersection curve
	// need to grow the cut-curve and tag all the traversed quads... 
	// to quadListMain[quadListMainNum...]
	// before that, need two things
	//		1. the list of the indices of quad containing each segment in edgeListMainQuad	done
	//		2. a Queue stucture																done

	int curEdgeDexInEdgeListMainQuad, curEdgeStartDex, curEdgeEndDex;

	// before everything, tag all the edges of the dilated-i-curve related quads to span
	int *quadEdgeToSpan= new int[mainFNumEstimate*4];
	memset(quadEdgeToSpan, 1, mainFNumEstimate*4);

	j=mCurveStartDex[cDex];
	for(i=0; i<mainInterPNum[cDex]; i++)
	{
		curEdgeDexInEdgeListMainQuad=mainInterCurveDexInEdgeListMainQuad[cDex][j];
		curEdgeStartDex=edgeListMainQuad[cDex][2*curEdgeDexInEdgeListMainQuad];
		curEdgeEndDex=edgeListMainQuad[cDex][2*curEdgeDexInEdgeListMainQuad+1];
		thisPoly=thisMesh->FaceList(edgeListMainQuadDex[cDex][curEdgeDexInEdgeListMainQuad]);
		pSides=thisPoly->numsides();
		for(k=0; k<pSides; k++)
		{
			// need to only span the quad from the edge that was not contributing to the mainInterCurve
			// in other words, only span the quad inward
			if(thisPoly->Pvertices[k]==curEdgeStartDex && thisPoly->Pvertices[(k+1)%pSides]==curEdgeEndDex)
			{
				quadEdgeToSpan[4*edgeListMainQuadDex[cDex][curEdgeDexInEdgeListMainQuad]+k]=0;
			}
		}
		j++;
	}

	int *visited=new int[mainFNumEstimate];
	memset(visited, 0, sizeof(int)*mainFNumEstimate);
	MyQueue que(mainVNumEstimate);

	// fill in the que for the first time
	j=mCurveStartDex[cDex];
	for(i=0; i<mainInterPNum[cDex]; i++)
	{
		curEdgeDexInEdgeListMainQuad=mainInterCurveDexInEdgeListMainQuad[cDex][j];
		curEdgeStartDex=edgeListMainQuad[cDex][2*curEdgeDexInEdgeListMainQuad];
		curEdgeEndDex=edgeListMainQuad[cDex][2*curEdgeDexInEdgeListMainQuad+1];
		thisPoly=thisMesh->FaceList(edgeListMainQuadDex[cDex][curEdgeDexInEdgeListMainQuad]);

	if(!visited[edgeListMainQuadDex[cDex][curEdgeDexInEdgeListMainQuad]])
	{
		visited[edgeListMainQuadDex[cDex][curEdgeDexInEdgeListMainQuad]]=1;

		pSides=thisPoly->numsides();
		for(k=0; k<pSides; k++)
		{
			// need to only span the quad from the edge that was not contributing to the mainInterCurve
			// in other words, only span the quad inward
			if(quadEdgeToSpan[4*edgeListMainQuadDex[cDex][curEdgeDexInEdgeListMainQuad]+k])
			//if(!(thisPoly->Pvertices[k]==curEdgeStartDex && thisPoly->Pvertices[(k+1)%pSides]==curEdgeEndDex))
			{
#ifdef QIONG_BLENDING_EXPANDED_MAIN_FIGURE_QUADS_TO_REMOVE
				if(thisPoly->Pneighbors[k]==15)
				{
					bool ttt=true;
				}
#endif
				if(!visited[thisPoly->Pneighbors[k]])
				{
					if(!polyDexList1[cDex][thisPoly->Pneighbors[k]])
					{
						polyDexList1[cDex][thisPoly->Pneighbors[k]]=3;

						quadListMain[cDex][quadListMainNum[cDex]++]=thisPoly->Pneighbors[k];
					}

					que.Enque(thisPoly->Pneighbors[k]);
				}
			}
		}
	}
		j++;
	}

	// by traversing the queue, expand quads related to the dilated intersection curve to all the quads
	// that should be removed from the main figure quads
	int curQuadDex;
	while(!que.Empty())
	{
		curQuadDex=que.Deque();
		if(curQuadDex==-1)
		{
			bool trap=true;
		}
		if(!visited[curQuadDex])
		{
			visited[curQuadDex]=1;
			thisPoly=thisMesh->FaceList(curQuadDex);
			pSides=thisPoly->numsides();
			for(k=0; k<pSides; k++)
			{
#ifdef QIONG_BLENDING_EXPANDED_MAIN_FIGURE_QUADS_TO_REMOVE
				if(thisPoly->Pneighbors[k]==15)
				{
					bool ttt=true;
				}
#endif
				if(!visited[thisPoly->Pneighbors[k]])
				{
					if(!polyDexList1[cDex][thisPoly->Pneighbors[k]])
					{
						polyDexList1[cDex][thisPoly->Pneighbors[k]]=3;

						quadListMain[cDex][quadListMainNum[cDex]++]=thisPoly->Pneighbors[k];
					}

					que.Enque(thisPoly->Pneighbors[k]);
				}
			}
		}
	}
	delete []visited;
	delete []quadEdgeToSpan;

#ifdef QIONG_BLENDING_EXPANDED_MAIN_FIGURE_QUADS_TO_REMOVE
	color[0]=0;
	color[1]=1;
	color[2]=0;
	//*
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "The list of main-figure quads that was found by spanning the dilated-i-curve quads...\n");
	#endif
if(attachmentsType[1] == INDENT) //if(IS_IT_AN_INDENTATION_MODEL)
{
	for(i=0; i<nTiles1; i++)
	{
		if(polyDexList1[i]==3)
		{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			if(i%8==0)
				fprintf(fDebug, "\n");
			fprintf(fDebug, "%d\t", i);
	#endif
			glRender(SUBDIV_LEVEL_TO_USE+1, i, rMode, color);
		}
	}
}
else
{
	for(i=0; i<nTiles1; i++)
	{
		if(polyDexList1[cDex][i]==3)
		{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
			if(i%8==0)
				fprintf(fDebug, "\n");
			fprintf(fDebug, "%d\t", i);
	#endif
			glRender(SUBDIV_LEVEL_TO_USE, i, rMode, color);
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////

	#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fprintf(fDebug, "\n");
	#endif
#endif

	// update the quad-tags indicating that they are in the cut region
	// also tag all the vertices that are entirely surrounded by these cut quads
	for(i=0; i<quadListMainNum[cDex]; i++)
	{
		thisPoly=thisMesh->FaceList(quadListMain[cDex][i]);
		pSides=thisPoly->numsides();
		for(j=0; j<pSides; j++)
		{
			// try to add curV into the vert list
			curV=thisPoly->Pvertices[j];
			k=0;
			// check whether curV is already in the list
			while(k<vertListMainQuadNum[cDex])
			{
				if(curV==vertListMainQuad[cDex][k])
					break;
				k++;
			}
			// if it's not, add it to the list
			//		and update the counter by 1 for the first time
			if(k==vertListMainQuadNum[cDex])
			{
				vertListMainQuad[cDex][vertListMainQuadNum[cDex]++]=curV;
				vertListMainQuadCnt[cDex][curV]++;
#ifndef QIONG_BLENDING_NO_DEBUG
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				if(curV>=mainVNumEstimate)
					fprintf(fDebug, "%d\t", curV);
#endif
#endif
			}
			// otherwise need to update the counter and find the shared vertices
			else
			{
				vertListMainQuadCnt[cDex][curV]++;
#ifndef QIONG_BLENDING_NO_DEBUG
#ifdef QIONG_DEBUG_FILE_TO_WRITE
				if(curV>=mainVNumEstimate)
					fprintf(fDebug, "%d\t", curV);
#endif
#endif
			}	
		}
	}

	// tag all the vertices on the main figure to remove
	for(i=0; i<vertListMainQuadNum[cDex]; i++)
	{
		thisVert=thisMesh->VertList(vertListMainQuad[cDex][i]);
		if(thisVert->v_valence==vertListMainQuadCnt[cDex][vertListMainQuad[cDex][i]])
		//if(thisVert->v_valence==vertListMainQuadCnt[i])
			vertListMain[cDex][vertListMainQuad[cDex][i]]=0;
			//vertListMain[i]=0;
#ifndef QIONG_BLENDING_NO_DEBUG
#ifdef QIONG_DEBUG_FILE_TO_WRITE
		if(vertListMainQuad[cDex][i]>=mainVNumEstimate)
			fprintf(fDebug, "%d\t", vertListMainQuad[cDex][i]);
#endif
#endif
	}

	// generate the cut curve on the main figure, stored in mainInterCurve[0...mainInterPNum]
	for(i=mCurveStartDex[cDex]; i<=mCurveEndDex[cDex]; i++)
	//for(i=0; i<=mainInterPNum; i++)
	{
		thisVert=thisMesh->VertList(mainInterCurveDex[cDex][i]);

		//		to assign a Bpoint object to a Bpoint2 object, only the values of
		//		pnt/norm/rad are tranferred so far!

		mainInterCurve[cDex][i-mCurveStartDex[cDex]].pnt[0]=thisVert->limitloc.x();
		mainInterCurve[cDex][i-mCurveStartDex[cDex]].pnt[1]=thisVert->limitloc.y();
		mainInterCurve[cDex][i-mCurveStartDex[cDex]].pnt[2]=thisVert->limitloc.z();
	}

#ifdef QIONG_BLENDING_MAIN_FIGURE_INTERSECTION_CURVE
	for(i=mCurveStartDex[cDex]; i<=mCurveEndDex[cDex]; i++) //mainInterPNum; i++)
	{
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "%d\t", mainInterCurveDex[cDex][i]);
	#else
		fprintf(stderr, "%d\t", mainInterCurveDex[cDex][i]);
	#endif
	}
	#ifdef QIONG_DEBUG_FILE_TO_WRITE
		fprintf(fDebug, "\n");
	#else
		fprintf(stderr, "\n");
	#endif

	color[0]=0;
	color[1]=1;
	color[2]=0;
	glRenderCurve(0, mainInterPNum[cDex], mainInterCurve[cDex], color);
#endif
}

	reduction=0;
	// reDex is the important array to re-index all the orginal main-figure vertices to the new array in vList
	for(i=0; i<MAX_INTERSECTION_NUMBER; i++) 
		reDex[0][i]=i;
//*
	// reindex the vertices for the main-figure first
	//		need to remove all the vertices according to vertListMain[]
	//		0 means to be removed
	for(vDex=0; vDex<vNums; vDex++)
	{
		if(!vertListMain[0][vDex] || !vertListMain[1][vDex])
		{
			reDex[0][vDex]=-1;
			reduction++;
			continue;
		}
		else
		// 110802 need to get rid of all the points inside the dilated intersection curve
		{
			thisVert=thisMesh->VertList(vDex);
			reDex[0][vDex]=vOffsetT+vOffset4;
			vList[vOffsetT+vOffset4].setvals(thisVert->limitloc.x(), thisVert->limitloc.y(), thisVert->limitloc.z(),
									thisVert->limitnorm.x(), thisVert->limitnorm.y(), thisVert->limitnorm.x(), 
									thisVert->rad, thisVert->uvt.u, thisVert->uvt.v, thisVert->uvt.t);
			vTag[vOffsetT+vOffset4]=3;
			vOffset4++;
		}
	}

	// keep all the remaining patches of the main-figure first
	for(fDex=0; fDex<fNums; fDex++)
	{
		CCPolygon *thisPoly=thisMesh->FaceList(fDex);
		minU=999;
		double maxU=-1;
		if(!polyDexList1[0][fDex] && !polyDexList1[1][fDex])
		{
			fList[fOffsetT+fOffset4].Create(thisPoly->numsides());
			for(fvdex=0; fvdex<thisPoly->numsides(); fvdex++)
			{
				fList[fOffsetT+fOffset4][fvdex]=reDex[0][thisPoly->Pvertices[fvdex]];
			}
			fTag[fOffsetT+fOffset4]=4;
			fOffset4++;
		}
	}

	vOffsetT+=vOffset4;
	fOffsetT+=fOffset4;

#ifdef QIONG_BLENDING_MAIN_FIGURE
	color[0]=1;
	color[1]=0;
	color[2]=1;
	glRender(vOffsetT, fOffset4, vList, &fList[fOffset+fOffset2+fOffset3], rMode, color);
#endif


//*
/// 5.0. connect the rest of main-figure to the intersection curve
	int vOffset5=0, fOffset5=0;

#ifdef QIONG_BLENDING_INTERSECTION_MAIN_FIGURE
#endif

for(cDex=0; cDex<2; cDex++)
{
	if(cDex==0)
		actualVOffset=vOffset;
	else
		actualVOffset=vOffset+interPNum[0];

	// connect the 2 curves in mainInterCurve[mainInterPNum] and interCurveD[interPNum]
	bool goWithMain;
	bool mDone=false, iDone=false;
	bool done=false;
	int mainDex, interDex;
	int mainNext, interNext;
	int mainPicked, interPicked;

	int interStartPos=0;

	MyFloat dis1, dis2, dis;
	MyFloat angleMain, angleInter;
	MyFloat angle, angle1, angle2;

	dis1=9999;
	for(i=0; i<interPNum[cDex]; i++)
	{
		dis2=(mainInterCurve[cDex][mCurveStartDex[cDex]].pnt[0]-interCurveD[cDex][i].pnt[0])* \
			 (mainInterCurve[cDex][mCurveStartDex[cDex]].pnt[0]-interCurveD[cDex][i].pnt[0])+ \
			 (mainInterCurve[cDex][mCurveStartDex[cDex]].pnt[1]-interCurveD[cDex][i].pnt[1])* \
			 (mainInterCurve[cDex][mCurveStartDex[cDex]].pnt[1]-interCurveD[cDex][i].pnt[1])+ \
			 (mainInterCurve[cDex][mCurveStartDex[cDex]].pnt[2]-interCurveD[cDex][i].pnt[2])* \
			 (mainInterCurve[cDex][mCurveStartDex[cDex]].pnt[2]-interCurveD[cDex][i].pnt[2]);
		if(dis1>dis2)
		{
			dis1=dis2;
			interStartPos=i;
		}
	}

#ifdef QIONG_BLENDING_INTERSECTION_MAIN_FIGURE
#define PI 3.1415926
#endif

	Bpoint interDir, interDir1, interDir2, interDir3;
	MyFloat dot;//, dot1, dot2;
	mainDex=0;
	interDex=interStartPos;
	mainPicked=0;
	interPicked=0;

	bool secondThoughtWithMain;

	// start the intersection curve from the closest point to the first poin on the main-figure curve
	while(!(mDone&&iDone))
	{
		//int mainPrev=INTMOD(mainDex-1, mainInterPNum);
		int interPrev=INTMOD(interDex-1, interPNum[cDex]);
		int c2Next2=INTMOD(interDex+2, interPNum[cDex]);

		mainNext=(mainDex+1)%mainInterPNum[cDex];
		interNext=(interDex+1)%interPNum[cDex];
		if(!mDone && !iDone)
		{
			#define TRIANGULATION_CRETERION .9
			for(k=0; k<3; k++)
				interDir1.pnt[k]=(interCurveD[cDex][interDex].pnt[k]+mainInterCurve[cDex][mainNext].pnt[k])/2;
			for(k=0; k<3; k++)
				interDir.pnt[k]=interCurveD[cDex][interDex].pnt[k]-interDir1.pnt[k];
			dis=sqrt(interDir.pnt[0]*interDir.pnt[0]+interDir.pnt[1]*interDir.pnt[1]+ \
				interDir.pnt[2]*interDir.pnt[2]);
			for(k=0; k<3; k++)
				interDir.pnt[k]=interCurveD[cDex][interNext].pnt[k]-interDir1.pnt[k];
			dis1=sqrt(interDir.pnt[0]*interDir.pnt[0]+interDir.pnt[1]*interDir.pnt[1]+ \
				interDir.pnt[2]*interDir.pnt[2]);
			for(k=0; k<3; k++)
				interDir.pnt[k]=mainInterCurve[cDex][mainDex].pnt[k]-interDir1.pnt[k];
			dis2=sqrt(interDir.pnt[0]*interDir.pnt[0]+interDir.pnt[1]*interDir.pnt[1]+ \
				interDir.pnt[2]*interDir.pnt[2]);
			if(dis1<dis2)
				angleMain=dis1/dis;
			else 
				angleMain=dis2/dis;

			for(k=0; k<3; k++)
				interDir2.pnt[k]=(interCurveD[cDex][interNext].pnt[k]+mainInterCurve[cDex][mainDex].pnt[k])/2;
			for(k=0; k<3; k++)
				interDir.pnt[k]=interCurveD[cDex][interNext].pnt[k]-interDir2.pnt[k];
			dis=sqrt(interDir.pnt[0]*interDir.pnt[0]+interDir.pnt[1]*interDir.pnt[1]+ \
				interDir.pnt[2]*interDir.pnt[2]);
			for(k=0; k<3; k++)
				interDir.pnt[k]=interCurveD[cDex][interDex].pnt[k]-interDir2.pnt[k];
			dis1=sqrt(interDir.pnt[0]*interDir.pnt[0]+interDir.pnt[1]*interDir.pnt[1]+ \
				interDir.pnt[2]*interDir.pnt[2]);
			for(k=0; k<3; k++)
				interDir.pnt[k]=mainInterCurve[cDex][mainNext].pnt[k]-interDir2.pnt[k];
			dis2=sqrt(interDir.pnt[0]*interDir.pnt[0]+interDir.pnt[1]*interDir.pnt[1]+ \
				interDir.pnt[2]*interDir.pnt[2]);
			if(dis1<dis2)
				angleInter=dis1/dis;
			else 
				angleInter=dis2/dis;

			if(angleMain<TRIANGULATION_CRETERION && angleInter<TRIANGULATION_CRETERION)
			{
			for(k=0; k<3; k++)
				interDir1.pnt[k]=(mainInterCurve[cDex][mainNext].pnt[k]-mainInterCurve[cDex][mainDex].pnt[k]);
			dot=sqrt(interDir1.pnt[0]*interDir1.pnt[0]+interDir1.pnt[1]*interDir1.pnt[1]+ \
				interDir1.pnt[2]*interDir1.pnt[2]);
			for(k=0; k<3; k++)
				interDir1.pnt[k]=interDir1.pnt[k]/dot;
			for(k=0; k<3; k++)
				interDir2.pnt[k]=(interCurveD[cDex][interDex].pnt[k]-mainInterCurve[cDex][mainNext].pnt[k]);
			dot=sqrt(interDir2.pnt[0]*interDir2.pnt[0]+interDir2.pnt[1]*interDir2.pnt[1]+ \
				interDir2.pnt[2]*interDir2.pnt[2]);
			for(k=0; k<3; k++)
				interDir2.pnt[k]=interDir2.pnt[k]/dot;
			for(k=0; k<3; k++)
				interDir3.pnt[k]=(mainInterCurve[cDex][mainDex].pnt[k]-interCurveD[cDex][interDex].pnt[k]);
			dot=sqrt(interDir3.pnt[0]*interDir3.pnt[0]+interDir3.pnt[1]*interDir3.pnt[1]+ \
				interDir3.pnt[2]*interDir3.pnt[2]);
			for(k=0; k<3; k++)
				interDir3.pnt[k]=interDir3.pnt[k]/dot;
			angle1=fabs(interDir1.pnt[0]*interDir2.pnt[0]+interDir1.pnt[1]*interDir2.pnt[1]+ \
				interDir1.pnt[2]*interDir2.pnt[2]);

			angle=fabs(interDir2.pnt[0]*interDir3.pnt[0]+interDir2.pnt[1]*interDir3.pnt[1]+ \
				interDir2.pnt[2]*interDir3.pnt[2]);
			if(angle1<angle) angle1=angle;

			angle=fabs(interDir3.pnt[0]*interDir1.pnt[0]+interDir3.pnt[1]*interDir1.pnt[1]+ \
				interDir3.pnt[2]*interDir1.pnt[2]);
			if(angle1<angle) angle1=angle;
			////////////////////////////
			for(k=0; k<3; k++)
				interDir1.pnt[k]=(mainInterCurve[cDex][mainNext].pnt[k]-interCurveD[cDex][interNext].pnt[k]);
			dot=sqrt(interDir1.pnt[0]*interDir1.pnt[0]+interDir1.pnt[1]*interDir1.pnt[1]+ \
				interDir1.pnt[2]*interDir1.pnt[2]);
			for(k=0; k<3; k++)
				interDir1.pnt[k]=interDir1.pnt[k]/dot;
			for(k=0; k<3; k++)
				interDir3.pnt[k]=(interCurveD[cDex][interNext].pnt[k]-interCurveD[cDex][interDex].pnt[k]);
			dot=sqrt(interDir3.pnt[0]*interDir3.pnt[0]+interDir3.pnt[1]*interDir3.pnt[1]+ \
				interDir3.pnt[2]*interDir3.pnt[2]);
			for(k=0; k<3; k++)
				interDir3.pnt[k]=interDir3.pnt[k]/dot;
			angle=fabs(interDir1.pnt[0]*interDir2.pnt[0]+interDir1.pnt[1]*interDir2.pnt[1]+ \
				interDir1.pnt[2]*interDir2.pnt[2]);
			if(angle1<angle) angle1=angle;
			angle=fabs(interDir2.pnt[0]*interDir3.pnt[0]+interDir2.pnt[1]*interDir3.pnt[1]+ \
				interDir2.pnt[2]*interDir3.pnt[2]);
			if(angle1<angle) angle1=angle;
			angle=fabs(interDir3.pnt[0]*interDir1.pnt[0]+interDir3.pnt[1]*interDir1.pnt[1]+ \
				interDir3.pnt[2]*interDir1.pnt[2]);
			if(angle1<angle) angle1=angle;

			for(k=0; k<3; k++)
				interDir1.pnt[k]=(interCurveD[cDex][interNext].pnt[k]-mainInterCurve[cDex][mainDex].pnt[k]);
			dot=sqrt(interDir1.pnt[0]*interDir1.pnt[0]+interDir1.pnt[1]*interDir1.pnt[1]+ \
				interDir1.pnt[2]*interDir1.pnt[2]);
			for(k=0; k<3; k++)
				interDir1.pnt[k]=interDir1.pnt[k]/dot;
			for(k=0; k<3; k++)
				interDir2.pnt[k]=(interCurveD[cDex][interDex].pnt[k]-interCurveD[cDex][interNext].pnt[k]);
			dot=sqrt(interDir2.pnt[0]*interDir2.pnt[0]+interDir2.pnt[1]*interDir2.pnt[1]+ \
				interDir2.pnt[2]*interDir2.pnt[2]);
			for(k=0; k<3; k++)
				interDir2.pnt[k]=interDir2.pnt[k]/dot;
			for(k=0; k<3; k++)
				interDir3.pnt[k]=(mainInterCurve[cDex][mainDex].pnt[k]-interCurveD[cDex][interDex].pnt[k]);
			dot=sqrt(interDir3.pnt[0]*interDir3.pnt[0]+interDir3.pnt[1]*interDir3.pnt[1]+ \
				interDir3.pnt[2]*interDir3.pnt[2]);
			for(k=0; k<3; k++)
				interDir3.pnt[k]=interDir3.pnt[k]/dot;
			angle2=fabs(interDir1.pnt[0]*interDir2.pnt[0]+interDir1.pnt[1]*interDir2.pnt[1]+ \
				interDir1.pnt[2]*interDir2.pnt[2]);
			angle=fabs(interDir2.pnt[0]*interDir3.pnt[0]+interDir2.pnt[1]*interDir3.pnt[1]+ \
				interDir2.pnt[2]*interDir3.pnt[2]);
			// 2003/06/01
			//if(angle<angle) angle2=angle;
			if(angle2<angle) angle2=angle;
			angle=fabs(interDir3.pnt[0]*interDir1.pnt[0]+interDir3.pnt[1]*interDir1.pnt[1]+ \
				interDir3.pnt[2]*interDir1.pnt[2]);
			if(angle2<angle) angle2=angle;
			////////////////////////////
			for(k=0; k<3; k++)
				interDir2.pnt[k]=(mainInterCurve[cDex][mainNext].pnt[k]-interCurveD[cDex][interNext].pnt[k]);
			dot=sqrt(interDir2.pnt[0]*interDir2.pnt[0]+interDir2.pnt[1]*interDir2.pnt[1]+ \
				interDir2.pnt[2]*interDir2.pnt[2]);
			for(k=0; k<3; k++)
				interDir2.pnt[k]=interDir2.pnt[k]/dot;
			for(k=0; k<3; k++)
				interDir3.pnt[k]=(mainInterCurve[cDex][mainDex].pnt[k]-mainInterCurve[cDex][mainNext].pnt[k]);
			dot=sqrt(interDir3.pnt[0]*interDir3.pnt[0]+interDir3.pnt[1]*interDir3.pnt[1]+ \
				interDir3.pnt[2]*interDir3.pnt[2]);
			for(k=0; k<3; k++)
				interDir3.pnt[k]=interDir3.pnt[k]/dot;
			angle=fabs(interDir1.pnt[0]*interDir2.pnt[0]+interDir1.pnt[1]*interDir2.pnt[1]+ \
				interDir1.pnt[2]*interDir2.pnt[2]);
			if(angle2<angle) angle2=angle;
			angle=fabs(interDir2.pnt[0]*interDir3.pnt[0]+interDir2.pnt[1]*interDir3.pnt[1]+ \
				interDir2.pnt[2]*interDir3.pnt[2]);
			if(angle2<angle) angle2=angle;
			angle=fabs(interDir3.pnt[0]*interDir1.pnt[0]+interDir3.pnt[1]*interDir1.pnt[1]+ \
				interDir3.pnt[2]*interDir1.pnt[2]);
			if(angle2<angle) angle2=angle;

			if(angle1<angle2)
				secondThoughtWithMain=true;
			else
				secondThoughtWithMain=false;
			}

			if(angleMain>angleInter)
			{
				if(angleMain<TRIANGULATION_CRETERION)
					goWithMain=secondThoughtWithMain;
				// this means both triangulation is invalid, need to pick the relatively better one
				// with the bigger minimal angle
				else
					goWithMain=true;
			}
			else
			{
				if(angleInter<TRIANGULATION_CRETERION)
					goWithMain=secondThoughtWithMain;
				// this means both triangulation is invalid, need to pick the relatively better one
				// with the bigger minimal angle
				else
					goWithMain=false;
			}
			if(goWithMain)
			{
				// this validity check might still be necessary for a real robust implementation
				// check the validity of the segment of Inter(dex)->Main(next) by Inter(dex)->Inter(next)
				// with Inter(dex)->Inter(prev)
				// pick m[mNext] and i[iDex]

				//if(mainPicked<4)
				{
					goWithMain=true;
					mainPicked++;
				}
			}
			else
			{
				// pick m[mDex] and i[iNext]
				//if(interPicked<4)
				{
					goWithMain=false;
					interPicked++;
				}
			}
		}
		else
			if(mDone)
			{
				// pick m[mNext] and i[iDex]
				goWithMain=false;
			}
			else
				if(iDone)
				{
					// pick m[mNext] and i[iDex]
					goWithMain=true;
				}

		if(goWithMain)
		{
			// pick m[mNext] and i[iDex]
			fList[fOffsetT+fOffset5].Create(3);
			fList[fOffsetT+fOffset5][0]=reDex[0][mainInterCurveDex[cDex][mainDex]];
			fList[fOffsetT+fOffset5][1]=actualVOffset+interDex%interPNum[cDex]; //reDex[subVertOldDex[0]];
			fList[fOffsetT+fOffset5][2]=reDex[0][mainInterCurveDex[cDex][mainNext]];

			fTag[fOffsetT+fOffset5]=3;

			fOffset5++;

			mainDex=mainNext;
			if(mainDex==0)
				mDone=true;
		}
		else
		{
			// pick m[mDex] and i[iNext]
			fList[fOffsetT+fOffset5].Create(3);
			fList[fOffsetT+fOffset5][0]=reDex[0][mainInterCurveDex[cDex][mainDex]];
			fList[fOffsetT+fOffset5][1]=actualVOffset+interDex%interPNum[cDex];//reDex[subVertOldDex[0]];
			fList[fOffsetT+fOffset5][2]=actualVOffset+interNext%interPNum[cDex];//reDex[subVertOldDex[0]];

			fTag[fOffsetT+fOffset5]=3;

			fOffset5++;

			interDex=interNext;
			if(interDex==interStartPos)
				iDone=true;
		}
	}
}
	// here is the bug, it is wrong that
	// the order of all the vertices in the triangles connecting the mainfigure and the dialated intersection curve
	for(i=fOffsetT; i<fOffsetT+fOffset5; i++)
	{
		int temp;
		temp=fList[i][0];
		fList[i][0]=fList[i][2];
		fList[i][2]=temp;
	}

	vOffsetT+=vOffset5;
	fOffsetT+=fOffset5;

#ifdef QIONG_BLENDING_INTERSECTION_MAIN_FIGURE
	color[0]=0;
	color[1]=1;
	color[2]=1;
	glRender(vOffsetT, fOffset5, vList, &fList[fOffset+fOffset2+fOffset3+fOffset4], rMode, color);
	//glRender(vOffsetT, 1, vList, &fList[fOffset+fOffset2+fOffset3+fOffset4+fOffset5/2+9], rMode, color);
#endif

#ifdef QIONG_DEBUG_20030510
	color[0]=1;
	color[1]=1;
	color[2]=1;
	//for(i=55; i<65; i++)
	//	glRenderPoint(vList[59], color);
	bool trap1=true;
#endif

//*
/// 7.0.
//	actually generate the subdivision surface for the overall blended object!
#ifdef QIONG_BLENDING_NO_DEBUG
	fprintf(stderr, "Time to calculate the blended OBJECT......");
	MyFloat color[3];
	color[0]=0;
	color[1]=0;
	color[2]=1;
	InitializeBlendedSubdivSurf_New(SUBDIV_LEVEL_TO_USE, vOffsetT, fOffsetT, vList, fList, vTag, fTag);
	fprintf(stderr, "Done!\n");
#endif
//*/

#ifdef QIONG_DEBUG_FILE_TO_WRITE
	fclose(fDebug);
#endif

#ifdef DEBUG
	cout << "Tree " << treeId << " yielded " << vertexLists[treeId].size() << " verticies" << endl;
#endif

	for(i=0; i<2; i++)
	{
		if(polyDexList1[i]!=NULL)
			delete polyDexList1[i];
		if(polyDexList2[i]!=NULL)
			delete polyDexList2[i];
		if(reDex[i]!=NULL)
			delete reDex[i];

		if(quadListMain[i]!=NULL)
			delete quadListMain[i];
		if(vertListMainQuad[i]!=NULL)
			delete vertListMainQuad[i];
		if(vertListMainQuadCnt[i]!=NULL)
			delete vertListMainQuadCnt[i];
		if(edgeListMainQuad[i]!=NULL)
			delete edgeListMainQuad[i];
		if(vertListMain[i]!=NULL)
			delete vertListMain[i];

		if(edgeListMainQuadDex[i]!=NULL)
			delete edgeListMainQuadDex[i];

		if(takenMainQuadEdge[i]!=NULL)
			delete takenMainQuadEdge[i];
	}

	if(polyDexList1!=NULL)
		delete []polyDexList1;
	if(polyDexList2!=NULL)
		delete []polyDexList2;

	if(reDex!=NULL)
		delete []reDex;

	if(vList!=NULL)
		delete []vList;
	if(fList!=NULL)
		delete []fList;
	if(vTag!=NULL)
		delete []vTag;
	if(fTag!=NULL)
		delete []fTag;

	if(quadListMain!=NULL)
		delete []quadListMain;
	if(vertListMainQuad!=NULL)
		delete []vertListMainQuad;
	if(vertListMainQuadCnt!=NULL)
		delete []vertListMainQuadCnt;
	if(edgeListMainQuad!=NULL)
		delete []edgeListMainQuad;
	if(vertListMain!=NULL)
		delete []vertListMain;

	if(edgeListMainQuadDex!=NULL)
		delete []edgeListMainQuadDex;

	if(takenMainQuadEdge!=NULL)
		delete []takenMainQuadEdge;

	return true;
#endif
	return false;
}


// Function to extend the lists of figural point lists and blending boxes
// to include the children of a node
double Pointlist_serverB::linearize(M3DFigureTreeNode * node)
{
    int childCount, c;
	double lmin, m;		// Used to compute cubeSize and bound
//	int parentId;

    lmin = 1000000.0;
	if ((childCount = node->getChildCount()) > 0) 
	{
		M3DFigureTreeNode * subNode;

//		parentId = tileLists.size() - 1;

		// Loop through all attached subfigures
		for (c = 0; c < childCount; c++) 
		{
			subNode = node->getChild(c);

			// Create the pointList and update the minimum length of the
			// object in model space
			m = create_plist(subNode);
			if (m < lmin)
				lmin = m;

			// Compute the bounding boxes of each attached child figure
			/*
			Intersection * box = new Intersection(*subNode, subObject, tileLists[parentId],
				tileLists[parentId + c + 1]);
			box->setMultiplier(1.1);		// AGG: amount here is arbitrary!
			blendRegions.push_back(box);
			*/

			// Get the blending parameters from the subfigures
//			blendExtents[parentId + c] = (float) subNode->getBlendExtent();
//			blendAmounts[parentId + c] = (float) subNode->getBlendAmount();
		}

		// Recurse to process every attached subfigure
		for (c = 0; c < childCount; c++) 
		{
			subNode = node->getChild(c);
			m = linearize(subNode);	// Extend in breadth-first order
			if (m < lmin)
				lmin = m;
		}
	}
	return lmin;
}

// Compute cutoff times the minimum R of atoms in figure, a minimum scale for
// the figure.  A point list could be used for this, but using the atomic radii
// is much faster, with almost the same accuracy.
/*
double Pointlist_serverB::computeBBoxSize(M3DQuadFigure * figure, double cutoff)
{
}
*/

double Pointlist_serverB::create_plist(M3DFigureTreeNode * node) 
{
	M3DQuadFigure * figure;
	ThallCode::Pointlist_server2 * plist;
	Xferlist * xferList;
	int figureId;

#ifdef DEBUG
	cout << "Pointlist_serverB::create_plist()" << endl;
#endif
	figureId = node->getFigureId();
	figure = dynamic_cast<M3DQuadFigure*>( subObject->getFigurePtr(figureId));
	xferList = convertM3DtoXfer(figure);
	plist = new ThallCode::Pointlist_server2;
	plist->init(xferList);

	// we only need level surface to calculate the intersection
	if(figureId==0 && IS_IT_AN_INDENTATION_MODEL)
	{
		plist->ComputeSubdivBoundaryTiles(SUBDIV_LEVEL_TO_USE+1);
	}
	else
		plist->ComputeSubdivBoundaryTiles(SUBDIV_LEVEL_TO_USE);

	delete [](xferList->atomlist);
	delete xferList;

//	tileLists.push_back(plist);

	return 0;
}

