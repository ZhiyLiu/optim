/**************************************************************************************************************
2003/05/25
	NEW FUNCTIONS
	int WithinPolygonAtSublvlUVT(CCPolygon *thisPoly, int subLvl, double u, double v, double t);
	int WithinPolyDexAtSublvlUVT(int pDex, int subLvl, double u, double v, double t);

	bool Pointlist_server2::InitializeBlendedSubdivSurf_1(int vNum, int fNum, Bpoint2 *vList, PindexList *fList, int *vTag, int *fTag, Pointlist_server2 **pListOld)
	bool Pointlist_server2::InitializeBlendedSubdivSurf_New(int vNum, int fNum, Bpoint2 *vList, PindexList *fList, int *vTag, int *fTag, Pointlist_server2 **pListOld)

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
	in bool Pointlist_server2::InitializeSubdivSurf()
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

	111302	Pointlist_server2::UpdateSubdivSurf(Xferlist *thisList) and 
			Pointlist_server2::UpdateSubdivSurf(Diatomgrid *thisGrid)
			1. added this new function call, keep the topology untouched!

	111202	Pointlist_server2::subdivBposition
			1. disabled the codes normalizing the interpolated normals!
			ask Guodong to update his subfigure transformation codes
*/

/********************************************************************************/
/*                                                                              */
/*      File    :  Pointlist_server2.cpp										*/
/*                                                                              */
/*      Description:  class functions for PointCloud pointlist server,			*/
/*			which takes an Xferlist or Subdivsurf, copies the Diatomgrid, and	*/
/*			and returns points, normals, and radii of the surface				*/
/*			given an interpolation of the medial mesh to generate a dense		*/
/*			pointset, but NOT generating a smooth subdivision surface.			*/
/*																				*/
/*		User can change individual Diatoms, given row and column and a single	*/
/*			element Xferlist with the new value.								*/
/*																				*/
/*		NOTE:  we'll need to manually generate points on arcs between v1 and b  */
/*			and v2 and b, to get a dense pointcloud for voxelization.			*/
/*			How about nu from theta to 0 from v1 to bvec, where we have			*/
/*			have r(phi) = r*eta^((theta-nu)/theta).								*/
/*			No, this doesn't work...makes sharp point...think about derivative. */
/*			INSTEAD:  using raised cosine.  Not great, but okay.				*/
/*                                                                              */
/*      Project :  Seurat                                                       */
/*                                                                              */
/*      Author  :  A. Thall, S. Joshi                                           */
/*                                                                              */
/*      Date    :  1. June 2000													*/
/*                                                                              */
/*      Modifications:                                                          */
/*			4. June -- moved Xferlist-reading code to Diatomgrid				*/
/*			17. July -- added CreateBoundaryTiles(), tileinfo() and				*/
/*				glRenderTiles().												*/
/*			17. July -- fixed bug in CreatePointCloud() and CreateBoundaryT()	*/
/*				gridsubdivedge should NOT ever = 0.								*/
/*			7. Aug -- (Joshi & Thall) added (u, v, t) parameter to Bpoint and	*/
/*				relevant routines.												*/
/*			10. Aug -- added Bposition(u, v, t) providing a look-up for 		*/
/*				boundary positions												*/
/*			11. Aug -- changed gridsubdivrows=0 to rowspacing=1 for no			*/
/*				subdivision, and likewise gridsubdivcols to colspacing=1		*/
/*				Changed meshsubdivrow= 0 to meshrowspacing=1, & likewise		*/	
/*			16. Aug -- fixed bug in ComputeBoundaryTiles() --					*/
/*				AB is along the	col0-colmax direction, not row0-rowmax.			*/
/*				AD is along the row0-rowmax direction, not col0-colmax.			*/	
/*			4. Dec 00 -- changed alphafun() to etafun()							*/
/*			7. Dec 00 -- also merged in some partial-rendering code from Graham */
/*				Gash's version of older Seurat in Pablo.						*/
/*			23. Jan 01 -- fixed bug in ComputeBoundaryTiles() -- wasn't			*/
/*				calculating last set of normals for tiles on crest				*/
/*			29. Jan 01 -- modified to use Interpolfig::init() function, called	*/
/*				in init().  Added update_mesh() to modify internal Diatomgrid	*/
/*				and update Interpolfig object.									*/
/*			4. Feb 01 -- fixed bug in Bposition() -- incorrect test of 4 edge-	*/
/*				cases caused incorrect point output.							*/
/*			5. Feb 01 -- fixed bug in ComputeBoundaryTiles() -- incorrect test	*/
/*				for 4 edge-cases caused errors for meshes 2 atoms wide.			*/
/*			3. May 01 --														*/
/*				Added functions ComputeSubdivBoundaryTiles()					*/
/*								ComputeSubdivPointCloud()						*/
/*                          and SubdivTileInfo()								*/
/*                          and SubdivBoundaryInfo()							*/
/*				to create and access these subdivision surface boundary lists	*/
/*				(Added subdivbpoints and subdivnumpoints for vertex lists.)		*/
/*			5. May -- added subdivtileneighbors() to return tile-indexed list	*/
/*				of neighboring tile indices (See Class Tileneighbor in			*/
/*				Subdivsurf.h, simply int array of indices)						*/
/*				Added Tileneighbors *tneighbor as class member					*/
/*				Modified ComputeSubdivBoundaryTiles to compute tneighbor list.	*/
/*			23. May -- added glRenderSubdiv{}() functions for rendering tiles	*/
/*				and pointclouds from subdivision surface routines.				*/
/*			24. May -- found bug in (u, v, t) for v1 surface (top) of M-rep		*/
/*				in ComputeBoundaryTiles() had a uD,uD rather than uD, vD.		*/
/*			5. June 01 -- modified ComputeSubdiv{} from range [1, 7] to			*/
/*			              0 <= subdivlevel < NUMSUBLEVELS						*/
/*			9. June 01 -- fixed huge bug in indexing of (u, v) in routine		*/
/*				ComputePointCloud(); not only misscaled, but mixed up row and	*/
/*				column-major ordering.											*/
/*			3. Aug 01 -- added subdivNearpoint() functions,						*/
/*				which return (u, v, t)  and (x, y, z) coords for the			*/
/*				nearest point on the subdivided surface (at given level of		*/
/*				interations) to an external (x0, y0, z0).						*/
/*				NOTE:  this should be able to return a vector of such values	*/
/*			8. Aug 01 -- added subdivNearSurface(), which finds the bilinearly	*/
/*				interpolated point on the nearest tile with a Phong-normal		*/
/*				nearest the point (x0, y0, z0) in angle							*/
/*          13. Aug 01 -- added num_ucoltiles, num_vrowtiles, num_t_tiles to	*/
/*				hold data on subdivision tileset needed for indexing of (u,v,t) */
/*			    values; these are set in ComputeSubdivBoundaryTiles()			*/
/*			14. Aug -- added subdivBposition() to return Bpoint on the tiled	*/
/*				subdivision generated boundary corresponding to a given (u,v,t) */
/*			9. Sept -- modified subdivBposition() to return un-normalized		*/
/*				normal vector---had been normalizing after Phong edge interp.,	*/
/*				but after face interp., would not be normalized anyway.			*/
/*				This way, get same normal as produced by subdivNearSurface().	*/
/*          23. Sept -- added subdivNearSurfaceSIDED(), as above but only		*/
/*				searches tiles on top (TOP_p) or bottom (!TOP_p) tilesheet,		*/
/*				and ignores crest tiles.  This used for finding corresponding	*/
/*              involutes and medial atoms for a given surface Bposition.		*/
/*			    Added subdivMedialStuff() to return medial atom and opposing	*/
/*				boundary position corresponding to an input Bposition.			*/
/*			 2. Oct -- added outsideFLAG parameter to subdivNearSurface()		*/
/*				which is set if passed a non-NULL value, else no action taken.	*/
/*			 18. Oct -- fixed bug in subdivBposition() that gave wrong values	*/
/*				on crest when on first tiles on u[0-1]v0 side.					*/
/*			 18. Nov -- mods to subdivvertexneighbors() to compute vneighbors[]	*/
/*				based on new code in Subdivsurf.								*/
/*			 Dec 01 -- mods to search entire tileset if necessary to find		*/
/*				nearpoint on surface.  Now use std::queue<> and std::set<>		*/
/*				templates, rather than home-fried lists.						*/
/*			 Dec 01 -- modified glRenderSubdivTiles() with a displacement		*/
/*				term to perturb boundary toward medial axis.					*/
/*			23. Jan 02 -- added subdivtileindexinfo()							*/
/*                              copySubdivTileIndexInfo(),						*/
/*					 returns list of quad-tiles									*/
/*					as sets of 4 indices into bptlist of vertices returned		*/
/*                  by subdivboundaryinfo()										*/
/*				copySubdivTileInfo() was changed to return ntiles, not number of*/
/*					of points in list (== ntiles*4).							*/
/*				added int *subdivtileindices list for above						*/
/*				Fixed bug, had vertlistptr=NULL, instead of *vertlistptr=NULL	*/
/*			9. May 02 -- changed to better search caches in SubdivNearsurface() */
/*			16. May 02 -- changed to PseudoSets and cleaned up code in			*/
/*				subdivNearSurface()												*/
/*			24. June 02 -- eliminated workgrid, all interpolfig code			*/
/*			24. June 02 -- removed all .plist Parentset code					*/
/********************************************************************************/


//#include <time.h>
//#include <sys/timeb.h>
#include <typeinfo> 
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <queue>
#include <set>

#include <assert.h>

#include "globalBuildParms.h"
#define D_XFERLIST
#define D_POINTLIST_SERVER2
#define D_ZEROFINDER
#include "Shapedepend.h"

// Uncomment to get rainbow coloration of figures
//#define COLOR_DEBUG

#define TOP_BOTTOM 0
#define CREST_U 1
#define CREST_V 2

//#define QIONG_BLENDING_PURE_SMOOTHING
#define QIONG_BLENDING_MAINFIGURE_INTERSECTION_CURVE
#define QIONG_DEBUG_20030425
#define QIONG_DEBUG_02232003
// to divide the quads into triangles
//#define QIONG_DEBUG_SUBDIVISION_1		// all triangles, each quad -> 2 triangles
//#define QIONG_DEBUG_SUBDIVISION_2		// all triangles, center of quad inserted, each quad -> 4 triangles
// 02162003
// to test, it proves to fail that if we insert the purturbed middle vertex of each quad
//#define DEBUG_QIONG_TESTING_DISTURBING_MIDDLE_VERTEX_INSERTED
// 02192003
// work in progress
// do the inverse step AFTER one step subdivision
#define QIONG_DEBUG_BLENDING_02192003
// 02032003
// to define the inverted control polygon
//#define DEBUG_QIONG_SUBDIVISION_02032003
//#define DEBUG_QIONG
//#define DEBUG_QIONG_BLENDING
//#define QIONG_DEBUG_SEURAT_OPTIMIZATION
//#define DEBUG_QIONG_SUBDIVISION

#ifdef SURFACE_VISUALIZER
#include "M3DObjectSurfaceVisualizer.h"
#include "SurfaceColorMap.h"
#endif

using namespace ThallCode;

// effective zero for testing of double values as differences of near-integers
//   (want a loose fit, here, not mach-eps).
const double TESTEPS = 1e-10;

/********************************************************************************/
/* init() -- initialize the Pointlist_server2, passing in either an Xferlist	*/
/*    or an already created Subdivsurf; this creates a local newgrid.			*/
/********************************************************************************/
void Pointlist_server2::init(Xferlist *mfig1)
{
	// Load values from (row-major-ordered) Xferlist mfig1 to a Diatomgrid
	//   (column-major-ordered)
	// Why there is no check on whether newGrid==NULL?
	if(newGrid==NULL)
		newGrid=new Diatomgrid();

	newGrid->readXferlist(mfig1);
	smoothness = mfig1->tolerance;
}

void Pointlist_server2::init(Subdivsurf *subfig1)
{
	if(newGrid==NULL)
		newGrid=new Diatomgrid();

	if (newGrid->dlist != NULL)
		delete []newGrid->dlist;

	// Load values subfig1 to newgrid
	Diatomgrid *dptr = subfig1->Diatomgridptr();
	int nrows = dptr->rows;
	int ncols = dptr->cols;

	newGrid->rows = nrows;
	newGrid->cols = ncols;
	newGrid->dlist = new Diatom[nrows*ncols];

	for (int dgrid_index = 0; dgrid_index < nrows*ncols; dgrid_index++) {
		newGrid->dlist[dgrid_index] = dptr->dlist[dgrid_index];
	}
}

/********************************************************************************/
/* update_mesh() -- update the newgrid with the given medial atom at the given	*/
/*		location, and update the interpolator for this changed location.		*/
/********************************************************************************/
void Pointlist_server2::update_mesh(XferAtom *thisXatom, int Xfer_col, int Xfer_row)
{
	newGrid->update_mesh(thisXatom, Xfer_col, Xfer_row);
}

void Pointlist_server2::update_mesh(Diatom *thisDiatom, int u_col, int v_row)
{
	newGrid->update_mesh(thisDiatom, u_col, v_row);
}

/********************************************************************************/
/* etafun() -- compute a multiplier for the etavec for a given theta and		*/
/*		a given rotation rotangle from v1 to v2 where angle(v1,v2) = 2*theta.	*/
/********************************************************************************/
inline double Pointlist_server2::etafun(double eta, double theta, double rotangle)
{
	// try a simple sinusoid---window a cosine between -pi and pi, shift and scale
	// So, rotangle goes from 0 to 2*theta, so subtract theta and multiply by pi/theta

	double sinusend = ((rotangle - theta)*M_PI/theta);
	double multiplier = (cos(sinusend) + 1.0) / 2.0;
	return 1.0 + multiplier*(eta - 1.0);
}

// as above, but only lerps rad and (u, v, t)
void lerp_bpoint2(Bpoint *Bpnt, double t, Bpoint *ptA, Bpoint *ptB)
{
	Bpnt->rad = (1.0 - t)*ptA->rad + t*ptB->rad;
	Bpnt->u = (1.0 - t)*ptA->u + t*ptB->u;
	Bpnt->v = (1.0 - t)*ptA->v + t*ptB->v;
	Bpnt->t = (1.0 - t)*ptA->t + t*ptB->t;
}

/********************************************************************************
/* interpolate_quad2() -- as above, but only interpolates rad					*
/*	 and (u, v, t),	not entire Bpoint.											*
/********************************************************************************/
//void Pointlist_server2::interpolate_quad2(Bpoint *Bpnt, double uval, double vval, 
void interpolate_quad2(Bpoint *Bpnt, double uval, double vval, 
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

		lerp_bpoint2(&workpointAC, u_t, ptA, ptC);
		lerp_bpoint2(&workpointAD, u_t, ptA, ptD);
		lerp_bpoint2(Bpnt, v_t, &workpointAD, &workpointAC);
	}
	else {
		// compute patch coordinates
		u_t = uval/vval;
		v_t = vval;

		lerp_bpoint2(&workpointAC, v_t, ptA, ptC);
		lerp_bpoint2(&workpointAB, v_t, ptA, ptB);
		lerp_bpoint2(Bpnt, u_t, &workpointAB, &workpointAC);
	}
}

/********************************************************************************/
/* subdivNearpoint() -- return the nearpoint on the subdivided tileset to the	*/
/*    input (x0, y0, z0)														*/
/********************************************************************************/
double Pointlist_server2::subdivNearpoint(double x0, double y0, double z0,
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
double Pointlist_server2::subdivNearQuads(double x0, double y0, double z0,
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
bool possible_phongneartile(const DbVector3& thispoint,
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
bool Pointlist_server2::subdivNearSurface(double x0, double y0, double z0, 
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

			if (possible_phongneartile(thispoint, A, B, C, D, nA, nB, nC, nD)) {

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

				if (possible_phongneartile(thispoint, A, B, C, D, nA, nB, nC, nD)) {

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

				if (possible_phongneartile(thispoint, A, B, C, D, nA, nB, nC, nD)) {

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

				if (possible_phongneartile(thispoint, A, B, C, D, nA, nB, nC, nD)) {

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

		if (possible_phongneartile(thispoint, A, B, C, D, nA, nB, nC, nD)) {
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
			// NOTE:  we undo this when we call interpolate_quad2()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			u *= v;

			// compute true uvt for this quad
			if (rootfound) {
				interpolate_quad2(&uvtvals, u, v, 
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
			// NOTE:  we undo this when we call interpolate_quad2()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			v *= u;

			// compute true uvt for this quad
			if (rootfound) {

				interpolate_quad2(&uvtvals, u, v, 
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
			// NOTE:  we undo this when we call interpolate_quad2()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			u *= v;

			// compute true uvt for this quad
			if (rootfound) {
				interpolate_quad2(&uvtvals, u, v, 
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
			// NOTE:  we undo this when we call interpolate_quad2()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			v *= u;

			// compute true uvt for this quad
			if (rootfound) {

				interpolate_quad2(&uvtvals, u, v, 
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
		fprintf(stderr, "ERR::Pointlist_server2::subdivNearSurface()::no solution found\n");
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
bool Pointlist_server2::subdivNearSurfaceSIDED(double x0, double y0, double z0, searchSIDE testside,
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

			if (possible_phongneartile(thispoint, A, B, C, D, nA, nB, nC, nD)) {
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

				if (possible_phongneartile(thispoint, A, B, C, D, nA, nB, nC, nD)) {

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

				if (possible_phongneartile(thispoint, A, B, C, D, nA, nB, nC, nD)) {

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

				if (possible_phongneartile(thispoint, A, B, C, D, nA, nB, nC, nD)) {

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

		if (possible_phongneartile(thispoint, A, B, C, D, nA, nB, nC, nD)) {
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
			// NOTE:  we undo this when we call interpolate_quad2()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			u *= v;

			// compute true uvt for this quad (not computed by findroot)
			if (rootfound) {

				interpolate_quad2(&uvtvals, u, v, 
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
			// NOTE:  we undo this when we call interpolate_quad2()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			v *= u;

			// compute true uvt for this quad (not computed by findroot)
			if (rootfound) {

				interpolate_quad2(&uvtvals, u, v, 
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
			// NOTE:  we undo this when we call interpolate_quad2()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			u *= v;

			// compute true uvt for this quad (not computed by findroot3)
			if (rootfound) {

				interpolate_quad2(&uvtvals, u, v, 
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
			// NOTE:  we undo this when we call interpolate_quad2()...slight timewaster, but keeps things
			//   less confusing.  Eliminate SPEEDBUMP if microoptimizing.
			v *= u;

			// compute true uvt for this quad (not computed by findroot3)
			if (rootfound) {

				interpolate_quad2(&uvtvals, u, v, 
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
bool Pointlist_server2::subdivMedialStuff(Bpoint *boundpnt, Bpoint *boundpnt_opp, Diatom *medialatom)
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
			subdivBposition(boundpnt_opp, boundpnt->u, boundpnt->v, - boundpnt->t);
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
				subdivBposition(boundpnt_opp, boundpnt->u, boundpnt->v, - boundpnt->t);
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
bool Pointlist_server2::subdivMedialStuff2(Bpoint *boundpnt, Bpoint *boundpnt_opp, Diatom *medialatom)
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
	subdivBposition(boundpnt_opp, boundpnt->u, boundpnt->v, -boundpnt->t);
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
void Pointlist_server2::subdivboundaryinfo(int *numpts, Bpoint **bptlistptr)
{
	*numpts = numsubdivpoints;
	*bptlistptr = subdivbpoints;
}

// Added by Tom, does same as above but returns a copy, not the internal
void Pointlist_server2::copySubdivBoundaryInfo(int * numpts, Bpoint **bptlistptr)
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
void Pointlist_server2::subdivtileinfo(int *ntiles, Bpoint **bptlistptr)
{
	*ntiles = numsubdivtiles;
	*bptlistptr = subdivtileset;
}

// As above, but returns tile vertices as indices into subdivbpoint list
void Pointlist_server2::subdivtileindexinfo(int *ntiles, int **vertlistptr)
{
	*ntiles = numsubdivtiles;
	*vertlistptr = subdivtileindices;
}

// Added by Tom, does same as above but returns a copy, not the internal
// Also returns number of points not tiles
void Pointlist_server2::copySubdivTileInfo(int *ntiles, Bpoint **bptlistptr)
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
void Pointlist_server2::copySubdivTileIndexInfo(int *ntiles, int **vertlistptr)
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
void Pointlist_server2::subdivtileneighbors(int *ntiles, Tileneighbors **neighbors)
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
void Pointlist_server2::subdivvertexneighbors(int *nverts, Tileneighbors **neighbors)
{
	*nverts = numsubdivpoints;
	*neighbors = vneighbors;
}

// returns pointer to list of neighboring vertices for area calculation 
//	in Mask class
void Pointlist_server2::subdivNeighborVertices(int subLevel, int *nverts, NeighborVertices **_neighboringVerts) {
	*nverts = numsubdivpoints;
	*_neighboringVerts = neighboringVerts;
}


/********************************************************************************/
// Distance squared between two centers
/********************************************************************************/
double d_sqr(double * ctr_0, double * ctr_1) {
	return (ctr_1[0] - ctr_0[0])*(ctr_1[0] - ctr_0[0])
		+ (ctr_1[1] - ctr_0[1])*(ctr_1[1] - ctr_0[1])
		+ (ctr_1[2] - ctr_0[2])*(ctr_1[2] - ctr_0[2]);
}

/********************************************************************************/
/* glRenderSubdivPoints() -- if subdivbpoints != NULL, render surface point cloud*/
/********************************************************************************/
void Pointlist_server2::glRenderSubdivPoints()
{
	glPushAttrib(GL_POINT_BIT);
#ifdef COLOR_DEBUG
	glColor3d(1.0, 1.0, 1.0);
#endif
	glPointSize(lineDotWidth);
	glBegin(GL_POINTS);
	for (int vcnt1 = 0; vcnt1 < numsubdivpoints; vcnt1++) {
		glVertex3d(subdivbpoints[vcnt1].pnt[0], 
				   subdivbpoints[vcnt1].pnt[1],
				   subdivbpoints[vcnt1].pnt[2]);
	}
	glEnd();
	glPopAttrib();
}

/********************************************************************************/
/* glRenderSubdivWires() -- if numsubdivtiles != 0, render boundary wire frame	*/
/********************************************************************************/
void Pointlist_server2::glRenderSubdivWires()
{
	glPushAttrib(GL_LINE_BIT);

#ifdef COLOR_DEBUG
	glColor3d(0.6, 0.6, 0.6);
#endif

	glLineWidth(lineDotWidth);
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

	glPopAttrib();
}

/********************************************************************************/
/* glRenderSubdivTiles() -- if numsubdivtiles != 0, render boundary tiles		*/
/********************************************************************************/
void Pointlist_server2::glRenderSubdivTiles(double displacement)
{
#ifdef COLOR_DEBUG
	double u, v, t;
#endif
	int vcnt;

	// this gives displacement between 1 and zero.
	displacement = (1.0 - cos(displacement * 0.3/M_PI))/2.0;
	displacement *= 0.9;

	glBegin(GL_QUADS);

//	int numtopbottomtiles = 2*num_ucoltiles*num_vrowtiles;
//	for (vcnt = 0; vcnt < 4*numtopbottomtiles;/*numsubdivtiles;*/ vcnt++) {
	for (vcnt = 0; vcnt < 4*numsubdivtiles; vcnt++) {
#ifdef COLOR_DEBUG
		u = subdivtileset[vcnt].u;
		v = subdivtileset[vcnt].v;
		t = subdivtileset[vcnt].t;

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

char * colorMapMessage =
	"Error: attempting to use the Object Surface Visualizer with an empty colormap";

/********************************************************************************/
// glRenderColoredSubdivPoints() -- Like glRenderSubdivPoints, but color the
// surface using a colormap.  Used by class M3DObjectSurfaceVisualizer.
/********************************************************************************/
void Pointlist_server2::glRenderColoredSubdivPoints(float * colors, int numColors,
													float min, float max)
{
#ifdef SURFACE_VISUALIZER
	float color;
	float scale;

	if (numColors != numsubdivpoints) {
		// This occurs when the surface level is changed after the surface is computed in Match
		std::cout << "Incorrect number of match surface colors; try reinitializing the optimizer"
			<< std::endl;
		glRenderSubdivWires();
		return;
	}

	glPushAttrib(GL_POINT_BIT);
#ifdef COLOR_DEBUG
	glColor3d(1.0, 1.0, 1.0);
#endif

	scale = 1.0f/(max - min);

	int n = surfaceColorMap.size();	// Colormap length
	if (n <= 0) {
		std::cout << colorMapMessage << std::endl;
		glRenderSubdivPoints();
		return;
	}
	std::vector<float> & r = *(surfaceColorMap.red());
	std::vector<float> & g = *(surfaceColorMap.green());
	std::vector<float> & b = *(surfaceColorMap.blue());

	glPointSize(lineDotWidth);
	glBegin(GL_POINTS);
	for (int vcnt1 = 0; vcnt1 < numsubdivpoints; vcnt1++) {

		// Convert tile point numbers to subdivision point numbers before looking up color
		color = colors[subdivtileindices[vcnt1]];

		// Threshold the color to the min/max range
		if (color < min)
			color = min;
		else if (color > max)
			color = max;
		color = scale*(color - min);
		// Here, 0 <= color <= 1

		// Use the thresholded color (single) value computed above, as an index
		// (interpolated) into the colormap for drawing the tile vertex.
		float fIdx = (n - 1)*color;
		int iIdx(fIdx);
		float fraction = fIdx - iIdx;
		if (fraction > 0.0f) {
			// Interpolate the color
			glColor3f(((1.0f - fraction) * r[iIdx]) + (fraction * r[iIdx+1]), 
				((1.0f - fraction) * g[iIdx]) + (fraction * g[iIdx+1]), 
				((1.0f - fraction) * b[iIdx]) + (fraction * b[iIdx+1]));
		}
		else	// Lookup the exact color
			glColor3f(r[iIdx], g[iIdx], b[iIdx]);

		glVertex3d(subdivbpoints[vcnt1].pnt[0], 
		   subdivbpoints[vcnt1].pnt[1],
		   subdivbpoints[vcnt1].pnt[2]);
	}
	glEnd();
	glPopAttrib();
#endif	/* SURFACE_VISUALIZER */
}

/********************************************************************************/
// glRenderColoredSubdivWires() -- Like glRenderSubdivWires, but color the
// surface using a colormap.  Used by class M3DObjectSurfaceVisualizer.
/********************************************************************************/
void Pointlist_server2::glRenderColoredSubdivWires(float * colors, int numColors,
												   float min, float max)
{
#ifdef SURFACE_VISUALIZER
	float color;
	float scale;

	if (numColors != numsubdivpoints) {
		// This occurs when the surface level is changed after the surface is computed in Match
		std::cout << "Incorrect number of match surface colors; try reinitializing the optimizer"
			<< std::endl;
		glRenderSubdivWires();
		return;
	}

	glPushAttrib(GL_LINE_BIT);

#ifdef COLOR_DEBUG
	glColor3d(0.6, 0.6, 0.6);
#endif

	scale = 1.0f/(max - min);

	int n = surfaceColorMap.size();	// Colormap length
	if (n <= 0) {
		std::cout << colorMapMessage << std::endl;
		glRenderSubdivWires();
		return;
	}
	std::vector<float> & r = *(surfaceColorMap.red());
	std::vector<float> & g = *(surfaceColorMap.green());
	std::vector<float> & b = *(surfaceColorMap.blue());

	glLineWidth(lineDotWidth);
	for (int vcnt = 0; vcnt < 4*numsubdivtiles; vcnt += 4) {

		// Convert tile point numbers to subdivision point numbers before looking up color
		color = colors[subdivtileindices[vcnt]];

		// Threshold the color to the min/max range
		if (color < min)
			color = min;
		else if (color > max)
			color = max;
		color = scale*(color - min);
		// Here, 0 <= color <= 1

		// Use the thresholded color (single) value computed above, as an index
		// (interpolated) into the colormap for drawing the tile vertex.
		float fIdx = (n - 1)*color;
		int iIdx(fIdx);
		float fraction = fIdx - iIdx;
		if (fraction > 0.0f) {
			// Interpolate the color
			glColor3f(((1.0f - fraction) * r[iIdx]) + (fraction * r[iIdx+1]), 
				((1.0f - fraction) * g[iIdx]) + (fraction * g[iIdx+1]), 
				((1.0f - fraction) * b[iIdx]) + (fraction * b[iIdx+1]));
		}
		else	// Lookup the exact color
			glColor3f(r[iIdx], g[iIdx], b[iIdx]);

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

	glPopAttrib();
#endif	/* SURFACE_VISUALIZER */
}

/********************************************************************************/
// glRenderColoredSubdivTiles() -- Like glRenderSubdivTiles, but color the
// surface using a colormap.  Used by class M3DObjectSurfaceVisualizer.
/********************************************************************************/
void Pointlist_server2::glRenderColoredSubdivTiles(float * colors, int numColors,
	float min, float max, bool outlineTiles, double displacement)
{
#ifdef SURFACE_VISUALIZER
	int vcnt;
	float color;
	float scale;

	if (numColors != numsubdivpoints) {
		// This occurs when the surface level is changed after the surface is computed in Match
		std::cout << "Incorrect number of match surface colors; try reinitializing the optimizer"
			<< std::endl;
		glRenderSubdivTiles(displacement);
		return;
	}

	// This gives displacement between 1 and zero.
	displacement = (1.0 - cos(displacement * 0.3/M_PI))/2.0;
	displacement *= 0.9;

	scale = 1.0f/(max - min);

	int n = surfaceColorMap.size();	// Colormap length
	if (n <= 0) {
		std::cout << colorMapMessage << std::endl;
		glRenderSubdivTiles(displacement);
		return;
	}
	std::vector<float> & r = *(surfaceColorMap.red());
	std::vector<float> & g = *(surfaceColorMap.green());
	std::vector<float> & b = *(surfaceColorMap.blue());

	glBegin(GL_QUADS);

	for (vcnt = 0; vcnt < 4*numsubdivtiles; vcnt++) {

		// Convert tile point numbers to subdivision point numbers before looking up color
		color = colors[subdivtileindices[vcnt]];

		// Threshold the color to the min/max range
		if (color < min)
			color = min;
		else if (color > max)
			color = max;
		color = scale*(color - min);
		// Here, 0 <= color <= 1

		// Use the thresholded color (single) value computed above, as an index
		// (interpolated) into the colormap for drawing the tile vertex.
		float fIdx = (n - 1)*color;
		int iIdx(fIdx);
		float fraction = fIdx - iIdx;
		if (fraction > 0.0f) {
			// Interpolate the color
			glColor3f(((1.0f - fraction) * r[iIdx]) + (fraction * r[iIdx+1]), 
				((1.0f - fraction) * g[iIdx]) + (fraction * g[iIdx+1]), 
				((1.0f - fraction) * b[iIdx]) + (fraction * b[iIdx+1]));
		}
		else	// Lookup the exact color
			glColor3f(r[iIdx], g[iIdx], b[iIdx]);

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

	if (! outlineTiles) return;

	// Draw thin gray lines to indicate tile shapes
	glColor3d(0.2, 0.2, 0.2);
	double oldLineDotWidth = lineDotWidth;
	lineDotWidth = 1;
	glRenderSubdivWires();
	lineDotWidth = oldLineDotWidth;
#endif	/* SURFACE_VISUALIZER */
}

/********************************************************************************/
/* glPartialRenderboundary() -- if bpoints != NULL, render points with a value	*/
/* on an array of distances less than a given cutoff distance, coloring border-	*/
/* status of points.															*/
/********************************************************************************/
void Pointlist_server2::glPartialRenderSubdivPoints(double * distances,
								double cutoffDistance, bool complete)
{
	glPushAttrib(GL_POINT_BIT);
	glPointSize(lineDotWidth);
	glBegin(GL_POINTS);
    glColor3fv(color1);
	bool cflag = false;
	for (int vcnt1 = 0; vcnt1 < numsubdivpoints; vcnt1++) {

		if (distances[vcnt1] < cutoffDistance) {
			// Render inside quadrilaterals
			if (! cflag) {
				glColor3fv(color2);
				cflag = true;
			}
			glVertex3f((float) subdivbpoints[vcnt1].pnt[0],
				   (float) subdivbpoints[vcnt1].pnt[1],
				   (float) subdivbpoints[vcnt1].pnt[2]);
		}
		else {
			if (complete) {
				// Render outside quadrilaterals
				if (cflag) {
					glColor3fv(color1);
					cflag = false;
				}
				glVertex3f((float) subdivbpoints[vcnt1].pnt[0],
					   (float) subdivbpoints[vcnt1].pnt[1],
					   (float) subdivbpoints[vcnt1].pnt[2]);
			}
		}
	}
	glEnd();
	glPopAttrib();
}

/********************************************************************************
    glRenderSubdivWires() -- if numtiles != NULL, render boundary tiles as a wire
    frame.  For each point, distance is provided in array 'distances'.  The
    'complete' flag determines whether or not outside points are actually
    rendered.  If it is true, then those points having a distance value less
    than a given cutoff distance are rendered in one color, while those at or
    beyond the distance receive another color.  If not, than points on the
    outside are simply not rendered.
********************************************************************************/
void Pointlist_server2::glPartialRenderSubdivWires(double * distances,
							double cutoffDistance, bool complete)
{
	int vcnt;
	glPushAttrib(GL_LINE_BIT);
	glLineWidth(lineDotWidth);
    glColor3fv(color1);
	bool cflag = false;
	for (int cnt = 0; cnt < numsubdivtiles; cnt++) {
		vcnt = 4*cnt;

//		if (distances[vcnt] < 0)
//			std::cout << "Negative distance: " << distances[vcnt] << std::endl;

		int inside = 0;
		if (distances[vcnt] < cutoffDistance)
			inside++;
		if (distances[vcnt + 1] < cutoffDistance)
			inside++;
		if (distances[vcnt + 2] < cutoffDistance)
			inside++;
		if (distances[vcnt + 3] < cutoffDistance)
			inside++;

		if (inside == 4)
		{
			// Render inside quadrilaterals
			if (! cflag) {
				glColor3fv(color2);
				cflag = true;
			}

			glBegin(GL_LINE_LOOP);
			for (int i = vcnt; i < vcnt + 4; i++) {
				glNormal3f((float) subdivtileset[i].norm[0],
					(float) subdivtileset[i].norm[1],
					(float) subdivtileset[i].norm[2]);
				glVertex3f((float) subdivtileset[i].pnt[0],
					(float) subdivtileset[i].pnt[1],
					(float) subdivtileset[i].pnt[2]);
//				if (i == vcnt) std::cout << "* ";
//				else std::cout << "  ";
//				std::cout << i << ". " << distances[i] << " < " << cutoffDistance
//					<< "  inside point: (" << tileset[i].pnt[0] << ", "
//					<< tileset[i].pnt[1] << ", " << tileset[i].pnt[2] << ')' << std::endl;
			}
			glEnd();
		}
		else if (inside == 0) {
			// Render outside quadrilaterals
			if (complete) {
				if (cflag) {
					glColor3fv(color1);
					cflag = false;
				}
				glBegin(GL_LINE_LOOP);
				for (int i = vcnt; i < vcnt + 4; i++) {
					glNormal3f((float) subdivtileset[i].norm[0],
						(float) subdivtileset[i].norm[1],
						(float) subdivtileset[i].norm[2]);
					glVertex3f((float) subdivtileset[i].pnt[0],
						(float) subdivtileset[i].pnt[1],
						(float) subdivtileset[i].pnt[2]);
				}
				glEnd();
			}
		}
		else {
			// Render quadrilaterals crossing the border
			if (complete) {
				glBegin(GL_LINE_LOOP);
				for (int i = vcnt; i < vcnt + 4; i++) {
					if (distances[i] < cutoffDistance) {
						if (! cflag) {
							glColor3fv(color2);
							cflag = true;
						}
						glNormal3f((float) subdivtileset[i].norm[0],
							(float) subdivtileset[i].norm[1],
							(float) subdivtileset[i].norm[2]);
						glVertex3f((float) subdivtileset[i].pnt[0],
							(float) subdivtileset[i].pnt[1],
							(float) subdivtileset[i].pnt[2]);
					}
					else {
						if (cflag) {
							glColor3fv(color1);
							cflag = false;
						}
						glNormal3f((float) subdivtileset[i].norm[0],
							(float) subdivtileset[i].norm[1],
							(float) subdivtileset[i].norm[2]);
						glVertex3f((float) subdivtileset[i].pnt[0],
							(float) subdivtileset[i].pnt[1],
							(float) subdivtileset[i].pnt[2]);
					}
				}
				glEnd();
			}
			else {
				glBegin(GL_LINE_LOOP);
				for (int i = vcnt; i < vcnt + 4; i++) {
					if (distances[i] < cutoffDistance) {
						if (! cflag) {
							glColor3fv(color2);
							cflag = true;
						}
						glNormal3f((float) subdivtileset[i].norm[0],
							(float) subdivtileset[i].norm[1],
							(float) subdivtileset[i].norm[2]);
						glVertex3f((float) subdivtileset[i].pnt[0],
							(float) subdivtileset[i].pnt[1],
							(float) subdivtileset[i].pnt[2]);
					}
					else {
						glColor3fv(color3);
						glNormal3f((float) subdivtileset[i].norm[0],
							(float) subdivtileset[i].norm[1],
							(float) subdivtileset[i].norm[2]);
						glVertex3f((float) subdivtileset[i].pnt[0],
							(float) subdivtileset[i].pnt[1],
							(float) subdivtileset[i].pnt[2]);
						glColor3fv(color2);
						cflag = true;
					}
					//					if (i == vcnt) std::cout << "* ";
					//					else std::cout << "  ";
					//					std::cout << i << ". " << distances[i] << " < " << cutoffDistance
					//						<< "  halfway point: (" << tileset[i].pnt[0] << ", "
					//						<< tileset[i].pnt[1] << ", " << tileset[i].pnt[2] << ')' << std::endl;
				}
				glEnd();
			}
		}
	}
	glPopAttrib();
}

/********************************************************************************
    glRenderTiles() -- if numtiles != NULL, render boundary tiles as a solid
    surface.  For each point, distance is provided in array 'distances'.  The
    'complete' flag determines whether or not outside points are actually
    rendered.  If it is true, then those points having a distance value less
    than a given cutoff distance are rendered in one color, while those at or
    beyond the distance receive another color.  If not, than points on the
    outside are simply not rendered.
********************************************************************************/
void Pointlist_server2::glPartialRenderSubdivTiles(double * distances,
							double cutoffDistance, bool complete)
{
	int vcnt;
	glBegin(GL_QUADS);
    glColor3fv(color1);
	bool cflag = false;

	for (int cnt = 0; cnt < numsubdivtiles; cnt++) {
		vcnt = 4*cnt;

//		if (distances[vcnt] < 0)
//			std::cout << "Negative distance: " << distances[vcnt] << std::endl;

		int inside = 0;
		if (distances[vcnt] < cutoffDistance)
			inside++;
		if (distances[vcnt + 1] < cutoffDistance)
			inside++;
		if (distances[vcnt + 2] < cutoffDistance)
			inside++;
		if (distances[vcnt + 3] < cutoffDistance)
			inside++;

		if (inside == 4)
		{
			// Render inside quadrilaterals
			if (! cflag) {
				glColor3fv(color2);
				cflag = true;
			}

			for (int i = vcnt; i < vcnt + 4; i++) {
				glNormal3f((float) subdivtileset[i].norm[0],
					(float) subdivtileset[i].norm[1],
					(float) subdivtileset[i].norm[2]);
				glVertex3f((float) subdivtileset[i].pnt[0],
					(float) subdivtileset[i].pnt[1],
					(float) subdivtileset[i].pnt[2]);
//				if (i == vcnt) std::cout << "* ";
//				else std::cout << "  ";
//				std::cout << i << ". " << distances[i] << " < " << cutoffDistance
//					<< "  inside point: (" << tileset[i].pnt[0] << ", "
//					<< tileset[i].pnt[1] << ", " << tileset[i].pnt[2] << ')' << std::endl;
			}
		}
		else if (inside == 0) {
			// Render outside quadrilaterals
			if (complete) {
				if (cflag) {
					glColor3fv(color1);
					cflag = false;
				}
				for (int i = vcnt; i < vcnt + 4; i++) {
					glNormal3f((float) subdivtileset[i].norm[0],
						(float) subdivtileset[i].norm[1],
						(float) subdivtileset[i].norm[2]);
					glVertex3f((float) subdivtileset[i].pnt[0],
						(float) subdivtileset[i].pnt[1],
						(float) subdivtileset[i].pnt[2]);
				}
			}
		}
		else {
			// Render quadrilaterals crossing the border
			if (complete) {
				for (int i = vcnt; i < vcnt + 4; i++) {
					if (distances[i] < cutoffDistance) {
						if (! cflag) {
							glColor3fv(color2);
							cflag = true;
						}
						glNormal3f((float) subdivtileset[i].norm[0],
							(float) subdivtileset[i].norm[1],
							(float) subdivtileset[i].norm[2]);
						glVertex3f((float) subdivtileset[i].pnt[0],
							(float) subdivtileset[i].pnt[1],
							(float) subdivtileset[i].pnt[2]);
					}
					else {
						if (cflag) {
							glColor3fv(color1);
							cflag = false;
						}
						glNormal3f((float) subdivtileset[i].norm[0],
							(float) subdivtileset[i].norm[1],
							(float) subdivtileset[i].norm[2]);
						glVertex3f((float) subdivtileset[i].pnt[0],
							(float) subdivtileset[i].pnt[1],
							(float) subdivtileset[i].pnt[2]);
					}
				}
			}
			else {
				for (int i = vcnt; i < vcnt + 4; i++) {
					if (distances[i] < cutoffDistance) {
						if (! cflag) {
							glColor3fv(color2);
							cflag = true;
						}
						glNormal3f((float) subdivtileset[i].norm[0],
							(float) subdivtileset[i].norm[1],
							(float) subdivtileset[i].norm[2]);
						glVertex3f((float) subdivtileset[i].pnt[0],
							(float) subdivtileset[i].pnt[1],
							(float) subdivtileset[i].pnt[2]);
					}
					else {
						glColor3fv(color3);
						glNormal3f((float) subdivtileset[i].norm[0],
							(float) subdivtileset[i].norm[1],
							(float) subdivtileset[i].norm[2]);
						glVertex3f((float) subdivtileset[i].pnt[0],
							(float) subdivtileset[i].pnt[1],
							(float) subdivtileset[i].pnt[2]);
						glColor3fv(color2);
						cflag = true;
					}
//					if (i == vcnt) std::cout << "* ";
//					else std::cout << "  ";
//					std::cout << i << ". " << distances[i] << " < " << cutoffDistance
//						<< "  halfway point: (" << tileset[i].pnt[0] << ", "
//						<< tileset[i].pnt[1] << ", " << tileset[i].pnt[2] << ')' << std::endl;
				}
			}
		}
	}
	glEnd();
}

//////////////////////////////////////////////////////////////////////////
// Added by CL 01/21/02.  Add displacements to boundary point list

void Pointlist_server2::applyDisplacements(int numPts, double * dispVals)
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
//  Added for Eli for Atom regions
bool Pointlist_server2::computeLvl0VertAtomList(int numRows, int numCols, int *list)
{
	return thisSurf->computeLvl0VertAtomList(numRows, numCols, list);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Pointlist_server2::ClearBoundaryTilesPointers()
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

void Pointlist_server2::ClearPointCloudPointers()
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
	if(neighboringVerts!=NULL) {
		delete []neighboringVerts;
		neighboringVerts = NULL;
	}
}

bool Pointlist_server2::UpdateSubdivSurf(int subLevel, Xferlist *thisList)
{
	int slevel;
	if(newGrid==NULL)
		return false;

	if( thisMesh != NULL ) {
		thisMesh->CopyXferList(thisList);
		thisMesh->UpdateVertexList();
		thisSurf->InitValue(thisMesh->numverts(), thisMesh->vertlist());
	}
	else {
		return false;
	}

	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
		thisSurf->LoadLevelZeroValue();
	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
		smoothness=thisList->tolerance;
		thisSurf->setSmoothness(smoothness);
	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////

		if(subLevel-1>highestInitializedSubLevel)
		{
			for(slevel=highestInitializedSubLevel+1; slevel<subLevel; slevel++)
				thisSurf->splitandaverage(slevel);
			highestInitializedSubLevel=subLevel-1;
		}

		for(slevel = 0; slevel < subLevel; slevel++)
			thisSurf->Interpolate(slevel);
	////////////////////////////////////////////////////////////////

	return true;
}

void Pointlist_server2::UpdateSubdivPointCloud(int subLevel, Xferlist *thisList)
{
	if(!subdivSurfInitialized) 
	{
		return;
	}

	if(subLevel!=currentPointCloudLevel)
	{
		currentPointCloudLevel=subLevel;
		ClearPointCloudPointers();
	}

	if(!UpdateSubdivSurf(subLevel, thisList))
		return;

	int numvertneighbors;
	thisSurf->SubdivBoundaryInfo(subLevel, &numsubdivpoints, &subdivbpoints);
	thisSurf->GetVertexNeighbors(subLevel, &numvertneighbors, &vneighbors);
	thisSurf->GetNeighboringVertices(subLevel, &neighboringVerts);

	if (numsubdivpoints!=numvertneighbors)
		fprintf(stderr, "Pointlistserver2::ComputeSubdivPointCloud:  ERROR, bad number of vertices\n");
}

void Pointlist_server2::UpdateSubdivBoundaryTiles(int subLevel, Xferlist *thisList)
{
	if(!subdivSurfInitialized) 
	{
		return;
	}

	if(subLevel!=currentBoundaryTilesLevel)
	{
		currentBoundaryTilesLevel=subLevel;
		ClearBoundaryTilesPointers();
	}

	if(!UpdateSubdivSurf(subLevel, thisList))
		return;

	// Deal with NULL pointers instead in thisSurf->SubdivTileInfo, GetTileNeighbors and GetTileIndices

	//int ntiles;
	thisSurf->SubdivTileInfo(subLevel, &numsubdivtiles, &subdivtileset);
	thisSurf->GetTileNeighbors(subLevel, &numsubdivtiles, &tneighbors);
	thisSurf->GetTileIndices(subLevel, &subdivtileindices);
}


bool Pointlist_server2::UpdateSubdivSurf(int subLevel, Diatomgrid *thisGrid)
{
	int slevel;
	if( newGrid == NULL || thisSurf == NULL ) {
		return false;
	}

	if( thisMesh != NULL ) {
		thisMesh->CopyDiatomGrid(thisGrid);
		thisMesh->UpdateVertexList();
		thisSurf->InitValue(thisMesh->numverts(), thisMesh->vertlist());
	}
	else {
		return false;
	}

	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
		thisSurf->LoadLevelZeroValue();
	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
		for(slevel = 0; slevel < subLevel; slevel++)
			thisSurf->Interpolate(slevel);
	////////////////////////////////////////////////////////////////

	return true;
}

void Pointlist_server2::UpdateSubdivPointCloud(int subLevel, Diatomgrid *thisGrid)
{
	if(!subdivSurfInitialized) 
	{
		return;
	}

	if(subLevel!=currentPointCloudLevel)
	{
		currentPointCloudLevel=subLevel;
		ClearPointCloudPointers();
	}
	// deal with the NULL pointers in thisSurf->SubdivBoundaryInfo and GetVertexNeighbors instead
	if(!UpdateSubdivSurf(subLevel, thisGrid))
		return;

	int numvertneighbors;
	thisSurf->SubdivBoundaryInfo(subLevel, &numsubdivpoints, &subdivbpoints);
	thisSurf->GetVertexNeighbors(subLevel, &numvertneighbors, &vneighbors);
	thisSurf->GetNeighboringVertices(subLevel, &neighboringVerts);

	if (numsubdivpoints!=numvertneighbors)
		fprintf(stderr, "Pointlistserver2::ComputeSubdivPointCloud:  ERROR, bad number of vertices\n");
}

void Pointlist_server2::UpdateSubdivBoundaryTiles(int subLevel, Diatomgrid *thisGrid)
{
	if(!subdivSurfInitialized) 
	{
		return;
	}


	if(subLevel!=currentBoundaryTilesLevel)
	{
		currentBoundaryTilesLevel=subLevel;
		ClearBoundaryTilesPointers();
	}

	if(!UpdateSubdivSurf(subLevel, thisGrid))
		return;

	//int ntiles;
	thisSurf->SubdivTileInfo(subLevel, &numsubdivtiles, &subdivtileset);
	thisSurf->GetTileNeighbors(subLevel, &numsubdivtiles, &tneighbors);
	thisSurf->GetTileIndices(subLevel, &subdivtileindices);
}

/********************************************************************************/
/* Compute bpointlist of subdivision surface boundary generated from initial	*/
/*    diatomgrid interpolated to mesh{row,col}spacing and then sent to Subdivsurf */
/*    for Subdivision 0 <= subdivlevel < NUMSUBLEVELS and then output of		*/
/*	  resulting tiles to subdivtileset											*/
/********************************************************************************/
void Pointlist_server2::ComputeSubdivBoundaryTiles(int subLevel)
{
	if(!subdivSurfInitialized) 
	{
		if(newGrid==NULL)
			return;
		if(!InitializeSubdivSurf(subLevel))
			return;
		subdivSurfInitialized=true;
	}

	if(subLevel!=currentBoundaryTilesLevel)
	{
		currentBoundaryTilesLevel=subLevel;
		ClearBoundaryTilesPointers();
	}
	// Deal with NULL pointers instead in thisSurf->SubdivTileInfo, GetTileNeighbors and GetTileIndices
	thisSurf->SubdivTileInfo(subLevel, &numsubdivtiles, &subdivtileset);
	thisSurf->GetTileNeighbors(subLevel, &numsubdivtiles, &tneighbors);
	thisSurf->GetTileIndices(subLevel, &subdivtileindices);
}

//	new function InitializeSubdivSurf() added
//	function ComputeSubdivPointCloud totally rewritten
/********************************************************************************/
/* ComputeSubdivPointCloud() -- initialize subdivbpoints[] with subdivision		*/
/*    surface boundary vertices generated from initial diatomgrid interpolated  */
/*    to mesh{row,col}spacing and then sent to Subdivsurf for subdivision		*/
/*    0 <= subdivlevel < NUMSUBLEVELS.	(#defined in Subdivsurf.h)				*/
/********************************************************************************/
void Pointlist_server2::ComputeSubdivPointCloud(int subLevel)
{
	if(!subdivSurfInitialized) 
	{
		if(newGrid==NULL)
			return;
		if(!InitializeSubdivSurf(subLevel))
			return;
		subdivSurfInitialized=true;
	}

	if(subLevel!=currentPointCloudLevel)
	{
		currentPointCloudLevel=subLevel;
		ClearPointCloudPointers();
	}
	// deal with the NULL pointers in thisSurf->SubdivBoundaryInfo and GetVertexNeighbors instead

	int numvertneighbors;
	thisSurf->SubdivBoundaryInfo(subLevel, &numsubdivpoints, &subdivbpoints);
	thisSurf->GetVertexNeighbors(subLevel, &numvertneighbors, &vneighbors);
	thisSurf->GetNeighboringVertices(subLevel, &neighboringVerts);

	if (numsubdivpoints!=numvertneighbors)
		fprintf(stderr, "Pointlistserver2::ComputeSubdivPointCloud:  ERROR, bad number of vertices\n");
}

bool Pointlist_server2::InitializeSubdivSurf(int sublvl)
{
	int slevel;

	if(newGrid==NULL)
		return false;

	if(thisSurf!=NULL)
		delete thisSurf;
	thisSurf=new CCSubdivsurf(NUMLEVELS);
	//QuadMesh thisquadmesh(newGrid);
	if(thisMesh!=NULL) {
		delete thisMesh;
		thisMesh	= NULL;
	}
	if(newGrid->cols>1 && newGrid->rows>1) {
		thisMesh=new QuadMesh(newGrid);
	}
	else {
		thisMesh=new TubeMesh(newGrid);
		thisSurf->setIsTubeMesh(true);
	}
	// initialize the size of the figure atom grid here
	thisSurf->init(thisMesh->numverts(), thisMesh->numpolys(),
				   thisMesh->vertlist(), thisMesh->polylist());

	thisSurf->loadlevelzero();
	thisSurf->setSmoothness(smoothness);

	// to calculate the TOPOLOGY once and for all
	for(slevel = 0; slevel < (sublvl<MIN_SUBDIV_LEVEL_TO_INTERPOLATE_BP?MIN_SUBDIV_LEVEL_TO_INTERPOLATE_BP:sublvl); slevel++)
		thisSurf->splitandaverage(slevel);

	highestInitializedSubLevel=slevel-1;

	return true;
}

void Pointlist_server2::InitializeSubdivSurf(Xferlist *xferList)
{
	if(subdivSurfInitialized) 
		return;

	init(xferList);

	if(!InitializeSubdivSurf(MIN_SUBDIV_LEVEL_TO_INTERPOLATE_BP)) //if(!InitializeSubdivSurf())
		return;
	else
		subdivSurfInitialized=true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Pointlist_server2::BuildLookUpTableForPDexByUVT()
{
	CCMesh *thisMesh=thisSurf->SubMeshPointer(0);
	int pDex, vDex;
	int quadLoc;							// quadLoc		1 top, -1 bottom, 0 crest
	Medcoord uvtAvg;
	//int flag;
	//MyFloat utotal, vTotal, tTotal;		// uvtAvg.t		1 top, -1 bottom, (-1, 1) crest
	int minuvtTotal, uvtTotal, cornerVDex;//, uTmp, vTmp, tTmp;

	CCPolygon *thisPoly;
	CCVertex *thisVert;

	int maxU = (int) max_u();
	int maxV = (int) max_v();

	int fNum=thisMesh->numfaces();
	for(pDex=0; pDex<fNum; pDex++)
	{
		uvtAvg.setvals(0, 0, 0);
		//uTotal=0
		//tTotal=0;
		//uTmp=99999;
		//vTmp=99999;
		//tTmp=99999;
		minuvtTotal=99999;
		thisPoly=thisMesh->FaceList(pDex);
		for(vDex=0; vDex<thisPoly->numsides(); vDex++)
		{
			thisVert=thisMesh->VertList(thisPoly->Pvertices[vDex]);
			//tTotal+=thisVert->uvt.t;
			uvtAvg=uvtAvg+thisVert->uvt;
			uvtTotal=(int) (thisVert->uvt.u+thisVert->uvt.v+thisVert->uvt.t);
			if(minuvtTotal>uvtTotal)
			{
				minuvtTotal=uvtTotal;
				cornerVDex=vDex;
			}
		}
		uvtAvg=uvtAvg/thisPoly->numsides();
		thisVert=thisMesh->VertList(thisPoly->Pvertices[cornerVDex]);
		if(uvtAvg.t==1 || uvtAvg.t==-1)
			quadLoc = (int) uvtAvg.t;
		else
			quadLoc = 0;

		switch(quadLoc)
		{
		case 1:
			// this quad is on the top of the model
			lookUpTable[(int)(thisVert->uvt.u*maxV+thisVert->uvt.v)]=pDex;
			break;
		case -1:
			// this quad is on the bottom of the model
			lookUpTable[(int)(maxU*maxV+thisVert->uvt.u*maxV+thisVert->uvt.v)]=pDex;
			break;
		default:
			// this quad is on the crest of the model
			if(uvtAvg.v==0 || uvtAvg.v==maxV)
			{
				// case #3, #7
				if(uvtAvg.v==0)
				{
					if(uvtAvg.t>0)
						lookUpTable[(int)(2*maxU*maxV+thisVert->uvt.u)]=pDex;
					else
						lookUpTable[(int)(2*(maxU*maxV+maxU+maxV)+thisVert->uvt.u)]=pDex;

				}
				// case #5, #9
				else
				{
					if(uvtAvg.t>0)
						lookUpTable[(int)(2*maxU*maxV+maxU+maxV+thisVert->uvt.u)]=pDex;
					else
						lookUpTable[(int)(2*maxU*maxV+3*(maxU+maxV)+thisVert->uvt.u)]=pDex;
				}
			}
			if(uvtAvg.u==0 || uvtAvg.u==maxU)
			{
				// case #4, #8
				if(uvtAvg.u==maxU)
				{
					if(uvtAvg.t>0)
						lookUpTable[(int)(2*maxU*maxV+maxU+thisVert->uvt.v)]=pDex;
					else
						lookUpTable[(int)(2*maxU*maxV+2*(maxU+maxV)+maxU+thisVert->uvt.v)]=pDex;
				}
				// case #6, #10
				else
				{
					if(uvtAvg.t>0)
						lookUpTable[(int)(2*maxU*maxV+maxU+maxV+maxU+thisVert->uvt.v)]=pDex;
					else
						lookUpTable[(int)(2*maxU*maxV+3*(maxU+maxV)+maxU+thisVert->uvt.v)]=pDex;
				}
			}
			break;
		}
	}
}

int Pointlist_server2::LookUpPDexByUVT(Bpoint2 pt)
{
	int uTmp, vTmp, tTmp;
	int maxU=(int)max_u(), maxV=(int)max_v();
	uTmp=(int)pt.u;
	vTmp=(int)pt.v;
	tTmp=(int)pt.t;
	if(pt.t==-1 || pt.t==1)
	{
		// case #1, #2
		if(pt.t==1)
			return lookUpTable[uTmp*maxV+vTmp];
		else
			return lookUpTable[maxU*maxV+uTmp*maxV+vTmp];
	}
	if(pt.v==0 || pt.v==maxV)
	{
		// case #3, #7
		if(pt.v==0)
		{
			if(pt.t>0)
				return lookUpTable[2*maxU*maxV+uTmp];
			else
				return lookUpTable[2*(maxU*maxV+maxU+maxV)+uTmp];

		}
		// case #5, #9
		else
		{
			if(pt.t>0)
				return lookUpTable[2*maxU*maxV+maxU+maxV+uTmp];
			else
				return lookUpTable[2*maxU*maxV+3*(maxU+maxV)+uTmp];

		}
	}
	if(pt.u==0 || pt.u==maxU)
	{
		// case #4, #8
		if(pt.u==maxU)
		{
			if(pt.t>0)
				return lookUpTable[2*maxU*maxV+maxU+vTmp];
			else
				return lookUpTable[2*maxU*maxV+2*(maxU+maxV)+maxU+vTmp];

		}
		// case #6, #10
		else
		{
			if(pt.t>0)
				return lookUpTable[2*maxU*maxV+maxU+maxV+maxU+vTmp];
			else
				return lookUpTable[2*maxU*maxV+3*(maxU+maxV)+maxU+vTmp];

		}
	}
	bool trap=false;
	trap=true;
	return -1;
}

void Pointlist_server2::glRender(int subLvl, RenderStyle renderMode)
{
#ifdef DEBUG_QIONG_SUBDIVISION_02032003
	thisSurf->glRender(-1, WIREFRAME);
#endif
	thisSurf->glRender(subLvl, renderMode);
}


bool Pointlist_server2::subdivBPositionControlMesh(int lvl, Bpoint *Bpnt, double u, double v, double t)
{
	// you can't call this function if the subdivision surfaces have not been calculated yet!
	if(!subdivSurfInitialized) 
		return false;

	// if (u, v, t) coords out of bounds, return false
	if (u < 0.0 || u > max_u() || v < 0.0 || v > max_v() || t < -1.0 || t > 1.0)
		return false;

	CCMesh *lvlMesh[NUMLEVELS];
	CCPolygon *thisPoly, *childPoly;
	//CCVertex *quadVert[4];

	int dex;
	for(dex=0; dex<NUMLEVELS; dex++)
		lvlMesh[dex]=thisSurf->SubMeshPointer(dex);

	// from level 0, find the root polygon for the given (u, v, t)
	bool found=false;
	dex=0;
	while(dex<lvlMesh[0]->numfaces() && !found)
	{
		thisPoly=lvlMesh[0]->FaceList(dex);
		//if(WithinPolygonUVT(lvlMesh[0], thisPoly, u, v, t))
		if(WithinPolygonUVT(0, thisPoly, u, v, t))
		{
			found=true;
		}
		dex++;
	}

	// if not found, this will be a potentially BIG problem!
#ifdef QIONG_DEBUG_20030425
	bool trap;
	if(!found)
	{
		trap=true;
		printf("In Pointlist_server2::subdivBPositionControlMesh, INVALID UVT...!\n");
	}
#endif

	// from now on, we check the child polygons of the current polygon
	//		after the loop, childPoly is the pointer to the found polygon
	//		containing the current given paramters of (u, v, t)
	int subLvl;
	if(lvl>=1)
	{
		for(subLvl=1; subLvl<=lvl; subLvl++)
		{
			//CCPolygon *thispoly = oldlevel->FaceList(polydex);
			int nsides=thisPoly->numsides();
			dex=0;
			found=false;
			while(dex<nsides && !found)
			{
				childPoly=lvlMesh[subLvl]->FaceList(thisPoly->childPolygons[dex]);
				//if(WithinPolygonUVT(lvlMesh[subLvl], childPoly, u, v, t))
				if(WithinPolygonUVT(subLvl, childPoly, u, v, t))
				{
					found=true;
				}
				dex++;
			}
			if(found)
				thisPoly=childPoly;
			else
			{
				trap=true;
			}
		}
		subLvl--;
	}
	else
		subLvl=0;

	// now it's time to actually do the interpolation
	// moved it to upstairs
	// subLvl=NUMLEVELS-1;
	CCVertex *BPnt=new CCVertex();
	if(!InterpolateQuadControlMesh(BPnt, u, v, t, thisPoly, subLvl))
		return false;
	BPnt->limitnorm.selfnormalize();
	Bpnt->pnt[0]=BPnt->loc.x();
	Bpnt->pnt[1]=BPnt->loc.y();
	Bpnt->pnt[2]=BPnt->loc.z();
	Bpnt->norm[0]=BPnt->limitnorm.x();
	Bpnt->norm[1]=BPnt->limitnorm.y();
	Bpnt->norm[2]=BPnt->limitnorm.z();
	Bpnt->rad=BPnt->rad;

	delete BPnt;
	return true;
}

bool Pointlist_server2::InterpolateQuadControlMesh(CCVertex *vert, double u, double v, double t, CCPolygon *poly, int subLvl)
{
	// you can't call this function if the subdivision surfaces have not been calculated yet!
	if(!subdivSurfInitialized) 
		return false;

	// if (u, v, t) coords out of bounds, return false
	if (u < 0.0 || u > max_u() || v < 0.0 || v > max_v() || t < -1.0 || t > 1.0)
		return false;

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

	thisVert=new CCVertex;
	thatVert=new CCVertex;

	InterpolateVertexControlMesh(thisVert, uu, quadVert[0], quadVert[3]);
	InterpolateVertexControlMesh(thatVert, uu, quadVert[1], quadVert[2]);
	InterpolateVertexControlMesh(vert, vv, thisVert, thatVert);

	delete thisVert;
	delete thatVert;

	return true;
}

bool Pointlist_server2::InterpolateVertexControlMesh(CCVertex *vert, double p, CCVertex *vert1, CCVertex *vert2)
{
	vert->loc=(1.0-p)*vert1->loc+p*vert2->loc;
	vert->limitnorm=(1.0-p)*vert1->limitnorm+p*vert2->limitnorm;
	//vert->uvt=(1.0-p)*vert1->uvt+p*vert2->uvt;
	vert->rad=(1.0-p)*vert1->rad+p*vert2->rad;
	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/********************************************************************************/
/* subdivBposition() -- given (u, v, t), return Bpoint on the tiled	subdivision */
/*    generated boundary corresponding to it.  u, v, and t are double in ranges	*/
/*	  [0, max_u()],[0, max_v()], [-1, 1].										*/
/*	If successfully computed point (if legal address), return true, else false	*/
/********************************************************************************/
//bool Pointlist_server2::subdivBPositionNormal(Bpoint *BPosition, Bpoint *BNormal, double u, double v, double t)
#define MAX 99999
bool Pointlist_server2::WithinPolygonUVT(CCMesh *thisMesh, CCPolygon *thisPoly, double u, double v, double t)
{
	int dex;
	CCVertex *thisVert;
	double uMin=MAX, uMax=-MAX, vMin=MAX, vMax=-MAX, tMin=MAX, tMax=-MAX;
	double vWrapMax = -MAX, vWrapMin = MAX;
        const bool isTubeMesh = ((typeid(*this->thisMesh) == typeid(TubeMesh)) == 0 ? false : true);
	int endVertex = -1;	// set to 0 to vertex index only in case of tube mesh
							// and one of the v's is -1.
	for(dex=0; dex<thisPoly->numsides(); dex++)
	{
		thisVert=thisMesh->VertList(thisPoly->Pvertices[dex]);
		if(uMin>thisVert->uvt.u) uMin=thisVert->uvt.u;
		if(uMax<thisVert->uvt.u) uMax=thisVert->uvt.u;
		if( isTubeMesh && thisVert->uvt.v == -1.0 ) {
			endVertex	= dex;
			assert( thisVert->uvt.t == 0.0 );
			// just to ensure medial coordinate was valid.
		}
		else {
			double v	= thisVert->uvt.v;
			if( vMin > v )	vMin	= v;
			if( vMax < v )	vMax	= v;
			if( v < 0.5 )	v	= v + 1.0;
			if( vWrapMin > v )	vWrapMin	= v;
			if( vWrapMax < v )	vWrapMax	= v;
		}
		if(tMin>thisVert->uvt.t) tMin=thisVert->uvt.t;
		if(tMax<thisVert->uvt.t) tMax=thisVert->uvt.t;
	}

	if( isTubeMesh ) {
		// special cases for tube mesh
		if( endVertex != -1 && v == -1.0 ) {
			// face contains end vertex and end vertex in query.
			assert (tMin == 0.0 && uMin == uMax);
			assert (t == 0.0);	// ensure medial coordinate was valid.
			return (u == uMin);
		}
		bool wrapsAround	= false;
		if( vMax - vMin > 0.5 ) {
			// face wraps around
			if( v < vWrapMin )	v	= v + 1.0;
			vMin	= vWrapMin;
			vMax	= vWrapMax;
			wrapsAround	= true;
		}
		if( tMin != tMax ) {
			// We are on a face outside the cylindrical region.
			// Here the (2pi*v,t) coordinate system is approximately a polar
			// coordinate system. Tthis is the closest approximation
			// we can have as tube end cap polygons are not planar.
			// However, we are searching in u,v,t co-ordinate space.
			// In this case, if one of the vertices is an end vertex, only
			// one of them will be, then we take it's v coord as the 
			// mean of the others.
			assert( uMin == uMax );
			if( u != uMin || t == 1.0 ) {
				return false;
			}
			// Do inside polygon test as polygons at end caps in a tube are not
			// parameter axis-aligned.
			double *pv, *pt;
			pv	= new double[thisPoly->numsides()];
			pt	= new double[thisPoly->numsides()];
			if( endVertex != -1 ) {
				pv[endVertex]	= 0.0;
			}
			for( dex = 0; dex<thisPoly->numsides(); dex++) {
				thisVert	= thisMesh->VertList(thisPoly->Pvertices[3]);
				pt[dex]	= thisVert->uvt.t;
				if( endVertex != dex ) {
					// If this is not an end vertex.
					pv[dex]	= thisVert->uvt.v;
					if( endVertex != -1) {
						pv[endVertex]	+=  pv[dex];
					}
					if(wrapsAround && pv[dex] < 0.5) {
						pv[dex]	= pv[dex] + 1.0;
					}
				}
			}
			if( endVertex != -1 ) {
				pv[endVertex] /= (thisPoly->numsides()-1);
			}
			bool side_l	= false;
			bool side_g	= false;
			for( dex = 0; dex < thisPoly->numsides(); dex++) {
				const int i	= dex;
				const int j	= (dex+1) % thisPoly->numsides();
				const double side	= (v - pv[i])*(pt[j]-pt[i]) - (t - pt[i])*(pv[j]-pv[i]);
				if( side < 0.0 )		side_l	= true;
				else if( side > 0.0 )	side_g	= false;
				if( side_l && side_g ) {
					// point lies on different sides of different lines.
					break;
				}
			}
			delete[] pv;
			delete[] pt;
			return !side_l || !side_g;
		}
		else {
			return (t == 1.0 && u>=uMin && u<=uMax && v>=vMin && v<=vMax);
		}
	}
	else {
		// slab mesh.
		return (u>=uMin && u<=uMax && v>=vMin && v<=vMax && t>=tMin && t<=tMax);
	}
}

bool Pointlist_server2::WithinPolygonUVT(int subLvl, CCPolygon *thisPoly, double u, double v, double t)
{
	CCMesh *thisMesh=thisSurf->SubMeshPointer(subLvl);
	return WithinPolygonUVT(thisMesh, thisPoly, u, v, t);
}

int Pointlist_server2::WithinPolyDexAtSublvlUVT(int pDex, int subLvl, double u, double v, double t)
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

int Pointlist_server2::WithinPolygonAtSublvlUVT(CCPolygon *thisPoly, int subLvl, double u, double v, double t)
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

bool Pointlist_server2::subdivBposition(Bpoint *Bpnt, double u, double v, double t)
{
	return subdivBPosition(Bpnt, u, v, t);
}

int Pointlist_server2::subdivNearestBPosition(Bpoint *Bpnt, double u, double v, double t) {

		// you can't call this function if the subdivision surfaces have not been calculated yet!
	if(!subdivSurfInitialized) 
		return -1;
	// if (u, v, t) coords out of bounds, return false
	if (!thisMesh->isValidMedcoord( Medcoord(u,v,t) ) )
		return -1;

	CCMesh *lvlMesh[NUMLEVELS];
	CCPolygon *thisPoly, *childPoly;
	//CCVertex *quadVert[4];

	int dex;
	int subLvl;
	bool found;

	for(dex=0; dex<=highestInitializedSubLevel; dex++) //for(dex=0; dex<NUMLEVELS; dex++)
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
	for(subLvl=1; subLvl<=highestInitializedSubLevel; subLvl++) //for(subLvl=1; subLvl<NUMLEVELS; subLvl++)
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
	subLvl=highestInitializedSubLevel; //subLvl=NUMLEVELS-1;

	//
	// FIXME: Check for tubes, and find the nearest rather than just returning
	// some arbitrary point on the quad, when none match.
	//
	CCMesh *mesh=thisSurf->SubMeshPointer(subLvl);
	CCVertex* InterpolatedBPnt;
	int pickedUpVertex = 0;
	for( ; pickedUpVertex != 4; ++pickedUpVertex ) {
		InterpolatedBPnt = mesh->VertList(thisPoly->Pvertices[pickedUpVertex]);
		if( InterpolatedBPnt->uvt.u == u &&
			InterpolatedBPnt->uvt.v == v &&
			InterpolatedBPnt->uvt.t == t ) {
			break;
		}
	}
	if( pickedUpVertex == 4 ) {
		pickedUpVertex--;
	}

	InterpolatedBPnt->limitnorm.selfnormalize();
	Bpnt->pnt[0]=InterpolatedBPnt->limitloc.x();
	Bpnt->pnt[1]=InterpolatedBPnt->limitloc.y();
	Bpnt->pnt[2]=InterpolatedBPnt->limitloc.z();
	Bpnt->norm[0]=InterpolatedBPnt->limitnorm.x();
	Bpnt->norm[1]=InterpolatedBPnt->limitnorm.y();
	Bpnt->norm[2]=InterpolatedBPnt->limitnorm.z();
	Bpnt->rad=InterpolatedBPnt->limitRad;

	return (thisPoly->Pvertices[pickedUpVertex]);

}




// need a new parameter for the lvl needed!
bool Pointlist_server2::subdivBPosition(Bpoint *Bpnt, double u, double v, double t)
{
	// you can't call this function if the subdivision surfaces have not been calculated yet!
	if(!subdivSurfInitialized) 
		return false;
	// if (u, v, t) coords out of bounds, return false
	if (!thisMesh->isValidMedcoord( Medcoord(u,v,t) ) )
		return false;

	CCMesh *lvlMesh[NUMLEVELS];
	CCPolygon *thisPoly, *childPoly;
	//CCVertex *quadVert[4];

	int dex;
	int subLvl;
	bool found;

	for(dex=0; dex<=highestInitializedSubLevel; dex++) //for(dex=0; dex<NUMLEVELS; dex++)
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

	if(!found) {
		return false;
	}

	// from now on, we check the child polygons of the current polygon
	//		after the loop, childPoly is the pointer to the found polygon
	//		containing the current given paramters of (u, v, t)
	for(subLvl=1; subLvl<=highestInitializedSubLevel; subLvl++) //for(subLvl=1; subLvl<NUMLEVELS; subLvl++)
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
		assert(found);	// found *should* be true,
		thisPoly=childPoly;
	}

	// now it's time to actually do the interpolation
	subLvl=highestInitializedSubLevel; //subLvl=NUMLEVELS-1;

	if( u == 0.5 && v == 0.0078125 && t == 1.0 ) {
		printf("");
	}
	CCVertex InterpolatedBPnt;
	if(!InterpolateQuad(&InterpolatedBPnt, u, v, t, thisPoly, subLvl))
		return false;
	InterpolatedBPnt.limitnorm.selfnormalize();
	Bpnt->pnt[0]=InterpolatedBPnt.limitloc.x();
	Bpnt->pnt[1]=InterpolatedBPnt.limitloc.y();
	Bpnt->pnt[2]=InterpolatedBPnt.limitloc.z();
	Bpnt->norm[0]=InterpolatedBPnt.limitnorm.x();
	Bpnt->norm[1]=InterpolatedBPnt.limitnorm.y();
	Bpnt->norm[2]=InterpolatedBPnt.limitnorm.z();
	Bpnt->rad=InterpolatedBPnt.limitRad;

	return true;
}

bool Pointlist_server2::InterpolateQuad(CCVertex *vert, double u, double v, double t, CCPolygon *poly, int subLvl)
{
	// you can't call this function if the subdivision surfaces have not been calculated yet!
	if (!subdivSurfInitialized) {
		printf("%s:%d Pointlist_server2::InterpolateQuad(...):\n"
			"subdifSurf has not been initialized.\n", __FILE__, __LINE__ );
		return false;
	}

	// if (u, v, t) coords out of bounds, return false
	Medcoord uvt(u,v,t);
	if (!thisMesh->isValidMedcoord(uvt) ) {
		printf("%s:%d Pointlist_server2::InterpolateQuad(...):\n"
			"An invalid medial coordinate (%f,%f,%f) was passed in.\n",
			__FILE__, __LINE__, u, v, t );
		return false;
	}

	int dex;
	CCMesh *mesh=thisSurf->SubMeshPointer(subLvl);
        const bool isTubeMesh = ((typeid(*this->thisMesh) == typeid(TubeMesh)) == 0 ? false : true);
	if (isTubeMesh) {
		//
		// Interpolate based on "distance" from end-points.
		//
		static double* w	= NULL;
		static double* dist	= NULL;
		static int numverts	= -1;
		if(poly->numsides() != numverts) {
			// Do not keep deleting and recreating weights array
			// and getting a performance penalty in the process.
			// NOTE: This will cause problems if someone calls this
			// function from two separate threads.
			if(w) {
				delete[] w;
				delete[] dist;
			}
			numverts	= poly->numsides();
			w		= new double[numverts];
			dist	= new double[numverts];
		}
		double sum_w	= 0.0;
		const CCVertex* poly_vertex;
		int i;
		for (dex = 0; dex != poly->numsides(); ++dex) {
			dist[dex]	= uvt.distance( mesh->VertList(poly->Pvertices[dex])->uvt,
				Medcoord::TubeMedCoord );
			w[dex]	= 1.0;
		}
		// This could be built in O(n log n) time easily instead of O(n2)
		// as we are dealing with small n, who cares ...
		for (dex = 0; dex != poly->numsides(); ++dex) {
			const double d	= dist[dex];
			for( i = 0; i < dex; ++i ) {
				w[i]	*= d;
			}
			for( i = dex+1; i < poly->numsides(); ++i ) {
				w[i]	*= d;
			}
		}
		for (dex = 0; dex != poly->numsides(); ++dex) {
			sum_w	+= w[dex];
		}
		for (dex = 0; dex != poly->numsides(); ++dex) {
			w[dex]	/= sum_w;
		}
		vert->limitloc.set(0.0,0.0,0.0);
		vert->limitnorm.set(0.0,0.0,0.0);
		vert->limitRad	= 0.0;
		for (dex = 0; dex != poly->numsides(); ++dex) {
			poly_vertex		= mesh->VertList(poly->Pvertices[dex]);
			vert->limitloc	+= w[dex] * poly_vertex->limitloc;
			vert->limitnorm	+= w[dex] * poly_vertex->limitnorm;
			vert->limitRad	+= w[dex] * poly_vertex->limitRad;
		}
	}
	else {
		CCVertex *quadVert[4];
		int startPos;
		CCVertex *thisVert, *thatVert;
		thisVert=mesh->VertList(poly->Pvertices[0]);
		thatVert=mesh->VertList(poly->Pvertices[1]);
		double uu, vv;

		int caseNum;
		startPos=0;
		assert( poly->numsides() == 4 );
		//
		// need to figure out how to interpolate different part of the patches
		// because 
		//		t is 1 or -1
		//			interpolate by (u, v)
		//		else
		//			interpolate by (v, t) or (u, t)
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
		CCVertex thisVert2;
		CCVertex thatVert2;
		InterpolateVertex(&thisVert2, uu, quadVert[0], quadVert[3]);
		InterpolateVertex(&thatVert2, uu, quadVert[1], quadVert[2]);
		InterpolateVertex(vert, vv, &thisVert2, &thatVert2);
	}
	return true;
}

bool Pointlist_server2::InterpolateVertex(CCVertex *vert, double p, CCVertex *vert1, CCVertex *vert2)
{
	vert->limitloc=(1.0-p)*vert1->limitloc+p*vert2->limitloc;
	vert->limitnorm=(1.0-p)*vert1->limitnorm+p*vert2->limitnorm;
	vert->limitRad=(1.0-p)*vert1->limitRad+p*vert2->limitRad;
	return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Initialize for rendering a surface directly from a list of points.
// See also M3DObjectSurfaceRenderer::renderFigureFromControlPoints().
//     Initialize the level-0 subdivision control mesh
// and subdivide it into finer meshes.
void Pointlist_server2::init(int numverts, int numfaces, Bpoint *pntlist, PindexList *faces, int sublvl)
{
   int slevel;

   // hard-coded default value, shouldn't matter for the visualization purpose
   smoothness = 90;

   if (thisSurf!=NULL)
       delete thisSurf;
   thisSurf = new CCSubdivsurf(NUMLEVELS);

   // initialize the level 0 control mesh
   thisSurf->init(numverts, numfaces, pntlist, faces);
   thisSurf->loadlevelzero();
   thisSurf->setSmoothness(smoothness);

   // to calculate the TOPOLOGY
   for (slevel = 0; slevel < (sublvl<MIN_SUBDIV_LEVEL_TO_INTERPOLATE_BP?MIN_SUBDIV_LEVEL_TO_INTERPOLATE_BP:sublvl); slevel++)
       thisSurf->splitandaverage(slevel);

   highestInitializedSubLevel = slevel-1;

   subdivSurfInitialized = true;
}


