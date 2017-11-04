/**********************************************************\
 * TubeMesh.cpp
 *		Rohit Saboo, Derek Merck, Qiong Han, Shawn Liu & A. Thall
 *      (Based on code original written by A. Thall, modified by Shawn Liu,
 *		 dug up and repaired by Derek Merck with help from Qiong Han and
 *		 Tom Fletcher. Further dug up and repaired by Rohit. Finally lost of track
 *       of changes here, look into svn log files.)
 *
 *      17. July 2002
 *
 *		A tubemesh class as interface from Pointlist_server2 to CCSubdivsurf.
 *		It takes a 1 x n Diatomgrid or Xferlist and outputs vertexlist and
 *		tile-index list.
 *
 *		Modifications:
 *      see svn log for further modifications.
 *      9/25/05: Rohit - Added a base Mesh class, did a complete
 *                overhaul and made sure it works.
 *      2/16/04:  Plugged back into code, added switch to PointList_server2
 *                  that checks for tube conditions before calling the level 0 mesh.
 *                Normals are inverted.
 *      2/23/04:  Fixed normal inversion (draw quads counter-clockwise)
 *                Fixed end-caps to be quads instead of tris so optimizer
 *					doesn't crash, but it still doesn't work properly.
 *      5/14/04:  Fixed parameterization, determined fixes required in Seurat's
 *                  sub-division method.
 *
\**********************************************************/
/*  Notes regarding parameterization: 
 *  
 *  Tube vertexes are parameterized slightly differently than those of regular slabs.
 *  
 *  Given:
 *  S = number of atoms along the tube
 *  N = number of samples around the tube
 *  (these variable names are used throughout the code)
 *  
 *
 *  ** TUBE BODY ** 
 *
 *  u runs the length of the tube in the range [0,S-1]
 *  v runs about the circumference of the tube such that v = n/N = n radians/(2*pi)
 *    in the range [0,1],
 *    also, v = 0 AND 1 at the 2*pi point.
 *  t is fixed at 1 on the tube body, except at the end caps.
 *  
 *     <----- u -----> 
 *  /--------\----------\    
 *  |         |v         |
 *  \--------/----------/
 *  
 *
 *  ** ENDCAPS ** 
 *
 *  u is fixed for a cross-section of the tube, such as the end-cap
 *  v runs about the circumference, but,
 *    v = -1 at the center.  It's value is actually entirely arbitrary, since it doesn't
 *        get used, but it's a useful key to make it -1
 *  t runs in the range [0,1], from 0 at the center of the disc cross-section to 1 at the edge.
 *    This is very much like the tao term.
 *  
 *            v=0.25
 *             t=1
 *          /---*---\
 *         /         \
 *        |      v=-1 |  
 *  v=0.5 *     *t=0  * v=0/1
 *   t=1  |           | t=1
 *         \         /
 *          \---*---/
 *           v=0.75
 *             t=1
 *
 *  Interpolating over this parameterization requires fixes to Pointlist_server2 (at least):
 *
 *  1. v_new = average(-1,v_old) = v_old
 *     that is, any interpolation between a v and the vertex at the center
 *     leaves v as is.
 *
 *  2. v_new = average(v0,v_new) = average(0,v_old) if v_old < 0.5 or
 *								   average(1,v_old) if v_old > 0.5
 *     that is, any interpolation using the v0 term (the starting v, where it is both 0
 *     and 2*pi) requires averaging with the correct form of v0.
 *
\********************************************************************************/

#define D_MESH
#define D_TUBEMESH
#define D_XFERLIST
#define D_PSEUDOSET
#define D_DIATOMGRID
#include "Shapedepend.h"
#include <stdio.h>
#include "M3DFigure.h"

using namespace ThallCode;


/**********************************************************\
 * isValidMedcoord
 *
 * returns true/false based on whether <u,v,t> is valid or not
\**********************************************************/
bool TubeMesh::isValidMedcoord( const Medcoord& uvt )
{
	return uvt.u >= 0.0 && uvt.u <= TubeLength - 1.0
			&& ((uvt.v >= 0.0 && uvt.v <= 1.0) ||
				((uvt.u == 0.0 || uvt.u == TubeLength - 1.0) && uvt.v == -1.0 ) )
			&& (uvt.t >= 0.0 && uvt.t <= 1.0 );
}

/**********************************************************\
 * isEndMedCoord
 * returns true false based on whether <u,v,t> is an end
 * co-ordinate or not.
 *
\**********************************************************/
bool TubeMesh::isEndMedcoord( const Medcoord& uvt )
{
	return 0.0 <= uvt.t && uvt.t < 1.0;
}


/**********************************************************\
 * average
 * returns the average of n medial coordinates.
\**********************************************************/
Medcoord TubeMesh::average(int n, const Medcoord* uvt )
{
	Medcoord newuvt;
	newuvt.average(n, uvt, Medcoord::TubeMedCoord);
	assert(isValidMedcoord(newuvt));
	return newuvt;
}

/**********************************************************\
 * init() -- initialize the Pointlist_server2,
 * passing in either an Xferlist or an already created
 * Subdivsurf; this creates a local newgrid.
\**********************************************************/
TubeMesh::TubeMesh(Xferlist *thislist) : Mesh()
{

	// Should probably fix this to take column-wise tubes, too.
	// Load values from (row-major-ordered) Xferlist mfig1 to a Diatomgrid
	//   (column-major-ordered)
	fig_grid = new Diatomgrid;
	fig_grid->readXferlist(thislist);

	TubeLength = fig_grid->rows;
	if (fig_grid->cols != 1)
		fprintf(stderr, "WARNING: TubeMesh::TubeMesh:  more than 1 col in thisgrid (ignored).\n");

	constructMesh0();
	computePolyList();
}

TubeMesh::TubeMesh(Diatomgrid *thisgrid) : Mesh()
{
    // Assign figuregrid
	if (thisgrid != NULL) {
		fig_grid = new Diatomgrid;
		fig_grid->rows = thisgrid->rows;
		fig_grid->cols = thisgrid->cols;
		fig_grid->figure = thisgrid->figure;
		// copy entries into internal fig_grid
		fig_grid->dlist = new Diatom[thisgrid->rows*thisgrid->cols];
		for (int cnt = 0; cnt < thisgrid->rows*thisgrid->cols; cnt++)
			fig_grid->dlist[cnt] = thisgrid->dlist[cnt];
	}
	else {
		fig_grid = new Diatomgrid;
		fig_grid->rows = 0;
		fig_grid->cols = 0;
		fig_grid->dlist = NULL;
		fig_grid->figure = NULL;
	}
	TubeLength = fig_grid->rows;
	if (fig_grid->cols != 1)
		fprintf(stderr, "WARNING: TubeMesh::TubeMesh:  more than 1 col in thisgrid (ignored).\n");

	constructMesh0();
	computePolyList();
}


/********************************************************************************/
/* TubeMesh::constructMesh0() -- loads level = 0 mesh from the Diatomgrid,	    */
/*   and construct vertex list upon them.      								    */
/* input n is number of vertex-samples per atom (not including endcaps)			*/
/********************************************************************************/
void TubeMesh::constructMesh0()
{
	//
	// FIXME: This function kind of works for now, but it needs to be completely overhauled
	// It's way toooo messy - rrs
	//

	const int N	= dynamic_cast<const M3DTubeFigure*>(fig_grid->figure)->getNumberOfSpokes();
	// N = about the ring, v values
	// S = along the tube, u values
	// t values all 1 except at center vertex in endcap

    Diatom *atom = fig_grid->dlist;
    Diatom workatom, nextatom, lastatom;

    DbVector3 p, v1, b, v1pnt, bpnt, v1norm;
	double uval, vval, tval, theta;
    double radius, tvalStart, tvalEnd;

	numdivs = N;
	const int nverts = N*TubeLength + 2;
	if( nverts != this->nverts ) {
		if( v_list != NULL ) {
			delete[] v_list;
			v_list	= NULL;
		}
		this->nverts	= nverts;
	}
	if( v_list == NULL ) {
		v_list = new Bpoint[nverts];
	}

	tvalStart = 0.0;
	tvalEnd = 1.0;

	Quat q;
	theta = 2 * M_PI/N;  // step size along arc = 2*pi/num steps

	vval = 0;
	int vertdex = 0;

	// s ranges along the number of atoms we have.
	for (int s = 0; s < TubeLength; s++) {
	// for each atom along the length of the tube

		workatom = atom[s];		  // this atom
		p = workatom.p_val();     // location of this atom
	    radius = workatom.r_val();// radius along perpendicular
		uval = (double) s;

		const M3DTubePrimitive* prim	= dynamic_cast<const M3DTubePrimitive*>(fig_grid->figure->getPrimitivePtr(s));

		// n ranges about samples
		for (int n = 0; n < N; n++) {
			// For each of the n spokes

			const Vector3D yN	= prim->getYN(n);
			v1.X()	= yN.getX();
			v1.Y()	= yN.getY();
			v1.Z()	= yN.getZ();
			v1pnt	= p + v1;

			vval = (double)n/((double)N);  // v/column is actually spin
			tval = 1.0;                    // Always on the boundary

			v1norm = v1.normalize();  // Normalize the spoke (do we have to do this over and over?)
			v_list[vertdex++].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1norm.x(), v1norm.y(), v1norm.z(), radius, 
			              uval, vval, tval );
		}
	}

	// Now assign endcap vertices in the last two spots
	workatom = atom[0];
	p = workatom.p_val();
	b = workatom.bvector();
	bpnt = p + b;
	radius = workatom.r_val();
	uval = 0.0;
	vval = -1;  // Negative 1 indicates that it's a center atom
	tval = 0.0; // 0 indicates that it's a center atom
	b.selfnormalize();

	v_list[vertdex++].setvals(bpnt.x(), bpnt.y(), bpnt.z(), b.x(), b.y(), b.z(), radius, uval, vval, tval);

	workatom = atom[TubeLength - 1];
	p = workatom.p_val();
	b = workatom.bvector();
	bpnt = p + b;
	radius = workatom.r_val();
	uval = (double) TubeLength - 1;		// These indices are screwy.  Think about what makes sense.
	vval = -1;  // Negative 1 indicates that it's a center atom
	tval = 0.0; // 0 indicates that it's a center atom
	b.selfnormalize();
	v_list[vertdex++].setvals(bpnt.x(), bpnt.y(), bpnt.z(), b.x(), b.y(), b.z(), radius, uval, vval, tval);

	if (nverts != vertdex)
		fprintf(stderr, "ERROR:  TubeMesh::constructMesh0() -- incorrect nvert count\n");
}

/********************************************************************************/
/* TubeMesh::computePolyList() -- compute the polygon list upon the vertex list */
/*available.                                                                    */
/********************************************************************************/
void TubeMesh::computePolyList()
{

	// N = about the ring
	// S = along the tube

	int N = numdivs;
	npolys = (TubeLength)*N;			// (TubeLength - 1)*N quads on sides,
										// and N/2 quads on each endcap = one more N
	pi_list = new PindexList[npolys];

	CCPindex pdex = 0;

	int s, n;
	for (s = 0; s < TubeLength - 1; s++) {
		// Create quads along the tube
		for (n = 0; n < N; n++) {
			pi_list[pdex].Create(4);	// create PindexList of length 4
			pi_list[pdex][3] = s*N + n;
            // This atom * numdivs (= sum to this ring) + where we are on this ring
			pi_list[pdex][2] = s*N + INTMOD(n + 1, N);
            // This atom * numdivs (= sum to this ring) + next vert around on this ring (mod pushes back to 0 at the end)
			pi_list[pdex][1] = (s + 1)*N + INTMOD(n + 1, N);
            // [1]'s neighbor one ring higher 
			pi_list[pdex++][0] = (s + 1)*N + n;
			// [0]'s neighbor one ring higher
		}
	}

	CCVindex end0 = nverts - 2;
	CCVindex endn = nverts - 1;

/*	// Can't be tri's!
	// Reversed the order of these tiles to get them to be counter-clockwise
	atomcnt = 0;
	for (ringcnt = 0; ringcnt < numdivs; ringcnt++) {
		// Tile the bottom cap, go around the ring
		pi_list[pdex].Create(3);	// create PindexList of length 3
		pi_list[pdex][2] = atomcnt*numdivs + ringcnt;
		pi_list[pdex][1] = end0;  // Connect it to the center point
		pi_list[pdex++][0] = atomcnt*numdivs + INTMOD(ringcnt + 1, numdivs);  // 
	}
	atomcnt = TubeLength - 1;
	for (ringcnt = 0; ringcnt < numdivs; ringcnt++) {

		pi_list[pdex].Create(3);	// create PindexList of length 3
		pi_list[pdex][2] = atomcnt*numdivs + INTMOD(ringcnt + 1, numdivs);
		pi_list[pdex][1] = endn;
		pi_list[pdex++][0] = atomcnt*numdivs + ringcnt;
	}
*/

	// Made quad end-caps because it doesn't like tris (perhaps it would with the
	// parameterization fixes?
	s = 0;
	for (n = 0; n < N; n+=2) {
		// Tile the bottom cap, go around the ring
		pi_list[pdex].Create(4);				// create PindexList of length 4 (no tris!)
		pi_list[pdex][3] = s*N + n;
		pi_list[pdex][2] = end0;				// Connect it to the center point
		pi_list[pdex][1] = s*N + INTMOD(n + 2, N);
		pi_list[pdex++][0] = s*N + INTMOD(n + 1, N); 
	}
	s = TubeLength - 1;
	for (n = 0; n < N; n+=2) {

		pi_list[pdex].Create(4);				// create PindexList of length 4 (no tris!)
		pi_list[pdex][3] = s*N + INTMOD(n + 1, N);
		pi_list[pdex][2] = s*N + INTMOD(n + 2, N); 
		pi_list[pdex][1] = endn;				// Connect it to the center point
		pi_list[pdex++][0] = s*N + n;
	}


	if (npolys != pdex)
		fprintf(stderr, "ERROR:  TubeMesh::constructMesh0() -- incorrect nvert count\n");
}

