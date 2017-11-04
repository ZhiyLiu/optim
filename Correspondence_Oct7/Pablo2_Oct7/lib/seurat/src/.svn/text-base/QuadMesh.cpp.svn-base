/**********************************************************\
 * QuadMesh.cpp
 * A. Thall
 * 18. July 2002
 *
 * A quadmesh class as interface from Pointlist_server2 to
 * CCSubdivsurf. It takes a nxm Diatomgrid or Xferlist and
 * outputs vertexlist and tile-index list.
 *
 * Modifications:
 *  2005/9/25: Rohit - Moved functionality into a base class
 *    to make provision for supporting different meshes.
 *
\**********************************************************/

#define D_MESH
#define D_QUADMESH
#define D_XFERLIST
#define D_PSEUDOSET
#define D_DIATOMGRID
#include "Shapedepend.h"
#include <stdio.h>

using namespace ThallCode;

/* The number of vertex-samples per atom
 * (not including endcaps)
 */
const int QuadMesh::NUM_SIDE_ROW	= 1;


/**********************************************************\
 * isValidMedcoord
 *
 * returns true/false based on whether <u,v,t> is valid or not
\**********************************************************/
bool QuadMesh::isValidMedcoord( const Medcoord& uvt )
{
	return !( uvt.u < 0.0 || uvt.u > numcols - 1
			|| uvt.v < 0.0 || uvt.v > numrows - 1
			|| uvt.t < -1.0 || uvt.t > 1.0 );
}

/**********************************************************\
 * isEndMedCoord
 * returns true false based on whether <u,v,t> is an end
 * co-ordinate or not.
 *
\**********************************************************/
bool QuadMesh::isEndMedcoord( const Medcoord& uvt )
{
	return -1.0 < uvt.t && uvt.t < +1.0 ;
}

/**********************************************************\
 * average
 * returns the average of n medial coordinates.
\**********************************************************/
Medcoord QuadMesh::average(int n, const Medcoord* uvt )
{
	Medcoord newuvt;
	newuvt.average(n, uvt, Medcoord::SlabMedCoord);
	assert(isValidMedcoord(newuvt));
	return newuvt;
}

/**********************************************************\
 * init() -- initialize the Pointlist_server2,
 * passing in either an Xferlist or an already created
 * Subdivsurf; this creates a local newgrid.
\**********************************************************/
QuadMesh::QuadMesh(Xferlist *thisList) : Mesh()
{
    std::cout << "QuadMesh" << std::endl;
	if(thisList!=NULL) 
	{
		fig_grid = new Diatomgrid;
		fig_grid->readXferlist(thisList);	// the dlist will be assigned in here!
		//fig_grid->rows = thisgrid->rows;
		//fig_grid->cols = thisgrid->cols;
		// copy entries into internal fig_grid
		//fig_grid->dlist = new Diatom[thisgrid->rows*thisgrid->cols];
		//CopyDiatomGrid(thisgrid);
		/*
		for(int cnt = 0; cnt < thisgrid->rows*thisgrid->cols; cnt++)
			fig_grid->dlist[cnt] = thisgrid->dlist[cnt];
		*/
	}
	else 
	{
		fig_grid = new Diatomgrid;
		fig_grid->rows = 0;
		fig_grid->cols = 0;
		fig_grid->dlist = NULL;
	}

	numdivs = NUM_SIDE_ROW;
	numrows = fig_grid->rows;
	numcols = fig_grid->cols;
	numsiderows = NUM_SIDE_ROW;
	numsidecols = 2*numrows + 2*numcols - 4;
	if (numcols < 2 || numrows < 2)
		fprintf(stderr, "WARNING: QuadMesh::QuadMesh:  Not a legitimate QuadMesh.\n");

	nverts = 2*numrows*numcols + numsiderows*numsidecols;
	v_list = new Bpoint[nverts];
	int numprows = numrows - 1;
	int numpcols = numcols - 1;
	int numsideprows = numsiderows + 1;
	int numsidepcols = numsidecols;		// because of wrap-around, same # of polygons
										//   as vertices around the side
	npolys = 2*numprows*numpcols + numsideprows*numsidepcols;
	pi_list = new PindexList[npolys];

	UpdateFaceList();
	UpdateVertexList();
}

// qiong
QuadMesh::QuadMesh(Diatomgrid *thisgrid) : Mesh()
{
    // Assign figuregrid
	if(thisgrid!=NULL) 
	{
		fig_grid = new Diatomgrid;
		fig_grid->rows = thisgrid->rows;
		fig_grid->cols = thisgrid->cols;
		// copy entries into internal fig_grid
		fig_grid->dlist = new Diatom[thisgrid->rows*thisgrid->cols];
		CopyDiatomGrid(thisgrid);
		/*
		for(int cnt = 0; cnt < thisgrid->rows*thisgrid->cols; cnt++)
			fig_grid->dlist[cnt] = thisgrid->dlist[cnt];
		*/
	}
	else 
	{
		fig_grid = new Diatomgrid;
		fig_grid->rows = 0;
		fig_grid->cols = 0;
		fig_grid->dlist = NULL;
	}

	if(v_list!=NULL)
		delete []v_list;
	v_list = NULL;
	if(pi_list!=NULL)
		delete []pi_list;
	pi_list = NULL;
	numdivs = NUM_SIDE_ROW;
	numrows = fig_grid->rows;
	numcols = fig_grid->cols;
	numsiderows = NUM_SIDE_ROW;
	numsidecols = 2*numrows + 2*numcols - 4;
	if (numcols < 2 || numrows < 2)
		fprintf(stderr, "WARNING: QuadMesh::QuadMesh:  Not a legitimate QuadMesh.\n");

	nverts = 2*numrows*numcols + numsiderows*numsidecols;
	v_list = new Bpoint[nverts];
	int numprows = numrows - 1;
	int numpcols = numcols - 1;
	int numsideprows = numsiderows + 1;
	int numsidepcols = numsidecols;		// because of wrap-around, same # of polygons
										//   as vertices around the side
	npolys = 2*numprows*numpcols + numsideprows*numsidepcols;
	pi_list = new PindexList[npolys];

	UpdateFaceList();
	UpdateVertexList();
}

/********************************************************************************/
/* QuadMesh::constructMesh0() -- loads level = 0 mesh from the Diatomgrid,	    */
/*   and construct vertex list upon them.      								    */
/* input n is number of vertex-samples per crest (not including v1 & v2)		*/
/*   eta will be linearly interpolated along crest								*/
/********************************************************************************/
void QuadMesh::constructMesh0()
{
	const int n	= NUM_SIDE_ROW;
    Diatom *atom = fig_grid->dlist;
    Diatom workatom, nextatom, lastatom;

    DbVector3 p, v1, v2, b, bperp, v1pnt, v2pnt, bpnt, v1norm, v2norm, vEnd, vEndpnt, uEnd;
	double uval, vval, tval, etaval, newetaval;
    double radius, t_decrement;
	int vertindex;

	int B = numrows*numcols;		// first bottom atom index
	int S = 2*numrows*numcols;		// first side index
	int E = numsidecols;			// rowlength along side

	// First, initialize top and bottom of mesh
	int col, row, atomindex;
	for (col = 0; col < numcols; col++) {
		for (row = 0; row < numrows; row++) {

			atomindex = col*numrows + row;
			workatom = atom[atomindex];
			p = workatom.p_val();
			v1 = workatom.v1_val();
			v2 = workatom.v2_val();

			v1pnt = p + v1;
			v2pnt = p + v2;
			radius = workatom.r_val();
			v1.selfnormalize();
			v2.selfnormalize();

			uval = (double) col;
			vval = (double) row;

			v_list[atomindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1.x(), v1.y(), v1.z(), workatom.r1_val(), 
					uval, vval, 1.0);
			v_list[atomindex + B].setvals(v2pnt.x(), v2pnt.y(), v2pnt.z(), v2.x(), v2.y(), v2.z(), workatom.r2_val(), 
					uval, vval, -1.0);
		}
	}

	// Now, initialize side of mesh, one edge at a time

	// column zero, all but last atom (row == numrows - 1)
	int siderow;
	Quat q;
	double theta, rottheta;

	col = 0;
	for (row = 0; row < numrows - 1; row++) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();
		radius = workatom.r_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);

		//q.rot(bperp, rottheta);
		q.setAxisAngle(Vector3D(bperp.X(), bperp.Y(), bperp.Z()), rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);                
                
		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + row;
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			//q.rotate_vec(v1);
			Vector3D v1Tmp(v1.X(), v1.Y(), v1.Z()); 
			q.rotateVector(v1Tmp);
			v1.set(v1Tmp.getX(), v1Tmp.getY(), v1Tmp.getZ());

			v1.selfnormalize();

			v1pnt = p + v1*radius*newetaval;
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			/*v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1norm.x(), v1norm.y(), v1norm.z(), radius, 
					uval, vval, tval);*/
                        p = workatom.p_val();
			uEnd = workatom.uEnd_val();
                        radius = workatom.rEnd_val();
                        vEnd = uEnd * radius;
                    

			vEndpnt = p + vEnd;
                        vEnd.selfnormalize();

			v_list[vertindex].setvals(vEndpnt.x(), vEndpnt.y(), vEndpnt.z(), vEnd.x(), vEnd.y(), vEnd.z(), radius,
                                        uval, vval, tval);
		}
	}

	// row (numrows - 1), all but last atom (col == numcols - 1)
	row = numrows - 1;
	for (col = 0; col < numcols - 1; col++) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();
		radius = workatom.r_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		//q.rot(bperp, rottheta);
		q.setAxisAngle(Vector3D(bperp.X(), bperp.Y(), bperp.Z()), rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + (numrows - 1) + col;
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			//q.rotate_vec(v1);
			Vector3D v1Tmp(v1.X(), v1.Y(), v1.Z()); 
			q.rotateVector(v1Tmp);
			v1.set(v1Tmp.getX(), v1Tmp.getY(), v1Tmp.getZ());

			v1.selfnormalize();

			v1pnt = p + v1*radius*newetaval;
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			/*v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1norm.x(), v1norm.y(), v1norm.z(), radius, 
					uval, vval, tval);*/
                        p = workatom.p_val();
			uEnd = workatom.uEnd_val();
                        radius = workatom.rEnd_val();
                        vEnd = uEnd * radius;
                    

			vEndpnt = p + vEnd;
                        vEnd.selfnormalize();

			v_list[vertindex].setvals(vEndpnt.x(), vEndpnt.y(), vEndpnt.z(), vEnd.x(), vEnd.y(), vEnd.z(), radius,
                                        uval, vval, tval);
		}
	}

	// col (numcols - 1), all but last atom (row == 0);
	col = numcols - 1;
	for (row = numrows - 1; row > 0; row--) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();
		radius = workatom.r_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		//q.rot(bperp, rottheta);
		q.setAxisAngle(Vector3D(bperp.X(), bperp.Y(), bperp.Z()), rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + (numrows - 1) + (numcols - 1) + (numrows - 1 - row);
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			//q.rotate_vec(v1);
			Vector3D v1Tmp(v1.X(), v1.Y(), v1.Z()); 
			q.rotateVector(v1Tmp);
			v1.set(v1Tmp.getX(), v1Tmp.getY(), v1Tmp.getZ());

			v1.selfnormalize();

			v1pnt = p + v1*radius*newetaval;
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			/*v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1norm.x(), v1norm.y(), v1norm.z(), radius, 
					uval, vval, tval);*/
                        p = workatom.p_val();
			uEnd = workatom.uEnd_val();
                        radius = workatom.rEnd_val();
                        vEnd = uEnd * radius;
                    

			vEndpnt = p + vEnd;
                        vEnd.selfnormalize();

			v_list[vertindex].setvals(vEndpnt.x(), vEndpnt.y(), vEndpnt.z(), vEnd.x(), vEnd.y(), vEnd.z(), radius,
                                        uval, vval, tval);
		}
	}

	// row 0, all but last atom (col == 0);
	row = 0;
	for (col = numcols - 1; col > 0; col--) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();
		radius = workatom.r_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		//q.rot(bperp, rottheta);
		q.setAxisAngle(Vector3D(bperp.X(), bperp.Y(), bperp.Z()), rottheta);


		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + 2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col);
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			//q.rotate_vec(v1);
			Vector3D v1Tmp(v1.X(), v1.Y(), v1.Z()); 
			q.rotateVector(v1Tmp);
			v1.set(v1Tmp.getX(), v1Tmp.getY(), v1Tmp.getZ());

			v1.selfnormalize();

			v1pnt = p + v1*radius*newetaval;
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			/*v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1norm.x(), v1norm.y(), v1norm.z(), radius, 
					uval, vval, tval);*/
                        p = workatom.p_val();
			uEnd = workatom.uEnd_val();
                        radius = workatom.rEnd_val();
                        vEnd = uEnd * radius;
                    

			vEndpnt = p + vEnd;
                        vEnd.selfnormalize();

			v_list[vertindex].setvals(vEndpnt.x(), vEndpnt.y(), vEndpnt.z(), vEnd.x(), vEnd.y(), vEnd.z(), radius,
                                        uval, vval, tval);
		}
	}
	/*
	// Now, initialize side of mesh, one edge at a time

	// column zero, all but last atom (row == numrows - 1)
	int siderow;
	Quat q;
	double theta, rottheta;

	col = 0;
	for (row = 0; row < numrows - 1; row++) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		b=workatom.bvector();

		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		q.rot(bperp, rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + row;
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			q.rotate_vec(v1);
			v1.selfnormalize();

			v1pnt = p + b;
			//b=v1*radius*newetaval;
			b.selfnormalize();
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)
			b=v1*radius*newetaval;

			v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), b.x(), b.y(), b.z(), radius, 
					uval, vval, tval);
		}
	}

	// row (numrows - 1), all but last atom (col == numcols - 1)
	row = numrows - 1;
	for (col = 0; col < numcols - 1; col++) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		b=workatom.bvector();

		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		q.rot(bperp, rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + (numrows - 1) + col;
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			q.rotate_vec(v1);
			v1.selfnormalize();

			v1pnt = p + b;
			b.selfnormalize();
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			b=v1*radius*newetaval;
			v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), b.x(), b.y(), b.z(), radius, 
					uval, vval, tval);
		}	}

	// col (numcols - 1), all but last atom (row == 0);
	col = numcols - 1;
	for (row = numrows - 1; row > 0; row--) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		b=workatom.bvector();

		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		q.rot(bperp, rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + (numrows - 1) + (numcols - 1) + (numrows - 1 - row);
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			q.rotate_vec(v1);
			v1.selfnormalize();

			v1pnt = p + b;
			b.selfnormalize();
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			b=v1*radius*newetaval;
			v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), b.x(), b.y(), b.z(), radius, 
					uval, vval, tval);
		}	}

	// row 0, all but last atom (col == 0);
	row = 0;
	for (col = numcols - 1; col > 0; col--) 
	{

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		b=workatom.bvector();

		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		q.rot(bperp, rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + 2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col);
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			q.rotate_vec(v1);
			v1.selfnormalize();

			v1pnt = p + b;
			b.selfnormalize();
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			b=v1*radius*newetaval;
			v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), b.x(), b.y(), b.z(), radius, 
					uval, vval, tval);
		}	
	}
	*/
}

/********************************************************************************/
/* QuadMesh::computePolyList() -- compute the polygon list upon the vertex list */
/*available.                                                                    */
/********************************************************************************/
void QuadMesh::computePolyList()
{
	int numprows = numrows - 1;
	int numpcols = numcols - 1;
	int numsideprows = numsiderows + 1;
	int numsidepcols = numsidecols;		// because of wrap-around, same # of polygons
										//   as vertices around the side

	int polyB = numprows*numpcols;
	int polyS = 2*numprows*numpcols;
	int polyE = numsidecols;

	int B = numrows*numcols;
	int S = 2*numrows*numcols;
	int E = numsidecols;

	CCPindex pdex;
	int vdexA, vdexB, vdexC, vdexD;

	int row, col;

	// Do top and bottom surface polygons first
	for (col = 0; col < numpcols; col++) 
	{
		for (row = 0; row < numprows; row++) 
		{

			vdexA = col*numrows + row;
			vdexB = col*numrows + row + 1;
			vdexC = (col + 1)*numrows + row + 1;
			vdexD = (col + 1)*numrows + row;

			pdex = col*numprows + row;

			// Add vertex indices of top and bottom corresponding faces
			pi_list[pdex].Create(4);	// create PindexList of length 4
			pi_list[pdex][0] = vdexA;
			pi_list[pdex][1] = vdexB;
			pi_list[pdex][2] = vdexC;
			pi_list[pdex][3] = vdexD;

			pi_list[pdex + polyB].Create(4);
			pi_list[pdex + polyB][0] = vdexA + B;
			pi_list[pdex + polyB][1] = vdexD + B;
			pi_list[pdex + polyB][2] = vdexC + B;
			pi_list[pdex + polyB][3] = vdexB + B;
		}
	}

	// Now, initialize side of mesh, one edge at a time
	// First, top row of polygons

	// pdex will be used as a counter index, now
	pdex = polyS;
	// column zero, all but last atom (row == numrows - 1)
	col = 0;
	for (row = 0; row < numrows - 1; row++) {

		vdexA = col*numrows + row;
		vdexB = S + row;
		vdexC = S + row + 1;
		vdexD = col*numrows + row + 1;

		// Add vertex indices of top corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;

		vdexA = B + col*numrows + row;
		vdexB = B + col*numrows + row + 1;
		vdexC = S + E*(numdivs - 1) + row + 1;
		vdexD = S + E*(numdivs - 1) + row;

		// Add vertex indices of bottom corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;
	}

	// row (numrows - 1), all but last atom (col == numcols - 1)
	row = numrows - 1;
	for (col = 0; col < numcols - 1; col++) {

		vdexA = col*numrows + row;
		vdexB = S + (numrows - 1) + col;
		vdexC = S + (numrows - 1) + col + 1;
		vdexD = (col + 1)*numrows + row;

		// Add vertex indices of top corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;

		vdexA = B + col*numrows + row;
		vdexB = B + (col + 1)*numrows + row;
		vdexC = S + E*(numdivs - 1) + (numrows - 1) + col + 1;
		vdexD = S + E*(numdivs - 1) + (numrows - 1) + col;

		// Add vertex indices of bottom corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;
	}

	// col (numcols - 1), all but last atom (row == 0);
	col = numcols - 1;
	for (row = numrows - 1; row > 0; row--) {

		vdexA = col*numrows + row;
		vdexB = S + (numrows - 1) + (numcols - 1) + (numrows - 1 - row);
		vdexC = S + (numrows - 1) + (numcols - 1) + (numrows - 1 - row) + 1;
		vdexD = col*numrows + row - 1;

		// Add vertex indices of top corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;

		vdexA = B + col*numrows + row;
		vdexB = B + col*numrows + row - 1;
		vdexC = S + E*(numdivs - 1) + (numrows - 1) + (numcols - 1) + (numrows - 1 - row) + 1;
		vdexD = S + E*(numdivs - 1) + (numrows - 1) + (numcols - 1) + (numrows - 1 - row);

		// Add vertex indices of bottom corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;
	}

	// row 0, all but last atom (col == 0);
	row = 0;
	for (col = numcols - 1; col > 0; col--) {

		vdexA = col*numrows + row;
		vdexB = S + 2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col);
		vdexC = S + INTMOD(2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col) + 1, numsidecols);
		vdexD = (col - 1)*numrows + row;

		// Add vertex indices of top corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;


		vdexA = B + col*numrows + row;
		vdexB = B + (col - 1)*numrows + row;
		vdexC = S + E*(numdivs - 1) + INTMOD(2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col) + 1, numsidecols);
		vdexD = S + E*(numdivs - 1) + 2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col);

		// Add vertex indices of bottom corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;
	}

	// Now, if numdivs > 1, compute side polygons

	if (numdivs > 1) {	
		for (col = 0; col < numsidecols; col++) {
			for (row = 1; row < numdivs; row++) {

				vdexA = S + E*(row - 1) + col;
				vdexB = S + E*row + col;
				vdexC = S + E*row + INTMOD(col + 1, numsidecols);
				vdexD = S + E*(row - 1) + INTMOD(col + 1, numsidecols);

				// Add vertex indices of bottom corresponding faces
				pi_list[pdex].Create(4);	// create PindexList of length 4
				pi_list[pdex][0] = vdexA;
				pi_list[pdex][1] = vdexB;
				pi_list[pdex][2] = vdexC;
				pi_list[pdex++][3] = vdexD;
			}
		}
	}

	if (npolys != pdex)
		fprintf(stderr, "ERROR:  QuadMesh::constructMesh0() -- incorrect nvert count\n");
}

/*
(Diatomgrid *thisgrid)
{
	// Assign figuregrid
	if(thisgrid!=NULL && fig_grid==NULL) 
	{
		fig_grid = new Diatomgrid;
		fig_grid->rows = thisgrid->rows;
		fig_grid->cols = thisgrid->cols;
		// copy entries into internal fig_grid
		fig_grid->dlist = new Diatom[thisgrid->rows*thisgrid->cols];
		for(int cnt = 0; cnt < thisgrid->rows*thisgrid->cols; cnt++)
			fig_grid->dlist[cnt] = thisgrid->dlist[cnt];
	}
	else 
	{
		if(fig_grid==NULL)
		{
			fig_grid = new Diatomgrid;
			fig_grid->rows = 0;
			fig_grid->cols = 0;
			fig_grid->dlist = NULL;
		}
	}

	if(v_list!=NULL) delete v_list;
	v_list = NULL;
	if(pi_list!=NULL) delete pi_list;
	pi_list = NULL;
	numrows = fig_grid->rows;
	numcols = fig_grid->cols;
	if (fig_grid->cols < 2 || fig_grid->rows < 2)
		fprintf(stderr, "WARNING: QuadMesh::QuadMesh:  Not a legitimate QuadMesh.\n");

    Diatom *atom = fig_grid->dlist;
    Diatom workatom, nextatom, lastatom;

    DbVector3 p, v1, v2, b, bperp, v1pnt, v2pnt, bpnt, v1norm, v2norm;
	double uval, vval, tval, etaval, newetaval;
    double radius, t_decrement;
	int vertindex;

	numdivs = n;
	numrows = fig_grid->rows;
	numcols = fig_grid->cols;

	numsiderows = n;
	numsidecols = 2*numrows + 2*numcols - 4;

	int B = numrows*numcols;		// first bottom atom index
	int S = 2*numrows*numcols;		// first side index
	int E = numsidecols;			// rowlength along side

	nverts = 2*numrows*numcols + numsiderows*numsidecols;

	v_list = new Bpoint[nverts];

	// First, initialize top and bottom of mesh
	int col, row, atomindex;
	for (col = 0; col < numcols; col++) 
	{
		for (row = 0; row < numrows; row++) {

			atomindex = col*numrows + row;
			workatom = atom[atomindex];
			p = workatom.p_val();
			v1 = workatom.v1_val();
			v2 = workatom.v2_val();

			v1pnt = p + v1;
			v2pnt = p + v2;
			radius = workatom.r_val();
			v1.selfnormalize();
			v2.selfnormalize();

			uval = (double) col;
			vval = (double) row;

			v_list[atomindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1.x(), v1.y(), v1.z(), radius, 
					uval, vval, 1.0);
			v_list[atomindex + B].setvals(v2pnt.x(), v2pnt.y(), v2pnt.z(), v2.x(), v2.y(), v2.z(), radius, 
					uval, vval, -1.0);
		}
	}
	* /

	// Now, initialize side of mesh, one edge at a time
	// column zero, all but last atom (row == numrows - 1)
	int siderow;
	Quat q;
	double theta, rottheta;

	col = 0;
	for (row = 0; row < numrows - 1; row++) 
	{

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();
		radius = workatom.r_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		q.rot(bperp, rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + row;
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			q.rotate_vec(v1);
			v1.selfnormalize();

			v1pnt = p + v1*radius*newetaval;
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1norm.x(), v1norm.y(), v1norm.z(), radius, 
					uval, vval, tval);
		}
	}

	// row (numrows - 1), all but last atom (col == numcols - 1)
	row = numrows - 1;
	for (col = 0; col < numcols - 1; col++) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();
		radius = workatom.r_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		q.rot(bperp, rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + (numrows - 1) + col;
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			q.rotate_vec(v1);
			v1.selfnormalize();

			v1pnt = p + v1*radius*newetaval;
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1norm.x(), v1norm.y(), v1norm.z(), radius, 
					uval, vval, tval);
		}
	}

	// col (numcols - 1), all but last atom (row == 0);
	col = numcols - 1;
	for (row = numrows - 1; row > 0; row--) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();
		radius = workatom.r_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		q.rot(bperp, rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + (numrows - 1) + (numcols - 1) + (numrows - 1 - row);
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			q.rotate_vec(v1);
			v1.selfnormalize();

			v1pnt = p + v1*radius*newetaval;
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1norm.x(), v1norm.y(), v1norm.z(), radius, 
					uval, vval, tval);
		}
	}

	// row 0, all but last atom (col == 0);
	row = 0;
	for (col = numcols - 1; col > 0; col--) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();
		radius = workatom.r_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		q.rot(bperp, rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + 2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col);
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			q.rotate_vec(v1);
			v1.selfnormalize();

			v1pnt = p + v1*radius*newetaval;
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1norm.x(), v1norm.y(), v1norm.z(), radius, 
					uval, vval, tval);
		}
	}

	int numprows = numrows - 1;
	int numpcols = numcols - 1;
	int numsideprows = numsiderows + 1;
	int numsidepcols = numsidecols;		// because of wrap-around, same # of polygons
									//   as vertices around the side
	npolys = 2*numprows*numpcols + numsideprows*numsidepcols;

	int polyB = numprows*numpcols;
	int polyS = 2*numprows*numpcols;
	int polyE = numsidecols;

	int B = numrows*numcols;
	int S = 2*numrows*numcols;
	int E = numsidecols;

	pi_list = new PindexList[npolys];

	CCPindex pdex;
	int vdexA, vdexB, vdexC, vdexD;

	int row, col;

	// Do top and bottom surface polygons first
	for (col = 0; col < numpcols; col++) {
		for (row = 0; row < numprows; row++) {

			vdexA = col*numrows + row;
			vdexB = col*numrows + row + 1;
			vdexC = (col + 1)*numrows + row + 1;
			vdexD = (col + 1)*numrows + row;

			pdex = col*numprows + row;

			// Add vertex indices of top and bottom corresponding faces
			pi_list[pdex].Create(4);	// create PindexList of length 4
			pi_list[pdex][0] = vdexA;
			pi_list[pdex][1] = vdexB;
			pi_list[pdex][2] = vdexC;
			pi_list[pdex][3] = vdexD;

			pi_list[pdex + polyB].Create(4);
			pi_list[pdex + polyB][0] = vdexA + B;
			pi_list[pdex + polyB][1] = vdexD + B;
			pi_list[pdex + polyB][2] = vdexC + B;
			pi_list[pdex + polyB][3] = vdexB + B;
		}
	}

	// Now, initialize side of mesh, one edge at a time
	// First, top row of polygons

	// pdex will be used as a counter index, now
	pdex = polyS;
	// column zero, all but last atom (row == numrows - 1)
	col = 0;
	for (row = 0; row < numrows - 1; row++) {

		vdexA = col*numrows + row;
		vdexB = S + row;
		vdexC = S + row + 1;
		vdexD = col*numrows + row + 1;

		// Add vertex indices of top corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;

		vdexA = B + col*numrows + row;
		vdexB = B + col*numrows + row + 1;
		vdexC = S + E*(numdivs - 1) + row + 1;
		vdexD = S + E*(numdivs - 1) + row;

		// Add vertex indices of bottom corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;
	}

	// row (numrows - 1), all but last atom (col == numcols - 1)
	row = numrows - 1;
	for (col = 0; col < numcols - 1; col++) {

		vdexA = col*numrows + row;
		vdexB = S + (numrows - 1) + col;
		vdexC = S + (numrows - 1) + col + 1;
		vdexD = (col + 1)*numrows + row;

		// Add vertex indices of top corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;

		vdexA = B + col*numrows + row;
		vdexB = B + (col + 1)*numrows + row;
		vdexC = S + E*(numdivs - 1) + (numrows - 1) + col + 1;
		vdexD = S + E*(numdivs - 1) + (numrows - 1) + col;

		// Add vertex indices of bottom corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;
	}

	// col (numcols - 1), all but last atom (row == 0);
	col = numcols - 1;
	for (row = numrows - 1; row > 0; row--) {

		vdexA = col*numrows + row;
		vdexB = S + (numrows - 1) + (numcols - 1) + (numrows - 1 - row);
		vdexC = S + (numrows - 1) + (numcols - 1) + (numrows - 1 - row) + 1;
		vdexD = col*numrows + row - 1;

		// Add vertex indices of top corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;

		vdexA = B + col*numrows + row;
		vdexB = B + col*numrows + row - 1;
		vdexC = S + E*(numdivs - 1) + (numrows - 1) + (numcols - 1) + (numrows - 1 - row) + 1;
		vdexD = S + E*(numdivs - 1) + (numrows - 1) + (numcols - 1) + (numrows - 1 - row);

		// Add vertex indices of bottom corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;
	}

	// row 0, all but last atom (col == 0);
	row = 0;
	for (col = numcols - 1; col > 0; col--) {

		vdexA = col*numrows + row;
		vdexB = S + 2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col);
		vdexC = S + INTMOD(2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col) + 1, numsidecols);
		vdexD = (col - 1)*numrows + row;

		// Add vertex indices of top corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;


		vdexA = B + col*numrows + row;
		vdexB = B + (col - 1)*numrows + row;
		vdexC = S + E*(numdivs - 1) + INTMOD(2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col) + 1, numsidecols);
		vdexD = S + E*(numdivs - 1) + 2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col);

		// Add vertex indices of bottom corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;
	}

	// Now, if numdivs > 1, compute side polygons

	if (numdivs > 1) {	
		for (col = 0; col < numsidecols; col++) {
			for (row = 1; row < numdivs; row++) {

				vdexA = S + E*(row - 1) + col;
				vdexB = S + E*row + col;
				vdexC = S + E*row + INTMOD(col + 1, numsidecols);
				vdexD = S + E*(row - 1) + INTMOD(col + 1, numsidecols);

				// Add vertex indices of bottom corresponding faces
				pi_list[pdex].Create(4);	// create PindexList of length 4
				pi_list[pdex][0] = vdexA;
				pi_list[pdex][1] = vdexB;
				pi_list[pdex][2] = vdexC;
				pi_list[pdex++][3] = vdexD;
			}
		}
	}

	if (npolys != pdex)
		fprintf(stderr, "ERROR:  QuadMesh::constructMesh0() -- incorrect nvert count\n");

}


(Diatomgrid *thisgrid)
{
	// Assign figuregrid
	if (thisgrid != NULL) 
	{
		fig_grid = new Diatomgrid;
		fig_grid->rows = thisgrid->rows;
		fig_grid->cols = thisgrid->cols;
		// copy entries into internal fig_grid
		fig_grid->dlist = new Diatom[thisgrid->rows*thisgrid->cols];
		for (int cnt = 0; cnt < thisgrid->rows*thisgrid->cols; cnt++)
			fig_grid->dlist[cnt] = thisgrid->dlist[cnt];
	}
	else 
	{
		fig_grid = new Diatomgrid;
		fig_grid->rows = 0;
		fig_grid->cols = 0;
		fig_grid->dlist = NULL;
	}
	v_list = NULL;
	pi_list = NULL;
	numrows = fig_grid->rows;
	numcols = fig_grid->cols;
	if (fig_grid->cols < 2 || fig_grid->rows < 2)
		fprintf(stderr, "WARNING: QuadMesh::QuadMesh:  Not a legitimate QuadMesh.\n");

    Diatom *atom = fig_grid->dlist;
    Diatom workatom, nextatom, lastatom;

    DbVector3 p, v1, v2, b, bperp, v1pnt, v2pnt, bpnt, v1norm, v2norm;
	double uval, vval, tval, etaval, newetaval;
    double radius, t_decrement;
	int vertindex;

	numdivs = n;
	numrows = fig_grid->rows;
	numcols = fig_grid->cols;

	numsiderows = n;
	numsidecols = 2*numrows + 2*numcols - 4;

	int B = numrows*numcols;		// first bottom atom index
	int S = 2*numrows*numcols;		// first side index
	int E = numsidecols;			// rowlength along side

	nverts = 2*numrows*numcols + numsiderows*numsidecols;

	v_list = new Bpoint[nverts];

	// First, initialize top and bottom of mesh
	int col, row, atomindex;
	for (col = 0; col < numcols; col++) {
		for (row = 0; row < numrows; row++) {

			atomindex = col*numrows + row;
			workatom = atom[atomindex];
			p = workatom.p_val();
			v1 = workatom.v1_val();
			v2 = workatom.v2_val();

			v1pnt = p + v1;
			v2pnt = p + v2;
			radius = workatom.r_val();
			v1.selfnormalize();
			v2.selfnormalize();

			uval = (double) col;
			vval = (double) row;

			v_list[atomindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1.x(), v1.y(), v1.z(), radius, 
					uval, vval, 1.0);
			v_list[atomindex + B].setvals(v2pnt.x(), v2pnt.y(), v2pnt.z(), v2.x(), v2.y(), v2.z(), radius, 
					uval, vval, -1.0);
		}
	}

	// Now, initialize side of mesh, one edge at a time

	// column zero, all but last atom (row == numrows - 1)
	int siderow;
	Quat q;
	double theta, rottheta;

	col = 0;
	for (row = 0; row < numrows - 1; row++) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();
		radius = workatom.r_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		q.rot(bperp, rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + row;
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			q.rotate_vec(v1);
			v1.selfnormalize();

			v1pnt = p + v1*radius*newetaval;
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1norm.x(), v1norm.y(), v1norm.z(), radius, 
					uval, vval, tval);
		}
	}

	// row (numrows - 1), all but last atom (col == numcols - 1)
	row = numrows - 1;
	for (col = 0; col < numcols - 1; col++) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();
		radius = workatom.r_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		q.rot(bperp, rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + (numrows - 1) + col;
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			q.rotate_vec(v1);
			v1.selfnormalize();

			v1pnt = p + v1*radius*newetaval;
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1norm.x(), v1norm.y(), v1norm.z(), radius, 
					uval, vval, tval);
		}
	}

	// col (numcols - 1), all but last atom (row == 0);
	col = numcols - 1;
	for (row = numrows - 1; row > 0; row--) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();
		radius = workatom.r_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		q.rot(bperp, rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + (numrows - 1) + (numcols - 1) + (numrows - 1 - row);
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			q.rotate_vec(v1);
			v1.selfnormalize();

			v1pnt = p + v1*radius*newetaval;
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1norm.x(), v1norm.y(), v1norm.z(), radius, 
					uval, vval, tval);
		}
	}

	// row 0, all but last atom (col == 0);
	row = 0;
	for (col = numcols - 1; col > 0; col--) {

		atomindex = col*numrows + row;
		workatom = atom[atomindex];
		p = workatom.p_val();
		v1 = workatom.v1_val();
		etaval = workatom.eta_val();
		theta = 2.0 * workatom.theta_val();	// Theta here is rotation from v1 to v2
		bperp = workatom.bperp_val();
		radius = workatom.r_val();

		v1.selfnormalize();
		bperp.selfnormalize();

		rottheta = theta/(numdivs + 1);
		q.rot(bperp, rottheta);

		uval = (double) col;
		vval = (double) row;

		t_decrement = 2.0/(numdivs + 1);

		for (siderow = 1; siderow < numdivs + 1; siderow++) {

			vertindex = S + E*(siderow - 1) + 2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col);
			tval = 1.0 - t_decrement*siderow;

			// eta linearly interpolates between 1 (at fabs(t) = 1) and eta (at fabs(t) = 0)
			newetaval = (1.0 - fabs(tval))*etaval + fabs(tval);

			q.rotate_vec(v1);
			v1.selfnormalize();

			v1pnt = p + v1*radius*newetaval;
			v1norm = v1;		// This is an incorrect normal if eta != 1, except if numdivs = 1
								//    (which is the usual case)

			v_list[vertindex].setvals(v1pnt.x(), v1pnt.y(), v1pnt.z(), v1norm.x(), v1norm.y(), v1norm.z(), radius, 
					uval, vval, tval);
		}
	}

	int numprows = numrows - 1;
	int numpcols = numcols - 1;
	int numsideprows = numsiderows + 1;
	int numsidepcols = numsidecols;		// because of wrap-around, same # of polygons
										//   as vertices around the side
	npolys = 2*numprows*numpcols + numsideprows*numsidepcols;

	int polyB = numprows*numpcols;
	int polyS = 2*numprows*numpcols;
	int polyE = numsidecols;

	int B = numrows*numcols;
	int S = 2*numrows*numcols;
	int E = numsidecols;

	pi_list = new PindexList[npolys];

	CCPindex pdex;
	int vdexA, vdexB, vdexC, vdexD;

	int row, col;

	// Do top and bottom surface polygons first
	for (col = 0; col < numpcols; col++) {
		for (row = 0; row < numprows; row++) {

			vdexA = col*numrows + row;
			vdexB = col*numrows + row + 1;
			vdexC = (col + 1)*numrows + row + 1;
			vdexD = (col + 1)*numrows + row;

			pdex = col*numprows + row;

			// Add vertex indices of top and bottom corresponding faces
			pi_list[pdex].Create(4);	// create PindexList of length 4
			pi_list[pdex][0] = vdexA;
			pi_list[pdex][1] = vdexB;
			pi_list[pdex][2] = vdexC;
			pi_list[pdex][3] = vdexD;

			pi_list[pdex + polyB].Create(4);
			pi_list[pdex + polyB][0] = vdexA + B;
			pi_list[pdex + polyB][1] = vdexD + B;
			pi_list[pdex + polyB][2] = vdexC + B;
			pi_list[pdex + polyB][3] = vdexB + B;
		}
	}

	// Now, initialize side of mesh, one edge at a time
	// First, top row of polygons

	// pdex will be used as a counter index, now
	pdex = polyS;
	// column zero, all but last atom (row == numrows - 1)
	col = 0;
	for (row = 0; row < numrows - 1; row++) {

		vdexA = col*numrows + row;
		vdexB = S + row;
		vdexC = S + row + 1;
		vdexD = col*numrows + row + 1;

		// Add vertex indices of top corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;

		vdexA = B + col*numrows + row;
		vdexB = B + col*numrows + row + 1;
		vdexC = S + E*(numdivs - 1) + row + 1;
		vdexD = S + E*(numdivs - 1) + row;

		// Add vertex indices of bottom corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;
	}

	// row (numrows - 1), all but last atom (col == numcols - 1)
	row = numrows - 1;
	for (col = 0; col < numcols - 1; col++) {

		vdexA = col*numrows + row;
		vdexB = S + (numrows - 1) + col;
		vdexC = S + (numrows - 1) + col + 1;
		vdexD = (col + 1)*numrows + row;

		// Add vertex indices of top corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;

		vdexA = B + col*numrows + row;
		vdexB = B + (col + 1)*numrows + row;
		vdexC = S + E*(numdivs - 1) + (numrows - 1) + col + 1;
		vdexD = S + E*(numdivs - 1) + (numrows - 1) + col;

		// Add vertex indices of bottom corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;
	}

	// col (numcols - 1), all but last atom (row == 0);
	col = numcols - 1;
	for (row = numrows - 1; row > 0; row--) {

		vdexA = col*numrows + row;
		vdexB = S + (numrows - 1) + (numcols - 1) + (numrows - 1 - row);
		vdexC = S + (numrows - 1) + (numcols - 1) + (numrows - 1 - row) + 1;
		vdexD = col*numrows + row - 1;

		// Add vertex indices of top corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;

		vdexA = B + col*numrows + row;
		vdexB = B + col*numrows + row - 1;
		vdexC = S + E*(numdivs - 1) + (numrows - 1) + (numcols - 1) + (numrows - 1 - row) + 1;
		vdexD = S + E*(numdivs - 1) + (numrows - 1) + (numcols - 1) + (numrows - 1 - row);

		// Add vertex indices of bottom corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;
	}

	// row 0, all but last atom (col == 0);
	row = 0;
	for (col = numcols - 1; col > 0; col--) {

		vdexA = col*numrows + row;
		vdexB = S + 2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col);
		vdexC = S + INTMOD(2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col) + 1, numsidecols);
		vdexD = (col - 1)*numrows + row;

		// Add vertex indices of top corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;


		vdexA = B + col*numrows + row;
		vdexB = B + (col - 1)*numrows + row;
		vdexC = S + E*(numdivs - 1) + INTMOD(2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col) + 1, numsidecols);
		vdexD = S + E*(numdivs - 1) + 2*(numrows - 1) + (numcols - 1) + (numcols - 1 - col);

		// Add vertex indices of bottom corresponding faces
		pi_list[pdex].Create(4);	// create PindexList of length 4
		pi_list[pdex][0] = vdexA;
		pi_list[pdex][1] = vdexB;
		pi_list[pdex][2] = vdexC;
		pi_list[pdex++][3] = vdexD;
	}

	// Now, if numdivs > 1, compute side polygons

	if (numdivs > 1) {	
		for (col = 0; col < numsidecols; col++) {
			for (row = 1; row < numdivs; row++) {

				vdexA = S + E*(row - 1) + col;
				vdexB = S + E*row + col;
				vdexC = S + E*row + INTMOD(col + 1, numsidecols);
				vdexD = S + E*(row - 1) + INTMOD(col + 1, numsidecols);

				// Add vertex indices of bottom corresponding faces
				pi_list[pdex].Create(4);	// create PindexList of length 4
				pi_list[pdex][0] = vdexA;
				pi_list[pdex][1] = vdexB;
				pi_list[pdex][2] = vdexC;
				pi_list[pdex++][3] = vdexD;
			}
		}
	}

	if (npolys != pdex)
		fprintf(stderr, "ERROR:  QuadMesh::constructMesh0() -- incorrect nvert count\n");

}
*/
