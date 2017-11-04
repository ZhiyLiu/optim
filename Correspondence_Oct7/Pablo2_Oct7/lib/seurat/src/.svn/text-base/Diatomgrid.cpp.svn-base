/****************************************************************************/
/*																			*/
/*  	File	:  Diatomgrid.cpp											*/
/*																			*/
/*	Description:  a few utility functions for operating on Diatom			*/
/*		Diatomgrids, which are an open struct that's been class				*/
/*		encapsulated for organizational reasons.							*/
/*																			*/
/*	Project :  Seurat														*/
/*	Author  :  A. Thall														*/
/*	Date	:  4. June 2000													*/
/*																			*/
/*	Modifications:  														*/
/*		4. June -- added true_badmesh() and true_badatom() functions		*/
/*		4. June -- made readXferlist() and readXferlist as constructor().	*/
/*		24. June 02 -- removed all .plist Parentset code					*/
/****************************************************************************/

#define D_DIATOMGRID
#define D_XFERLIST
#include "Shapedepend.h"

using namespace ThallCode;

static const double EPSILON = 0.001;

bool nearzero(double blength) {
	if ((blength > -EPSILON) && (blength < EPSILON))
		return true;
	else
		return false;
}

/*
	Andrew changes: 
	1. deleted a plist of Parantset type
	2. deleted last line of
		'
		// add initialization of parent-grid ID
		plist[ndex].addmember(ndex);
		'
*/

/********************************************************************************/
/* constructor() -- read in Xferlist to initialize Diatomgrid					*/
/********************************************************************************/
Diatomgrid::Diatomgrid(Xferlist* mfig1)
{
    std::cout << "constructor" << std::endl;
	// Load values from (row-major-ordered) Xferlist mfig1 to a Diatomgrid
	//   (column-major-ordered)
	rows = mfig1->numrows;
	cols = mfig1->numcols;

	figure	= mfig1->figure;

	dlist = new Diatom[rows*cols];

	int Xfer_index, ndex;
	for (int row = 0; row < rows; row++) {
		for (int col = 0; col < cols; col++) {
			Xfer_index = row*cols + col;
			ndex = col*rows + row;

			dlist[ndex].set_p(
				DbVector3(mfig1->atomlist[Xfer_index].X_p[0],
						  mfig1->atomlist[Xfer_index].X_p[1],
						  mfig1->atomlist[Xfer_index].X_p[2]));
			dlist[ndex].set_q(
				Quat(mfig1->atomlist[Xfer_index].X_q[0],
					 mfig1->atomlist[Xfer_index].X_q[1],
					 mfig1->atomlist[Xfer_index].X_q[2],
					 mfig1->atomlist[Xfer_index].X_q[3]));
			dlist[ndex].set_r(
				mfig1->atomlist[Xfer_index].X_r);
			dlist[ndex].set_eta(
				mfig1->atomlist[Xfer_index].X_eta);
			dlist[ndex].set_phi(
				M_PI/2.0 - mfig1->atomlist[Xfer_index].X_theta);
			dlist[ndex].set_rho(
				mfig1->atomlist[Xfer_index].X_rho);
			dlist[ndex].set_r1(
				mfig1->atomlist[Xfer_index].X_r0);
			dlist[ndex].set_r2(
				mfig1->atomlist[Xfer_index].X_r1);
                        dlist[ndex].set_rEnd(
                                mfig1->atomlist[Xfer_index].X_rEnd);
                        
                        dlist[ndex].set_u1(
                                DbVector3(mfig1->atomlist[Xfer_index].X_u0[0],
						  mfig1->atomlist[Xfer_index].X_u0[1],
						  mfig1->atomlist[Xfer_index].X_u0[2]));
                        dlist[ndex].set_u2(
                                DbVector3(mfig1->atomlist[Xfer_index].X_u1[0],
						  mfig1->atomlist[Xfer_index].X_u1[1],
						  mfig1->atomlist[Xfer_index].X_u1[2]));
                        dlist[ndex].set_uEnd(
                                DbVector3(mfig1->atomlist[Xfer_index].X_uEnd[0],
						  mfig1->atomlist[Xfer_index].X_uEnd[1],
						  mfig1->atomlist[Xfer_index].X_uEnd[2]));
                        

			//std::cout << "Diatomgrid" << std::endl;

			// Now, set correct DiatomTYPE (INTERNAL_M, CORNER_M, EDGE_M)

			if (row == 0) {
				if (col == 0)
					dlist[ndex].Diatomtype(CORNER_M);
				else if (col == cols - 1)
					dlist[ndex].Diatomtype(CORNER_M);
				else 
					dlist[ndex].Diatomtype(EDGE_M);
			}
			else if (row == rows - 1) {
				if (col == 0)
					dlist[ndex].Diatomtype(CORNER_M);
				else if (col == cols - 1)
					dlist[ndex].Diatomtype(CORNER_M);
				else 
					dlist[ndex].Diatomtype(EDGE_M);
			}
			else if (col == 0) 
				dlist[ndex].Diatomtype(EDGE_M);
			else if (col == cols - 1)
				dlist[ndex].Diatomtype(EDGE_M);
			else 
				dlist[ndex].Diatomtype(INTERNAL_M);
		}
	}
}

// qiong added 111702
void Diatomgrid::CopyXferlist(Xferlist *mfig1)
{
    //std::cout << "copy" << std::endl;
	int Xfer_index, ndex;
	figure	= mfig1->figure;
	for (int row = 0; row < rows; row++) {
		for (int col = 0; col < cols; col++) {
			Xfer_index = row*cols + col;
			ndex = col*rows + row;

			dlist[ndex].set_p(
				DbVector3(mfig1->atomlist[Xfer_index].X_p[0],
						  mfig1->atomlist[Xfer_index].X_p[1],
						  mfig1->atomlist[Xfer_index].X_p[2]));
			dlist[ndex].set_q(
				Quat(mfig1->atomlist[Xfer_index].X_q[0],
					 mfig1->atomlist[Xfer_index].X_q[1],
					 mfig1->atomlist[Xfer_index].X_q[2],
					 mfig1->atomlist[Xfer_index].X_q[3]));
			dlist[ndex].set_r(
				mfig1->atomlist[Xfer_index].X_r);
			dlist[ndex].set_eta(
				mfig1->atomlist[Xfer_index].X_eta);
			dlist[ndex].set_phi(
				M_PI/2.0 - mfig1->atomlist[Xfer_index].X_theta);
			dlist[ndex].set_rho(
				mfig1->atomlist[Xfer_index].X_rho);
			dlist[ndex].set_r1(
				mfig1->atomlist[Xfer_index].X_r0);
			dlist[ndex].set_r2(
				mfig1->atomlist[Xfer_index].X_r1);
                        dlist[ndex].set_rEnd(
                                mfig1->atomlist[Xfer_index].X_rEnd);                        
                        
                        dlist[ndex].set_u1(
                                DbVector3(mfig1->atomlist[Xfer_index].X_u0[0],
						  mfig1->atomlist[Xfer_index].X_u0[1],
						  mfig1->atomlist[Xfer_index].X_u0[2]));
                        dlist[ndex].set_u2(
                                DbVector3(mfig1->atomlist[Xfer_index].X_u1[0],
						  mfig1->atomlist[Xfer_index].X_u1[1],
						  mfig1->atomlist[Xfer_index].X_u1[2]));
                        dlist[ndex].set_uEnd(
                                DbVector3(mfig1->atomlist[Xfer_index].X_uEnd[0],
						  mfig1->atomlist[Xfer_index].X_uEnd[1],
						  mfig1->atomlist[Xfer_index].X_uEnd[2]));

			//std::cout << "Diatomrid::CopyXferlist" << std::endl;
			//std::cout << "r: " << dlist[ndex].r_val() << std::endl;
			//std::cout << "r0: " << dlist[ndex].r1_val() << std::endl;
			//std::cout << "r1: " << dlist[ndex].r2_val() << std::endl;

/*			Quat rot90;
			DbVector3 one(1.0, 0.0, 0.0);
			rot90.rot(one, -M_PI/2.0);

			dlist[ndex].set_q(dlist[ndex].q_val() * rot90);
 */
			// Now, set correct DiatomTYPE (INTERNAL_M, CORNER_M, EDGE_M)

			if (row == 0) {
				if (col == 0)
					dlist[ndex].Diatomtype(CORNER_M);
				else if (col == cols - 1)
					dlist[ndex].Diatomtype(CORNER_M);
				else 
					dlist[ndex].Diatomtype(EDGE_M);
			}
			else if (row == rows - 1) {
				if (col == 0)
					dlist[ndex].Diatomtype(CORNER_M);
				else if (col == cols - 1)
					dlist[ndex].Diatomtype(CORNER_M);
				else 
					dlist[ndex].Diatomtype(EDGE_M);
			}
			else if (col == 0) 
				dlist[ndex].Diatomtype(EDGE_M);
			else if (col == cols - 1)
				dlist[ndex].Diatomtype(EDGE_M);
			else 
				dlist[ndex].Diatomtype(INTERNAL_M);
		}
	}
/*	DbVector3 av_position = DbVector3(0.0, 0.0, 0.0);
	for (int dex1 = 0; dex1 < rows*cols; dex1++)
		av_position += dlist[dex1].p_val();
	av_position /= rows*cols;

	for (int dex2 = 0; dex2 < rows*cols; dex2++)
		dlist[dex2].set_p(dlist[dex2].p_val() - av_position);
*/		
}

/*
	Andrew changes:
		again, the plist is gone!
*/

/********************************************************************************/
/* readXlist() -- read in Xferlist to initialize Diatomgrid						*/
/********************************************************************************/
void Diatomgrid::readXferlist(Xferlist* mfig1)
{
    std::cout << "read" << std::endl;
	// Load values from (row-major-ordered) Xferlist mfig1 to a Diatomgrid
	//   (column-major-ordered)
	rows = mfig1->numrows;
	cols = mfig1->numcols;
	figure	= mfig1->figure;
	if (dlist != NULL)
		delete []dlist;
	dlist = new Diatom[rows*cols];
	CopyXferlist(mfig1);
}

/********************************************************************************/
/* update_mesh() -- update (row, col) element of Diatomgrid with the values in	*/
/*     the Xferatom																*/
/********************************************************************************/
void Diatomgrid::update_mesh(XferAtom *thisatom, int modrow, int modcol)
{

    std::cout << "update_mesh1" << std::endl;
	int ndex = modrow*cols + modcol;

	dlist[ndex].set_p(
		DbVector3(thisatom->X_p[0], thisatom->X_p[1], thisatom->X_p[2]));
	dlist[ndex].set_q(
		Quat(thisatom->X_q[0], thisatom->X_q[1], thisatom->X_q[2], thisatom->X_q[3]));
	dlist[ndex].set_r(thisatom->X_r);
	dlist[ndex].set_eta(thisatom->X_eta);
	dlist[ndex].set_phi(M_PI/2.0 - thisatom->X_theta);
	dlist[ndex].set_rho(thisatom->X_rho);
	dlist[ndex].set_r1(thisatom->X_r0);
	dlist[ndex].set_r2(thisatom->X_r1);

	//std::cout << "DiatomGrid::update_mesh" << std::endl;

	Quat rot90;
	DbVector3 one(1.0, 0.0, 0.0);

	//rot90.rot(one, -M_PI/2.0);
	rot90.setAxisAngle(Vector3D(one.X(), one.Y(), one.Z()), -M_PI/2.0);

	dlist[ndex].set_q(dlist[ndex].q_val() * rot90);
}

/********************************************************************************/
/* update_mesh() -- update (row, col) element of Diatomgrid with the values in	*/
/*     the Xferatom																*/
/********************************************************************************/
void Diatomgrid::update_mesh(Diatom *thisatom, int modrow, int modcol)
{
    std::cout << "update_mesh2" << std::endl;

	int ndex = modrow*cols + modcol;

	dlist[ndex] = *thisatom;
}

/********************************************************************************/
/* true_badmesh() -- This trues ALL boundary (EDGE_M and CORNER_M) Diatoms		*/
/*		in the mesh																*/
/********************************************************************************/
void Diatomgrid::true_badmesh()
{
	int row, col;

	row = 0;
	for (col = 0; col < cols; col++)
		true_badatom(row, col);
	row = rows - 1;
	for (col = 0; col < cols; col++)
		true_badatom(row, col);
	col = 0;
	for (row = 0; row < rows; row++)
		true_badatom(row, col);
	col = cols - 1;
	for (row = 0; row < rows; row++)
		true_badatom(row, col);
}
/********************************************************************************/
/* This rotates EDGE_M and CORNER_M Diatoms to set bvec (q frame) to point		*/
/*	 outward from the Grid, based on mesh vectors between the bad Diatom and	*
/*   its two neighbors on the edge.												*/
/********************************************************************************/
/********************************************************************************/
/* We'll project the vectors to the neighbors onto the b-bperp plane, and then	*/
/*	 rotate b until it's midway between the projected vectors.					*/
/********************************************************************************/
void Diatomgrid::true_badatom(int row, int col)
{
	DbVector3 neighbor_r, neighbor_l;	// left and right neighbors as seen from top

	Diatom *thisatom = &dlist[idx(row, col)];

	DbVector3 thispoint = thisatom->p_val();

	if (row == 0 && col == 0) {
		neighbor_r = dlist[idx(row + 1, col    )].p_val() - thispoint;
		neighbor_l = dlist[idx(row    , col + 1)].p_val() - thispoint;
	}
	else if (row == 0 && col == cols - 1) {
		neighbor_r = dlist[idx(row    , col - 1)].p_val() - thispoint;
		neighbor_l = dlist[idx(row + 1, col    )].p_val() - thispoint;
	}
	else if (row == rows - 1 && col == cols - 1) {
		neighbor_r = dlist[idx(row - 1, col    )].p_val() - thispoint;
		neighbor_l = dlist[idx(row    , col - 1)].p_val() - thispoint;
	}
	else if (row == rows - 1 && col == 0) {
		neighbor_r = dlist[idx(row    , col + 1)].p_val() - thispoint;
		neighbor_l = dlist[idx(row - 1, col    )].p_val() - thispoint;
	}
	else if (row == 0) {
		neighbor_r = dlist[idx(row    , col - 1)].p_val() - thispoint;
		neighbor_l = dlist[idx(row    , col + 1)].p_val() - thispoint;
	}
	else if (row == rows - 1) {
		neighbor_r = dlist[idx(row    , col + 1)].p_val() - thispoint;
		neighbor_l = dlist[idx(row    , col - 1)].p_val() - thispoint;
	}
	else if (col == 0) {
		neighbor_r = dlist[idx(row + 1, col    )].p_val() - thispoint;
		neighbor_l = dlist[idx(row - 1, col    )].p_val() - thispoint;
	}
	else if (col == cols - 1) {
		neighbor_r = dlist[idx(row - 1, col    )].p_val() - thispoint;
		neighbor_l = dlist[idx(row + 1, col    )].p_val() - thispoint;
	}

	neighbor_r.selfnormalize();
	neighbor_l.selfnormalize();

	DbVector3 b, bperp, n, new_b, new_bperp, nr_proj, nl_proj;

	b = thisatom->b_val();
	bperp = thisatom->bperp_val();
	n = thisatom->n_val();

	nr_proj = b*(neighbor_r.dot(b)) + bperp*(neighbor_r.dot(bperp));

    nl_proj = b*(neighbor_l.dot(b)) + bperp*(neighbor_l.dot(bperp));

	nr_proj.selfnormalize();
	nl_proj.selfnormalize();

	b = nr_proj + nl_proj;

	if (nearzero(b.length()))
		b = nr_proj.cross(n);
	else {
		if (n.dot(nr_proj.cross(b)) > 0)
			b = -b;
	}
	b.selfnormalize();

	bperp = n.cross(b);
	bperp.selfnormalize();

	Quat new_q;
	//new_q.q_from_framevecs(b, bperp, n);
	new_q.buildFromFrame(Vector3D(b.X(), b.Y(), b.Z()), Vector3D(n.X(), n.Y(), n.Z()));

	thisatom->set_q(new_q);
}

