/********************************************************************************/
/*                                                                              */
/*      File    :  Xferlist.cpp													*/
/*                                                                              */
/*      Description:  function to load from M3DQuadfigure into Xferlist to		*/
/*			transfer data to Diatomgrids in ThallCode							*/
/*																				*/
/*      Project :  Seurat                                                       */
/*                                                                              */
/*      Author  :  A. Thall                                                     */
/*                                                                              */
/*      Date    :  28. May 2000													*/
/*                                                                              */
/*      Modifications:                                                          */
/********************************************************************************/
#define D_XFERLIST
#include "Shapedepend.h"
#include <typeinfo>
#include "M3DQuadFigure.h"
#include "M3DTubeFigure.h"

Xferlist *convertM3DtoXfer(const M3DFigure *figure)
{
	Xferlist* xferlist	= NULL;
	if( typeid(*figure) == typeid(M3DQuadFigure) ) {
		xferlist	= convertM3DtoXfer(dynamic_cast<const M3DQuadFigure*>(figure));
	}
	else if( typeid(*figure) == typeid(M3DTubeFigure) ) {
		xferlist	= convertM3DtoXfer(dynamic_cast<const M3DTubeFigure*>(figure));
	}
	else {
		// ERROR
	}
	return xferlist;
}

//dibyendu - modified this function to enable s-rep display

Xferlist *convertM3DtoXfer(const M3DQuadFigure *thisfigure)
{
	// Read values from M3DQuadfigure
	Xferlist *newlist = new Xferlist;
	newlist->figure  = thisfigure;
	newlist->numrows = thisfigure->getRowCount();
	newlist->numcols = thisfigure->getColumnCount();
	newlist->tolerance = thisfigure->getTolerance();
	newlist->atomlist = new XferAtom[newlist->numrows * newlist->numcols];

	XferAtom *thisatom = newlist->atomlist;
	M3DPrimitive * currPrim;

	Quat q, rotateQ;

	rotateQ.setAxisAngle(Vector3D(1.0, 0.0, 0.0), R_HALF_PI);

	for (int row = 0; row < newlist->numrows; row++) {
		for (int col = 0; col < newlist->numcols; col++) {

			currPrim = thisfigure->getPrimitivePtr(row, col);

			q = currPrim->getQ();	// * rotateQ;

			thisatom->X_p[0] = currPrim->getX().getX();
			thisatom->X_p[1] = currPrim->getX().getY();
			thisatom->X_p[2] = currPrim->getX().getZ();
			thisatom->X_q[0] = q.getX();
			thisatom->X_q[1] = q.getY();
			thisatom->X_q[2] = q.getZ();
			thisatom->X_q[3] = q.getW();
			thisatom->X_r = currPrim->getR();

			if(currPrim->type() == M3D_END_PRIMITIVE)
				thisatom->X_eta = (dynamic_cast<M3DEndPrimitive*>(currPrim))->getElongation();
			else
				thisatom->X_eta = 1.0;

			thisatom->X_theta = currPrim->getTheta();
			thisatom->X_rho = 0.25;

			// dibyendu - on 06/17/2011, added the following code to facilitate s-reps

			thisatom->X_r0 = (dynamic_cast<M3DQuadPrimitive*>(currPrim))->getR0() ;
			thisatom->X_r1 = (dynamic_cast<M3DQuadPrimitive*>(currPrim))->getR1() ;

			thisatom->X_u0[0] = currPrim->getU0().getX() ;
			thisatom->X_u0[1] = currPrim->getU0().getY() ;
			thisatom->X_u0[2] = currPrim->getU0().getZ() ;

			thisatom->X_u1[0] = currPrim->getU1().getX() ;
			thisatom->X_u1[1] = currPrim->getU1().getY() ;
			thisatom->X_u1[2] = currPrim->getU1().getZ() ;

			//// this is the mid-spoke information for the end atoms

			if( (dynamic_cast<M3DEndPrimitive*>(currPrim)) != NULL )
				thisatom->X_rEnd = (dynamic_cast<M3DEndPrimitive*>(currPrim))->getREnd() ;
			else
				thisatom->X_rEnd = currPrim->getR() ;

			thisatom->X_uEnd[0] = (dynamic_cast<M3DQuadPrimitive*>(currPrim))->getUEnd().getX() ;
			thisatom->X_uEnd[1] = (dynamic_cast<M3DQuadPrimitive*>(currPrim))->getUEnd().getY() ;
			thisatom->X_uEnd[2] = (dynamic_cast<M3DQuadPrimitive*>(currPrim))->getUEnd().getZ() ;
			
			// -----------------------------------------------------------------------

			thisatom++;
		}
	}

	return newlist;
}

Xferlist *convertM3DtoXfer(const M3DTubeFigure *_thisfigure)
{
	const M3DTubeFigure* thisfigure	= dynamic_cast<const M3DTubeFigure*>(_thisfigure);
//	M3DTubeFigure* thisfigure	= dynamic_cast<M3DTubeFigure*>(_thisfigure->clone());
//	thisfigure->subdivide();
//	thisfigure->subdivide();
//	thisfigure->subdivide();
	// Read values from M3DTubefigure
	Xferlist *newlist = new Xferlist;
	newlist->figure  = thisfigure;
	newlist->numcols = 1;	// TubeMesh thinks the reverse way :(
	newlist->numrows = thisfigure->getColumnCount();
	newlist->tolerance = thisfigure->getTolerance();
	newlist->atomlist = new XferAtom[newlist->numrows * newlist->numcols];

	XferAtom *thisatom = newlist->atomlist;
	M3DPrimitive * currPrim;

	Quat q, rotateQ;

	rotateQ.setAxisAngle(Vector3D(1.0, 0.0, 0.0), R_HALF_PI);

	for (int row = 0; row < newlist->numrows; row++) {
		const int col	= 0;
		currPrim = thisfigure->getPrimitivePtr(row);
		//
		// The first primitive for a tube atom has theta > 90 degrees and tangent
		// pointing in the positive direction. Fix this before converting to Xferlist.
		//
		if( row == 0 ) {
			//
			// Check the condition mentioned above
			// sometimes we would want to create tubes with no end atoms or those that have theta
			// <= 90 degrees for the first atom (note that theta should not be 90 degrees, it
			// causes some math error due to precision problems later.)
			//
			const M3DPrimitive* nextPrim = thisfigure->getPrimitivePtr(row+1);
			if( currPrim->getB() * (nextPrim->getX() - currPrim->getX()) > 0.0 ) {
			//if( currPrim->getTheta() >= R_PI/2.0 ) {
				q.setAxisAngle( thisfigure->getPrimitivePtr(row)->getN(), R_PI );
				q	= q * currPrim->getQ() * rotateQ;
				thisatom->X_theta = R_PI - currPrim->getTheta();
			}
		}
		else {
			q = currPrim->getQ() * rotateQ;
			thisatom->X_theta = currPrim->getTheta();
		}

		thisatom->X_p[0] = currPrim->getX().getX();
		thisatom->X_p[1] = currPrim->getX().getY();
		thisatom->X_p[2] = currPrim->getX().getZ();
		thisatom->X_q[0] = q.getX();
		thisatom->X_q[1] = q.getY();
		thisatom->X_q[2] = q.getZ();
		thisatom->X_q[3] = q.getW();
		thisatom->X_r = currPrim->getR();

		if(currPrim->type() == M3D_END_PRIMITIVE)
			thisatom->X_eta = (dynamic_cast<M3DEndPrimitive*>(currPrim))->getElongation();
		else
			thisatom->X_eta = 0.00001;	// some non-zero small value.

		thisatom->X_rho = 0.25;
		thisatom++;
	}

// FIXME: memory leak, identify when it needs to be deleted.
//	delete thisfigure;
	return newlist;
}
