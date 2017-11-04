/*-----------------------------------------------------------
	To predict the rest of a figure by a group of link atoms
	by Qiong Han, 06/04/2004
-----------------------------------------------------------*/

#include <math.h>
#include <stdio.h>

#define D_XFERLIST
#include "Shapedepend.h"

#include "matrix.h"
#include "M3DFigure.h"
#include "M3DFigurePredictor.h"

// if enabled, the predictor will apply a bunch of similarity
// transformations to the data and try to recover the exact
// transformations, in order to test predictor
// #define DEBUG_PREDICTOR


M3DFigurePredictor::M3DFigurePredictor
(int linkCount, int *linkIds, M3DFigure *figure, M3DFigure *newFigure)
{
	M3DPrimitive **newLinkAtoms=new M3DPrimitive*[linkCount];
	for(int i=0; i<linkCount; i++)
	{
		if(newFigure->getPrimitivePtr(linkIds[i])->type()==M3D_STANDARD_PRIMITIVE)
			newLinkAtoms[i]=new M3DQuadPrimitive( *(dynamic_cast<M3DQuadPrimitive *>(newFigure->getPrimitivePtr(linkIds[i]))) );
		else
			newLinkAtoms[i]=new M3DQuadEndPrimitive( *(dynamic_cast<M3DQuadEndPrimitive *>(newFigure->getPrimitivePtr(linkIds[i]))) );

	}
	predictFigureByLinkAtoms(linkCount, linkIds, figure, newLinkAtoms);
	delete []newLinkAtoms;
}

M3DFigurePredictor::M3DFigurePredictor
(int linkCount, int *linkIds, M3DFigure *figure, M3DPrimitive **newPrimitives)
{
	predictFigureByLinkAtoms(linkCount, linkIds, figure, newPrimitives);
}

// to estimate the best (in the term of least squared distance) 
// similarity transformation by two given set of points []x and []y
void M3DFigurePredictor::estimateSimilarityTransformation
(int n, Vector3D *x, Vector3D *y, SimilarityTransform *sTrans, Vector3D *COR) {
	int i;
	double weight, weights;
	Vector3D xAvg, yAvg;
	double xDeltaSquare, yDeltaSquare;

	xAvg.set(0, 0, 0);
	yAvg.set(0, 0, 0);
	weights=0;
	for(i=0; i<n; i++) {
		weight=1;
		weights+=weight;
		xAvg+=weight*x[i];
		yAvg+=weight*y[i];
	}

	sTrans->scale=0;
	sTrans->cov.setSize(3, 3);
	xAvg/=weights;
	yAvg/=weights;

	// the center of roration is defaulted to be the center of the gravity
	// of the points x[], unless it is specified by the user
	if(COR==NULL)
		sTrans->COR=xAvg;
	else
		sTrans->COR=*COR;

	sTrans->transTemp.set(0, 0, 0);

	xDeltaSquare=0;
	yDeltaSquare=0;
	weights=0;
	for(i=0; i<n; i++) {
		weight=1;
		weights++;
		xDeltaSquare+=xAvg.distSquare(x[i]);
		yDeltaSquare+=yAvg.distSquare(y[i]);
	}
	xDeltaSquare/=weights;
	yDeltaSquare/=weights;

	// calculate the covariance matrix by the two given set of points
	weights=0;
	for(i=0; i<n; i++) {
		Vector3D xDiff=(x[i]-xAvg);
		Vector3D yDiff=(y[i]-yAvg);
		Matrix m1(3, 1, yDiff.getX(), yDiff.getY(), yDiff.getZ());
		Matrix m2(1, 3, xDiff.getX(), xDiff.getY(), xDiff.getZ());
		weight=1;
		weights+=weight;
		Matrix mT=weight*m1*m2;
		sTrans->cov+=mT;
	}
	sTrans->cov/=weights;

	// get the rotation and scale by a SVD of the covariance matrix
	Vector sv;
	Matrix U, Vt;
	sTrans->cov.factorSVD(U, Vt, sv);
	sTrans->rotM=U*Vt;

#ifdef DEBUG_PREDICTOR
	double uDet=U.det();
	double vDet=Vt.det();
	if(uDet*vDet<0)
	{
		printf("determinant of U is %f, determinant of V is %f.\n\n", uDet, vDet);
		printf("after SVD:\n\tU is\n");
		U.print();
		printf("\tV' is\n");
		Vt.print();
		printf("\tD is\n");
		sv.print();
		printf("\trotM is' is\n");
		sTrans->rotM.print();
	}
#endif

	double traceOfD;
	traceOfD=0;
	for(i=0; i<sv.size(); i++)
		traceOfD+=sv(i);

	// then estimate the translation according to the recovered
	// rotation, scale, and the given COR
	sTrans->scale=traceOfD/xDeltaSquare;
	xAvg=xAvg-sTrans->COR;
	Matrix temp(3, 1, xAvg.getX(), xAvg.getY(), xAvg.getZ());
	Vector temp2(3);
	temp2 = (sTrans->rotM*temp);
	sTrans->trans = yAvg-(sTrans->COR+sTrans->scale*Vector3D(temp2(0), temp2(1), temp2(2)));

	// convert the rotation matrix into a quaternion, to be used to rotate the atom frame later
	Vector3D vecX, vecY, vecZ;
	vecX.set((sTrans->rotM.getColumn(0))(0), (sTrans->rotM.getColumn(0))(1), (sTrans->rotM.getColumn(0))(2));
	vecY.set((sTrans->rotM.getColumn(1))(0), (sTrans->rotM.getColumn(1))(1), (sTrans->rotM.getColumn(1))(2));
	vecZ.set((sTrans->rotM.getColumn(2))(0), (sTrans->rotM.getColumn(2))(1), (sTrans->rotM.getColumn(2))(2));
	sTrans->rotQ=Quat(vecX, vecY, vecZ);
	sTrans->rotQ.normalize();
}


// the main prediction function goes here
// linkCount/linkIds	the # and indices of the "link" atoms
// figure				the entire figure, with the atoms before transformation
// newLinkAtoms			all the "link" atoms after the transformation
// the predicted figure will be stored back to *figure
void M3DFigurePredictor::predictFigureByLinkAtoms
(int linkCount, int *linkIds, M3DFigure *figure, M3DPrimitive **newLinkAtoms) {
	if(linkCount<1 || linkCount>figure->getPrimitiveCount() || figure==NULL || newLinkAtoms==NULL)
		return;

	if (linkCount == figure->getPrimitiveCount())
	{
		int dex;
		M3DPrimitive *prim;
		// no need to do anything else but copying if all the atoms are predictors
		for (dex=0; dex<linkCount; dex++)
		{
			prim = figure->getPrimitivePtr(dex);
			*prim = *newLinkAtoms[dex];
			if (prim->type() == M3D_END_PRIMITIVE)
				(dynamic_cast<M3DQuadEndPrimitive *>(prim))->setElongation((dynamic_cast<M3DQuadEndPrimitive *>(newLinkAtoms[dex]))->getElongation());
		}
		return;
	}

	int primitiveId;
	M3DPrimitive **primitive=NULL;

	int i, j;
	int atomNum=figure->getPrimitiveCount();

	primitive=new M3DPrimitive*[linkCount];
	for(j=0; j<linkCount; j++) {
		primitive[j]=figure->getPrimitivePtr(linkIds[j]);
	}

	// use the two spokeends, and the atom position as the given set of points
	Vector3D *x, *y;
	x=new Vector3D[3*linkCount];
	y=new Vector3D[3*linkCount];
	for(j=0; j<linkCount; j++) {
		x[3*j]=primitive[j]->getX()+primitive[j]->getY0();
		x[3*j+1]=primitive[j]->getX()+primitive[j]->getY1();
		x[3*j+2]=primitive[j]->getX();
		y[3*j]=newLinkAtoms[j]->getX()+newLinkAtoms[j]->getY0();
		y[3*j+1]=newLinkAtoms[j]->getX()+newLinkAtoms[j]->getY1();
		y[3*j+2]=newLinkAtoms[j]->getX();
	}

	// the center of gravity of the ENTIRE figure is used as the 
	// center of rotation
	// the default COR is the center of the set of points only
	Vector3D COR=figure->getCOG(true);
	estimateSimilarityTransformation(3*linkCount, x, y, &bestSim, &COR);

#ifdef DEBUG_PREDICTOR
	printf("the estimated rotation, scale and translation are:\n");
	bestSim.rotM.print();
	printf("scale is: %2.12lf\n", bestSim.scale);
	bestSim.trans.print();
	printf("\n");

	// to recover the estimated transformation after applying it to []x
	SimilarityTransform rTrans;
	Matrix rot(3, 3);
	Vector3D TRANS;
	double SCALE;
	// use the same COR
	// apply the just found rotation, scale and translation to []x
	rot=bestSim.rotM;
	SCALE=bestSim.scale;
	TRANS=bestSim.trans;
	for(j=0; j<linkCount; j++) {
		Vector3D diff;
		Vector temp1(3), temp2(3);

		diff=x[3*j]-COR;
		temp1(0)=diff.getX();
		temp1(1)=diff.getY();
		temp1(2)=diff.getZ();
		temp2=rot*temp1;
		diff.set(temp2(0), temp2(1), temp2(2));
		diff*=SCALE;
		diff+=COR;
		diff+=TRANS;
		y[3*j]=diff;

		diff=x[3*j+1]-COR;
		temp1(0)=diff.getX();
		temp1(1)=diff.getY();
		temp1(2)=diff.getZ();
		temp2=rot*temp1;
		diff.set(temp2(0), temp2(1), temp2(2));
		diff*=SCALE;
		diff+=COR;
		diff+=TRANS;
		y[3*j+1]=diff;

		diff=x[3*j+2]-COR;
		temp1(0)=diff.getX();
		temp1(1)=diff.getY();
		temp1(2)=diff.getZ();
		temp2=rot*temp1;
		diff.set(temp2(0), temp2(1), temp2(2));
		diff*=SCALE;
		diff+=COR;
		diff+=TRANS;
		y[3*j+2]=diff;
	}
	estimateSimilarityTransformation(3*linkCount, x, y, &rTrans, &COR);

	printf("the recovered rotation, scale and translation are:\n");
	rTrans.rotM.print();
	printf("scale is: %2.12lf\n", rTrans.scale);
	rTrans.trans.print();
	printf("\n");
#endif

	//printf("bestSim = %f\n", bestSim.scale);
	//bestSim.rotQ.print();
	//bestSim.trans.print();


	// apply the recovered transformations to the rest of the atoms
	M3DPrimitive *prim;
	for (i=0; i<atomNum; i++)
	{
		j=0;
		while(j<linkCount && i!=linkIds[j])
			j++;
		if(j>=linkCount) {
			primitiveId=i;
			prim=figure->getPrimitivePtr(primitiveId);

			Vector3D diff=prim->getX()-bestSim.COR;
			Vector temp(3, diff.getX(), diff.getY(), diff.getZ());
			temp=bestSim.scale*(bestSim.rotM*temp);
			Vector3D temp2(temp(0), temp(1), temp(2));
			prim->setX(temp2+bestSim.COR+bestSim.trans);
			prim->setR(prim->getR()*bestSim.scale);

			prim->rotateBy(bestSim.rotQ);

			if (prim->type() == M3D_END_PRIMITIVE)
				(dynamic_cast<M3DQuadEndPrimitive *>(prim))->setElongation((dynamic_cast<M3DQuadEndPrimitive *>(prim))->getElongation()*bestSim.scale);
		}
		else {
			primitiveId=linkIds[j];
            prim=figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;
			prim->setX(newLinkAtoms[j]->getX());
			prim->setQ(newLinkAtoms[j]->getQ());
			prim->setR(newLinkAtoms[j]->getR());
			prim->setTheta(newLinkAtoms[j]->getTheta());
			if(prim->type()==M3D_END_PRIMITIVE)
				(dynamic_cast<M3DQuadEndPrimitive *>(prim))->setElongation((dynamic_cast<M3DQuadEndPrimitive *>(newLinkAtoms[j]))->getElongation());
		}
	}

	delete []primitive;
	delete []x;
	delete []y;

	return;
}


// the following codes are to predict the end of the bones that were
// arbitrarily cut in the image, by the rest of the atoms in the 
// mean model
// this works exactly the same as the previous predictor, however only
// a subset of the "link atoms" are used as the predictors, which are
// the close neighboring atoms to the atoms at the cut-end
M3DFigurePredictor::M3DFigurePredictor
(int linkCount, int *linkIds, int realLinkCount, int *realLinkIds, 
 M3DFigure *figure, M3DFigure *newFigure) {
	int i;

	M3DPrimitive **newLinkAtoms=new M3DPrimitive*[linkCount];
	for(i=0; i<linkCount; i++)
	{
		if(newFigure->getPrimitivePtr(linkIds[i])->type()==M3D_STANDARD_PRIMITIVE)
			newLinkAtoms[i]=new M3DQuadPrimitive( *(dynamic_cast<M3DQuadPrimitive *>(newFigure->getPrimitivePtr(linkIds[i]))) );
		else
			newLinkAtoms[i]=new M3DQuadEndPrimitive( *(dynamic_cast<M3DQuadEndPrimitive *>(newFigure->getPrimitivePtr(linkIds[i]))) );

	}
	predictFigureByLinkAtoms(linkCount, linkIds, realLinkCount, realLinkIds, figure, newLinkAtoms);
	for(i=0; i<linkCount; i++)
	{
		delete newLinkAtoms[i];
	}
	delete []newLinkAtoms;
}

void M3DFigurePredictor::predictFigureByLinkAtoms
(int linkCount, int *linkIds, int realLinkCount, int *realLinkIds, 
 M3DFigure *figure, M3DPrimitive **newLinkAtoms) {

	if(linkCount<1 || linkCount>figure->getPrimitiveCount() || figure==NULL || newLinkAtoms==NULL)
		return;

	if(linkCount==figure->getPrimitiveCount())
	{
		int dex;
		for(dex=0; dex<linkCount; dex++)
			*(figure->getPrimitivePtr(dex))=*newLinkAtoms[dex];
	}

	int primitiveId;
	M3DPrimitive **primitive;

	int i, j;
	int atomNum=figure->getPrimitiveCount();

	// calculate the initial moment of initia of the
	//	starting link atoms

	primitive=new M3DPrimitive*[realLinkCount];
	for(j=0; j<realLinkCount; j++) {
		primitive[j]=figure->getPrimitivePtr(realLinkIds[j]);
	}

	Vector3D *x, *y;
	x=new Vector3D[3*realLinkCount];
	y=new Vector3D[3*realLinkCount];
	for(j=0; j<realLinkCount; j++) {
		x[3*j]=primitive[j]->getX()+primitive[j]->getY0();
		x[3*j+1]=primitive[j]->getX()+primitive[j]->getY1();
		x[3*j+2]=primitive[j]->getX();
		y[3*j]=newLinkAtoms[j]->getX()+newLinkAtoms[j]->getY0();
		y[3*j+1]=newLinkAtoms[j]->getX()+newLinkAtoms[j]->getY1();
		y[3*j+2]=newLinkAtoms[j]->getX();
	}

	// the center of gravity of the ENTIRE figure is used as the 
	// center of rotation
	Vector3D COR=figure->getCOG(true);
	estimateSimilarityTransformation(3*linkCount, x, y, &bestSim, &COR);

	M3DPrimitive *prim;
	for(i=0; i<atomNum; i++)
	{
		j=0;
		while(j<linkCount && i!=linkIds[j])
			j++;
		if(j>=linkCount) {
			primitiveId=i;
			prim=figure->getPrimitivePtr(primitiveId);

			Vector3D diff=prim->getX()-bestSim.COR;
			Vector temp(3, diff.getX(), diff.getY(), diff.getZ());
			temp=bestSim.scale*(bestSim.rotM*temp);
			Vector3D temp2(temp(0), temp(1), temp(2));
			prim->setX(temp2+bestSim.COR+bestSim.trans);
			prim->setR(prim->getR()*bestSim.scale);

			prim->rotateBy(bestSim.rotQ);
		}
		else {
			primitiveId=linkIds[j];
            prim=figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;
			prim->setX(newLinkAtoms[j]->getX());
			prim->setQ(newLinkAtoms[j]->getQ());
			prim->setR(newLinkAtoms[j]->getR());
			prim->setTheta(newLinkAtoms[j]->getTheta());
			if(prim->type()==M3D_END_PRIMITIVE)
				(dynamic_cast<M3DQuadEndPrimitive *>(prim))->setElongation((dynamic_cast<M3DQuadEndPrimitive *>(newLinkAtoms[j]))->getElongation());
		}
	}

	delete []primitive;
	delete []x;
	delete []y;

	return;
}



