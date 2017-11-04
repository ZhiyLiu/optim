#include <math.h>
#include <stdio.h>

#define D_XFERLIST
#define D_POINTLIST_SERVER2
#include "Shapedepend.h"

#include "matrix.h"
#include "M3DSubfigureTransformation.h"
#include "M3DQuadFigure.h"

/*
#include "M3DSubfigureTransformation.h"
#include "M3DQuadFigure.h"
#include <math.h>
#include <stdio.h>

#define D_XFERLIST
#include "Shapedepend.h"
*/
//#define QIONG_DEBUG_LINALG
const double MAX_CREST_T = 1.0 - R_SMALL_TOLERANCE;
const double MIN_CREST_T = -1.0 + R_SMALL_TOLERANCE;

const double LOG_BASE = 1.2;
const double INVERSED_LOG_BASE = 1.0/LOG_BASE;

M3DSubfigureTransformation::M3DSubfigureTransformation(M3DObject * _object,
                                M3DFigureTreeNode * _figureTreeNode)
{
    referenceFigure = NULL;

    init(_object, _figureTreeNode);
}

void M3DSubfigureTransformation::init(M3DObject * _object, M3DFigureTreeNode * _figureTreeNode)
{
    int figureId;
    M3DFigure * figure;

    object= _object;
    figureTreeNode = _figureTreeNode;

    if(figureTreeNode != NULL && object != NULL)
    {
		M3DFigureTreeNode *parentNode=figureTreeNode->getParent();
		M3DQuadFigure *parentFigure=dynamic_cast<M3DQuadFigure*>((object->getFigurePtr(parentNode->getFigureId())));
		maxParentU=parentFigure->getRowCount()-1;
		maxParentV=parentFigure->getColumnCount()-1;
		factorT2U=1.0;
		factorT2V=1.0;

        figureId = figureTreeNode->getFigureId();
        figure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(figureId));
        if(figure != NULL)
            referenceFigure = figure->assign();
			referenceFigureOld = referenceFigure->assign();
    }

    invalidLinks = false;
    checkForValidLinks();

    lastScaleValue = 1.0;
    lastElongateValue = 1.0;
}

void M3DSubfigureTransformation::initStep()
{
    int figureId;
    M3DFigure * figure;

    if(figureTreeNode != NULL && object != NULL)
    {
        figureId = figureTreeNode->getFigureId();
        figure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(figureId));
//		cout << "FIGURE: " << endl;
//		dynamic_cast<M3DQuadFigure*>(figure)->print(true);
//		cout << endl << "REFERENCE:" << endl;
//		dynamic_cast<M3DQuadFigure*>(referenceFigure)->print(true);
        if(figure != NULL)
            (* dynamic_cast<M3DQuadFigure*>(referenceFigure)) = (* dynamic_cast<M3DQuadFigure*>(figure));
//		cout << "FIGURE: " << endl;
//		dynamic_cast<M3DQuadFigure*>(figure)->print(true);
//		cout << endl << "REFERENCE:" << endl;
//		dynamic_cast<M3DQuadFigure*>(referenceFigure)->print(true);

    }
}

double AngleBetweenVectors(Vector3D v1, Vector3D v2)
{
	double dotP=v1*v2;
	if(dotP>1||dotP<-1)
	{
#ifdef QIONG_DEBUG_LINALG
		printf("Error: M3DSubfigureTransformation.cpp/AngelBetweenVectors, v1*v2 is out of range [-1, 1]!\n");
#endif
		if(dotP>1)
			return 0;
		else
			return M_PI;
	}
	else
	{
		return acos(dotP);
	}
}

void M3DSubfigureTransformation::reInit()
{
    int figureId;
    M3DFigure * figure;

    if(figureTreeNode != NULL && object != NULL)
    {
        figureId = figureTreeNode->getFigureId();
        figure = object->getFigurePtr(figureId);
        if(figure != NULL)
            (* referenceFigure) = (* referenceFigureOld);
    }

	// these two lies are actually crucial to use the previous transformed sub-figure 
	// as the reference figure in the next step!
    lastScaleValue = 1.0;
    lastElongateValue = 1.0;
}

M3DSubfigureTransformation::~M3DSubfigureTransformation()
{
    if(referenceFigure != NULL)
        delete referenceFigure;
}


void M3DSubfigureTransformation::updateSubfigure(M3DQuadFigure * oldParentFigure)
{
    M3DFigureTreeNode * parentNode;
    M3DQuadFigure * figure;
    M3DQuadFigure * parentFigure;
    M3DPrimitive * primitive;
    M3DQuadPrimitive startPrimitive;
    M3DQuadPrimitive finishPrimitive;
    Vector3D oldNormal, newNormal;
    Vector3D startPosition, endPosition;
    Vector3D deltaX;
    Vector3D axis;
    Quat deltaQ;
	bool doRotations;
    Real norm;
	double angle;
    M3DPrimitiveLinkInfo * link;
	int deltaIndex, startIndex, middleIndex;

    double maxU, maxV;
	int theFlag=0;

    int linkCount,
        rowCount, columnCount,
        i, j,
        primitiveId,
        primitiveIndexDelta,
        neighborIndexDelta,
        numPrimitivesPerLink;

    int prevPrimitiveId,
        neighborPrimitiveId;




    if(figureTreeNode == NULL || object == NULL || referenceFigure == NULL)
        return;

    figure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(figureTreeNode->getFigureId()));
    if(figure == NULL)
        return;

    parentNode = figureTreeNode->getParent();
    if(parentNode == NULL)
        return;

    parentFigure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(parentNode->getFigureId()));

    rowCount = figure->getRowCount();
    columnCount = figure->getColumnCount();

	maxU = (double) (parentFigure->getColumnCount() - 1);
    maxV = (double) (parentFigure->getRowCount() - 1);

    if(rowOrientedLinks)
    {
        primitiveIndexDelta = columnCount;
        numPrimitivesPerLink = rowCount;
        neighborIndexDelta = 1;
    }
    else
    {
        primitiveIndexDelta = 1;
        numPrimitivesPerLink = columnCount;
        neighborIndexDelta = columnCount;
    }

	//		integrate the '(u, v, t)->BPosition and BNormal' into this function
	ThallCode::Pointlist_server2 *pList1=new ThallCode::Pointlist_server2();
	Xferlist *xferList = convertM3DtoXfer(oldParentFigure);
	pList1->InitializeSubdivSurf(xferList);
	delete []xferList->atomlist;
	delete xferList;

	ThallCode::Pointlist_server2 *pList2=new ThallCode::Pointlist_server2();
	xferList = convertM3DtoXfer(parentFigure);
	pList2->InitializeSubdivSurf(xferList);
	delete []xferList->atomlist;
	delete xferList;

    linkCount = figureTreeNode->getLinkCount();
	middleIndex = linkCount / 2;
	startIndex = middleIndex;

	link = figureTreeNode->getLink(middleIndex);

	for(deltaIndex = -1; deltaIndex <= 1; deltaIndex += 2)
    {
        for(i = startIndex; i < linkCount && i >= 0; i += deltaIndex)
        {				
            link = figureTreeNode->getLink(i);

			getSurfacePointAndNormal(pList1, startPosition, oldNormal,
                                     link->u, link->v, link->t);
			getSurfacePointAndNormal(pList2, endPosition, newNormal,
                                     link->u, link->v, link->t);       

			primitiveId = link->primitiveId;
            primitive = figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;

            deltaX = endPosition - startPosition;
            axis = oldNormal.cross(newNormal);

			angle=AngleBetweenVectors(oldNormal, newNormal);
            //angle = acos(oldNormal * newNormal);

			deltaQ.setAxisAngle(axis, angle);

            // Make sure rotation is not zero
            norm = deltaQ.norm();
            if(norm > 1.0 - R_SMALL_TOLERANCE &&
               norm < 1.0 + R_SMALL_TOLERANCE)
                doRotations = true;
            else
                doRotations = false;

            primitive = figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;

            if(doRotations)
            {
                primitive->rotateBy(deltaQ);

                axis = primitive->getX() - startPosition;
                deltaQ.rotateVector(axis);
                primitive->setX(axis + startPosition);
            }

            primitive->translateBy(deltaX);
			//primitive->setX(endPosition);

            if(i == middleIndex)
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;

                    setPrimitiveToPreviousPrediction(figure, primitiveId,
                        prevPrimitiveId, lastElongateValue);
                }
            }

            else
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;
                    neighborPrimitiveId = primitiveId -
                        deltaIndex * neighborIndexDelta;

                    setPrimitiveToAveragePrediction(figure, primitiveId,
                        prevPrimitiveId, neighborPrimitiveId, lastScaleValue,
                        lastElongateValue);
                }
            }
			theFlag = 0;
        }

        startIndex = middleIndex + 1;
    }

	delete pList1;
	delete pList2;
}




void M3DSubfigureTransformation::translateStep(double theDeltaU, double theDeltaV)
{
    M3DFigureTreeNode * parentNode;
    M3DQuadFigure * figure;
    M3DQuadFigure * parentFigure;
    M3DPrimitive * primitive;
    M3DQuadPrimitive startPrimitive;
    M3DQuadPrimitive finishPrimitive;
    Vector3D oldNormal, newNormal;
    Vector3D startPosition, endPosition;
    Vector3D deltaX;
    Vector3D axis;
    double angle, deltaU, deltaV, deltaU1, deltaU2, deltaV1, deltaV2, deltaT;
	double centerU, centerV, centerT;
	double alphaU=factorT2U;
	double alphaV=factorT2V;
    Quat deltaQ;
    M3DPrimitiveLinkInfo * link;
	int deltaIndex, startIndex, middleIndex;

    double maxU, maxV;
    double u, v, t;
	int flag, theFlag=0;
    bool doRotations;
    Real norm;

    int linkCount,
        rowCount, columnCount,
        i, j,
        primitiveId,
        primitiveIndexDelta,
        neighborIndexDelta,
        numPrimitivesPerLink;

    int prevPrimitiveId,
        neighborPrimitiveId;




    if(figureTreeNode == NULL || object == NULL || referenceFigure == NULL)
        return;

    figure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(figureTreeNode->getFigureId()));
    if(figure == NULL)
        return;

    parentNode = figureTreeNode->getParent();
    if(parentNode == NULL)
        return;

    parentFigure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(parentNode->getFigureId()));

    rowCount = figure->getRowCount();
    columnCount = figure->getColumnCount();

    //maxU = (double) (parentFigure->getRowCount() - 1);
    //maxV = (double) (parentFigure->getColumnCount() - 1);
	maxU = (double) (parentFigure->getColumnCount() - 1);
    maxV = (double) (parentFigure->getRowCount() - 1);

    if(rowOrientedLinks)
    {
        primitiveIndexDelta = columnCount;
        numPrimitivesPerLink = rowCount;
        neighborIndexDelta = 1;
    }
    else
    {
        primitiveIndexDelta = 1;
        numPrimitivesPerLink = columnCount;
        neighborIndexDelta = columnCount;
    }

	//		integrate the '(u, v, t)->BPosition and BNormal' into this function
	ThallCode::Pointlist_server2 *pList2=new ThallCode::Pointlist_server2();
	Xferlist *xferList = convertM3DtoXfer(parentFigure);
	pList2->InitializeSubdivSurf(xferList);
	delete []xferList->atomlist;
	delete xferList;

    linkCount = figureTreeNode->getLinkCount();
	middleIndex = linkCount / 2;
	startIndex = middleIndex;


	for (i = 0; i < linkCount; i++)
	{
		link = figureTreeNode->getLink(i);
		if (link->t > -1)
			if (link->u == 0 || link->u == maxU)
				flag  = 1;
			else if (link->v == 0 || link->v == maxU)
				flag = 2;
			else
				flag = 0;
	}


	link = figureTreeNode->getLink(middleIndex);
	centerU = link->u;
	centerV = link->v;
	centerT = link->t;

	for(deltaIndex = -1; deltaIndex <= 1; deltaIndex += 2)
    {
        for(i = startIndex; i < linkCount && i >= 0; i += deltaIndex)
        {				
            link = figureTreeNode->getLink(i);
			if ((link->u == 0 || link->u == maxU) && (centerV == 0 || centerV == maxV))
				theFlag = 1;
			if ((link->v == 0 || link->v == maxU) && (centerU == 0 || centerU == maxV))
				theFlag = 1;

            if(link == NULL)
                continue;

			if (link->t >-1 ){
				deltaU=theDeltaU;
				deltaV=theDeltaV;
			}
			else{// on the bottom face
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
			u=link->u; v=link->v; t=link->t;
			//printf("i= %d, before: u= %f, v= %f t= %f deltaU= %f deltaV= %f\n", i, u, v, t, deltaU, deltaV);
			if (link->t == -1 || link->t == 1){ 
			// Neither u nor v is in the crest region.
				if (link->u + deltaU <= maxU && v + deltaV <= maxV && 
						link->u + deltaU >= 0 && v + deltaV >= 0){
					//Neither u+deltaU nor v+deltaV is in the crest region.
					u=link->u+deltaU; v=link->v+deltaV;
				}
				else{
				    if (link->u + deltaU > maxU){
						// u+deltaU is in the crest region.
						deltaU1=maxU-link->u;
						deltaU2=deltaU-deltaU1;
						u=link->u+deltaU1; 
						deltaT=deltaU2*alphaU;	//deltaT=deltaU2/alpha;
						//vDelta=link->v/maxV*(maxV+2*alpha)-link->v-alpha;
						//deltaT=deltaU2/sqrt(alpha*alpha+vDelta*vDelta);
						if (link->t == 1)
							t=link->t-fabs(deltaT);
						else
							t=link->t+fabs(deltaT);
						}
					else if (link->u + deltaU < 0){
						// u+deltaU is in the crest region.
						deltaU1=0-link->u;
						deltaU2=deltaU-deltaU1;
						u=link->u+deltaU1; 
						deltaT=deltaU2*alphaU;	//deltaT=deltaU2/alpha;
						//vDelta=link->v/maxV*(maxV+2*alpha)-link->v-alpha;
						//deltaT=deltaU2/sqrt(alpha*alpha+vDelta*vDelta);
						if (link->t == 1)
							t=link->t-fabs(deltaT);
						else
							t=link->t+fabs(deltaT);
					}
					else
						u=link->u+deltaU;

					if ( link->v + deltaV > maxV){ 
						//link->v+deltaV is in the crest region.
						deltaV1=maxV-link->v;
						deltaV2=deltaV-deltaV1;
						v=link->v+deltaV1;
						deltaT=deltaV2*alphaV;	//deltaT=deltaV2/alpha;
						//uDelta=link->u/maxU*(maxU+2*alpha)-link->u-alpha;
						//deltaT=deltaV2/sqrt(alpha*alpha+uDelta*uDelta);
						if (link->t == 1)
							t=link->t-fabs(deltaT);
						else
							t=link->t+fabs(deltaT); 
					}
					else if (link->v + deltaV < 0){
						//link->v + deltaV < 0) 
						//link->v+deltaV is in the crest region.
						deltaV1=0-link->v;
						deltaV2=deltaV-deltaV1;
						v=link->v+deltaV1;
						deltaT=deltaV2*alphaV;	//deltaT=deltaV2/alpha;
						//uDelta=link->u/maxU*(maxU+2*alpha)-link->u-alpha;
						//deltaT=deltaV2/sqrt(alpha*alpha+uDelta*uDelta);
						if (link->t == 1)
							t=link->t-fabs(deltaT);
						else
							t=link->t+fabs(deltaT); 
					}
					else
						v=link->v+deltaV;
				}
			}
			else if ((link->u == 0) || (link->u == maxU)){
				// u is in the crest region.				
				if (deltaU != 0){ // moving in u direction.	
					if (theFlag == 1){
						// actually moving in v direction
						if ( link->v == 0)
							v=fabs(deltaU);
						else
							v=link->v-fabs(deltaU);
					}				
					else{
						//vDelta=link->v/maxV*(maxV+2*alpha)-link->v-alpha;
						//if (link->t >= 0)
						if (link->u == maxU)
							deltaT=-deltaU*alphaU;	//deltaT=-deltaU/alpha;
							//deltaT=-deltaU/sqrt(alpha*alpha+vDelta*vDelta);
						else 
							deltaT=deltaU*alphaU;	//deltaT=deltaU/alpha;
							//deltaT=deltaU/sqrt(alpha*alpha+vDelta*vDelta);
						if (link->t+deltaT >= -1 && link->t+deltaT <= 1){
							// still ends in the crest region.
							t=link->t+deltaT; 
						}
						else{
							// ends in the non-crest region.
							if (link->t+deltaT < -1){
								deltaT=-1-link->t;
								t=-1;
								deltaU1=deltaT/alphaU;	//deltaU1=deltaT*alpha;
								//deltaU1=deltaT*sqrt(alpha*alpha+vDelta*vDelta);
								deltaU2=deltaU-deltaU1;
								if (link->u == 0)
									u=link->u+fabs(deltaU-deltaU1);
								else
									u=link->u-fabs(deltaU+deltaU1);
							}
							else{
								deltaT=1-link->t;
								t=1;
								deltaU1=deltaT/alphaU;	//deltaU1=deltaT*alpha;
								//deltaU1=deltaT*sqrt(alpha*alpha+vDelta*vDelta);
								//deltaU2=deltaU-deltaU1;
								if (link->u == 0)
									u=link->u+fabs(deltaU-deltaU1);
								else
									u=link->u-fabs(deltaU+deltaU1);						
							}
						}
					}
				}
				else { // moving in v direction.
					if (link->v+deltaV >= 0 && link->v+deltaV <= maxV){
						// still in the range of V, normal movement
						v=link->v+deltaV;
					}
					else{
						if (link->v+deltaV < 0){
							deltaV1=0-link->v;
							v=0;
							deltaU2=deltaV-deltaV1;
							if ( link->u == 0)
								u=fabs(deltaU2);
							else
								u=link->u-fabs(deltaU2);
						}
						else{
							deltaV1=maxV-link->v;
							v=maxV;
							deltaU2=deltaV-deltaV1;
							if ( link->u == 0)
								u=fabs(deltaU2);
							else
								u=link->u-fabs(deltaU2);
						}
					}
				}
			}
			else{
				// v is in the crest region.
				if (deltaV != 0){ // moving in v direction.
				if (theFlag == 1){
					// actually moving in u direction
					if ( link->u == 0)
							u=fabs(deltaV);
						else
							u=link->u-fabs(deltaV);
				}				
				else{
					deltaT=deltaV*alphaV;	//deltaT=deltaV/alpha;
					//uDelta=link->u/maxU*(maxU+2*alpha)-link->u-alpha;
					//deltaT=deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					//if (link->t >= 0)
					if (link->v == maxV)
						deltaT=-deltaV*alphaV;	//deltaT=-deltaV/alpha;
						//deltaT=-deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					else 
						deltaT=deltaV*alphaV;	//deltaT=deltaV/alpha;
						//deltaT=deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					if (link->t+deltaT >= -1 && link->t+deltaT <= 1){
						// still ends in the crest region.
						t=link->t+deltaT;
					}
					else{
						// ends in the non-crest region.
						if (link->t+deltaT < -1){
							deltaT=-1-link->t;
							t=-1;
							deltaV1=deltaT/alphaV;	//deltaV1=deltaT*alpha;
							//deltaV1=deltaT*sqrt(alpha*alpha+uDelta*uDelta);
							//deltaV2=deltaV-deltaV1;
							if (link->v == 0)
								v=link->v+fabs(deltaV-deltaV1);
							else
								v=link->v-fabs(deltaV+deltaV1);
						}
						else{
							deltaT=1-link->t;
							t=1;
							deltaV1=deltaT/alphaV;	//deltaV1=deltaT*alpha;
							//deltaV1=deltaT*sqrt(alpha*alpha+uDelta*uDelta);
							//deltaV2=deltaV-deltaV1;
							if (link->v == 0)
								v=link->v+fabs(deltaV-deltaV1);
							else
								v=link->v-fabs(deltaV+deltaV1);						
						}
					}
				}
				}
				else { // moving in u direction.
					if (link->u+deltaU <= maxU && link->u+deltaU >= 0){
						// still in the range of U, normal movement
						u=link->u+deltaU;
					}
					else{
						if (link->u+deltaU < 0){
							deltaU1=0-link->u;
							u=0;
							deltaV2=deltaU-deltaU1;
							if ( link->v == 0)
								v=fabs(deltaV2);
							else
								v=link->v-fabs(deltaV2);
						}
						else{ // link->u+deltaU <= maxU
							deltaU1=maxU-link->u;
							u=maxU;
							deltaV2=deltaU-deltaU1;
							if ( link->v == 0)
								v=fabs(deltaU2);
							else
								v=link->v-fabs(deltaV2);
						}
					}
				}
			}


            getSurfacePointAndNormal(pList2, startPosition, oldNormal,
                                     link->u, link->v, link->t);

            getSurfacePointAndNormal(pList2, endPosition, newNormal,
                                     u, v, t);


            //getSurfacePointAndNormal(startPosition, oldNormal, parentFigure,
            //                         link->u, link->v, link->t);

            //getSurfacePointAndNormal(endPosition, newNormal, parentFigure,
            //                         u, v, t);

            link->u = u;
            link->v = v;
            link->t = t;
			//printf("i= %d, after: u= %f, v= %f t= %f\n", i, u, v, t);

            primitiveId = link->primitiveId;
            primitive = figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;

            deltaX = endPosition - startPosition;
            axis = oldNormal.cross(newNormal);

			angle=AngleBetweenVectors(oldNormal, newNormal);
            //angle = acos(oldNormal * newNormal);

			deltaQ.setAxisAngle(axis, angle);

            // Make sure rotation is not zero
            norm = deltaQ.norm();
            if(norm > 1.0 - R_SMALL_TOLERANCE &&
               norm < 1.0 + R_SMALL_TOLERANCE)
                doRotations = true;
            else
                doRotations = false;

            primitive = figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;

            if(doRotations)
            {
                primitive->rotateBy(deltaQ);

                axis = primitive->getX() - startPosition;
                deltaQ.rotateVector(axis);
                primitive->setX(axis + startPosition);
            }

            primitive->translateBy(deltaX);

            if(i == middleIndex)
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;

                    setPrimitiveToPreviousPrediction(figure, primitiveId,
                        prevPrimitiveId, lastElongateValue);
                }
            }

            else
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;
                    neighborPrimitiveId = primitiveId -
                        deltaIndex * neighborIndexDelta;

                    setPrimitiveToAveragePrediction(figure, primitiveId,
                        prevPrimitiveId, neighborPrimitiveId, lastScaleValue,
                        lastElongateValue);
                }
            }
			theFlag = 0;
        }

        startIndex = middleIndex + 1;
    }

	delete pList2;
}


void M3DSubfigureTransformation::translate(double deltaU, double deltaV)
{
	int i;
	int uIterations, vIterations;
	double delta, lastDeltaU, lastDeltaV;

	delta = .25;

	// regularize U by the circular u space
	while(deltaU>=2*(maxParentU+2*factorT2U))
		deltaU-=2*(maxParentU+2*factorT2U);
	while(deltaU<=-2*(maxParentU+2*factorT2U))
		deltaU+=2*(maxParentU+2*factorT2U);
	if(deltaU>maxParentU+2*factorT2U)
	{
		deltaU=-(2*(maxParentU+2*factorT2U)-deltaU);
	}
	else 
		if(deltaU<-maxParentU-2*factorT2U)
		{
			deltaU=2*(maxParentU+2*factorT2U)+deltaU;
		}

	// regularize V by the circular v space
	while(deltaV>=2*(maxParentV+2*factorT2V))
		deltaV-=2*(maxParentV+2*factorT2V);
	while(deltaV<=-2*(maxParentV+2*factorT2V))
		deltaV+=2*(maxParentV+2*factorT2V);
	if(deltaV>maxParentV+2*factorT2V)
	{
		deltaV=-(2*(maxParentV+2*factorT2V)-deltaV);
	}
	else 
		if(deltaV<-maxParentV-2*factorT2V)
		{
			deltaV=2*(maxParentV+2*factorT2V)+deltaV;
		}

#ifdef DEBUG
	fprintf(stderr, "dU: %f, dV: %f in ACTUAL translation.\n", deltaU, deltaV);
#endif

	uIterations = (int) fabs(deltaU / delta);
	if (deltaU > 0.0)
		lastDeltaU = deltaU - delta * uIterations;
	else 
		lastDeltaU = deltaU + delta * uIterations;

	for ( i=0; i < uIterations; i++)
	{
		if (deltaU > 0.0) 
			translateStep(delta, 0);
		else 
			translateStep(-delta, 0);
	}
	translateStep(lastDeltaU, 0);

	vIterations = (int) fabs(deltaV / delta);
	if (deltaV > 0.0)
		lastDeltaV = deltaV - delta * vIterations;
	else 
		lastDeltaV = deltaV + delta * vIterations;

	for ( i=0; i < vIterations; i++)
	{
		if (deltaV > 0.0) 
			translateStep(0, delta);
		else 
			translateStep(0, -delta);
		initStep();
	}
	translateStep(0, lastDeltaV);
//	translateStep(deltaU, deltaV);
	reInit();
}


/*
void M3DSubfigureTransformation::translateHingeAtom(double theDeltaU, double theDeltaV, int hingeIndex)
{
    M3DFigureTreeNode * parentNode;
    M3DQuadFigure * figure;
    M3DQuadFigure * parentFigure;
    M3DPrimitive * primitive;
    M3DPrimitive startPrimitive;
    M3DPrimitive finishPrimitive;
    Vector3D oldNormal, newNormal;
    Vector3D startPosition, endPosition;
    Vector3D deltaX;
    Vector3D axis;
    double angle, deltaU, deltaV, deltaU1, deltaU2, deltaV1, deltaV2, deltaT;
	double centerU, centerV, centerT;
	double alpha=1.0;
    Quat deltaQ;
    M3DPrimitiveLinkInfo * link;

    double maxU, maxV;
    double u, v, t;
	int theFlag=0;
    bool doRotations;
    Real norm;

    int primitiveId;


	if(figureTreeNode == NULL || object == NULL || referenceFigure == NULL)
        return;

    figure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(figureTreeNode->getFigureId()));
    if(figure == NULL)
        return;

    parentNode = figureTreeNode->getParent();
    if(parentNode == NULL)
        return;

    parentFigure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(parentNode->getFigureId()));

	maxU = (double) (parentFigure->getColumnCount() - 1);
    maxV = (double) (parentFigure->getRowCount() - 1);



	//		integrate the '(u, v, t)->BPosition and BNormal' into this function
	ThallCode::Pointlist_server2 *pList2=new ThallCode::Pointlist_server2();
	Xferlist *xferList = convertM3DtoXfer(parentFigure);
	pList2->InitializeSubdivSurf(xferList);
	delete []xferList->atomlist;
	delete xferList;



	link = figureTreeNode->getLink(hingeIndex);
	centerU = link->u;
	centerV = link->v;
	centerT = link->t;



			if ((link->u == 0 || link->u == maxU) && (centerV == 0 || centerV == maxV))
				theFlag = 1;
			if ((link->v == 0 || link->v == maxU) && (centerU == 0 || centerU == maxV))
				theFlag = 1;

			if (link->t >-1 ){
				deltaU=theDeltaU;
				deltaV=theDeltaV;
			}
			else{
				deltaU=-theDeltaU;
				deltaV=-theDeltaV;
			}
			u=link->u; v=link->v; t=link->t;
			//printf("i= %d, before: u= %f, v= %f t= %f deltaU= %f deltaV= %f\n", i, u, v, t, deltaU, deltaV);
			if (link->t == -1 || link->t == 1){ 
			// Neither u nor v is in the crest region.
				if (link->u + deltaU <= maxU && v + deltaV <= maxV && 
						link->u + deltaU >= 0 && v + deltaV >= 0){
					//Neither u+deltaU nor v+deltaV is in the crest region.
					u=link->u+deltaU; v=link->v+deltaV;
				}
				else{
				    if (link->u + deltaU > maxU){
						// u+deltaU is in the crest region.
						deltaU1=maxU-link->u;
						deltaU2=deltaU-deltaU1;
						u=link->u+deltaU1; 
						deltaT=deltaU2/alpha;
						//vDelta=link->v/maxV*(maxV+2*alpha)-link->v-alpha;
						//deltaT=deltaU2/sqrt(alpha*alpha+vDelta*vDelta);
						if (link->t == 1)
							t=link->t-fabs(deltaT);
						else
							t=link->t+fabs(deltaT);
						}
					else if (link->u + deltaU < 0){
						// u+deltaU is in the crest region.
						deltaU1=0-link->u;
						deltaU2=deltaU-deltaU1;
						u=link->u+deltaU1; 
						deltaT=deltaU2/alpha;
						//vDelta=link->v/maxV*(maxV+2*alpha)-link->v-alpha;
						//deltaT=deltaU2/sqrt(alpha*alpha+vDelta*vDelta);
						if (link->t == 1)
							t=link->t-fabs(deltaT);
						else
							t=link->t+fabs(deltaT);
					}
					else
						u=link->u+deltaU;

					if ( link->v + deltaV > maxV){ 
						//link->v+deltaV is in the crest region.
						deltaV1=maxV-link->v;
						deltaV2=deltaV-deltaV1;
						v=link->v+deltaV1;
						deltaT=deltaV2/alpha;
						//uDelta=link->u/maxU*(maxU+2*alpha)-link->u-alpha;
						//deltaT=deltaV2/sqrt(alpha*alpha+uDelta*uDelta);
						if (link->t == 1)
							t=link->t-fabs(deltaT);
						else
							t=link->t+fabs(deltaT); 
					}
					else if (link->v + deltaV < 0){
						//link->v + deltaV < 0) 
						//link->v+deltaV is in the crest region.
						deltaV1=0-link->v;
						deltaV2=deltaV-deltaV1;
						v=link->v+deltaV1;
						deltaT=deltaV2/alpha;
						//uDelta=link->u/maxU*(maxU+2*alpha)-link->u-alpha;
						//deltaT=deltaV2/sqrt(alpha*alpha+uDelta*uDelta);
						if (link->t == 1)
							t=link->t-fabs(deltaT);
						else
							t=link->t+fabs(deltaT); 
					}
					else
						v=link->v+deltaV;
				}
			}
			else if ((link->u == 0) || (link->u == maxU)){
				// u is in the crest region.				
				if (deltaU != 0){ // moving in u direction.	
					if (theFlag == 1){
						// actually moving in v direction
						if ( link->v == 0)
							v=fabs(deltaU);
						else
							v=link->v-fabs(deltaU);
					}				
					else{
						//vDelta=link->v/maxV*(maxV+2*alpha)-link->v-alpha;
						//if (link->t >= 0)
						if (link->u == maxU)
							deltaT=-deltaU/alpha;
							//deltaT=-deltaU/sqrt(alpha*alpha+vDelta*vDelta);
						else 
							deltaT=deltaU/alpha;
							//deltaT=deltaU/sqrt(alpha*alpha+vDelta*vDelta);
						if (link->t+deltaT >= -1 && link->t+deltaT <= 1){
							// still ends in the crest region.
							t=link->t+deltaT; 
						}
						else{
							// ends in the non-crest region.
							if (link->t+deltaT < -1){
								deltaT=-1-link->t;
								t=-1;
								deltaU1=deltaT*alpha;
								//deltaU1=deltaT*sqrt(alpha*alpha+vDelta*vDelta);
								deltaU2=deltaU-deltaU1;
								if (link->u == 0)
									u=link->u+fabs(deltaU-deltaU1);
								else
									u=link->u-fabs(deltaU+deltaU1);
							}
							else{
								deltaT=1-link->t;
								t=1;
								deltaU1=deltaT*alpha;
								//deltaU1=deltaT*sqrt(alpha*alpha+vDelta*vDelta);
								//deltaU2=deltaU-deltaU1;
								if (link->u == 0)
									u=link->u+fabs(deltaU-deltaU1);
								else
									u=link->u-fabs(deltaU+deltaU1);						
							}
						}
					}
				}
				else { // moving in v direction.
					if (link->v+deltaV >= 0 && link->v+deltaV <= maxV){
						// still in the range of V, normal movement
						v=link->v+deltaV;
					}
					else{
						if (link->v+deltaV < 0){
							deltaV1=0-link->v;
							v=0;
							deltaU2=deltaV-deltaV1;
							if ( link->u == 0)
								u=fabs(deltaU2);
							else
								u=link->u-fabs(deltaU2);
						}
						else{
							deltaV1=maxV-link->v;
							v=maxV;
							deltaU2=deltaV-deltaV1;
							if ( link->u == 0)
								u=fabs(deltaU2);
							else
								u=link->u-fabs(deltaU2);
						}
					}
				}
			}
			else{
				// v is in the crest region.
				if (deltaV != 0){ // moving in v direction.
				if (theFlag == 1){
					// actually moving in u direction
					if ( link->u == 0)
							u=fabs(deltaV);
						else
							u=link->u-fabs(deltaV);
				}				
				else{
					deltaT=deltaV/alpha;
					//uDelta=link->u/maxU*(maxU+2*alpha)-link->u-alpha;
					//deltaT=deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					//if (link->t >= 0)
					if (link->v == maxV)
						deltaT=-deltaV/alpha;
						//deltaT=-deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					else 
						deltaT=deltaV/alpha;
						//deltaT=deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					if (link->t+deltaT >= -1 && link->t+deltaT <= 1){
						// still ends in the crest region.
						t=link->t+deltaT;
					}
					else{
						// ends in the non-crest region.
						if (link->t+deltaT < -1){
							deltaT=-1-link->t;
							t=-1;
							deltaV1=deltaT*alpha;
							//deltaV1=deltaT*sqrt(alpha*alpha+uDelta*uDelta);
							//deltaV2=deltaV-deltaV1;
							if (link->v == 0)
								v=link->v+fabs(deltaV-deltaV1);
							else
								v=link->v-fabs(deltaV+deltaV1);
						}
						else{
							deltaT=1-link->t;
							t=1;
							deltaV1=deltaT*alpha;
							//deltaV1=deltaT*sqrt(alpha*alpha+uDelta*uDelta);
							//deltaV2=deltaV-deltaV1;
							if (link->v == 0)
								v=link->v+fabs(deltaV-deltaV1);
							else
								v=link->v-fabs(deltaV+deltaV1);						
						}
					}
				}
				}
				else { // moving in u direction.
					if (link->u+deltaU <= maxU && link->u+deltaU >= 0){
						// still in the range of U, normal movement
						u=link->u+deltaU;
					}
					else{
						if (link->u+deltaU < 0){
							deltaU1=0-link->u;
							u=0;
							deltaV2=deltaU-deltaU1;
							if ( link->v == 0)
								v=fabs(deltaV2);
							else
								v=link->v-fabs(deltaV2);
						}
						else{ // link->u+deltaU <= maxU
							deltaU1=maxU-link->u;
							u=maxU;
							deltaV2=deltaU-deltaU1;
							if ( link->v == 0)
								v=fabs(deltaV2);
							else
								v=link->v-fabs(deltaV2);
						}
					}
				}
			}


            getSurfacePointAndNormal(pList2, startPosition, oldNormal,
                                     link->u, link->v, link->t);

            getSurfacePointAndNormal(pList2, endPosition, newNormal,
                                     u, v, t);


            //getSurfacePointAndNormal(startPosition, oldNormal, parentFigure,
            //                         link->u, link->v, link->t);

            //getSurfacePointAndNormal(endPosition, newNormal, parentFigure,
            //                         u, v, t);

            link->u = u;
            link->v = v;
            link->t = t;
			//printf("i= %d, after: u= %f, v= %f t= %f\n", i, u, v, t);

            primitiveId = link->primitiveId;
            primitive = figure->getPrimitivePtr(primitiveId);

            deltaX = endPosition - startPosition;
            axis = oldNormal.cross(newNormal);

			angle=AngleBetweenVectors(oldNormal, newNormal);
            //angle = acos(oldNormal * newNormal);

			deltaQ.setAxisAngle(axis, angle);

            // Make sure rotation is not zero
            norm = deltaQ.norm();
            if(norm > 1.0 - R_SMALL_TOLERANCE &&
               norm < 1.0 + R_SMALL_TOLERANCE)
                doRotations = true;
            else
                doRotations = false;

            primitive = figure->getPrimitivePtr(primitiveId);

            if(doRotations)
            {
                primitive->rotateBy(deltaQ);

                axis = primitive->getX() - startPosition;
                deltaQ.rotateVector(axis);
                primitive->setX(axis + startPosition);
            }

            primitive->translateBy(deltaX);                    

	delete pList2;
}
*/



void M3DSubfigureTransformation::rotateStep(double theta)
{
    M3DFigureTreeNode * parentNode;
    M3DQuadFigure * figure;
    M3DQuadFigure * parentFigure;
    M3DPrimitive * primitive;
    M3DQuadPrimitive startPrimitive;
    M3DQuadPrimitive finishPrimitive;
    M3DQuadPrimitive centerPrimitive;
    Vector3D oldNormal, newNormal;
    Vector3D startPosition, endPosition;
    Vector3D deltaX;
    Vector3D axis;
    Quat deltaQ;
    M3DPrimitiveLinkInfo * link;

    double maxU, maxV;
    double centerU, centerV, centerT, tempU, tempV, tempUV, tempT1, tempT2;
    double cosTheta,
           sinTheta;
    double angle;
    double u, v, t,
           uDiff, vDiff,
           uRot, vRot;

	double  deltaU, deltaV, deltaU1, deltaU2, deltaV1, deltaV2, deltaT;
	double alpha=1.0;

    int linkCount,
        rowCount, columnCount,
        i, j,
        primitiveId,
        primitiveIndexDelta,
        neighborIndexDelta,
        numPrimitivesPerLink;

    int prevPrimitiveId,
        neighborPrimitiveId;


    if(figureTreeNode == NULL || object == NULL)
        return;

    figure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(figureTreeNode->getFigureId()));
    if(figure == NULL)
        return;

    parentNode = figureTreeNode->getParent();
    if(parentNode == NULL)
        return;

    parentFigure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(parentNode->getFigureId()));

    rowCount = figure->getRowCount();
    columnCount = figure->getColumnCount();

    maxU = (double) (parentFigure->getColumnCount() - 1);
    maxV = (double) (parentFigure->getRowCount() - 1);

    if(rowOrientedLinks)
    {
        primitiveIndexDelta = columnCount;
        numPrimitivesPerLink = rowCount;
        neighborIndexDelta = 1;
    }
    else
    {
        primitiveIndexDelta = 1;
        numPrimitivesPerLink = columnCount;
        neighborIndexDelta = columnCount;
    }

	//		integrate the '(u, v, t)->BPosition and BNormal' into this function
	ThallCode::Pointlist_server2 *pList2=new ThallCode::Pointlist_server2();
	Xferlist *xferList = convertM3DtoXfer(parentFigure);
	pList2->InitializeSubdivSurf(xferList);
	delete []xferList->atomlist;
	delete xferList;

	link = figureTreeNode->getLink(0);
    centerU = link->u;
    centerV = link->v;
	centerT = link->t;

    linkCount = figureTreeNode->getLinkCount();
    for(i = 1; i < linkCount; i++)
    {
        link = figureTreeNode->getLink(i);
        if(link == NULL)
            continue;


		if (centerT <= -1 || centerT >= 1) // center is not in the crest region
		{
			if (link->t == -1 && centerT <= -1 || link->t == 1 && centerT >= 1)
			{
				// both points are on the same side but not in the crest region
				centerU = (centerU * i / (i+1)) + link->u * 1 / (i + 1);
				centerV = (centerV * i / (i+1)) + link->v * 1 / (i + 1);
			}

			else // current point in in the crest region
			{
				if (link->u < maxU && link->u > 0) // not in the crest region of U
					centerU = (centerU * i / (i+1)) + link->u * 1 / (i + 1);
				else // in the crest region of U
				{
					if (centerT == 1)
					{
						if (link->u == 0)
						{
							tempU = centerU * i / (i+1) + (link->t - 1) * 1 / (i + 1);
							if (tempU < 0){
								centerU = 0;
								tempT1 = tempU + 1;
							}
							else
								centerU = tempU;
						}
						else // link->u == maxU
						{
							tempU = centerU * i / (i+1) + (maxU + 1 - link->t) * 1 / (i + 1);
							if (tempU > maxU){
								centerU = maxU;
								tempT1 = maxU - tempU + 1;
							}
							else
								centerU = tempU;
						}
					}
					else // centerT == -1
					{
						if (link->u == 0)
						{
							tempU = centerU * i / (i+1) - (link->t + 1) * 1 / (i + 1);
							if (tempU < 0){
								centerU = 0;
								tempT2 = -tempU - 1;
							}
							else
								centerU = tempU;
						}
						else // link->u == maxU
						{
							tempU = centerU * i / (i+1) + ( maxU + link->t + 1) * 1 / (i + 1);
							if (tempU > maxU){
								centerU = maxU;
								tempT1 = -maxU + tempU - 1;
							}
							else
								centerU = tempU;
						}
					}
				}
				if (link->v < maxV && link->v > 0) // not in the crest region of V
					centerV = (centerV * i / (i + 1)) + link->v * 1 / (i + 1);
				else // in the crest region of V
				{
					if (centerT == 1)
					{
						if (link->v == 0)
						{
							tempV = centerV * i / (i + 1) + (link->t - 1) * 1 / (i + 1);
							if (tempV < 0){
								centerV = 0;
								tempT2 = tempV + 1;
							}
							else
								centerV = tempV;
						}
						else // link->v == maxV
						{
							tempV = centerV * i / (i + 1) + (maxV + 1 - link->t) * 1 / (i + 1);
							if (tempV > maxV){
								centerV = maxV;
								tempT1 = maxV - tempV + 1;
							}
							else
								centerV = tempV;
						}
					}
					else // centerT == -1
					{
						if (link->v == 0)
						{
							tempV = centerV * i / (i + 1) - (link->t + 1) * 1 / (i + 1);
							if (tempV < 0){
								centerV = 0;
								tempT2 = -tempV - 1;
							}
							else
								centerV = tempV;
						}
						else // link->v == maxV
						{
							tempV = centerV * i / (i + 1) + ( maxV + link->t + 1) * 1 / (i + 1);
							if (tempV > maxV){
								centerV = maxV;
								tempT1 = -maxV + tempV - 1;
							}
							else
								centerV = tempV;
						}
					}
				}
				if ((link->u >= maxU || link->u <= 0) && (link->v >= maxV || link->v <= 0)) 
					centerT = (tempT1 + tempT2)/2;
			}
		}
		else // the center is in the crest region
		{
			if (fabs(link->t) >= 1) // current point is not in the crest region 
			{
				if (centerU < maxU && centerU > 0) // not in the crest region of U
					centerU = (centerU * i / (i + 1)) + link->u * 1 / (i + 1);
				else // in the crest region of U
				{
					if (link->t == 1)
					{
						if (centerU == 0)
						{
							tempU = (centerT - 1) * i / (i + 1) + (link->u + link->t - 1) * 1 / (i + 1);
							if (tempU < 0){
								centerU = 0;
								tempT1 = tempU + 1;
							}
							else
								centerU = tempU;
						}
						else // CenterU == maxU
						{
							tempU = (centerU + 1 - centerT) * i / (i + 1) + (link->u + 1 - link->t) * 1 / (i + 1);
							if (tempU > maxU){
								centerU = maxU;
								tempT1 = maxU - tempU + 1;
							}
							else
								centerU = tempU;
						}
					}
					else // link->T == -1
					{
						if (centerU == 0)
						{
							tempU = -(centerT + 1) * i / (i + 1) + (link->u - link->t - 1) * 1 / (i + 1);
							if (tempU < 0){
								centerU = 0;
								tempT1 = -tempU - 1;
							}
							else
								centerU = tempU;
						}
						else // centerU == maxU
						{
							tempU = (centerU + centerT + 1) * i / (i + 1) + ( link->u + link->t + 1) * 1 / (i + 1);
							if (tempU > maxU){
								centerU = maxU;
								tempT1 = -maxU + tempU - 1;
							}
							else
								centerU = tempU;
						}
					}
				}
				if (centerV < maxV && centerV > 0) // not in the crest region of V
					centerV = (centerV * i / (i + 1)) + link->v * 1 / (i + 1);
				else // in the crest region of V
				{
					if (link->t == 1)
					{
						if (centerV == 0)
						{
							tempV = (centerT - 1) * i / (i + 1) + (link->v + link->t - 1) * 1 / (i + 1);
							if (tempV < 0){
								centerV = 0;
								tempT2 = tempV + 1;
							}
							else
								centerV = tempV;
						}
						else // centerV == maxV
						{
							tempV = (centerV + 1 - centerT) * i / (i + 1) + (link->v + 1 - link->t) * 1 / (i + 1);
							if (tempV > maxV){
								centerV = maxV;
								tempT2 = maxV - tempV + 1;
							}
							else
								centerV = tempV;
						}
					}
					else // link->t == -1
					{
						if (centerV == 0)
						{
							tempV = -(centerT + 1) * i / (i + 1) + (link->v - link->t - 1) * 1 / (i + 1);
							if (tempV < 0){
								centerV = 0;
								tempT2 = -tempV - 1;
							}
							else
								centerV = tempV;
						}
						else // centerV == maxV
						{
							tempV = (centerV + centerT + 1) * i / (i + 1) + ( link->v + link->t + 1) * 1 / (i + 1);
							if (tempV > maxV){
								centerV = maxV;
								tempT2 = -maxV + tempV - 1;
							}
							else
								centerV = tempV;
						}
					}
				}
				if ((centerU >= maxU || centerU <= 0) && (centerV >= maxV || centerV <= 0)) 
					centerT = (tempT1 + tempT2)/2;
			}


			else  
			{
				if (fabs(link->t) < 1) // current point is also in the crest region
				{
					if ((centerU > 0 && centerU < maxU) && (link->v > 0 && link->v < maxV))   
					{
						// center in the V crest and current point in the U crest
						if (centerV == 0)
						{
							if (link->u == 0)
							{
								tempUV = centerU * i / (i + 1) - link->v * 1 / (i + 1);
								if (tempUV > 0)
									centerU = tempUV;
								else
								{
									centerU = link->u;
									centerV = -tempUV;
								}
							}
							else // link->u = maxU
							{
								tempUV = centerU * i / (i + 1) + (maxV+link->v) * 1 / (i + 1);
								if (tempUV < maxU)
									centerU = tempUV;
								else
								{
									centerU = link->u;
									centerV = -tempUV;
								}
							}
						}
						else // centerV == maxV)
						{
							if (link->u == 0)
							{
								tempUV = centerU * i / (i + 1) - (maxV -link->v) * 1 / (i + 1);
								if (tempUV > 0)
									centerU = tempUV;
								else
								{
									centerU = link->u;
									centerV = maxV + tempUV;
								}
							}
							else // link->u = maxU
							{
								tempUV = centerU * i / (i + 1) + (maxU+maxV-link->v) * 1 / (i + 1);
								if (tempUV < maxU)
									centerU = tempUV;
								else
								{
									centerU = link->u;
									centerV = maxU + maxV - tempUV;
								}
							}
						}
					}

					else if ((centerV > 0 && centerV < maxV) && (link->u > 0 && link->u < maxU))   
					{
						// center in the U crest and current point in the V crest
						if (centerU == 0)
						{
							if (link->v == 0)
							{
								tempUV = centerV * i / (i + 1) - link->u * 1 / (i + 1);
								if (tempUV > 0)
									centerV = tempUV;									
								else
								{
									centerU = -tempUV;
									centerV = link->v;
								}
							}
							else // link->v = maxV
							{
								tempUV = centerV * i / (i + 1) + (maxV+link->u) * 1 / (i + 1);
								if (tempUV < maxV)
									centerV = tempUV;
								else
								{
									centerU = maxV-tempUV;
									centerV = link->v;
								}
							}
						}
						else // centerU == maxU)
						{
							if (link->v == 0)
							{
								tempUV = centerV * i / (i + 1) - (maxU-link->u) * 1 / (i + 1);
								if (tempUV > 0)
									centerV = tempUV;
								else
								{
									centerU = -tempUV;
									centerV = link->v;
								}
							}
							else // link->v = maxV
							{
								tempUV = centerV * i / (i + 1) + (maxV+maxU-link->u) * 1 / (i + 1);
								if (tempUV < maxV)
									centerV = tempUV;
								else
								{
									centerU = maxU - tempUV + maxV;
									centerV = link->v;
								}
							}
						}
						centerT = centerT * i / (i + 1) + link->t * 1 / (i + 1);

					}
					else // both points are in the same crest region 
					{
						centerU = (centerU * i / (i+1)) + link->u * 1 / (i + 1);
						centerV = (centerV * i / (i+1)) + link->v * 1 / (i + 1);
						centerT = centerT * i / (i + 1) + link->t * 1 / (i + 1);
					}
				}			
			}			        
		}
	}

	//theta = 3.14159 / 6;
    cosTheta = cos(theta);
    sinTheta = sin(theta);

//	getInterpolatedPrimitive(centerPrimitive, parentFigure, centerU, centerV, t);

    int middleIndex = linkCount / 2;
    int deltaIndex;

    int startIndex = middleIndex;
    for(deltaIndex = -1; deltaIndex <= 1; deltaIndex += 2)
    {
        for(i = startIndex; i < linkCount && i >= 0; i += deltaIndex)
        {
            link = figureTreeNode->getLink(i);
            if(link == NULL)
                continue;
			if (centerT <= -1 || centerT >= 1) // center is not in the crest region
			{
				if (link->t == -1 && centerT <= -1 || link->t == 1 && centerT >= 1)
				{
					// both points are on the same side but not in the crest region
					uDiff = link->u - centerU;
					vDiff = link->v - centerV;
				}			
				else // current point is in the crest region
				{
					if (link->u < maxU && link->u > 0) // not in the crest region of U
						uDiff = link->u - centerU;
					else // in the crest region of U
					{
						if (centerT == 1)
						{
							if (link->u == 0)							
								uDiff = -centerU+link->t-1;
							else // link->u == maxU
								uDiff = maxU-centerU+1-link->t;
						}
						else // centerT == -1
						{
							if (link->u == 0)							
								uDiff = -centerU-link->t-1;
							else // link->u == maxU
								uDiff = maxU-centerU+1+link->t;
						}
					}
					if (link->v < maxV && link->v > 0) // not in the crest region of V
						vDiff = link->v - centerV;
					else // in the crest region of V
					{
						if (centerT == 1)
						{
							if (link->v == 0)							
								vDiff = -centerV+link->t-1;
							else // link->u == maxU
								vDiff = maxV-centerV+1-link->t;
						}
						else // centerT == -1
						{
							if (link->v == 0)							
								vDiff = -centerV-link->t-1;
							else // link->v == maxV
								vDiff = maxV-centerV+1+link->t;
						}
					}
				}
			}
			else // (centerT > -1 || centerT < 1) // center is in the crest region
			{
				if (link->t == -1  || link->t == 1)
				{
					// current point is not in the crest region
					if (centerU < maxU && centerU > 0) // not in the crest region of U
						uDiff = link->u - centerU;
					else // in the crest region of U
					{
						if (link->t == 1)
						{
							if (centerU == 0)							
								uDiff = centerU-link->t+1;
							else // centerU == maxU
								uDiff = -maxU+centerU-1+link->t;
						}
						else // link->t == -1
						{
							if (centerU == 0)							
								uDiff = centerU+link->t+1;
							else // centerU == maxU
								uDiff = -maxU+centerU-1-link->t;
						}
					}
					if (link->v < maxV && link->v > 0) // not in the crest region of V
						vDiff = link->v - centerV;
					else // in the crest region of V
					{
						if (link->t == 1)
						{
							if (centerV == 0)							
								vDiff = centerV-link->t+1;
							else // centerV == maxV
								vDiff = -maxV+centerV-1+link->t;
						}
						else // link->t == -1
						{
							if (centerV == 0)							
								vDiff = centerV+link->t+1;
							else // centerV == maxV
								vDiff = -maxV+centerV-1-link->t;
						}
					}
				}
				else // (link->t == -1  || link->t == 1)
				{
					// current point is in the crest region, thus both pts in the region
					if (centerU == 0 && link->u == 0) // region A
					{
						uDiff = (link->t - centerT) * alpha;
						/*if (centerT * link->t > 0 || fabs(centerT) > fabs(link->t))						
							uDiff = -fabs(link->t-centerT);						
						else						
							uDiff = fabs(link->t-centerT);*/
						vDiff = link->v - centerV;
					}
					else if (centerU == maxU && link->u == maxU) // region B
					{
						uDiff = (centerT - link->t) * alpha;
						/*if (centerT * link->t > 0 || fabs(centerT) > fabs(link->t))						
							uDiff = fabs(link->t-centerT);						
						else						
							uDiff = -fabs(link->t-centerT);*/
						vDiff = link->v - centerV;
					}
					else if (centerV == 0 && link->v == 0) // region C
					{
						vDiff = (link->t - centerT) * alpha; 
						/*if (centerT * link->t > 0 || fabs(centerT) > fabs(link->t))						
							vDiff = -fabs(link->t-centerT);						
						else						
							vDiff = fabs(link->t-centerT);*/
						uDiff = link->u - centerU;
					}
					else // if (centerV == maxV && link->v == maxV) // region D
					{
						vDiff = (centerT - link->t) * alpha;
						/*if (centerT * link->t > 0 || fabs(centerT) > fabs(link->t))						
							vDiff = fabs(link->t-centerT);						
						else						
							vDiff = -fabs(link->t-centerT);*/
						uDiff = link->u - centerU;
					}
				}
			}					

            uRot = cosTheta * uDiff - sinTheta * vDiff;
            vRot = sinTheta * uDiff + cosTheta * vDiff;

			if (link->t >-1 ){
				deltaU=uRot-uDiff;
				deltaV=vRot-vDiff;
			}
			else{
				deltaU=uDiff-uRot;
				deltaV=vDiff-vRot;
			}
			u=link->u; v=link->v; t=link->t;
			if (link->t == -1 || link->t == 1){ 
			// Neither u nor v is in the crest region.
				if (link->u + deltaU <= maxU && v + deltaV <= maxV && 
						link->u + deltaU >= 0 && v + deltaV >= 0){
					//Neither u+deltaU nor v+deltaV is in the crest region.
					u=link->u+deltaU; v=link->v+deltaV;
				}
				else if (link->u + deltaU > maxU){
					// u+deltaU is in the crest region.
					deltaU1=maxU-link->u;
					deltaU2=deltaU-deltaU1;
					u=maxU;
					//u=link->u+deltaU1; 
					v=link->v+deltaV;
					deltaT=deltaU2/alpha;
					//vDelta=link->v/maxV*(maxV+2*alpha)-link->v-alpha;
					//deltaT=deltaU2/sqrt(alpha*alpha+vDelta*vDelta);
					if (link->t == 1)
						t=link->t-fabs(deltaT);
					else
						t=link->t+fabs(deltaT);
				}
				else if (link->u + deltaU < 0){
					// u+deltaU is in the crest region.
					deltaU1=0-link->u;
					deltaU2=deltaU-deltaU1;
					u=0;
					//u=link->u+deltaU1; 
					v=link->v+deltaV;
					deltaT=deltaU2/alpha;
					//vDelta=link->v/maxV*(maxV+2*alpha)-link->v-alpha;
					//deltaT=deltaU2/sqrt(alpha*alpha+vDelta*vDelta);
					if (link->t == 1)
						t=link->t-fabs(deltaT);
					else
						t=link->t+fabs(deltaT);
				}
				else if ( link->v + deltaV > maxV){ 
					//link->v+deltaV is in the crest region.
					deltaV1=maxV-link->v;
					deltaV2=deltaV-deltaV1;
					u=link->u+deltaU; 
					v=maxV;
					//v=link->v+deltaV1;
					deltaT=deltaV2/alpha;
					//uDelta=link->u/maxU*(maxU+2*alpha)-link->u-alpha;
					//deltaT=deltaV2/sqrt(alpha*alpha+uDelta*uDelta);
					if (link->t == 1)
						t=link->t-fabs(deltaT);
					else
						t=link->t+fabs(deltaT); 
				}
				else {
					//link->v + deltaV < 0) 
					//link->v+deltaV is in the crest region.
					deltaV1=0-link->v;
					deltaV2=deltaV-deltaV1;
					u=link->u+deltaU; 
					v=0;
					//v=link->v+deltaV1;
					deltaT=deltaV2/alpha;
					//uDelta=link->u/maxU*(maxU+2*alpha)-link->u-alpha;
					//deltaT=deltaV2/sqrt(alpha*alpha+uDelta*uDelta);
					if (link->t == 1)
						t=link->t-fabs(deltaT);
					else
						t=link->t+fabs(deltaT); 
				}
			}
			else if ((link->u == 0) || (link->u == maxU)){
				// u is in the crest region.				
				if (deltaU != 0){ 
					if (link->u == maxU)
						deltaT=-deltaU/alpha;
						//deltaT=-deltaU/sqrt(alpha*alpha+vDelta*vDelta);
					else 
						deltaT=deltaU/alpha;
						//deltaT=deltaU/sqrt(alpha*alpha+vDelta*vDelta);
					if (link->t+deltaT >= -1 && link->t+deltaT <= 1){
						// still ends in the crest region.
						t=link->t+deltaT; 
					}
					else{
						// ends in the non-crest region.
						if (link->t+deltaT < -1){
							deltaT=-1-link->t;
							t=-1;
							deltaU1=deltaT*alpha;
							//deltaU1=deltaT*sqrt(alpha*alpha+vDelta*vDelta);
							deltaU2=deltaU-deltaU1;
							if (link->u == 0)
								u=link->u+fabs(deltaU2);
							else
								u=link->u-fabs(deltaU2);
						}
						else{
							deltaT=1-link->t;
							t=1;
							deltaU1=deltaT*alpha;
							//deltaU1=deltaT*sqrt(alpha*alpha+vDelta*vDelta);
							deltaU2=deltaU-deltaU1;
							if (link->u == 0)
								u=link->u+fabs(deltaU2);
							else
								u=link->u-fabs(deltaU2);						
						}
					}
				}
				if (deltaV != 0) { // moving in v direction.
					if (link->v+deltaV >= 0 && link->v+deltaV <= maxV){
						// still in the range of V, normal movement
						v=link->v+deltaV;
					}
					else{
						if (link->v+deltaV < 0){
							deltaV1=0-link->v;
							v=0;
							deltaU2=deltaV-deltaV1;
							if ( link->u == 0)
								u=fabs(deltaU2);
							else
								u=link->u-fabs(deltaU2);
						}
						else{
							deltaV1=maxV-link->v;
							v=maxV;
							deltaU2=deltaV-deltaV1;
							if ( link->u == 0)
								u=fabs(deltaU2);
							else
								u=link->u-fabs(deltaU2);
						}
					}
				}
			}
			else{
				// v is in the crest region.
				if (deltaV != 0){ // moving in v direction.
					//uDelta=link->u/maxU*(maxU+2*alpha)-link->u-alpha;
					//deltaT=deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					//if (link->t >= 0)
					if (link->v == maxV)
						deltaT=-deltaV/alpha;
						//deltaT=-deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					else 
						deltaT=deltaV/alpha;
						//deltaT=deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					if (link->t+deltaT >= -1 && link->t+deltaT <= 1){
						// still ends in the crest region.
						t=link->t+deltaT;
					}
					else{
						// ends in the non-crest region.
						if (link->t+deltaT < -1){
							deltaT=-1-link->t;
							t=-1;
							deltaV1=deltaT*alpha;
							//deltaV1=deltaT*sqrt(alpha*alpha+uDelta*uDelta);
							deltaV2=deltaV-deltaV1;
							if (link->v == 0)
								v=link->v+fabs(deltaV2);
							else
								v=link->v-fabs(deltaV2);
						}
						else{
							deltaT=1-link->t;
							t=1;
							deltaV1=deltaT*alpha;
							//deltaV1=deltaT*sqrt(alpha*alpha+uDelta*uDelta);
							deltaV2=deltaV-deltaV1;
							if (link->v == 0)
							v=link->v+fabs(deltaV2);
							else
								v=link->v-fabs(deltaV2);						
						}
					}
				}
				if (deltaU != 0) { // moving in u direction.
					if (fabs(link->u+deltaU) <= maxU){
						// still in the range of U, normal movement
						u=link->u+deltaU;
					}
					else{
						if (link->u+deltaU < 0){
							deltaU1=0-link->u;
							u=0;
							deltaV2=deltaU-deltaU1;
							if ( link->v == 0)
								v=fabs(deltaV2);
							else
								v=link->v-fabs(deltaV2);
						}
						else{
							deltaU1=maxU-link->u;
							u=maxU;
							deltaV2=deltaU-deltaU1;
							if ( link->v == 0)
								v=fabs(deltaU2);
							else
								v=link->v-fabs(deltaV2);
						}
					}
				}
			}



            getSurfacePointAndNormal(pList2, startPosition, oldNormal,
                                     link->u, link->v, link->t);
            getSurfacePointAndNormal(pList2, endPosition, newNormal,
                                     u, v, t);
			/*
            getSurfacePointAndNormal(startPosition, oldNormal, parentFigure,
                                     link->u, link->v, link->t);
            getSurfacePointAndNormal(endPosition, newNormal, parentFigure,
                                     u, v, t);
			*/

            link->u = u;
            link->v = v;
			link->t = t;

            primitiveId = link->primitiveId;
            primitive = figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;

            deltaX = endPosition - startPosition;
            axis = oldNormal.cross(newNormal);

			angle=AngleBetweenVectors(oldNormal, newNormal);
            //angle = acos(oldNormal * newNormal);

			deltaQ.setAxisAngle(axis, angle);

            // Find the primitive rotation due to position change
            Vector3D initDiff;
            Vector3D finalDiff;
            Vector3D centerPos;
            Quat primQ;

            primQ.setAxisAngle(newNormal, -theta);

            // Make sure rotation is not zero
            bool doRotations;
            bool doPrimRotations;
            Real norm = deltaQ.norm();
            if(norm > 1.0 - R_SMALL_TOLERANCE &&
               norm < 1.0 + R_SMALL_TOLERANCE)
                doRotations = true;
            else
                doRotations = false;

            norm = primQ.norm();
            if(norm > 1.0 - R_SMALL_TOLERANCE &&
               norm < 1.0 + R_SMALL_TOLERANCE)
                doPrimRotations = true;
            else
                doPrimRotations = false;

            primitive = figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;

            if(doRotations)
            {
                primitive->rotateBy(deltaQ);

                axis = primitive->getX() - startPosition;
                deltaQ.rotateVector(axis);
                primitive->setX(axis + startPosition);
            }

            if(doPrimRotations)
                primitive->rotateBy(primQ);

            primitive->translateBy(deltaX);

            if(i == middleIndex)
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;

                    setPrimitiveToPreviousPrediction(figure, primitiveId,
                        prevPrimitiveId, lastElongateValue);
                }
            }

            else
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;
                    neighborPrimitiveId = primitiveId -
                        deltaIndex * neighborIndexDelta;

                    setPrimitiveToAveragePrediction(figure, primitiveId,
                        prevPrimitiveId, neighborPrimitiveId, lastScaleValue,
                        lastElongateValue);
                }
            }
        }

        startIndex = middleIndex + 1;
    }

	delete pList2;
}

#define TWO_PI 2*M_PI
void M3DSubfigureTransformation::rotate(double theta)
{
	int i;
	double delta;
	double lastTheta;

	delta=.5;	//delta = 1.0;

	while(theta>TWO_PI)
		theta-=TWO_PI;
	while(theta<-TWO_PI)
		theta+=TWO_PI;

#ifdef DEBUG
	fprintf(stderr, "dTheta: %f degree(s) in ACTUAL rotation.\n", theta/M_PI*180);
#endif

	int thetaIterations = (int) fabs(theta / delta);
	if (theta > 0.0)
		lastTheta = theta - delta * thetaIterations;
	else 
		lastTheta = theta + delta * thetaIterations;

	for ( i=0; i < thetaIterations; i++)
	{
		//if (thetaIterations > 0.0) 
		if (theta > 0.0)
			rotateStep(delta);
		else 
			rotateStep(-delta);
		initStep();
	}
	rotateStep(lastTheta);
	//rotateStep(theta);
	reInit();
}

void M3DSubfigureTransformation::widenStep(double value)
{
    M3DFigureTreeNode * parentNode;
    M3DQuadFigure * figure;
    M3DQuadFigure * parentFigure;
    M3DPrimitive * primitive;
    M3DQuadPrimitive startPrimitive;
    M3DQuadPrimitive finishPrimitive;
    M3DQuadPrimitive centerPrimitive;
    Vector3D startPosition, endPosition;
    Vector3D oldNormal, newNormal;
    Vector3D deltaX;
    Vector3D axis;
    Quat deltaQ;
    M3DPrimitiveLinkInfo * link;


	double maxU, maxV, angle;
    double centerU, centerV, centerT, tempU, tempV, tempUV, tempT1, tempT2;
    double u, v, t, uDiff, vDiff;

	double  deltaU, deltaV, deltaU1, deltaU2, deltaV1, deltaV2, deltaT;
	double alpha=1.0;

    int linkCount,
        rowCount, columnCount,
        i, j,
        primitiveId,
        primitiveIndexDelta,
        neighborIndexDelta,
        numPrimitivesPerLink;

    int prevPrimitiveId,
        neighborPrimitiveId;


    if(figureTreeNode == NULL || object == NULL)
        return;

    figure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(figureTreeNode->getFigureId()));
    if(figure == NULL)
        return;

    parentNode = figureTreeNode->getParent();
    if(parentNode == NULL)
        return;

    parentFigure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(parentNode->getFigureId()));

    rowCount = figure->getRowCount();
    columnCount = figure->getColumnCount();

	maxU = (double) (parentFigure->getColumnCount() - 1);
    maxV = (double) (parentFigure->getRowCount() - 1);

    //maxU = (double) (parentFigure->getRowCount() - 1);
    //maxV = (double) (parentFigure->getColumnCount() - 1);

    if(rowOrientedLinks)
    {
        primitiveIndexDelta = columnCount;
        numPrimitivesPerLink = rowCount;
        neighborIndexDelta = 1;
    }
    else
    {
        primitiveIndexDelta = 1;
        numPrimitivesPerLink = columnCount;
        neighborIndexDelta = columnCount;
    }

	//		integrate the '(u, v, t)->BPosition and BNormal' into this function
	ThallCode::Pointlist_server2 *pList2=new ThallCode::Pointlist_server2();
	Xferlist *xferList = convertM3DtoXfer(parentFigure);
	pList2->InitializeSubdivSurf(xferList);
	delete []xferList->atomlist;
	delete xferList;

	link = figureTreeNode->getLink(0);
    centerU = link->u;
    centerV = link->v;
	centerT = link->t;

    linkCount = figureTreeNode->getLinkCount();
    for(i = 1; i < linkCount; i++)
    {
        link = figureTreeNode->getLink(i);
        if(link == NULL)
            continue;


		if (centerT <= -1 || centerT >= 1) // center is not in the crest region
		{
			if (link->t == -1 && centerT <= -1 || link->t == 1 && centerT >= 1)
			{
				// both points are on the same side but not in the crest region
				centerU = (centerU * i / (i+1)) + link->u * 1 / (i + 1);
				centerV = (centerV * i / (i+1)) + link->v * 1 / (i + 1);
			}

			else // current point in in the crest region
			{
				if (link->u < maxU && link->u > 0) // not in the crest region of U
					centerU = (centerU * i / (i+1)) + link->u * 1 / (i + 1);
				else // in the crest region of U
				{
					if (centerT == 1)
					{
						if (link->u == 0)
						{
							tempU = centerU * i / (i+1) + (link->t - 1) * 1 / (i + 1);
							if (tempU < 0){
								centerU = 0;
								tempT1 = tempU + 1;
							}
							else
								centerU = tempU;
						}
						else // link->u == maxU
						{
							tempU = centerU * i / (i+1) + (maxU + 1 - link->t) * 1 / (i + 1);
							if (tempU > maxU){
								centerU = maxU;
								tempT1 = maxU - tempU + 1;
							}
							else
								centerU = tempU;
						}
					}
					else // centerT == -1
					{
						if (link->u == 0)
						{
							tempU = centerU * i / (i+1) - (link->t + 1) * 1 / (i + 1);
							if (tempU < 0){
								centerU = 0;
								tempT2 = -tempU - 1;
							}
							else
								centerU = tempU;
						}
						else // link->u == maxU
						{
							tempU = centerU * i / (i+1) + ( maxU + link->t + 1) * 1 / (i + 1);
							if (tempU > maxU){
								centerU = maxU;
								tempT1 = -maxU + tempU - 1;
							}
							else
								centerU = tempU;
						}
					}
				}
				if (link->v < maxV && link->v > 0) // not in the crest region of V
					centerV = (centerV * i / (i + 1)) + link->v * 1 / (i + 1);
				else // in the crest region of V
				{
					if (centerT == 1)
					{
						if (link->v == 0)
						{
							tempV = centerV * i / (i + 1) + (link->t - 1) * 1 / (i + 1);
							if (tempV < 0){
								centerV = 0;
								tempT2 = tempV + 1;
							}
							else
								centerV = tempV;
						}
						else // link->v == maxV
						{
							tempV = centerV * i / (i + 1) + (maxV + 1 - link->t) * 1 / (i + 1);
							if (tempV > maxV){
								centerV = maxV;
								tempT1 = maxV - tempV + 1;
							}
							else
								centerV = tempV;
						}
					}
					else // centerT == -1
					{
						if (link->v == 0)
						{
							tempV = centerV * i / (i + 1) - (link->t + 1) * 1 / (i + 1);
							if (tempV < 0){
								centerV = 0;
								tempT2 = -tempV - 1;
							}
							else
								centerV = tempV;
						}
						else // link->v == maxV
						{
							tempV = centerV * i / (i + 1) + ( maxV + link->t + 1) * 1 / (i + 1);
							if (tempV > maxV){
								centerV = maxV;
								tempT1 = -maxV + tempV - 1;
							}
							else
								centerV = tempV;
						}
					}
				}
				if ((link->u >= maxU || link->u <= 0) && (link->v >= maxV || link->v <= 0)) 
					centerT = (tempT1 + tempT2)/2;
			}
		}
		else // the center is in the crest region
		{
			if (fabs(link->t) >= 1) // current point is not in the crest region 
			{
				if (centerU < maxU && centerU > 0) // not in the crest region of U
					centerU = (centerU * i / (i + 1)) + link->u * 1 / (i + 1);
				else // in the crest region of U
				{
					if (link->t == 1)
					{
						if (centerU == 0)
						{
							tempU = (centerT - 1) * i / (i + 1) + (link->u + link->t - 1) * 1 / (i + 1);
							if (tempU < 0){
								centerU = 0;
								tempT1 = tempU + 1;
							}
							else
								centerU = tempU;
						}
						else // CenterU == maxU
						{
							tempU = (centerU + 1 - centerT) * i / (i + 1) + (link->u + 1 - link->t) * 1 / (i + 1);
							if (tempU > maxU){
								centerU = maxU;
								tempT1 = maxU - tempU + 1;
							}
							else
								centerU = tempU;
						}
					}
					else // link->T == -1
					{
						if (centerU == 0)
						{
							tempU = -(centerT + 1) * i / (i + 1) + (link->u - link->t - 1) * 1 / (i + 1);
							if (tempU < 0){
								centerU = 0;
								tempT1 = -tempU - 1;
							}
							else
								centerU = tempU;
						}
						else // centerU == maxU
						{
							tempU = (centerU + centerT + 1) * i / (i + 1) + ( link->u + link->t + 1) * 1 / (i + 1);
							if (tempU > maxU){
								centerU = maxU;
								tempT1 = -maxU + tempU - 1;
							}
							else
								centerU = tempU;
						}
					}
				}
				if (centerV < maxV && centerV > 0) // not in the crest region of V
					centerV = (centerV * i / (i + 1)) + link->v * 1 / (i + 1);
				else // in the crest region of V
				{
					if (link->t == 1)
					{
						if (centerV == 0)
						{
							tempV = (centerT - 1) * i / (i + 1) + (link->v + link->t - 1) * 1 / (i + 1);
							if (tempV < 0){
								centerV = 0;
								tempT2 = tempV + 1;
							}
							else
								centerV = tempV;
						}
						else // centerV == maxV
						{
							tempV = (centerV + 1 - centerT) * i / (i + 1) + (link->v + 1 - link->t) * 1 / (i + 1);
							if (tempV > maxV){
								centerV = maxV;
								tempT2 = maxV - tempV + 1;
							}
							else
								centerV = tempV;
						}
					}
					else // link->t == -1
					{
						if (centerV == 0)
						{
							tempV = -(centerT + 1) * i / (i + 1) + (link->v - link->t - 1) * 1 / (i + 1);
							if (tempV < 0){
								centerV = 0;
								tempT2 = -tempV - 1;
							}
							else
								centerV = tempV;
						}
						else // centerV == maxV
						{
							tempV = (centerV + centerT + 1) * i / (i + 1) + ( link->v + link->t + 1) * 1 / (i + 1);
							if (tempV > maxV){
								centerV = maxV;
								tempT2 = -maxV + tempV - 1;
							}
							else
								centerV = tempV;
						}
					}
				}
				if ((centerU >= maxU || centerU <= 0) && (centerV >= maxV || centerV <= 0)) 
					centerT = (tempT1 + tempT2)/2;
			}


			else  
			{
				if (fabs(link->t) < 1) // current point is also in the crest region
				{
					if ((centerU > 0 && centerU < maxU) && (link->v > 0 && link->v < maxV))   
					{
						// center in the V crest and current point in the U crest
						if (centerV == 0)
						{
							if (link->u == 0)
							{
								tempUV = centerU * i / (i + 1) - link->v * 1 / (i + 1);
								if (tempUV > 0)
									centerU = tempUV;
								else
								{
									centerU = link->u;
									centerV = -tempUV;
								}
							}
							else // link->u = maxU
							{
								tempUV = centerU * i / (i + 1) + (maxV+link->v) * 1 / (i + 1);
								if (tempUV < maxU)
									centerU = tempUV;
								else
								{
									centerU = link->u;
									centerV = -tempUV;
								}
							}
						}
						else // centerV == maxV)
						{
							if (link->u == 0)
							{
								tempUV = centerU * i / (i + 1) - (maxV -link->v) * 1 / (i + 1);
								if (tempUV > 0)
									centerU = tempUV;
								else
								{
									centerU = link->u;
									centerV = maxV + tempUV;
								}
							}
							else // link->u = maxU
							{
								tempUV = centerU * i / (i + 1) + (maxU+maxV-link->v) * 1 / (i + 1);
								if (tempUV < maxU)
									centerU = tempUV;
								else
								{
									centerU = link->u;
									centerV = maxU + maxV - tempUV;
								}
							}
						}
					}

					else if ((centerV > 0 && centerV < maxV) && (link->u > 0 && link->u < maxU))   
					{
						// center in the U crest and current point in the V crest
						if (centerU == 0)
						{
							if (link->v == 0)
							{
								tempUV = centerV * i / (i + 1) - link->u * 1 / (i + 1);
								if (tempUV > 0)
									centerV = tempUV;									
								else
								{
									centerU = -tempUV;
									centerV = link->v;
								}
							}
							else // link->v = maxV
							{
								tempUV = centerV * i / (i + 1) + (maxV+link->u) * 1 / (i + 1);
								if (tempUV < maxV)
									centerV = tempUV;
								else
								{
									centerU = maxV-tempUV;
									centerV = link->v;
								}
							}
						}
						else // centerU == maxU)
						{
							if (link->v == 0)
							{
								tempUV = centerV * i / (i + 1) - (maxU-link->u) * 1 / (i + 1);
								if (tempUV > 0)
									centerV = tempUV;
								else
								{
									centerU = -tempUV;
									centerV = link->v;
								}
							}
							else // link->v = maxV
							{
								tempUV = centerV * i / (i + 1) + (maxV+maxU-link->u) * 1 / (i + 1);
								if (tempUV < maxV)
									centerV = tempUV;
								else
								{
									centerU = maxU - tempUV + maxV;
									centerV = link->v;
								}
							}
						}
						centerT = centerT * i / (i + 1) + link->t * 1 / (i + 1);

					}
					else // both points are in the same crest region 
					{
						centerU = (centerU * i / (i+1)) + link->u * 1 / (i + 1);
						centerV = (centerV * i / (i+1)) + link->v * 1 / (i + 1);
						centerT = centerT * i / (i + 1) + link->t * 1 / (i + 1);
					}
				}			
			}			        
		}
	}

	int middleIndex = linkCount / 2;
    int deltaIndex;

    int startIndex = middleIndex;
    for(deltaIndex = -1; deltaIndex <= 1; deltaIndex += 2)
    {
        for(i = startIndex; i < linkCount && i >= 0; i += deltaIndex)
        {
            link = figureTreeNode->getLink(i);
            if(link == NULL)
                continue;
			if (centerT <= -1 || centerT >= 1) // center is not in the crest region
			{
				if (link->t == -1 && centerT <= -1 || link->t == 1 && centerT >= 1)
				{
					// both points are on the same side but not in the crest region
					uDiff = link->u - centerU;
					vDiff = link->v - centerV;
				}			
				else // current point is in the crest region
				{
					if (link->u < maxU && link->u > 0) // not in the crest region of U
						uDiff = link->u - centerU;
					else // in the crest region of U
					{
						if (centerT == 1)
						{
							if (link->u == 0)							
								uDiff = -centerU+link->t-1;
							else // link->u == maxU
								uDiff = maxU-centerU+1-link->t;
						}
						else // centerT == -1
						{
							if (link->u == 0)							
								uDiff = -centerU-link->t-1;
							else // link->u == maxU
								uDiff = maxU-centerU+1+link->t;
						}
					}
					if (link->v < maxV && link->v > 0) // not in the crest region of V
						vDiff = link->v - centerV;
					else // in the crest region of V
					{
						if (centerT == 1)
						{
							if (link->v == 0)							
								vDiff = -centerV+link->t-1;
							else // link->u == maxU
								vDiff = maxV-centerV+1-link->t;
						}
						else // centerT == -1
						{
							if (link->v == 0)							
								vDiff = -centerV-link->t-1;
							else // link->v == maxV
								vDiff = maxV-centerV+1+link->t;
						}
					}
				}
			}
			else // (centerT > -1 || centerT < 1) // center is in the crest region
			{
				if (link->t == -1  || link->t == 1)
				{
					// current point is not in the crest region
					if (centerU < maxU && centerU > 0) // not in the crest region of U
						uDiff = link->u - centerU;
					else // in the crest region of U
					{
						if (link->t == 1)
						{
							if (centerU == 0)							
								uDiff = centerU-link->t+1;
							else // centerU == maxU
								uDiff = -maxU+centerU-1+link->t;
						}
						else // link->t == -1
						{
							if (centerU == 0)							
								uDiff = centerU+link->t+1;
							else // centerU == maxU
								uDiff = -maxU+centerU-1-link->t;
						}
					}
					if (link->v < maxV && link->v > 0) // not in the crest region of V
						vDiff = link->v - centerV;
					else // in the crest region of V
					{
						if (link->t == 1)
						{
							if (centerV == 0)							
								vDiff = centerV-link->t+1;
							else // centerV == maxV
								vDiff = -maxV+centerV-1+link->t;
						}
						else // link->t == -1
						{
							if (centerV == 0)							
								vDiff = centerV+link->t+1;
							else // centerV == maxV
								vDiff = -maxV+centerV-1-link->t;
						}
					}
				}
				else // (link->t == -1  || link->t == 1)
				{
					// current point is in the crest region, thus both pts in the region
					if (centerU == 0 && link->u == 0) // region A
					{
						uDiff = (link->t - centerT) * alpha;
						/*if (centerT * link->t > 0 || fabs(centerT) > fabs(link->t))						
							uDiff = -fabs(link->t-centerT);						
						else						
							uDiff = fabs(link->t-centerT);*/
						vDiff = link->v - centerV;
					}
					else if (centerU == maxU && link->u == maxU) // region B
					{
						uDiff = (centerT - link->t) * alpha;
						/*if (centerT * link->t > 0 || fabs(centerT) > fabs(link->t))						
							uDiff = fabs(link->t-centerT);						
						else						
							uDiff = -fabs(link->t-centerT);*/
						vDiff = link->v - centerV;
					}
					else if (centerV == 0 && link->v == 0) // region C
					{
						vDiff = (link->t - centerT) * alpha; 
						/*if (centerT * link->t > 0 || fabs(centerT) > fabs(link->t))						
							vDiff = -fabs(link->t-centerT);						
						else						
							vDiff = fabs(link->t-centerT);*/
						uDiff = link->v - centerV;
					}
					else // if (centerV == maxV && link->v == maxV) // region D
					{
						vDiff = (centerT - link->t) * alpha;
						/*if (centerT * link->t > 0 || fabs(centerT) > fabs(link->t))						
							vDiff = fabs(link->t-centerT);						
						else						
							vDiff = -fabs(link->t-centerT);*/
						uDiff = link->v - centerV;
					}
				}
			}

			if (centerT >-1 ){
				deltaU = uDiff * value;
				deltaV = vDiff * value;
			}
			else{
				deltaU= -uDiff * value;
				deltaV= vDiff * value;;
			}
			u=centerU; v=centerV; t=centerT;
			if (centerT == -1 || centerT == 1){ 
			// Neither u nor v is in the crest region.
				if (centerU + deltaU <= maxU && v + deltaV <= maxV && 
						centerU + deltaU >= 0 && v + deltaV >= 0){
					//Neither u+deltaU nor v+deltaV is in the crest region.
					u=centerU+deltaU; v=centerV+deltaV;
				}
				else if (centerU + deltaU > maxU){
					// u+deltaU is in the crest region.
					deltaU1=maxU-centerU;
					deltaU2=deltaU-deltaU1;
					u=centerU+deltaU1; v=centerV+deltaV;
					deltaT=deltaU2/alpha;
					//vDelta=centerV/maxV*(maxV+2*alpha)-centerV-alpha;
					//deltaT=deltaU2/sqrt(alpha*alpha+vDelta*vDelta);
					if (centerT == 1)
						t=centerT-fabs(deltaT);
					else
						t=centerT+fabs(deltaT);
				}
				else if (centerU + deltaU < 0){
					// u+deltaU is in the crest region.
					deltaU1=0-centerU;
					deltaU2=deltaU-deltaU1;
					u=centerU+deltaU1; v=centerV+deltaV;
					deltaT=deltaU2/alpha;
					//vDelta=centerV/maxV*(maxV+2*alpha)-centerV-alpha;
					//deltaT=deltaU2/sqrt(alpha*alpha+vDelta*vDelta);
					if (centerT == 1)
						t=centerT-fabs(deltaT);
					else
						t=centerT+fabs(deltaT);
				}
				else if ( centerV + deltaV > maxV){ 
					//centerV+deltaV is in the crest region.
					deltaV1=maxV-centerV;
					deltaV2=deltaV-deltaV1;
					u=centerU+deltaU; v=centerV+deltaV1;
					deltaT=deltaV2/alpha;
					//uDelta=centerU/maxU*(maxU+2*alpha)-centerU-alpha;
					//deltaT=deltaV2/sqrt(alpha*alpha+uDelta*uDelta);
					if (centerT == 1)
						t=centerT-fabs(deltaT);
					else
						t=centerT+fabs(deltaT); 
				}
				else {
					//centerV + deltaV < 0) 
					//centerV+deltaV is in the crest region.
					deltaV1=0-centerV;
					deltaV2=deltaV-deltaV1;
					u=centerU+deltaU; v=centerV+deltaV1;
					deltaT=deltaV2/alpha;
					//uDelta=centerU/maxU*(maxU+2*alpha)-centerU-alpha;
					//deltaT=deltaV2/sqrt(alpha*alpha+uDelta*uDelta);
					if (centerT == 1)
						t=centerT-fabs(deltaT);
					else
						t=centerT+fabs(deltaT); 
				}
			}
			else if ((centerU == 0) || (centerU == maxU)){
				// u is in the crest region.				
				if (deltaU != 0){ // moving in u direction.	
					if (centerU == maxU)
						deltaT=-deltaU/alpha;
						//deltaT=-deltaU/sqrt(alpha*alpha+vDelta*vDelta);
					else 
						deltaT=deltaU/alpha;
						//deltaT=deltaU/sqrt(alpha*alpha+vDelta*vDelta);
					if (centerT+deltaT >= -1 && centerT+deltaT <= 1){
						// still ends in the crest region.
						t=centerT+deltaT; 
					}
					else{
						// ends in the non-crest region.
						if (centerT+deltaT < -1){
							deltaT=-1-centerT;
							t=-1;
							deltaU1=deltaT*alpha;
							//deltaU1=deltaT*sqrt(alpha*alpha+vDelta*vDelta);
							deltaU2=deltaU-deltaU1;
							if (centerU == 0)
								u=centerU+fabs(deltaU2);
							else
								u=centerU-fabs(deltaU2);
						}
						else{
							deltaT=1-centerT;
							t=1;
							deltaU1=deltaT*alpha;
							//deltaU1=deltaT*sqrt(alpha*alpha+vDelta*vDelta);
							deltaU2=deltaU-deltaU1;
							if (centerU == 0)
								u=centerU+fabs(deltaU2);
							else
								u=centerU-fabs(deltaU2);						
						}
					}
				}
				if (deltaV != 0) { // moving in v direction.
					if (centerV+deltaV >= 0 && centerV+deltaV <= maxV){
						// still in the range of V, normal movement
						v=centerV+deltaV;
					}
					else{
						if (centerV+deltaV < 0){
							deltaV1=0-centerV;
							v=0;
							deltaU2=deltaV-deltaV1;
							if ( centerU == 0)
								u=fabs(deltaU2);
							else
								u=centerU-fabs(deltaU2);
						}
						else{
							deltaV1=maxV-centerV;
							v=maxV;
							deltaU2=deltaV-deltaV1;
							if ( centerU == 0)
								u=fabs(deltaU2);
							else
								u=centerU-fabs(deltaU2);
						}
					}
				}
			}
			else{
				// v is in the crest region.
				if (deltaV != 0){ // moving in v direction.
					if (centerV == maxV)
						deltaT=-deltaV/alpha;
						//deltaT=-deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					else 
						deltaT=deltaV/alpha;
						//deltaT=deltaV/sqrt(alpha*alpha+uDelta*uDelta);
					if (centerT+deltaT >= -1 && centerT+deltaT <= 1){
						// still ends in the crest region.
						t=centerT+deltaT;
					}
					else{
						// ends in the non-crest region.
						if (centerT+deltaT < -1){
							deltaT=-1-centerT;
							t=-1;
							deltaV1=deltaT*alpha;
							//deltaV1=deltaT*sqrt(alpha*alpha+uDelta*uDelta);
							deltaV2=deltaV-deltaV1;
							if (centerV == 0)
							v=centerV+fabs(deltaV2);
							else
								v=centerV-fabs(deltaV2);
						}
						else{
							deltaT=1-centerT;
							t=1;
							deltaV1=deltaT*alpha;
							//deltaV1=deltaT*sqrt(alpha*alpha+uDelta*uDelta);
							deltaV2=deltaV-deltaV1;
							if (centerV == 0)
							v=centerV+fabs(deltaV2);
							else
								v=centerV-fabs(deltaV2);						
						}
					}
				}
				if (deltaU != 0) { // moving in u direction.
					if (fabs(centerU+deltaU) <= maxU){
						// still in the range of U, normal movement
						u=centerU+deltaU;
					}
					else{
						if (centerU+deltaU < 0){
							deltaU1=0-centerU;
							u=0;
							deltaV2=deltaU-deltaU1;
							if ( centerV == 0)
								v=fabs(deltaV2);
							else
								v=centerV-fabs(deltaV2);
						}
						else{
							deltaU1=maxU-centerU;
							u=maxU;
							deltaV2=deltaU-deltaU1;
							if ( centerV == 0)
								v=fabs(deltaU2);
							else
								v=centerV-fabs(deltaV2);
						}
					}
				}
			}

			lastScaleValue = value; //* lastScaleValue;  

            getSurfacePointAndNormal(pList2, startPosition, oldNormal,
                                     link->u, link->v, link->t);

            getSurfacePointAndNormal(pList2, endPosition, newNormal,
                                     u, v, t);
			/*
            getSurfacePointAndNormal(startPosition, oldNormal, parentFigure,
                                     link->u, link->v, link->t);

            getSurfacePointAndNormal(endPosition, newNormal, parentFigure,
                                     u, v, t);
			*/

            link->u = u;
            link->v = v;
			link->t = t;

            primitiveId = link->primitiveId;
            primitive = figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;

            deltaX = endPosition - startPosition;
            axis = oldNormal.cross(newNormal);

			angle=AngleBetweenVectors(oldNormal, newNormal);
            //angle = acos(oldNormal * newNormal);

			deltaQ.setAxisAngle(axis, angle);

            // Make sure rotation is not zero
            bool doRotations;
            Real norm = deltaQ.norm();
            if(norm > 1.0 - R_SMALL_TOLERANCE &&
               norm < 1.0 + R_SMALL_TOLERANCE)
                doRotations = true;
            else
                doRotations = false;

            primitive = figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;

            if(doRotations)
            {
                primitive->rotateBy(deltaQ);

                axis = primitive->getX() - startPosition;
                deltaQ.rotateVector(axis);
                primitive->setX(axis + startPosition);
            }

            primitive->translateBy(deltaX);
			primitive->scaleBy(lastScaleValue);

            if(i == middleIndex)
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;

                    setPrimitiveToPreviousPrediction(figure, primitiveId,
                        prevPrimitiveId, lastElongateValue);

					figure->getPrimitivePtr(primitiveId)->scaleBy(lastScaleValue);
                }
            }

            else
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;
                    neighborPrimitiveId = primitiveId -
                        deltaIndex * neighborIndexDelta;

                    setPrimitiveToAveragePrediction(figure, primitiveId,
                        prevPrimitiveId, neighborPrimitiveId,
                        lastScaleValue, lastElongateValue);

					figure->getPrimitivePtr(primitiveId)->scaleBy(lastScaleValue);
                }
            }
        }

        startIndex = middleIndex + 1;
    }
	delete pList2;
}


#define WIDEN_UPPER_LIMIT		10
#define WIDEN_LOWER_LIMIT		.1
#define ELONGATE_UPPER_LIMIT	10
#define ELONGATE_LOWER_LIMIT	.1

#define EPSILON 1E-4
bool IsNumericallyEqual(double v1, double v2)
{
	if(fabs(v1-v2)<EPSILON)
		return true;
	else 
		return false;
}

void M3DSubfigureTransformation::widen(double value)
{
	int i;

	if(value>WIDEN_UPPER_LIMIT)
		value=WIDEN_UPPER_LIMIT;
	if(value<WIDEN_LOWER_LIMIT)
		value=WIDEN_LOWER_LIMIT;

#ifdef DEBUG
	fprintf(stderr, "dWiden: %f times in ACTUAL widening.\n", value);
#endif

	double loggedValue=log(value)/log(LOG_BASE);
	double lastValue;
	int widenIterations=(int)loggedValue;
	if(loggedValue>0)
		lastValue=pow(LOG_BASE, loggedValue-widenIterations);
	else
		lastValue=pow(LOG_BASE, loggedValue-widenIterations);

	if(loggedValue>0)
	{
		for ( i=0; i < widenIterations; i++)
		{
			widenStep(LOG_BASE);
			initStep();
		}
	}
	else
	{
		for ( i=0; i < -widenIterations; i++)
		{
			widenStep(INVERSED_LOG_BASE);
			initStep();
		}
	}

	if(!IsNumericallyEqual(lastValue, 1.0))
		widenStep(lastValue);
	reInit();
}
// the old version does not really do the job it was supposed to do
/*
void M3DSubfigureTransformation::widen(double value)
{
	int i;
	double loggedValue=ln(value)/ln(1.1);

	double delta = 0.1;
	double lastValue;


	int widenIterations = (int) fabs(value / delta);

	if (value > 0)
	{
		lastValue = value - delta * widenIterations;
		for ( i=0; i < widenIterations; i++)
		{
			widenStep(delta + 1.0);
			initStep();
		}
	}
	else
	{
		lastValue = value + delta * widenIterations;
		for ( i=0; i < widenIterations; i++)
		{
			widenStep(delta + 1.0);
			initStep();
		}
	}
	widenStep(lastValue + 1.0);
	reInit();
}
*/

void M3DSubfigureTransformation::elongateStep(double value)
{
//    M3DFigureTreeNode * parentNode;
    M3DQuadFigure * figure;
//    M3DQuadFigure * parentFigure;
    M3DPrimitive * primitive;
//    M3DPrimitive startPrimitive;
//    M3DPrimitive finishPrimitive;
//    M3DPrimitive centerPrimitive;
//    Vector3D startPosition, endPosition;
//    Vector3D oldNormal, newNormal;
//    Vector3D deltaX;
//    Vector3D axis;
//    Quat deltaQ;
    M3DPrimitiveLinkInfo * link;

    int linkCount,
        rowCount, columnCount,
        i, j,
        primitiveId,
        primitiveIndexDelta,
        neighborIndexDelta,
        numPrimitivesPerLink;

    int prevPrimitiveId,
        neighborPrimitiveId;


    if(figureTreeNode == NULL || object == NULL)
        return;

    figure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(figureTreeNode->getFigureId()));
    if(figure == NULL)
        return;

	//figure->print();

//    parentNode = figureTreeNode->getParent();
//    if(parentNode == NULL)
//        return;

//    parentFigure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(parentNode->getFigureId()));

    rowCount = figure->getRowCount();
    columnCount = figure->getColumnCount();

    if(rowOrientedLinks)
    {
        primitiveIndexDelta = columnCount;
        numPrimitivesPerLink = rowCount;
        neighborIndexDelta = 1;
    }
    else
    {
        primitiveIndexDelta = 1;
        numPrimitivesPerLink = columnCount;
        neighborIndexDelta = columnCount;
    }

    linkCount = figureTreeNode->getLinkCount();

    int middleIndex = linkCount / 2;
    int deltaIndex;

    lastElongateValue = value; //* lastElongateValue;

    int startIndex = middleIndex;
    for(deltaIndex = -1; deltaIndex <= 1; deltaIndex += 2)
    {
        for(i = startIndex; i < linkCount && i >= 0; i += deltaIndex)
        {
            link = figureTreeNode->getLink(i);
            if(link == NULL)
                continue;

            primitiveId = link->primitiveId;
            primitive = figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;

            if(i == middleIndex)
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;

                    setPrimitiveToPreviousPrediction(figure, primitiveId,
                        prevPrimitiveId, lastElongateValue);
                }
				/*
				printf("if(i == middleIndex) and (%d, %d)\n", i, middleIndex);
				figure->print();
				*/
            }
            else
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;
                    neighborPrimitiveId = primitiveId -
                        deltaIndex * neighborIndexDelta;

                    setPrimitiveToAveragePrediction(figure, primitiveId,
                        prevPrimitiveId, neighborPrimitiveId,
                        lastScaleValue, lastElongateValue);

					/*
					printf("for(j = 1; j < numPrimitivesPerLink; j++) and (%d, %d)\n", j, numPrimitivesPerLink);
					figure->print();
					*/
                }
				//printf("for(j = 1; j < numPrimitivesPerLink; j++) and (%d, %d)\n", j, numPrimitivesPerLink);
				//figure->print();
            }
        }

        startIndex = middleIndex + 1;
    }
}

void M3DSubfigureTransformation::elongate(double value)
{
	int i;

	if(value>ELONGATE_UPPER_LIMIT)
		value=ELONGATE_UPPER_LIMIT;
	if(value<ELONGATE_LOWER_LIMIT)
		value=ELONGATE_LOWER_LIMIT;

#ifdef DEBUG
	fprintf(stderr, "dElongate: %f times in ACTUAL elongation.\n", value);
#endif

	double loggedValue=log(value)/log(LOG_BASE);
	double lastValue;
	int elongateIterations=(int)loggedValue;
	if(loggedValue>0)
		lastValue=pow(LOG_BASE, loggedValue-elongateIterations);
	else
		lastValue=pow(LOG_BASE, loggedValue-elongateIterations);

	if(loggedValue>0)
	{
		for ( i=0; i < elongateIterations; i++)
		{
			elongateStep(LOG_BASE);
			initStep();
		}
	}
	else
	{
		for ( i=0; i < -elongateIterations; i++)
		{
			elongateStep(INVERSED_LOG_BASE);
			initStep();
		}
	}

	if(!IsNumericallyEqual(lastValue, 1.0))
		elongateStep(lastValue);
	reInit();
}
/*
void M3DSubfigureTransformation::elongate(double value)
{
	int i;
	double delta = 0.1;
	double lastValue;

	int elongateIterations = (int) fabs(value / delta);

	if (value > 0)
	{
		lastValue = value - delta * elongateIterations;
		for ( i=0; i < elongateIterations; i++)
		{
			elongateStep(delta + 1.0);
			initStep();
		}
	}
	else
	{
		lastValue = value + delta * elongateIterations;
		for ( i=0; i < elongateIterations; i++)
		{
			elongateStep(delta + 1.0);
			initStep();
		}
	}
	elongateStep(lastValue + 1.0);
//	elongateStep(value);
	reInit();
}
*/


void M3DSubfigureTransformation::hingeStep(double angle)
{
    M3DFigureTreeNode * parentNode;
    M3DQuadFigure * figure;
    M3DQuadFigure * parentFigure;
    M3DPrimitive * primitive, * primitivePost, * primitivePre;
    M3DQuadPrimitive startPrimitive;
    M3DQuadPrimitive finishPrimitive;
    Vector3D oldNormal, newNormal;
    Vector3D startPosition, endPosition;
    Vector3D deltaX;
    Vector3D axis, vector1, vector2;
    Quat deltaQ;
    M3DPrimitiveLinkInfo * link, * linkPost, * linkPre;

    bool doRotations;
    Real norm;

    int linkCount,
        rowCount, columnCount,
        i, j,
        primitiveId,
		primitivePostId,
		primitivePreId,
        primitiveIndexDelta,
        neighborIndexDelta,
        numPrimitivesPerLink;

    int prevPrimitiveId,
        neighborPrimitiveId;


    if(figureTreeNode == NULL || object == NULL || referenceFigure == NULL)
        return;

    figure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(figureTreeNode->getFigureId()));
    if(figure == NULL)
        return;

    parentNode = figureTreeNode->getParent();
    if(parentNode == NULL)
        return;

    parentFigure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(parentNode->getFigureId()));

    rowCount = figure->getRowCount();
    columnCount = figure->getColumnCount();


    if(rowOrientedLinks)
    {
        primitiveIndexDelta = columnCount;
        numPrimitivesPerLink = rowCount;
        neighborIndexDelta = 1;
    }
    else
    {
        primitiveIndexDelta = 1;
        numPrimitivesPerLink = columnCount;
        neighborIndexDelta = columnCount;
    }

    linkCount = figureTreeNode->getLinkCount();

    int middleIndex = linkCount / 2;
    int deltaIndex;

    int startIndex = middleIndex;
    for(deltaIndex = -1; deltaIndex <= 1; deltaIndex += 2)
    {
        for(i = startIndex; i < linkCount && i >= 0; i += deltaIndex)
        {
            link = figureTreeNode->getLink(i);
            if(link == NULL)
                continue;            

            primitiveId = link->primitiveId;
            primitive = figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;

            primitive = figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;

			linkPost = figureTreeNode->getLink(i+1);
            if(linkPost != NULL){
                primitivePostId = linkPost->primitiveId;
				primitivePost = figure->getPrimitivePtr(primitivePostId);
				if(primitivePost != NULL)
					vector1=primitivePost->getX()-primitive->getX();
			}
			linkPre = figureTreeNode->getLink(i-1);
            if(linkPre != NULL){
                primitivePreId = linkPre->primitiveId;
				primitivePre = figure->getPrimitivePtr(primitivePreId);
				if(primitivePre != NULL)
					vector2=primitive->getX()-primitivePre->getX();
			}

			if (linkPost == NULL || primitivePost == NULL)
				if (linkPre == NULL || primitivePre == NULL)
					continue;
				else
					axis=vector2;
			else
				if (linkPre == NULL || primitivePre == NULL)
					axis=vector1;
				else
					axis=(vector1+vector2)*0.5;

			axis.normalize();
			deltaQ.setAxisAngle(axis, angle);

            // Make sure rotation is not zero
            norm = deltaQ.norm();
            if(norm > 1.0 - R_SMALL_TOLERANCE &&
               norm < 1.0 + R_SMALL_TOLERANCE)
                doRotations = true;
            else
                doRotations = false;

            primitive = figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;

            if(doRotations)
            {
                primitive->rotateBy(deltaQ);
            }


            if(i == middleIndex)
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;

                    setPrimitiveToPreviousPrediction(figure, primitiveId,
                        prevPrimitiveId, lastElongateValue);
                }
            }

            else
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;
                    neighborPrimitiveId = primitiveId -
                        deltaIndex * neighborIndexDelta;

                    setPrimitiveToAveragePrediction(figure, primitiveId,
                        prevPrimitiveId, neighborPrimitiveId, lastScaleValue,
                        lastElongateValue);
                }
            }
        }

        startIndex = middleIndex + 1;
    }
}

void M3DSubfigureTransformation::hinge(double angle)
{
	int i;
	double delta = 0.1;
	int hingeIterations;
	double lastAngle;


	hingeIterations = (int) fabs(angle / delta);
	if (angle > 0.0)
		lastAngle = angle - delta * hingeIterations;
	else 
		lastAngle = angle + delta * hingeIterations;

	for ( i=0; i < hingeIterations; i++)
	{
		if (angle > 0.0)
			hingeStep(delta);
		else 
			hingeStep(-delta);
		initStep();
	}
	hingeStep(lastAngle);

//	hingeStep(angle);
	reInit();
}









// Make sure links are aligned along a row or a column
// Save whether they are row or column oriented
void M3DSubfigureTransformation::checkForValidLinks()
{
    int linkCount;
    int rowCount,
        columnCount;
    int primitiveId;
    M3DQuadFigure * figure;

    M3DPrimitiveLinkInfo * link;

    if(figureTreeNode == NULL || object == NULL)
        return;

    figure = dynamic_cast<M3DQuadFigure*>( object->getFigurePtr(figureTreeNode->getFigureId()));
    if(figure == NULL)
        return;

    linkCount = figureTreeNode->getLinkCount();

    rowCount = figure->getRowCount();
    columnCount = figure->getColumnCount();

    // No links
    if(linkCount < 1)
    {
        printf("Error: No links defined on subfigure.\n");
        invalidLinks = true;
        return;
    }

    // Single link (should be a single row (tube) figure)
    if(linkCount == 1)
    {
        // Invalid links
        if(rowCount != 1 && columnCount != 1)
        {
            printf("Error: Links are not an entire row or column\n");
            invalidLinks = true;
            return;
        }

        link = figureTreeNode->getLink(0);
        primitiveId = link->primitiveId;

        // Invalid links
        if((rowCount == 1 && (primitiveId != 0 && primitiveId != columnCount - 1)) ||
           (columnCount == 1 && (primitiveId != 0 && primitiveId != rowCount - 1)))
        {
            printf("Error: Links are not an entire row or column\n");
            invalidLinks = true;
            return;
        }

        // Links are oriented in column direction
        rowOrientedLinks = false;
        return;
    }

    link = figureTreeNode->getLink(0);
    primitiveId = link->primitiveId;

    link = figureTreeNode->getLink(1);
    if((link->primitiveId % columnCount == 0 &&
        primitiveId % columnCount == 0) ||
       (link->primitiveId % columnCount == columnCount-1 &&
        primitiveId % columnCount == columnCount-1))
        rowOrientedLinks = false;
    else
        rowOrientedLinks = true;

    // TODO : test the validity of the links
 /*
    for(i = 0; i < linkCount; i++)
    {
        link = figureTreeNode->getLink(i);
        if(link == NULL)
            continue;
    }
*/
}

bool M3DSubfigureTransformation::getInterpolatedPrimitive(
    M3DPrimitive & prim, M3DQuadFigure * figure, double u, double v, double t)
{
//    Xferlist * xfer;
//    ThallCode::Diatomgrid diatomGrid;
//    ThallCode::DiatomInterp diatomInterp;
//    ThallCode::Diatom diatom;
    Quat q, rotateQ;
	double umax, vmax;

	umax = (double) (figure->getColumnCount() - 1);
	vmax = (double) (figure->getRowCount() - 1);
    if(u > umax || u < 0.0 || v > vmax || v < 0.0)
        return false;

    rotateQ.setAxisAngle(Vector3D(1.0, 0.0, 0.0), -R_HALF_PI);

	fprintf(stderr, "ERROR: M3DSubfigureTransformation::getInterpolatedPrimitive()\n");
	fprintf(stderr, "   This code should not be used unless rewritten.\n");
/*
    xfer = convertM3DtoXfer(figure);


    diatomGrid.readXferlist(xfer);
    diatomInterp.init(&diatomGrid);

    diatom = diatomInterp.interpolatePrim(u, v);

    prim.setX(diatom.p_val().X(), diatom.p_val().Y(), diatom.p_val().Z());
    prim.setR(diatom.r_val());
    prim.setTheta(R_HALF_PI - diatom.phi_val());

    q.setX(diatom.q_val()(ThallCode::X));
    q.setY(diatom.q_val()(ThallCode::Y));
    q.setZ(diatom.q_val()(ThallCode::Z));
    q.setW(diatom.q_val()(ThallCode::W));

    q = q * rotateQ;
    prim.setQ(q);

    delete (xfer->atomlist);
    delete xfer;
*/
    return true;
}

/*
void M3DSubfigureTransformation::getSurfacePointAndNormal(
    Vector3D & point, Vector3D & normal, M3DQuadFigure * figure,
    double u, double v, double t)
{
//    Xferlist * xfer;
//    ThallCode::Diatomgrid diatomGrid;
 //   ThallCode::DiatomInterp diatomInterp;
//    ThallCode::Diatom diatom;

	fprintf(stderr, "ERROR: M3DSubfigureTransformation::getSurfacePointAndNormal()\n");
	fprintf(stderr, "   This code should not be used unless rewritten.\n");
//*

    xfer = convertM3DtoXfer(figure);

    diatomGrid.readXferlist(xfer);
    diatomInterp.init(&diatomGrid);

    diatom = diatomInterp.interpolatePrim(u, v);

    ThallCode::DbVector3 dbPoint, dbNormal;
    diatom.edgeboundary_t(t, dbPoint, dbNormal);

    point.set(dbPoint.X(), dbPoint.Y(), dbPoint.Z());
    normal.set(dbNormal.X(), dbNormal.Y(), dbNormal.Z());

    delete (xfer->atomlist);
    delete xfer;
// * /
}
*/

void M3DSubfigureTransformation::getSurfacePointAndNormal(
	ThallCode::Pointlist_server2* pList2, Vector3D & point, Vector3D & normal, double u, double v, double t)
{
    //ThallCode::DbVector3 dbPoint, dbNormal;
	Bpoint bpoint;
	pList2->subdivBPosition(&bpoint, u, v, t);
    point.set(bpoint.pnt);
	normal.set(bpoint.norm);
	normal.normalize();
}

void M3DSubfigureTransformation::setPrimitiveToPreviousPrediction(
    M3DQuadFigure * figure, int primitiveId, int prevPrimitiveId, double scale)
{
    M3DPrimitive * primitive,
                 * prevPrimitive,
                 * refPrimitive,
                 * refPrevPrimitive;

    Vector3D prevDiff;

    Quat q, dq;


    primitive = figure->getPrimitivePtr(primitiveId);
    prevPrimitive = figure->getPrimitivePtr(prevPrimitiveId);

    refPrimitive = referenceFigure->getPrimitivePtr(primitiveId);
    refPrevPrimitive = referenceFigure->getPrimitivePtr(prevPrimitiveId);

    if(primitive == NULL || prevPrimitive == NULL ||
       refPrimitive == NULL || refPrevPrimitive == NULL)
        return;

    prevDiff = scale * (refPrimitive->getX() - refPrevPrimitive->getX());

    dq = prevPrimitive->getQ() * refPrevPrimitive->getQ().conj();
    dq.rotateVector(prevDiff);

    q = refPrevPrimitive->getQ().conj() * refPrimitive->getQ();

    primitive->setX(prevPrimitive->getX() + prevDiff);
    primitive->setQ(prevPrimitive->getQ() * q);
}

void M3DSubfigureTransformation::setPrimitiveToAveragePrediction(
    M3DQuadFigure * figure, int primitiveId, int prevPrimitiveId,
    int neighborPrimitiveId, double scale, double elongation)
{
    M3DPrimitive * primitive,
                 * prevPrimitive,
                 * neighborPrimitive,
                 * refPrimitive,
                 * refPrevPrimitive,
                 * refNeighborPrimitive;

    M3DQuadPrimitive prevPrediction,
                 neighborPrediction;

    Vector3D prevDiff,
             neighborDiff;
    double prevDist,
           neighborDist,
           weight;

    Quat q, dq;

//    double angle;

    primitive = figure->getPrimitivePtr(primitiveId);
    prevPrimitive = figure->getPrimitivePtr(prevPrimitiveId);
    neighborPrimitive = figure->getPrimitivePtr(neighborPrimitiveId);

    refPrimitive = referenceFigure->getPrimitivePtr(primitiveId);
    refPrevPrimitive = referenceFigure->getPrimitivePtr(prevPrimitiveId);
    refNeighborPrimitive = referenceFigure->getPrimitivePtr(neighborPrimitiveId);

    if(primitive == NULL || prevPrimitive == NULL ||
       neighborPrimitive == NULL || refPrimitive == NULL ||
       refPrevPrimitive == NULL || refNeighborPrimitive == NULL)
        return;

    // Calculate the averaging weight based on prev/neighbor distances
    prevDiff = elongation * (refPrimitive->getX() - refPrevPrimitive->getX());
    neighborDiff = scale * (refPrimitive->getX() - refNeighborPrimitive->getX());
//    prevDist = prevDiff.norm();
//    neighborDist = neighborDiff.norm();
//    weight = neighborDist / (prevDist + neighborDist);

    dq = prevPrimitive->getQ() * refPrevPrimitive->getQ().conj();

    // Take quaternion to the "elongate" power
//    angle = 2.0 * acos(dq.getW()) * elongation;
//    dq.setAxisAngle(dq.getVector(), angle);

    dq.rotateVector(prevDiff);

    q = refPrevPrimitive->getQ().conj() * refPrimitive->getQ();

    prevPrediction.setX(prevPrimitive->getX() + prevDiff);
    prevPrediction.setQ(prevPrimitive->getQ() * q);
    prevPrediction.setR(primitive->getR());
    prevPrediction.setTheta(primitive->getTheta());

    dq = neighborPrimitive->getQ() * refNeighborPrimitive->getQ().conj();

    // Take quaternion to the "scale" power
//    angle = 2.0 * acos(dq.getW()) * scale;
//    dq.setAxisAngle(dq.getVector(), angle);

    dq.rotateVector(neighborDiff);

    q = refNeighborPrimitive->getQ().conj() * refPrimitive->getQ();

    neighborPrediction.setX(neighborPrimitive->getX() + neighborDiff);
    neighborPrediction.setQ(neighborPrimitive->getQ() * q);
    neighborPrediction.setR(primitive->getR());
    neighborPrediction.setTheta(primitive->getTheta());

    prevDiff = primitive->getX() - prevPrimitive->getX();
    neighborDiff = primitive->getX() - neighborPrimitive->getX();
    prevDist = prevDiff.norm();
    neighborDist = neighborDiff.norm();

    weight = elongation * neighborDist / (scale * prevDist + elongation * neighborDist);

    averagePrimitives(*primitive, prevPrediction, neighborPrediction, weight);
    //(*primitive) = neighborPrediction;
}

void M3DSubfigureTransformation::averagePrimitives(M3DPrimitive & avePrim,
    const M3DPrimitive & prim1, const M3DPrimitive & prim2, double weight)
{
    Quat q;

    double weight2 = 1.0 - weight;

    avePrim.setX(weight2 * prim1.getX() + weight * prim2.getX());
    avePrim.setQ(q.slerp(weight, prim1.getQ(), prim2.getQ()));
    avePrim.setR(weight2 * prim1.getR() + weight * prim2.getR());
    avePrim.setTheta(weight2 * prim1.getTheta() + weight * prim2.getTheta());
}

