/* The purpose of this class is to optimize spoke length by NEWUOA
 *
 * Zhiyuan Liu
 * 2018.01
 */
#include "M3DSpokeAngleOptimizer.h"

#include <iostream>
#include <exception>
#include "newuoa.h"
#include "toolsfunc.h"
#include "M3DQuadFigure.h"
#include "movespokes.h"
#include "M3DAtomPredictorQuad.h"
#include "NormalMatchComputer.h"
#include "SimilarityComputer.h"
#include "M3DObjectFile.h"
//#include "quaternion.h"

 M3DSpokeAngleOptimizer::M3DSpokeAngleOptimizer()
     :mSreps(NULL){

}

M3DSpokeAngleOptimizer::M3DSpokeAngleOptimizer(M3DObject* sreps)
{
    if(sreps != NULL)
    {
        mSreps = sreps;
    }

}

double M3DSpokeAngleOptimizer::computeSradPenalty()
{
    double sradPenalty = 0.0;
    // assume figure is quad figure
    M3DAtomPredictorQuad *atomPredictor = new M3DAtomPredictorQuad();
    M3DFigure *figure = mSreps->getFigurePtr(mFigureIndex);
    sradPenalty = atomPredictor->getFigureRSradPenalty(figure, -1, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold));
    return sradPenalty;
}

int iterNum = 0;
void M3DSpokeAngleOptimizer::updateSpokeDir(const double *coeff, int spokeIndex, Vector3D& dir)
{
    int baseAddr = spokeIndex * 4;
    double qw = coeff[baseAddr];
    double qx = coeff[baseAddr + 1];
    double qy = coeff[baseAddr + 2];
    double qz = coeff[baseAddr + 3];

    // update direction vector here
    rotateByQuaternion(qw, qx, qy, qz, dir);

    dir.normalize();
}

double a, b, c,d; // debug code
/* Update spoke properties after each optimization */
void M3DSpokeAngleOptimizer::updateFigure(const double *coeff, int figureId)
{
    try
    {
        std::cout << "Updating figure id: " << figureId << ", after the " << iterNum++ << "th iteration" << std::endl;

        M3DFigure* figure = mSreps->getFigurePtr(figureId);
        int spokeCount = figure->getSpokeCount();
        for(int iSpoke = 0; iSpoke < spokeCount; ++iSpoke)
        {
            std::cout << "The angle of the rotation of " << iSpoke + 1 << "th spoke is:" << 114.6 * acos(coeff[4*iSpoke]) << " degree" << std::endl;
            std::cout << "around the vector [" << coeff[4*iSpoke + 1] << " " << coeff[4*iSpoke + 2] << " " << coeff[4*iSpoke + 3] << "]" << std::endl; // keep sin(theta/2) factor in this vector
        }
        int primitiveCount = figure->getPrimitiveCount();
        int spokeIndex = 0;
        for(int i = 0; i < primitiveCount; ++i)
        {
            //debug code
            if(i == 25)
            {
                if(coeff[spokeIndex * 4] != a || coeff[spokeIndex * 4 + 1] != b || coeff[spokeIndex* 4 +2] != c || coeff[spokeIndex*4 + 3] != d)
                {
                    bool isHere = true;
                }
            }
            // The i-th primitive
            M3DPrimitive* currentPrimitive = figure->getPrimitivePtr(i);
            M3DQuadEndPrimitive* endPrimitive = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);

            if(endPrimitive == NULL)
            {
                // it is an standard primitive
                M3DQuadPrimitive* quadAtom = dynamic_cast<M3DQuadPrimitive*>(currentPrimitive);

                Vector3D U0 = quadAtom->getU0();
                Vector3D U1 = quadAtom->getU1();

                // compute new U0 and U1
                updateSpokeDir(coeff, spokeIndex, U0);
                updateSpokeDir(coeff, spokeIndex+1, U1);

                // update U0 and U1 in quadAtom
                quadAtom->setU0(U0);
                quadAtom->setU1(U1);

                // update spokeIndex
                spokeIndex += 2;
            }
            else
            {
                M3DQuadPrimitive* standardVersion = dynamic_cast<M3DQuadPrimitive*> (currentPrimitive);

                Vector3D U0 = standardVersion->getU0();
                Vector3D U1 = standardVersion->getU1();

                // compute new U0 and U1
                updateSpokeDir(coeff, spokeIndex, U0);
                updateSpokeDir(coeff, spokeIndex+1, U1);

                // update U0 and U1 in quadAtom
                standardVersion->setU0(U0);
                standardVersion->setU1(U1);

                //update spokeIndex
                spokeIndex += 2;
                
                M3DQuadEndPrimitive* endVersion = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);
                Vector3D UEnd = endVersion->getUEnd();

                // update UEnd
                updateSpokeDir(coeff, spokeIndex, UEnd);
                endVersion->setUEnd(UEnd);
                spokeIndex += 1;
            }
        }
    }
    catch(std::exception& e)
    {
        std::cout << "[Error]updateFigure has an exception:"<< e.what() << std::endl;
    }

}
double minImageDist = 10000.0;

/* Compute the total entropy. */
double M3DSpokeAngleOptimizer::getObjectiveFunctionValue(const double *coeff, double w1, double w2)
{
    // input coeff from newuoa is the coeff of length. it could be positive or negative

    // 0. Update new length to each spoke
    updateFigure(coeff, mFigureIndex);

    double objFunctionValue = 0.0;

    // 1. spoke normal match
    double w_NormalPenalty = 9.0;
/*
    NormalMatchComputer normalMatchComputer;
    normalMatchComputer.setSrepModel(mSreps);
    normalMatchComputer.setTargetFigureIndex(mFigureIndex);
    normalMatchComputer.setTargetImage(mSignedDistanceImage);

    double normalMatch = 0.0;
    if(normalMatchComputer.compute(&normalMatch) == false)
    {
        std::cout << "[Error]Error encountered when compute normal match measure" << std::endl;
        return -999.0;
    }
*/
    // 2. measure sum square of distance from implied boundary to expected image boundary
    SimilarityComputer similarityCompter;
    similarityCompter.setSrepModel(mSreps);
    similarityCompter.setTargetFigureIndex(mFigureIndex);
    similarityCompter.setTargetImage(mSignedDistanceImage);

    double imageDist = 0.0;
    double normalMatch = 0.0;
    if(similarityCompter.compute(&imageDist, &normalMatch) == false)
    {
        std::cout << "[Error]Error encountered when compute similarity measure" << std::endl;
        return -999.0;
    }

    // debug code
    if(imageDist < minImageDist)
    {
        minImageDist = imageDist;
    }
    double w_ImageMatch = 5.0;
    objFunctionValue = w_ImageMatch * imageDist /*+ w_NormalPenalty * normalMatch*/;

    return objFunctionValue;
}

double M3DSpokeAngleOptimizer::operator() (double *coeff)
{
    double cost = 0.0;
    cost = this->getObjectiveFunctionValue(coeff, 15, 1);
    recoverSrep();
    return cost;
}

void M3DSpokeAngleOptimizer::setImage(ImageDistanceMap* imageDistanceMap)
{
    mSignedDistanceImage = imageDistanceMap;
}

void M3DSpokeAngleOptimizer::setObject(M3DObject* sreps)
{
    mSreps = sreps;
}

void M3DSpokeAngleOptimizer::rotateByQuaternion(double a, double b, double c, double d, Vector3D& target)
{
    // a is cos(theta/2)
    // refer to wiki for construction of rotation matrix from quaternion 
    double rotArray[9];
    // normalize the quaternion
    double length = sqrt(a*a + b*b + c*c + d*d);
    a /= length; b/= length; c/=length; d/= length;
//    rotArray[0] = a * a + b * b - c * c - d * d;
    rotArray[1] = 2 * (b * c - a * d);
    rotArray[2] = 2 * (b * d + a * c);
    rotArray[3] = 2 * (b * c + a * d);
//    rotArray[4] = a * a - b * b + c * c - d * d;
    rotArray[5] = 2 * (c * d - a * b);
    rotArray[6] = 2 * (b * d - a * c);
    rotArray[7] = 2 * (c * d + a * b);
//    rotArray[8] = a * a - b * b - c * c + d * d;


    rotArray[4] = 1- 2 * (b * b + d * d);
    rotArray[0] = 1 - 2 * (c * c + d * d);
    rotArray[8] = 1 - 2 * (b * b + c * c);
    double oldX = target.getX();
    double oldY = target.getY();
    double oldZ = target.getZ();
    target.setX(rotArray[0] * oldX + rotArray[1] * oldY + rotArray[2] * oldZ);
    target.setY(rotArray[3] * oldX + rotArray[4] * oldY + rotArray[5] * oldZ);
    target.setZ(rotArray[6] * oldX + rotArray[7] * oldY + rotArray[8] * oldZ);

}

void M3DSpokeAngleOptimizer::recoverSrep()
{
    if(mSreps == NULL || mSignedDistanceImage == NULL)
    {
        std::cout << "[Error]object or signed distance image is empty in recover function" << std::endl;
        return;
    }

    M3DFigure* figure = mSreps->getFigurePtr(mFigureIndex);
    int spokeIndex = 0;
    int primitiveCount = figure->getPrimitiveCount();
    for(int i = 0; i < primitiveCount; ++i)
    {
        M3DPrimitive* currentPrimitive = figure->getPrimitivePtr(i);
        M3DQuadPrimitive* quadAtom = dynamic_cast<M3DQuadPrimitive*>(currentPrimitive);
        quadAtom->setU0(mOrigDirs[spokeIndex]);
        quadAtom->setU1(mOrigDirs[spokeIndex+1]);

        M3DQuadEndPrimitive* endPrimitive = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);

        if(endPrimitive != NULL)
        {
            endPrimitive->setUEnd(mOrigDirs[spokeIndex+2]);
            spokeIndex += 1;
        }
        spokeIndex += 2;
    }

}

/* Main entry of optimizer */
int M3DSpokeAngleOptimizer::perform(M3DObject* outputModel)
{
//    using namespace quaternion;
//    typedef rotation_matrix<double> RM;
    if(mSreps == NULL || mSignedDistanceImage == NULL)
    {
        std::cout << "[Error]object or signed distance image is empty in perform function" << std::endl;
        return -1;
    }

    // optimization
    M3DFigure* figure = mSreps->getFigurePtr(mFigureIndex);
    // initialize quaternions
    int primitiveCount = figure->getPrimitiveCount();
    for(int i = 0; i < primitiveCount; ++i)
    {
    // The i-th primitive
        M3DPrimitive* currentPrimitive = figure->getPrimitivePtr(i);
        M3DQuadPrimitive* quadAtom = dynamic_cast<M3DQuadPrimitive*>(currentPrimitive);

        Vector3D U0 = quadAtom->getU0();
        Vector3D U1 = quadAtom->getU1();

        U0.normalize();
        U1.normalize();
        mOrigDirs.push_back(U0);
        mOrigDirs.push_back(U1);

/****************************test code for quaternion****************************/
/*
        Vector3D hinge = U0.cross(U1);
        hinge.normalize();
        double theta = 0.5 * acos(U0 * U1);
        double sinTheta = sin(theta);
        rotateByQuaternion(cos(theta), sinTheta * hinge.getX(), sinTheta * hinge.getY(), sinTheta * hinge.getZ(), U0);
        U0.normalize();
        quadAtom->setU0(U0);
//        rotateByQuaternion(0.7071, 0.0, 0.7071, 0.0, p); // rotate 90 degree about y axis
*/
/****************************test code for quaternion****************************/

        //Rotation q0, q1;
        // rotate 0 radians about any vector
        //q0.a = 1;        q0.b = 0;  q0.c = 0; q0.d = 0;
        //q1.a = 1;        q1.b = 0;  q1.c = 0; q1.d = 0;
        //mRotations.push_back(q0);
        //mRotations.push_back(q1);
        mQuaternionSet.push_back(1); mQuaternionSet.push_back(0); mQuaternionSet.push_back(0); mQuaternionSet.push_back(0); // for top spoke
        mQuaternionSet.push_back(1); mQuaternionSet.push_back(0); mQuaternionSet.push_back(0); mQuaternionSet.push_back(0); // for bottom spoke
//        mQuaternionSet.push_back(0.866); mQuaternionSet.push_back(0.5); mQuaternionSet.push_back(0); mQuaternionSet.push_back(0); // for top spoke


        M3DQuadEndPrimitive* endPrimitive = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);

        if(endPrimitive != NULL)
        {
            mQuaternionSet.push_back(1); mQuaternionSet.push_back(0); mQuaternionSet.push_back(0); mQuaternionSet.push_back(0); // for end spoke
            Vector3D UEnd = endPrimitive->getUEnd();
            mOrigDirs.push_back(UEnd);
        }
    }

    // for each rotation structure, there are 4 elements
/*    double coeffArray[mQuaternionSet.size()];
    for(int i = 0; i < mRotations.size(); ++i)
    {
        coeffArray[i]   = mRotations[i].a;
        coeffArray[i+1] = mRotations[i].b;
        coeffArray[i+2] = mRotations[i].c;
        coeffArray[i+3] = mRotations[i].d;
    }
*/
    double *coeffArray = &mQuaternionSet[0];
    a = coeffArray[100]; b = coeffArray[101]; c = coeffArray[102]; d = coeffArray[103];//debug code 
    min_newuoa(mQuaternionSet.size(), coeffArray, *this, 1.0, 0.001, 300);

    M3DObjectFile objectFile;
    objectFile.write("out.m3d", *mSreps);

    return 0;

}

