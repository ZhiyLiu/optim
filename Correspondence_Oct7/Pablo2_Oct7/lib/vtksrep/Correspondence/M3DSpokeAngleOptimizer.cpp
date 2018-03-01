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
    int baseAddr = spokeIndex * 3;
    double x = coeff[baseAddr];
    double y = coeff[baseAddr+1];
    double z = coeff[baseAddr+2];
    
    Vector3D newDir(x, y, z);
    newDir.normalize();
    dir.normalize();
    double dotProduct = newDir * dir;
    if(acos(dotProduct) > 30)
    {
        return;
    }
    dir.setX(x);
    dir.setY(y);
    dir.setZ(z);
//    double qw = coeff[baseAddr];
//    double qx = coeff[baseAddr + 1];
//    double qy = coeff[baseAddr + 2];
//    double qz = coeff[baseAddr + 3];

  
// update direction vector here
//    rotateByQuaternion(qw, qx, qy, qz, dir);

    dir.normalize();
}

/* Update spoke properties after each optimization */
void M3DSpokeAngleOptimizer::updateFigure(const double *coeff, int figureId)
{
    try
    {
        std::cout << "Updating figure id: " << figureId << ", after the " << iterNum++ << "th iteration" << std::endl;

        M3DFigure* figure = mSreps->getFigurePtr(figureId);

        int primitiveCount = figure->getPrimitiveCount();
        int spokeIndex = 0;
        for(int i = 0; i < primitiveCount; ++i)
        {
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

                Vector3D U0;
                Vector3D U1;

                // compute new U0 and U1
                updateSpokeDir(coeff, spokeIndex, U0);
                updateSpokeDir(coeff, spokeIndex+1, U1);

                Vector3D oldU0 = standardVersion->getU0();
                double theta = acos(oldU0 * U0);
                if(theta > 0.01)
                {
                    bool isHere = true;
                }
                // update U0 and U1 in quadAtom
                standardVersion->setU0(U0);
                standardVersion->setU1(U1);

                M3DQuadEndPrimitive* endVersion = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);
                Vector3D UEnd;

                // update UEnd
                updateSpokeDir(coeff, spokeIndex+2, UEnd);
                endVersion->setUEnd(UEnd);
                spokeIndex += 3;
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
    
    // 2. measure how far from regular srep model
    double sradPenalty = computeSradPenalty();

    // 3. spoke model penalty
    DistanceType distType = (enum DistanceType) (int) tuningWt(BpSpokeDistanceType);
    double spokeModelpenalty = mSreps->dist2FromObject(mSreps->loadedObject(),
                                                       mFigureIndex, distType);

    double w_sradPenalty = 1.0;
    double w_srepModelPenalty = 9.0;  // may be too large, was 10

    // normal match: the larger, the better
    objFunctionValue = 1/normalMatch;// + w_ImageMatch * imageDist + w_sradPenalty * sradPenalty + w_srepModelPenalty * spokeModelpenalty;

    recoverSrep();
    // finally, we will get the minimum of objective function value
    return objFunctionValue;
}

double M3DSpokeAngleOptimizer::operator() (double *coeff)
{
    double cost = 0.0;
    cost = this->getObjectiveFunctionValue(coeff, 15, 1);

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

    //limit the theta/2 in the range [-45, 45]
    if(fabs(a) < 0.07071)
    {
        return;
    }
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

    if(target.getX() != oldX || target.getY() != oldY || target.getZ() != oldZ)
    {
        bool imHere = true;
    }
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
        int baseAddr = spokeIndex * 3;
        quadAtom->setU0(mOrigDirs[baseAddr]);
        quadAtom->setU1(mOrigDirs[baseAddr+1]);

        M3DQuadEndPrimitive* endPrimitive = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);

        if(endPrimitive != NULL)
        {
            endPrimitive->setUEnd(mOrigDirs[baseAddr+2]);
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
    int spokeCount = figure->getSpokeCount();
    double coeffArray[spokeCount * 3];
    int spokeIndex = 0;
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

        int baseAddr = spokeIndex * 3;
        coeffArray[baseAddr] = U0.getX();
        coeffArray[baseAddr+1] = U0.getY();
        coeffArray[baseAddr+2] = U0.getZ();
        coeffArray[baseAddr+3] = U1.getX();
        coeffArray[baseAddr+4] = U1.getY();
        coeffArray[baseAddr+5] = U1.getZ();

        M3DQuadEndPrimitive* endPrimitive = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);

        if(endPrimitive != NULL)
        {
            Vector3D UEnd = endPrimitive->getUEnd();
            mOrigDirs.push_back(UEnd);
            coeffArray[baseAddr+6] = UEnd.getX();
            coeffArray[baseAddr+7] = UEnd.getY();
            coeffArray[baseAddr+8] = UEnd.getZ();
            spokeIndex += 1;
        }
        spokeIndex += 2;
    }

    min_newuoa(spokeCount*3, coeffArray, *this, 0.8, 0.001, 5000);

    M3DObjectFile objectFile;
    objectFile.write("out.m3d", *mSreps);

    return 0;

}

