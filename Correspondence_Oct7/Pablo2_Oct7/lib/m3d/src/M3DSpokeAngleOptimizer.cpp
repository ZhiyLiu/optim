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
#define PI 3.1415926
#define RATIO 180/PI
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

int iterNum2 = 0;

void M3DSpokeAngleOptimizer::updateSpokeDir(const double *coeff, int spokeIndex, Vector3D& dir)
{
    int baseAddr = spokeIndex * 3;
    double x = coeff[baseAddr];
    double y = coeff[baseAddr+1];
    double z = coeff[baseAddr+2];

    Vector3D newDir(x, y, z);
    newDir.normalize();
    dir.setX(x);
    dir.setY(y);
    dir.setZ(z);
    dir.normalize();
}

/* Update spoke properties after each optimization */
void M3DSpokeAngleOptimizer::updateFigure(const double *coeff, int figureId)
{
    try
    {
        std::cout << "Updating figure id: " << figureId << ", after the " << iterNum2++ << "th iteration" << std::endl;

        M3DFigure* figure = mSreps->getFigurePtr(figureId);
//        M3DQuadPrimitive* quadAtom = dynamic_cast<M3DQuadPrimitive*>(figure->getPrimitivePtr(25));

//        Vector3D U0 = quadAtom->getU0();
//        updateSpokeDir(coeff, 0, U0);
//        quadAtom->setU0(U0);


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
double normalMatch = 0.0;
double sradPenalty = 0.0;
double imageDist = 0.0;

/* Compute the total entropy. */
double M3DSpokeAngleOptimizer::getObjectiveFunctionValue(const double *coeff, double w1, double w2)
{
    // input coeff from newuoa is the coeff of length. it could be positive or negative

    // 0. Update new length to each spoke
    updateFigure(coeff, mFigureIndex);

    double objFunctionValue = 0.0;


    // 2. measure sum square of distance from implied boundary to expected image boundary
    SimilarityComputer similarityCompter;
    similarityCompter.setSrepModel(mSreps);
    similarityCompter.setTargetFigureIndex(mFigureIndex);
    similarityCompter.setTargetImage(mSignedDistanceImage);


    if(similarityCompter.compute(&imageDist, &normalMatch) == false)
    {
        std::cout << "[Error]Error encountered when compute similarity measure" << std::endl;
        return -999.0;
    }

    // These setting are best for refining modify_angle.m3d
    // 1. spoke normal match
//    double w_NormalPenalty = 1000.0;
//    double w_ImageMatch = 0.005;
//    double w_sradPenalty = 20.0;
//    double w_srepModelPenalty = 9.0;
    // These are for real data
//    double w_ImageMatch = 0.00001;
//    double w_sradPenalty = 2000.0;
//    double w_NormalPenalty = 50.0;
    double w_ImageMatch = mWtImageMatch;
    double w_sradPenalty = mWtSradPenalty;
    double w_NormalPenalty = mWtNormalPenalty;
    
    // 2. measure how far from regular srep model
    sradPenalty = computeSradPenalty();

    // 3. spoke model penalty
    DistanceType distType = (enum DistanceType) (int) tuningWt(BpSpokeDistanceType);
    double spokeModelpenalty = mSreps->dist2FromObject(mSreps->loadedObject(),
                                                       mFigureIndex, distType);



    // normal match: the larger, the better
    objFunctionValue =  w_ImageMatch * imageDist + w_NormalPenalty * normalMatch + w_sradPenalty * sradPenalty ;

    if(iterNum2 == 1)
    {
        // the first iteration, output the object function terms for tuning weights
        std::cout << "For the input, the image distance is:" << imageDist << ", the normal penalty is:" << normalMatch << ", sradPenalty is:" << sradPenalty << std::endl;
    }
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
        int baseAddr = spokeIndex; // mOrigDirs store direction vector3D for each spoke
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

    // read parameters from config file
    double newuoaParams[12];
    std::ifstream infile("/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/newuoa.config");
    if(!infile.is_open())
    {
        std::cerr << "Failed to open config file: newuoa.config." << std::endl;
        return -1;
    }

    std::string line;
    std::getline(infile, line);
    int i = 0;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        newuoaParams[i] = atof(line.c_str());
        i++;
        std::getline(infile,line);
    }

    infile.close();
    mWtImageMatch = newuoaParams[9];
    mWtSradPenalty = newuoaParams[10];
    mWtNormalPenalty = newuoaParams[11];

//    min_newuoa(spokeCount*3, coeffArray, *this, 0.01, 0.001, 1000);
    min_newuoa(spokeCount*3, coeffArray, *this, newuoaParams[6], newuoaParams[7], newuoaParams[8]);
    // update figure last time
    updateFigure(coeffArray, mFigureIndex);

    // output file
    M3DObjectFile objectFile;
    objectFile.write("out.m3d", *mSreps);
    return 0;

}

