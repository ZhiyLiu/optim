#include "M3DNEWUOAOptimizer.h"
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

M3DNEWUOAOptimizer::M3DNEWUOAOptimizer()
{

}
void M3DNEWUOAOptimizer::setImage(ImageDistanceMap* imageDistanceMap)
{
    mSignedDistanceImage = imageDistanceMap;
}

void M3DNEWUOAOptimizer::setObject(M3DObject* sreps)
{
    mSreps = sreps;
}
double M3DNEWUOAOptimizer::computeSradPenalty()
{
    double sradPenalty = 0.0;
    // assume figure is quad figure
    M3DAtomPredictorQuad *atomPredictor = new M3DAtomPredictorQuad();
    M3DFigure *figure = mSreps->getFigurePtr(mFigureIndex);
    sradPenalty = atomPredictor->getFigureRSradPenalty(figure, -1, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold));
    return sradPenalty;
}
double M3DNEWUOAOptimizer::operator() (double *coeff)
{
    double cost = 0.0;
    cost = this->getObjectiveFunctionValue(coeff, 15, 1);

    return cost;
}
void M3DNEWUOAOptimizer::updateSpokeDir(const double *coeff, int spokeIndex, Vector3D& dir)
{
    int baseAddr = spokeIndex * 4;
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
void M3DNEWUOAOptimizer::updateFigure(const double *coeff, int figureId)
{
    try
    {
        //std::cout << "Updating figure id: " << figureId << ", after the " << iterNum2++ << "th iteration" << std::endl;

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

                // update R0 and R1
                double scale0 = exp(coeff[spokeIndex * 4 + 3]);
                double oldR0 = quadAtom->getR0();
                quadAtom->setR0(oldR0 * scale0);

                double scale1 = exp(coeff[(spokeIndex+1)*4 + 3]);
                double oldR1 = quadAtom->getR1();
                quadAtom->setR1(oldR1 * scale1);

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

                // update U0 and U1 in quadAtom
                standardVersion->setU0(U0);
                standardVersion->setU1(U1);
                // update R0 and R1
                double scale0 = exp(coeff[spokeIndex * 4 + 3]);
                double oldR0 = standardVersion->getR0();
                standardVersion->setR0(oldR0 * scale0);

                double scale1 = exp(coeff[(spokeIndex+1)*4 + 3]);
                double oldR1 = standardVersion->getR1();
                standardVersion->setR1(oldR1 * scale1);

                M3DQuadEndPrimitive* endVersion = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);
                Vector3D UEnd;

                // update UEnd
                updateSpokeDir(coeff, spokeIndex+2, UEnd);
                endVersion->setUEnd(UEnd);
                // update REnd
                double scaleEnd = exp(coeff[(spokeIndex+2) * 4 + 3]);
                double oldREnd = endVersion->getREnd();
                endVersion->setREnd(oldREnd * scaleEnd);

                spokeIndex += 3;
            }
        }

    }
    catch(std::exception& e)
    {
        std::cout << "[Error]updateFigure has an exception:"<< e.what() << std::endl;
    }

}

int iterNum3 = 0;
/* Compute the total entropy. */
double M3DNEWUOAOptimizer::getObjectiveFunctionValue(const double *coeff, double w1, double w2)
{
    // input coeff from newuoa is the coeff of length. it could be positive or negative
    iterNum3++;
    std::cout << "The " << iterNum3 << " th iteration..." << std::endl;
    double normalMatch = 0.0;
    double sradPenalty = 0.0;
    double imageDist = 0.0;
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

    if(iterNum3 == 1)
    {
        std::cout << "Adjust weights to balance between imageDist:" << imageDist << ", normalMatch:" << normalMatch << ", and sradPenalty:" << sradPenalty << std::endl;
    }
    recoverSrep();
    // finally, we will get the minimum of objective function value
    return objFunctionValue;
}
void M3DNEWUOAOptimizer::recoverSrep()
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

        quadAtom->setR0(mOrigLength[spokeIndex]);
        quadAtom->setR1(mOrigLength[spokeIndex+1]);

        M3DQuadEndPrimitive* endPrimitive = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);

        if(endPrimitive != NULL)
        {
            endPrimitive->setUEnd(mOrigDirs[baseAddr+2]);
            endPrimitive->setREnd(mOrigLength[baseAddr+2]);
            spokeIndex += 1;
        }
        spokeIndex += 2;
    }

}

/* Main entry of optimizer */
int M3DNEWUOAOptimizer::perform(M3DObject* outputModel)
{
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
    double coeffArray[spokeCount * 4];
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
        mOrigLength.push_back(quadAtom->getR0());
        mOrigLength.push_back(quadAtom->getR1());

        int baseAddr = spokeIndex * 4;
        coeffArray[baseAddr] = U0.getX();
        coeffArray[baseAddr+1] = U0.getY();
        coeffArray[baseAddr+2] = U0.getZ();
        coeffArray[baseAddr+3] = 0; // initially, the power of the exponential is 0

        coeffArray[baseAddr+4] = U1.getX();
        coeffArray[baseAddr+5] = U1.getY();
        coeffArray[baseAddr+6] = U1.getZ();
        coeffArray[baseAddr+7] = 0;

        M3DQuadEndPrimitive* endPrimitive = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);

        if(endPrimitive != NULL)
        {
            Vector3D UEnd = endPrimitive->getUEnd();
            mOrigDirs.push_back(UEnd);
            mOrigLength.push_back(endPrimitive->getREnd());
            coeffArray[baseAddr+8] = UEnd.getX();
            coeffArray[baseAddr+9] = UEnd.getY();
            coeffArray[baseAddr+10] = UEnd.getZ();
            coeffArray[baseAddr+11] = 0;
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

//    min_newuoa(spokeCount*4, coeffArray, *this, 0.01, 0.0001, 2000);
    min_newuoa(spokeCount*4, coeffArray, *this, newuoaParams[6], newuoaParams[7], newuoaParams[8]);
    // update figure last time
    updateFigure(coeffArray, mFigureIndex);

    // output file
    M3DObjectFile objectFile;
    objectFile.write("out.m3d", *mSreps);
    return 0;

}
