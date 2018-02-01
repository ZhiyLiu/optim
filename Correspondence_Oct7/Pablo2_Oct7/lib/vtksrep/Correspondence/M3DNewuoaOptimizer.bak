/* This is an wrapper of optimizer using NEWUOA
 *
 * Zhiyuan Liu
 * 2017.11
 */
#include "M3DSpokeLengthOptimizer.h"
#include <ctime>
#include <iostream>
#include <exception>
#include "newuoa.h"
#include "toolsfunc.h"
#include "M3DQuadFigure.h"
#include "movespokes.h"
#include "M3DAtomPredictorQuad.h"
#include "SimilarityComputer.h"
#include "M3DObjectFile.h"

M3DSpokeLengthOptimizer::M3DSpokeLengthOptimizer()
    :mSreps(NULL){}

M3DSpokeLengthOptimizer::M3DSpokeLengthOptimizer(M3DObject* sreps)
{
    if(sreps != NULL)
    {
        mSreps = sreps;
    }
}

double M3DSpokeLengthOptimizer::computeSradPenalty()
{
    double sradPenalty = 0.0;
    // assume figure is quad figure
    M3DAtomPredictorQuad *atomPredictor = new M3DAtomPredictorQuad();
    M3DFigure *figure = mSreps->getFigurePtr(mFigureIndex);
    sradPenalty = atomPredictor->getFigureRSradPenalty(figure, -1, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold));
    return sradPenalty;
}

int iterNum = 0;
/* Update spoke properties after each optimization */
void M3DSpokeLengthOptimizer::updateFigure(const double *coeff, int figureId)
{
    try
    {
        std::cout << "Updating figure id: " << figureId << ", after the " << iterNum++ << "th iteration" << std::endl;
        
        M3DFigure* figure = mSreps->getFigurePtr(figureId);
        int spokeCount = figure->getSpokeCount();
        for(int iSpoke = 0; iSpoke < spokeCount; ++iSpoke)
        {
            std::cout << "The length of the " << iSpoke + 1 << "th spoke is:" << exp(coeff[iSpoke]) << std::endl;
        }
        int primitiveCount = figure->getPrimitiveCount();
        int spokeIndex = 0;
        for(int i = 0; i < primitiveCount; ++i)
        {
            // The i-th primitive
            M3DPrimitive* currentPrimitive = figure->getPrimitivePtr(i);
            M3DQuadEndPrimitive* endPrimitive = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);
            double episilon = 1e-5;
            if(endPrimitive == NULL)
            {
                // it is an standard primitive
                M3DQuadPrimitive* primitive = dynamic_cast<M3DQuadPrimitive*>(currentPrimitive);
                double r0 = primitive->getR0();
                double tempR0 = r0 * exp(coeff[spokeIndex]);
                // if they are different significantly
                if(abs(tempR0 - r0) > episilon)
                {
                    primitive->setR0(tempR0);
                }
                double r1 = primitive->getR1();
                double tempR1 = r1 * exp(coeff[spokeIndex+1]);
                // if they are different significantly
                if(abs(tempR1 - r1) > episilon)
                {
                    primitive->setR1(tempR1);
                }

                spokeIndex += 2;

            }
            else
            {
                M3DQuadPrimitive* standardVersion = dynamic_cast<M3DQuadPrimitive*> (currentPrimitive);
                double r0 = standardVersion->getR0();

                double tempR0 = r0 * exp(coeff[spokeIndex]);
                // if they are different significantly
                if(abs(tempR0 - r0) > episilon)
                {
                    standardVersion->setR0(tempR0);
                }

                double r1 = standardVersion->getR1();
                double tempR1 = r1 * exp(coeff[spokeIndex+1]);
                // if they are different significantly
                if(abs(tempR1 - r1) > episilon)
                {
                    standardVersion->setR1(tempR1);
                }

                M3DQuadEndPrimitive* endVersion = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);
                double r = endVersion->getREnd();
                double tempR = r * exp(coeff[spokeIndex+2]);
                // if they are different significantly
                if(abs(tempR - r) > episilon)
                {
                    endVersion->setREnd(tempR);
                }

                spokeIndex += 3;
            }
        }
    }
    catch(std::exception& e)
    {
        std::cout << "[Error]updateFigure has an exception:"<< e.what() << std::endl;
    }
}
/* Compute the total entropy. */
double M3DSpokeLengthOptimizer::getObjectiveFunctionValue(const double *coeff, double w1, double w2)
{
    // input coeff from newuoa is the coeff of length. it could be positive or negative

    // 0. Update new length to each spoke
    updateFigure(coeff, mFigureIndex);

    double objFunctionValue = 0.0;
    // coeff are now lengths of spokes

    // 1. Image match
    // TODO: should read from config file
    double w_ImageMatch = 9.0;
    double w_sradPenalty = 1.0;
    double w_srepModelPenalty = 9.0;  // may be too large, was 10

    // 1.1 Interpolate spokes

    // 1.2 measure sum square of distance from implied boundary to expected image boundary
    SimilarityComputer similarityCompter;
    similarityCompter.setSrepModel(mSreps);
    similarityCompter.setTargetFigureIndex(mFigureIndex);
    similarityCompter.setTargetImage(mSignedDistanceImage);

    double imageMatch = 0.0;
    if(similarityCompter.compute(&imageMatch) == false)
    {
        std::cout << "[Error]Error encountered when compute similarity measure" << std::endl;
        return -999.0;
    }

    // 2. measure how far from regular srep model
    double sradPenalty = computeSradPenalty();

    // 3. spoke model penalty
    DistanceType distType = (enum DistanceType) (int) tuningWt(BpSpokeDistanceType);
    double spokeModelpenalty = mSreps->dist2FromObject(mSreps->loadedObject(),
                                                       mFigureIndex, distType);

    // 4. image normal match
//    double imageNormalMatch = 
    objFunctionValue = w_ImageMatch * imageMatch + w_sradPenalty * sradPenalty
        + w_srepModelPenalty * spokeModelpenalty;
    
    return objFunctionValue;
}
double M3DSpokeLengthOptimizer::operator() (double *coeff)
{
    double cost = 0.0;
    cost = this->getObjectiveFunctionValue(coeff, 15, 1);
    return cost;
}

void M3DSpokeLengthOptimizer::setImage(ImageDistanceMap* imageDistanceMap)
{
    mSignedDistanceImage = imageDistanceMap;
}

void M3DSpokeLengthOptimizer::setObject(M3DObject* sreps)
{
    mSreps = sreps;
}

/* Main entry of optimizer */
int M3DSpokeLengthOptimizer::perform(M3DObject* outputModel)
{
    if(mSreps == NULL || mSignedDistanceImage == NULL)
    {
        std::cout << "[Error]object or signed distance image is empty in perform function" << std::endl;
        return -1;
    }
    
    // optimization
    M3DFigure* figure = mSreps->getFigurePtr(mFigureIndex);
    int spokeCount = figure->getSpokeCount();
    double coeffOfLength[spokeCount];
    for(unsigned int i =0;i<spokeCount;i++){
        coeffOfLength[i] = 0.0;
    }
    const clock_t begin_time = clock();
    min_newuoa(spokeCount,coeffOfLength,*this,0.1, 1e-9, 3000);

    std::cout << "Finish optimization, use time: " << float(clock() - begin_time)/ CLOCKS_PER_SEC << std::endl;
    // coeffOfLength here is variable at optimum
    std::cout << "The optimum objective function value is: " << getObjectiveFunctionValue(coeffOfLength, 1, 1) << std::endl;
    M3DObjectFile objectFile;
    objectFile.write("out.m3d", *mSreps);

    delete[] coeffOfLength;

    return 0;

}
