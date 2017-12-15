/* This is an wrapper of optimizer using NEWUOA
 *
 * Zhiyuan Liu
 * 2017.11
 */
#include "M3DNewuoaOptimizer.h"

#include <iostream>
#include <exception>
#include "newuoa.h"
#include "toolsfunc.h"
#include "M3DQuadFigure.h"
#include "movespokes.h"
#include "M3DAtomPredictorQuad.h"

M3DNewuoaOptimizer::M3DNewuoaOptimizer(char* mVarFileDir, const std::vector< std::string > &inputSreps, int interpolationLevel, int spokeType) :
    mInterpolationLevel(interpolationLevel),
    mSpokeType(spokeType),
    mInputSreps(inputSreps),
    mVarFileDir(mVarFileDir){
}

void M3DNewuoaOptimizer::initialize() {
    // step 1. compute total number of atoms
    toolsfunc tools;
    M3DQuadFigure* mShiftingQuadFig; // as a template, storing a s-rep for shifting.
    // readin a arbtrary srep as a template, holding the current srep's info, the optimization iteratively modifing this srep.
    mShiftingQuadFig = tools.GetQuadFigure(mInputSreps[0].c_str());

    int rowNum = mShiftingQuadFig->getRowCount();
    int colNum = mShiftingQuadFig->getColumnCount();

    int exteriorAtomNums = (rowNum+colNum)*2 - 4;
    int interiorAtomNums = rowNum*colNum - exteriorAtomNums;

    int mTotalDimensions;
    int dimensionNum;

    if(mSpokeType==2){ // crest regularity
        dimensionNum = exteriorAtomNums - 4;
        mTotalDimensions = exteriorAtomNums;
    }
    else { // up or down regularity variable number.
        dimensionNum = exteriorAtomNums*1 - 4 + interiorAtomNums*2; // 46. coeff for standard spokes (up or down spokes).
        mTotalDimensions = exteriorAtomNums + interiorAtomNums; // the total atoms on each srep, one atom corresponding to one spoke.
    }

    std::cout<<"Message from optimizationusingnewuoa: Exterior atom number: "<< exteriorAtomNums<<std::endl;
    std::cout<<"Message from optimizationusingnewuoa: Interior atom number: "<< interiorAtomNums<<std::endl;
    std::cout<<"Message from optimizationusingnewuoa: Variable number: "<< dimensionNum<<std::endl;
    std::cout<<"Message from optimizationusingnewuoa: The input srep number: "<< mInputSreps.size()<<std::endl;

    mTotalDimensions = mInputSreps.size() * dimensionNum;

    //step 2. Initialized the srep quadFig list
    for(unsigned int i = 0; i < mInputSreps.size(); i++){
        //Initialize mQuadFigList, which holding the original sreps, the changing during each iteration will based on this.
        std::cout<<"---------Reading:"<<mInputSreps[i]<<std::endl;
        M3DQuadFigure* qFig = tools.GetQuadFigure(mInputSreps[i].c_str());
        this->mQuadFigList.push_back(qFig);

        // Initialize srepfig.
        vtkSmartPointer<vtkReadSRep> readsrep = vtkSmartPointer<vtkReadSRep>::New();
        readsrep->SetFileName(mInputSreps[i].c_str());
        readsrep->Update();
        vtkSmartPointer<vtkSRep> srepfig = readsrep->GetOutput();
        this->mSrepFigList.push_back(srepfig);
    }

    // step 3. Compute subdivision of quad
    uvmap mp;
    // holding the 0~1 subdivision of the line.
    this->mSubdivisions = mp.getUVSubdivision(mInterpolationLevel);

}

// --------------------------------------------Move all moveable atoms----------------------------------------------
/* The cost function was iterated automaticly by the optimizer, each time it change the values of coeff, we use the new coeff
 * to update the srep, and got entropy as the objective function. And then from the obj function result, the optimizer will reset one
 * element in coeff to a new value and repeat the process, untill find a better minimum solution.
 * coeff is a N*46 dimesion arrary. Store the 46 spoke's u or v of all sreps. We will split this array to get each 46 for one srep.
 * coeff_changing: is the coeff generate by optimizer. For some purpose, we only move specific spokes, so the coeff generated is less than
 * the coeff when shifting all the spokes. So here we set the non-shifting spokes coeff to 0.
*/
double M3DNewuoaOptimizer::getCost(const double *coeff) {
    double w1 = 15;
    double w2 = 1;

    return getObjectiveFunctionValue(coeff, w1, w2);
}

double M3DNewuoaOptimizer::computeSradPenalty()
{
    double sradPenalty = 0.0;
    M3DAtomPredictorQuad *atomPredictor = new M3DAtomPredictorQuad();
//    sradPenalty = atomPredictor->getRSradPenalty();
    return sradPenalty;
}

double M3DNewuoaOptimizer::computeSRepImageMatch(double dilationFactor)
{
    double match = 0.0;
    toolsfunc tools;
    // step 1: find relevant spokes
    std::vector<M3DSpoke> allSpokes;
    interpolateSRep(&allSpokes);
    
   // step 2: compute the sum of square distances as the match value
    std::vector<Vector3D> surfacePoints;
    for(int i = 0; i < allSpokes.size(); ++i)
    {
        surfacePoints.push_back(allSpokes[i].getB());
    }

    int pointsNum = surfacePoints.size();
    if(pointsNum == 0)
    {
        // ERROR
        std::cout << "[ERROR]Implied boundary after interpolation is empty in method: M3DNewuoaOptimizer::computeSRepImageMatch!" << std::endl;
        return -999.0;
    }

    double max_dist = 0.0, min_dist = +1e20, avg_dist = 0.0;


    for(int i = 0; i < pointsNum; ++i)
    {
        // TODO: DistanceMap should be singleton class that provide utilities
        //double distance = fabs(binaryDistanceMap->getDistance(surfacePoints[i]) - dilationFactor);

        double distance = 0.0;
        if(distance > max_dist)
        {
            max_dist = distance;
        }

        if(distance < min_dist)
        {
            min_dist = distance;
        }

        match += distance * distance;
        avg_dist += distance;
    }

    
    avg_dist /= pointsNum;
    match /= pointsNum;

    MatchUtility::max_dist = max_dist;
	MatchUtility::min_dist = min_dist;
	MatchUtility::nPoints = pointsNum;
	MatchUtility::avg_dist = avg_dist;
	MatchUtility::avg_dist_squared = match;
    return match;
}

void M3DNewuoaOptimizer::interpolateSRep(std::vector<M3DSpoke> *outputSpokes)
{
    // interpolation level from config file
    int level = 2;
    // TODO: interpolator should be singleton tool manufactured by a factory to provide utilities
}

/* Update spoke properties after each optimization */
void M3DNewuoaOptimizer::updateSpokes(const double *coeff)
{
    try
    {
        std::vector<M3DSpoke> allSpokes;
        for(int i = 0; i < allSpokes.size(); ++i)
        {
            allSpokes[i].setR(exp(coeff[i]));
        }
    }
    catch(std::exception& e)
    {
        std::cout << "[Error]updateSpokes has an exception" << std::endl;
    }
}
/* Compute the total entropy. */
double M3DNewuoaOptimizer::getObjectiveFunctionValue(const double *coeff, double w1, double w2)
{
    // input coeff from newuoa is the coeff of length. it could be positive or negative
    updateSpokes(coeff);
    
    double objFunctionValue = 0.0;
    // coeff are now lengths of spokes
    // 0. Update new length to each spoke
    
    // 1. Image match
    // TODO: should read from config file
    double w_ImageMatch = 9999;
    double w_sradPenalty = 1.0;
    double dilationFactor = 0.0;

    // measure sum square of distance from implied boundary to expected image boundary
    double imageMatch = computeSRepImageMatch(dilationFactor);

    // measure how far from regular srep model
    double sradPenalty = computeSradPenalty();

    objFunctionValue = w_ImageMatch * imageMatch + w_sradPenalty * sradPenalty;
    return objFunctionValue;
}
double M3DNewuoaOptimizer::operator() (double *coeff) {
        this->mIterationCounter++;
        std::cout<<"The "<<this->mIterationCounter<<"th iteration. "<<std::endl;

        double cost = 0.0;
        cost = this->getCost(coeff);
        // output the coeff every 1000 iteration
        if(this->mIterationCounter%500 == 0){
            std::cout<<"Variables: ";
            for(unsigned int i = 0; i < this->mTotalDimensions; i++){
                std::cout << coeff[i] << "  ";
            }
        }
        std::cout << std::endl;

        return cost;
}

/* Main entry of optimizer */
int M3DNewuoaOptimizer::perform(const std::string& srep) {
    // initialize
    
    //std::cout<<"this->mTotalDimensions---------"<<this->mTotalDimensions<<std::endl;
    double *coeffOfLength = new double[this->mTotalDimensions];  // optimize length's coefficients instead of length itself, to avoid negative length output by newuoa

    for(unsigned int i =0;i<this->mTotalDimensions;i++){
        coeffOfLength[i] = 1.0;
    }
    min_newuoa(this->mTotalDimensions,coeffOfLength,*this,0.1, 0.000001, 30000);

    delete[] coeffOfLength;
    return 0;

}
