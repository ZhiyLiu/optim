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

/* Compute the total entropy. */
double M3DNewuoaOptimizer::getObjectiveFunctionValue(const double *coeff, double w1, double w2) const {
        //Step 1: Move spokes, compute the regularity entropy and save geometry matrix to .txt.
    movespokes mSpokes(mVarFileDir, mTotalDimensions, this->mQuadFigList.size(), mSpokeType, this->mInterpolationLevel);
    double regEntropy = mSpokes.calculateRegEntropy(coeff,this->mShiftingQuadFig, mSubdivisions, mTotalDimensions, mQuadFigList, mSrepFigList);

    toolsfunc tools;

    //Step 2: Call matlab script to compute geometry entropy.
    //matlab -nodesktop -nodisplay -nojvm -nosplash -r "dataDir='/work/ltu/WorkSpace/Sep_3_newFeature/downSide/temp_data/';geoEntropy"
    string s = string("matlab -nodesktop -nodisplay -nojvm -nosplash -singleCompThread -r ") + '"' + string("dataDir='") + tools.connectStr(mVarFileDir, "/temp_data/");

    string script_str;
    if(this->mSpokeType == 0)
        script_str = s + string("';upSide") + '"';
    else if(this->mSpokeType == 1)
        script_str = s + string("';downSide") + '"';
    else
        script_str = s + string("';crestSide") + '"';

    int returnValue1 = system(script_str.c_str());
    if (returnValue1 != 0 ) std::cout << "something has happened in Matlab, following value returned: " << returnValue1 << std::endl;


    //Step 3: Read the geo entropy from .txt
    string geoEntropyResultPath = tools.connectStr(tools.connectStr(mVarFileDir, "/temp_data/"), "geoEntropyResult.txt");
    ifstream myfile(geoEntropyResultPath.c_str());
    string line;
    double geoEntropy = 0;
    if (myfile.is_open()) {
        while ( getline (myfile,line) )
        {
            geoEntropy = atof(line.c_str());
            std::cout <<"----geoEntropy: "<<geoEntropy<<std::endl;
        }
        myfile.close();
    }

    else std::cout << "Unable to open file: " << geoEntropyResultPath << std::endl;

    double objectFuc = w1*geoEntropy - w2*regEntropy;
    std::cout <<"----objective function: "<<objectFuc<<std::endl;

    return objectFuc;

}
double M3DNewuoaOptimizer::operator() (double *coeff) {
    
        this->mIterationCounter++;
        std::cout<<"The "<<this->mIterationCounter<<"th iteration. "<<std::endl;

        double cost = 0.0;
        if(isCorrectMove(coeff, this->mTotalDimensions, 1)){
            // call the cost function defined in optimizationusingnewuoa
            cost = this->getCost(coeff);
        }
        else {
            cost = 10000.1;
        }

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

int M3DNewuoaOptimizer::perform(const char* logFileName, bool initialOpt) {
    using namespace std;
    //std::cout<<"this->mTotalDimensions---------"<<this->mTotalDimensions<<std::endl;
    double *changingUV = new double [this->mTotalDimensions];

    // for the first optimization, all the coeff set to 0.
    if(initialOpt){
        for(unsigned int i =0;i<this->mTotalDimensions;i++){
            changingUV[i] = 0.0;
        }
    }
    else {
        string varFilePath = string(this->mVarFileDir) + string("/vars.txt");
        ifstream ss(varFilePath.c_str());

        if (! ss.is_open()){
            std::cerr << "Msg from optimizationusingnewuoa::callower: Unable to open the vars file \"" << varFilePath << "\"!" << std::endl;
            return 1;
        }

        string token;

        int i = 0;
        while (ss >> token && i <this->mTotalDimensions) {
            changingUV[i] = atof(token.c_str());
            i++;
        }
    }

    std::ofstream out(logFileName, std::ios_base::out | std::ios_base::app);
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

    min_newuoa(this->mTotalDimensions,changingUV,*this,0.1, 0.000001, 30000);

    //std::cout<<", The objectFunction is: "<<cost(ary)<<std::endl;

    delete changingUV;

    //reset to standard output
    std::cout.rdbuf(coutbuf);

    return 0;

}

bool M3DNewuoaOptimizer::isCorrectMove( double *coeff, int length, double moveDis) {
    // If one of the coeff generated is bigger than 1 or smaller than -1, throw away this tuple.
    for(unsigned int i = 0; i < length; i++){
        if(coeff[i] > moveDis || coeff[i] < -moveDis){
            std::cout << "Step bigger than 0.5, throw away!" << std::endl;
            return false;
        }
    }

    return true;
}