/* This is an wrapper of optimizer using NEWUOA
 *
 * Zhiyuan Liu
 * 2017.11
 */
#include "M3DNewuoaOptimizer.h"

#include <iostream>
#include <exception>
#include "newoa.h"
#include "toolsfunc.h"
#include "M3DQuadFigure.h"
#include "movespokes.h"

M3DNewuaOptimizer::M3DNewuaOptimizer(const char* rootDir, const std::vector< std::string > &inputSreps, int interpolationLevel, int spokeType) :
    mInterpolationLevel(interpolationLevel),
    mSpokeType(spokeType),
    mInputSreps(inputSreps),
    mRootDir(rootDir){
    }

void M3DNewuaOptimizer::initialize() {
    toolsfunc tools;
    M3DQuadFigure* shiftingQuadFig; // as a template, storing a s-rep for shifting.
    // readin a arbtrary srep as a template, holding the current srep's info, the optimization iteratively modifing this srep.
    shiftingQuadFig = tools.GetQuadFigure(inputSreps[0].c_str());

    int rowNum = shiftingQuadFig->getRowCount();
    int colNum = shiftingQuadFig->getColumnCount();

    int exteriorAtomNums = (rowNum+colNum)*2 - 4;
    int interiorAtomNums = rowNum*colNum - exteriorAtomNums;

    int spokeNum;
    int dimensionNum;

    if(side==2){ // crest regularity
        dimensionNum = exteriorAtomNums - 4;
        spokeNum = exteriorAtomNums;
    }
    else { // up or down regularity variable number.
        dimensionNum = exteriorAtomNums*1 - 4 + interiorAtomNums*2; // 46. variables for standard spokes (up or down spokes).
        spokeNum = exteriorAtomNums + interiorAtomNums; // the total atoms on each srep, one atom corresponding to one spoke.
    }

    std::cout<<"Message from optimizationusingnewuoa: Exterior atom number: "<< exteriorAtomNums<<std::endl;
    std::cout<<"Message from optimizationusingnewuoa: Interior atom number: "<< interiorAtomNums<<std::endl;
    std::cout<<"Message from optimizationusingnewuoa: Variable number: "<< varsNum<<std::endl;
    std::cout<<"Message from optimizationusingnewuoa: The input srep number: "<< mInputSreps.size()<<std::endl;

    mTotalDimensions = mInputSreps.size() * dimensionNum;

}

// --------------------------------------------Move all moveable atoms----------------------------------------------
/* The cost function was iterated automaticly by the optimizer, each time it change the values of coeff, we use the new coeff
 * to update the srep, and got entropy as the objective function. And then from the obj function result, the optimizer will reset one
 * element in coeff to a new value and repeat the process, untill find a better minimum solution.
 * coeff is a N*46 dimesion arrary. Store the 46 spoke's u or v of all sreps. We will split this array to get each 46 for one srep.
 * coeff_changing: is the variables generate by optimizer. For some purpose, we only move specific spokes, so the variables generated is less than
 * the variables when shifting all the spokes. So here we set the non-shifting spokes variables to 0.
*/
double M3DNewuaOptimizer::getCost(const double coeff) {
    double w1 = 15;
    double w2 = 1;

    return getObjectiveFunctionValue(coeff, w1, w2);
}

/* Compute the total entropy. */
double M3DNewuaOptimizer::getObjectiveFunctionValue(const double wholeCoeff, double w1, double w2) const {
        //Step 1: Move spokes, compute the regularity entropy and save geometry matrix to .txt.
    movespokes mSpokes(rootDir, spokeNum, this->quadFigList.size(), side, this->interpolationLevel);
    double regEntropy = mSpokes.calculateRegEntropy(wholeCoeff,this->shiftingQuadFig, subVs, varsNum, quadFigList, srepfigList);

    toolsfunc tools;

    //Step 2: Call matlab script to compute geometry entropy.
    //matlab -nodesktop -nodisplay -nojvm -nosplash -r "dataDir='/work/ltu/WorkSpace/Sep_3_newFeature/downSide/temp_data/';geoEntropy"
    string s = string("matlab -nodesktop -nodisplay -nojvm -nosplash -singleCompThread -r ") + '"' + string("dataDir='") + tools.connectStr(rootDir, "/temp_data/");

    string script_str;
    if(this->side == 0)
        script_str = s + string("';upSide") + '"';
    else if(this->side == 1)
        script_str = s + string("';downSide") + '"';
    else
        script_str = s + string("';crestSide") + '"';

    int returnValue1 = system(script_str.c_str());
    if (returnValue1 != 0 ) std::cout << "something has happened in Matlab, following value returned: " << returnValue1 << std::endl;


    //Step 3: Read the geo entropy from .txt
    string geoEntropyResultPath = tools.connectStr(tools.connectStr(rootDir, "/temp_data/"), "geoEntropyResult.txt");
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
double M3DNewuaOptimizer::operator() (double stepSize) {
    
        this->mIterationCounter++;
        std::cout<<"The "<<this->mIterationCounter<<"th iteration. "<<std::endl;

        double cost = 0.0;
        if(isCorrectMove(stepSize, this->mTotalDimensions, 1)){
            // call the cost function defined in optimizationusingnewuoa
            cost = this->getCost(stepSize);
        }
        else {
            cost = 10000.1;
        }

        // output the variables every 1000 iteration
        if(this->mIterationCounter%500 == 0){
            std::cout<<"Variables: ";
            for(unsigned int i = 0; i < this->mTotalDimensions; i++){
                std::cout << arg[i] << "  ";
            }
        }
        std::cout << std::endl;

        return cost;
}

int M3DNewuaOptimizer::perform(const char* logFileName, bool initialOpt) {
    using namespace std;
    //std::cout<<"this->mTotalDiimensions---------"<<this->mTotalDiimensions<<std::endl;
    double *changingUV = new double [this->mTotalDimensions];

    // for the first optimization, all the variables set to 0.
    if(initialOpt){
        for(unsigned int i =0;i<this->mTotalDiimensions;i++){
            changingUV[i] = 0.0;
        }
    }
    else {
        string varFilePath = string(this->mRootDir) + string("/vars.txt");
        ifstream ss(varFilePath.c_str());

        if (! ss.is_open()){
            std::cerr << "Msg from optimizationusingnewuoa::callNEWUOA: Unable to open the vars file \"" << varFilePath << "\"!" << std::endl;
            return 1;
        }

        string token;

        int i = 0;
        while (ss >> token && i <this->mTotalDiimensions) {
            changingUV[i] = atof(token.c_str());
            i++;
        }
    }

    std::ofstream out(logFileName, std::ios_base::out | std::ios_base::app);
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

    min_newuoa(this->mTotalDiimensions,changingUV,this,0.1, 0.000001, 30000);

    //std::cout<<", The objectFunction is: "<<cost(ary)<<std::endl;

    delete changingUV;

    //reset to standard output
    std::cout.rdbuf(coutbuf);

    return 0;

}