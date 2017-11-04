/* This is for shifting the crest spokes only.
 *
 * Liyun Tu
 * Apr 22, 2014
*/


#include "optimizationusingnewuoa.h"
#include "cost_function.h"
#include "newuoa.h"
#include <iostream>

//#include <sstream>

#include <exception>


//Initialize
optimizationusingnewuoa::optimizationusingnewuoa()
{    
}




/* wirteSrep: if set to true, will write out the moved s-reps.
 *            default false, not wirte the shifting sreps.
*/
optimizationusingnewuoa::optimizationusingnewuoa(const char* rootDir, vector< std::string > inputSreps, int interpolationLevel, int side){
    this->rootDir = rootDir;
    this->srepNum = inputSreps.size();

    toolsfunc tools;
    M3DQuadFigure* shiftingQuadFig; // as a template, storing a s-rep for shifting.
    // readin a arbtrary srep as a template, holding the current srep's info, the optimization iteratively modifing this srep.    
    shiftingQuadFig = tools.GetQuadFigure(inputSreps[0].c_str());

    int rowNum = shiftingQuadFig->getRowCount();
    int colNum = shiftingQuadFig->getColumnCount();

    int exteriorAtomNums = (rowNum+colNum)*2 - 4;
    int interiorAtomNums = rowNum*colNum - exteriorAtomNums;

    int spokeNum;

    if(side==2){ // crest regularity
        varsNum = exteriorAtomNums - 4;
        spokeNum = exteriorAtomNums;
    }
    else { // up or down regularity variable number.
        varsNum = exteriorAtomNums*1 - 4 + interiorAtomNums*2; // 46. variables for standard spokes (up or down spokes).
        spokeNum = exteriorAtomNums + interiorAtomNums; // the total atoms on each srep, one atom corresponding to one spoke.
    }

    std::cout<<"Message from optimizationusingnewuoa: Exterior atom number: "<< exteriorAtomNums<<std::endl;
    std::cout<<"Message from optimizationusingnewuoa: Interior atom number: "<< interiorAtomNums<<std::endl;
    std::cout<<"Message from optimizationusingnewuoa: Variable number: "<< varsNum<<std::endl;
    std::cout<<"Message from optimizationusingnewuoa: The input srep number: "<< inputSreps.size()<<std::endl;

    this->totalVars = inputSreps.size() * varsNum;    

    objectivefunction obj_instance(rootDir, side, varsNum, spokeNum, shiftingQuadFig, interpolationLevel, inputSreps);
    this->objFunc = obj_instance;
}



// --------------------------------------------Move all moveable atoms----------------------------------------------
/* The cost function was iterated automaticly by the optimizer, each time it change the values of coeff, we use the new coeff
 * to update the srep, and got entropy as the objective function. And then from the obj function result, the optimizer will reset one
 * element in coeff to a new value and repeat the process, untill find a better minimum solution.
 * coeff is a N*46 dimesion arrary. Store the 46 spoke's u or v of all sreps. We will split this array to get each 46 for one srep.
 * coeff_changing: is the variables generate by optimizer. For some purpose, we only move specific spokes, so the variables generated is less than
 * the variables when shifting all the spokes. So here we set the non-shifting spokes variables to 0.
*/
double optimizationusingnewuoa::cost(double *coeff){

    double w1 = 15; // this value can be observed by run this optimizer once and evaluate a value to balance geo and reg.
    double w2 = 1;

    return objFunc.objectiveFunction(coeff, w1, w2);
}


int optimizationusingnewuoa::callNEWUOA(const char * logFileName, bool initialOpt){

    //std::cout<<"this->totalVars---------"<<this->totalVars<<std::endl;
    double *changingUV = new double [this->totalVars];

    // for the first optimization, all the variables set to 0.
    if(initialOpt){
        for(unsigned int i =0;i<this->totalVars;i++){
            changingUV[i] = 0.0;
        }
    }
    else {
        string varFilePath = string(this->rootDir) + string("/vars.txt");
        ifstream ss(varFilePath.c_str());

        if (! ss.is_open()){
            std::cerr << "Msg from optimizationusingnewuoa::callNEWUOA: Unable to open the vars file \"" << varFilePath << "\"!" << std::endl;
            EXIT_FAILURE;
        }

        string token;

        int i = 0;
        while (ss >> token && i <this->totalVars) {
            changingUV[i] = atof(token.c_str());
            i++;
        }
    }

    std::ofstream out(logFileName, std::ios_base::out | std::ios_base::app);
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

    cost_function costFunc(this, this->totalVars);
    min_newuoa(this->totalVars,changingUV,costFunc,0.1, 0.000001, 30000);

    //std::cout<<", The objectFunction is: "<<cost(ary)<<std::endl;

    delete changingUV;

    //reset to standard output
    std::cout.rdbuf(coutbuf);

    return 0;
}

// -----------------------End------------------Move all moveable atoms--------------------------End-----------------




//// --------------------------------------------Move the given one atoms----------------------------------------------
//double optimizationusingnewuoa::cost(double *coeff){

//    double w1 = 4; // this value can be observed by run this optimizer once and evaluate a value to balance geo and reg.
//    double w2 = 1;

//    return objFunc.objectiveFunction(coeff, w1, w2);
//}


//int optimizationusingnewuoa::callNEWUOA(const char * logFileName, bool initialOpt){

//    //std::cout<<"this->totalVars---------"<<this->totalVars<<std::endl;
//    double *changingUV = new double [this->totalVars];

//    // for the first optimization, all the variables set to 0.
//    if(initialOpt){
//        for(unsigned int i =0;i<this->totalVars;i++){
//            changingUV[i] = 0.0;
//        }
//    }
//    else {
//        string varFilePath = string(this->rootDir) + string("/vars.txt");
//        ifstream ss(varFilePath.c_str());

//        if (! ss.is_open()){
//            std::cerr << "Msg from optimizationusingnewuoa::callNEWUOA: Unable to open the vars file \"" << varFilePath << "\"!" << std::endl;
//            EXIT_FAILURE;
//        }

//        string token;

//        int i = 0;
//        while (ss >> token && i <this->totalVars) {
//            changingUV[i] = atof(token.c_str());
//            i++;
//        }
//    }

//    std::ofstream out(logFileName, std::ios_base::out | std::ios_base::app);
//    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
//    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

//    cost_function costFunc(this, this->totalVars);
//    min_newuoa(this->totalVars,changingUV,costFunc,0.1, 0.00001, 1);

//    //std::cout<<", The objectFunction is: "<<cost(ary)<<std::endl;

//    delete changingUV;

//    //reset to standard output
//    std::cout.rdbuf(coutbuf);

//    return 0;
//}




// *********************************************************************************
// ---- Move one up spoke ----
// *********************************************************************************
//double optimizationusingnewuoa::cost(double *coeff){

//    double w1 = 20; // this value can be observed by run this optimizer once and evaluate a value to balance geo and reg.
//    double w2 = 1;

//    // set coeff_changing  the 17th crest spokes for demo.
//    double * wholeCoeff = new double [this->totalVars];

//    //set the value generated by the optimizer to the 17th crest spokes (its varable index is 14) for demo.
//    for(int i = 0; i < this->srepNum; i++) {
//        wholeCoeff[this->varsNum * i + 40] = coeff[i*2]; // 17th crest spoke only v change
//        wholeCoeff[this->varsNum * i + 41] = coeff[i*2+1]; // 17th crest spoke only v change
//    }

//    // Output all the variables for all the sreps in this movement.
//    /*cout<<"---------------------the wholeCoeff is : "<<endl;
//    for(unsigned int j =0; j < wholeCoeffLength; j++){
//        cout<<wholeCoeff[j]<<"  ";
//    }
//    cout<<endl;
//    cout<<"----wholeCoeff's size is: "<<wholeCoeffLength<<endl;*/

//    double geoEntropy = objFunc.objectiveFunction(wholeCoeff, w1, w2);

//    delete wholeCoeff;

//    return geoEntropy;
//}


//int optimizationusingnewuoa::callNEWUOA(const char * logFileName, bool initialOpt){
//    int changingVarNum = this->srepNum * 2; // Each srep only has 2 variables.

//    double *changingUV = new double [changingVarNum];

//    // for the first optimization, all the variables set to 0.
//    if(initialOpt){
//        for(unsigned int i = 0; i < changingVarNum; i++){
//            changingUV[i] = 0.0;
//        }
//    }
//    else {
//        string varFilePath = string(this->rootDir) + string("/vars.txt");
//        ifstream ss(varFilePath.c_str());

//        if (! ss.is_open()){
//            std::cerr << "Msg from optimizationusingnewuoa::callNEWUOA: Unable to open the vars file \"" << varFilePath << "\"!" << std::endl;
//            EXIT_FAILURE;
//        }

//        string token;

//        int i = 0;
//        while (ss >> token && i < changingVarNum) {
//            changingUV[i] = atof(token.c_str());
//            i++;
//        }
//    }

//    std::ofstream out(logFileName, std::ios_base::out | std::ios_base::app);
//    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
//    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

//    cost_function costFunc(this, changingVarNum);
//    min_newuoa(changingVarNum,changingUV,costFunc,0.1, 0.000001, 30000);

//    //std::cout<<", The objectFunction is: "<<cost(ary)<<std::endl;

//    delete changingUV;

//    //reset to standard output
//    std::cout.rdbuf(coutbuf);

//    return 0;
//}


