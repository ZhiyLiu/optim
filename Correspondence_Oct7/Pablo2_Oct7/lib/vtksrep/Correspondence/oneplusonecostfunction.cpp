/* This is the costfunction of the itkOnePlusOneEvolutionaryOptimizer.
 * GetValue is the iteration part, its read only, derived from a virtual function. We cannot change any value in this function.
 * We even can't define objectivefunction objFunc as a member variable, set value in constructor and call its member function in GetValue, why???
 *
*/


#include "oneplusonecostfunction.h"


//using namespace std;

// Define a galoble object, to make the fixed parameters only initialize once, before the iteration.


int iterationCounter;


oneplusonecostfunction::oneplusonecostfunction()
{
    iterationCounter = 0;
}


oneplusonecostfunction::oneplusonecostfunction(const char* rootDir, vector< std::string > inputSreps, int interpolationLevel, int side){
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

    std::cout<<"Message from movespokes::movespokes: Exterior atom number: "<< exteriorAtomNums<<std::endl;
    std::cout<<"Message from movespokes::movespokes: Interior atom number: "<< interiorAtomNums<<std::endl;
    std::cout<<"Message from movespokes::movespokes: Variable number: "<< varsNum<<std::endl;
    std::cout<<"Message from movespokes::movespokes: The input srep number: "<< inputSreps.size()<<std::endl;

    this->totalVars = inputSreps.size() * varsNum;

    objectivefunction obj_instance(rootDir, side, varsNum, spokeNum, shiftingQuadFig, interpolationLevel, inputSreps);
    this->objFunc = obj_instance;

    iterationCounter = 0;
}



// *********************************************************************************
// ---- Move all moveable atoms ----
// *********************************************************************************
oneplusonecostfunction::MeasureType oneplusonecostfunction::GetValue( const ParametersType & parameters ) const  {
    iterationCounter++;
    std::cout << "The " << iterationCounter << "th iteration. " << std::endl;

    double *changingUV = new double [this->totalVars];

    // If one of the coeff generated is bigger than 1 or smaller than -1, throw away this tuple.
    for(unsigned int i = 0; i < this->totalVars; i++){
        if(parameters[i] > 1 || parameters[i] < -1){
            std::cout<<"The delta u or v shouldn't bigger than 1, this tuple of coeffes were throw away!"<<std::endl;
            return 10000.1;
        }

        changingUV[i] = parameters[i];
    }

    // output the variables every 100 iteration
    //if(iterationCounter%200 == 0){
        std::cout<<"Variables: ";
        for(unsigned int i = 0; i < this->totalVars; i++){
            std::cout << changingUV[i] << "  ";
        }
    //}
    std::cout << std::endl;

    double w1 = 1; // this value can be observed by run this optimizer one time and evaluate a value to balance geo and reg.
    double w2 = 1;
    MeasureType measure = this->objFunc.objectiveFunction(changingUV, w1, w2);

    delete changingUV;

    return measure;
}



// *********************************************************************************
// ---- Move one crest spoke ----
// *********************************************************************************
///* parameters[i] are the changing variables gotten from the optimizer. */
//oneplusonecostfunction::MeasureType oneplusonecostfunction::GetValue( const ParametersType & parameters ) const  {
//    iterationCounter++;
//    std::cout << "The " << iterationCounter << "th iteration. " << std::endl;

//    int changingVarNum = this->srepNum;

//    // set coeff_changing  the 17th crest spokes for demo.
//    double * wholeCoeff = new double [this->totalVars];

//    // If one of the coeff generated is bigger than 1 or smaller than -1, throw away this tuple.
//    for(unsigned int i = 0; i < this->srepNum; i++){
//        if(parameters[i] > 0.5 || parameters[i] < -0.5){
//            std::cout<<"The delta u or v shouldn't bigger than 0.5, this tuple of coeffes were throw away!"<<std::endl;
//            return 10000.1;
//        }

//        wholeCoeff[this->varsNum * i + 14] = parameters[i]; // 17th crest spoke only v change, its variable index is 14.
//    }

//    // output the variables every 1000 iteration
//    if(iterationCounter%1000 == 0){
//        std::cout<<"Variables: ";
//        for(unsigned int i = 0; i < changingVarNum; i++){
//            std::cout << parameters[i] << "  ";
//        }
//    }
//    std::cout << std::endl;

//    double w1 = 30; // this value can be observed by run this optimizer one time and evaluate a value to balance geo and reg.
//    double w2 = 1;

//    // Output all the variables for all the sreps in this movement.
//    /*cout<<"---------------------the wholeCoeff is : "<<endl;
//    for(unsigned int j =0; j < wholeCoeffLength; j++){
//        cout<<wholeCoeff[j]<<"  ";
//    }
//    cout<<endl;
//    cout<<"----wholeCoeff's size is: "<<wholeCoeffLength<<endl;*/

//    MeasureType measure = objFunc.objectiveFunction(wholeCoeff, w1, w2);

//    delete wholeCoeff;

//    return measure;
//}




// *********************************************************************************
// ---- Move one up spoke ----
// *********************************************************************************
/* parameters[i] are the changing variables gotten from the optimizer. */
//oneplusonecostfunction::MeasureType oneplusonecostfunction::GetValue( const ParametersType & parameters ) const  {
//    iterationCounter++;
//    std::cout << "The " << iterationCounter << "th iteration. " << std::endl;

//    int changingVarNum = this->srepNum * 2;

//    // set coeff_changing  the 17th crest spokes for demo.
//    double * wholeCoeff = new double [this->totalVars];

//    // If one of the coeff generated is bigger than 1 or smaller than -1, throw away this tuple.
//    for(unsigned int i = 0; i < this->srepNum; i++){
//        if(parameters[i] > 0.5 || parameters[i] < -0.5){
//            std::cout<<"The delta u or v shouldn't bigger than 0.5, this tuple of coeffes were throw away!"<<std::endl;
//            return 10000.1;
//        }

//        wholeCoeff[this->varsNum * i + 40] = parameters[i*2];
//        wholeCoeff[this->varsNum * i + 41] = parameters[i*2+1];
//    }

//    // output the variables every 1000 iteration

//    std::cout<<"Variables: ";
//    for(unsigned int i = 0; i < changingVarNum; i++){
//        std::cout << parameters[i] << "  ";
//    }

//    std::cout << std::endl;

//    double w1 = 20; // this value can be observed by run this optimizer one time and evaluate a value to balance geo and reg.
//    double w2 = 1;

//    // Output all the variables for all the sreps in this movement.
//    /*cout<<"---------------------the wholeCoeff is : "<<endl;
//    for(unsigned int j =0; j < wholeCoeffLength; j++){
//        cout<<wholeCoeff[j]<<"  ";
//    }
//    cout<<endl;
//    cout<<"----wholeCoeff's size is: "<<wholeCoeffLength<<endl;*/

//    MeasureType measure = objFunc.objectiveFunction(wholeCoeff, w1, w2);

//    delete wholeCoeff;

//    return measure;
//}

