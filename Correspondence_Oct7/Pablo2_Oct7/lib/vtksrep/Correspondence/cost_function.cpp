#include "cost_function.h"


cost_function::cost_function()
{
    this->iterationCounter = 0;
}



bool cost_function::correctMove(double *arg, int length, double moveDis){
    // If one of the coeff generated is bigger than 1 or smaller than -1, throw away this tuple.
    for(unsigned int i = 0; i < length; i++){
        if(arg[i] > moveDis || arg[i] < -moveDis){
            std::cout << "Step bigger than 0.5, throw away!" << std::endl;
            return false;
        }
    }
    return true;
}


