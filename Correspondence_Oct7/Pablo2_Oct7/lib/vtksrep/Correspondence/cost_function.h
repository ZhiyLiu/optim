#ifndef COST_FUNCTION_H
#define COST_FUNCTION_H

#include "optimizationusingnewuoa.h"

class cost_function
{
public:
    cost_function();

    cost_function(optimizationusingnewuoa *instance, int totalVars){
        m_instance = instance;

        this->iterationCounter = 0;

        this->argLength = totalVars;
    }

    bool correctMove(double *arg, int length, double moveDis);

    double operator () (double *arg){
        this->iterationCounter++;
        std::cout<<"The "<<this->iterationCounter<<"th iteration. "<<std::endl;

        double cost = 0.0;
        if(correctMove(arg, this->argLength, 1)){
            // call the cost function defined in optimizationusingnewuoa
            cost = m_instance->cost(arg);
        }
        else {
            cost = 10000.1;
        }

        // output the variables every 1000 iteration
        if(this->iterationCounter%500 == 0){
            std::cout<<"Variables: ";
            for(unsigned int i = 0; i < this->argLength; i++){
                std::cout << arg[i] << "  ";
            }
        }
        std::cout << std::endl;

        return cost;
    }


private:
    optimizationusingnewuoa *m_instance;

    int iterationCounter;
    int argLength; // the coefficients number
};

#endif // COST_FUNCTION_H
