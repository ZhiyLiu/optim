#include "newuoa.h"


class ObjectiveFunction
{
public:
    double operator() (double *coeff)
        {
            double objFunctionValue = 200.0;
            if(coeff[70] <= 1) objFunctionValue += 5-coeff[70];
            else objFunctionValue += coeff[70] + 3;
            return objFunctionValue;
        }
};
int main()
{
    double x[74];
    for(int i = 0; i < 74; ++i)
    {
        x[i] = 2;
    }
    ObjectiveFunction f;
    min_newuoa(74, x, f, 0.5, 0.0001, 300);
    return 0;
}
