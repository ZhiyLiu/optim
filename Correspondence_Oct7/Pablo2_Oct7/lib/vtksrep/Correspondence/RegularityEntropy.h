/* The purpose of this class is to provide
 * calculator of regularity entropy
 * Zhiyuan Liu
 * 2017.11
 */
#ifndef REGULARITYENTROPY_H
#define REGULARITYENTROPY_H

class RegularityEntropy {
public:
    RegularityEntropy();
    double calculate();
    void setCoefficients(double *coeff);

private:
    double *mCoeff;
};


#endif