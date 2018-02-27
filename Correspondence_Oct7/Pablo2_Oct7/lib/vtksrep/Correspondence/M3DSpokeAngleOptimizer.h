/* The purpose of this class is to optimize spoke angle by NEWUOA
 *
 * Zhiyuan Liu
 * 2018.01
 */
#ifndef _M3DSPOKEANGLEOPTIMIZER_H_
#define _M3DSPOKEANGLEOPTIMIZER_H_

#include <vector>
#include <string>
#include "OptimizerBase.h"
#include "vtksrep.h"
#include <vtkSmartPointer.h>

class M3DQuadFigure;
class ImageDistanceMap;

// rotation struct.
// In euler angle: a-angle, b/c/d-vector
// In quaternion: a-scalar, b/c/d-coeff of imagenary part
struct Rotation
{
    double a;
    double b;
    double c;
    double d;
};
class M3DSpokeAngleOptimizer
{
public:
    M3DSpokeAngleOptimizer();

    M3DSpokeAngleOptimizer(M3DObject* sreps);

    // Description: Return a cost value computed according to pre-defined entropy minimization objective function
    // Parameter: @coeff: input coefficients in each iterariton
    // @w1/w2: weight for geometric entropy and regularity
    double getObjectiveFunctionValue(const double *wholeCoeff, double w1, double w2);

    // Description: Override operator (). Required by min_newuoa.
    // Parameter: @coeff: the pointer to coefficients
    double operator () (double *coeff);

    // major entry to invoke this optimizer
    // Each object might contain multiple figures
    // Each figure is an srep. But currently each object has only one srep
    // return: error code
    // output: s-rep model at optimum
    int perform(M3DObject* outputModel);

    // Set object want to optimize
    // Each object might contain multiple figures
    // Each figure is an srep. But currently each object has only one srep
    void setObject(M3DObject* sreps);

    // Description: Set signed distance image related to srep model
    // Input: signed distance image
    void setImage(ImageDistanceMap* imageDistanceMap);
private:
    double computeSradPenalty();

     // update a figure after each optimization 
    void updateFigure(const double *coeff, int figureId);

    void rotateByQuaternion(double a, double b, double c, double d, Vector3D& target);

    void updateSpokeDir(const double *coeff, int spokeIndex, Vector3D& dir);

    void recoverSrep();

private:
    M3DObject*                          mSreps;
    std::vector<Vector3D>               mOrigDirs;
    //  all the rotation relating to each spoke. indexed by primitive index as well as the side.
    // For different sides: 0 - top; 1- bottom; 2 - end(if any)
    std::vector<Rotation>               mRotations;
    // for each spoke ,there are 4 elements stored in this set;
    std::vector<double>                 mQuaternionSet;
    ImageDistanceMap*                   mSignedDistanceImage;
    int                                 mFigureIndex = 0; // The figure user want to optimize, currently have only one
};

#endif
