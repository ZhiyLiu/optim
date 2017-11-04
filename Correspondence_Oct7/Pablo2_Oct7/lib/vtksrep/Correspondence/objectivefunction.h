#ifndef OBJECTIVEFUNCTION_H
#define OBJECTIVEFUNCTION_H


#include <vector>
#include "M3DQuadFigure.h"
#include "P3DControl.h"
#include "movespokes.h"
#include "toolsfunc.h"

#include <vtkSmartPointer.h>


class objectivefunction
{
public:
    objectivefunction();

    objectivefunction(const char* rootDir, int side, int varsNum, int spokeNum, M3DQuadFigure* shiftingQuadFig, int interpolationLevel,
                      vector< std::string > inputSreps);

    double objectiveFunction(const double * wholeCoeff, double w1, double w2) const;


//    vector<double> generateShiftingVariables(double * changingCoeff, int index1, int index2) const;

    int getIterationCount(int iteCounter) const;


private:
    int varsNum; // variables number on each srep. For up or down entropy, 46; for crest entropy 28.

    vector<M3DQuadFigure *> quadFigList; //holding the input sreps. Reading in before iteration.
    vector<vtkSmartPointer<vtkSRep> > srepfigList; // hoding srepfig, used for crest spoke interpolate.

//    int iteratorCount; // record the times the objective function was called (count the optimizer loop times).

    int side;

    const char* rootDir;

    DoubleVec subVs;

    int interpolationLevel;
    M3DQuadFigure* shiftingQuadFig;

    int spokeNum;

};

#endif // OBJECTIVEFUNCTION_H
