#ifndef SLIDECRESTSPOKES_H
#define SLIDECRESTSPOKES_H


#include "regularityentropy.h"
#include "calcrestregularityfeatures.h"



class slidecrestspokes
{
public:
    slidecrestspokes();

    slidecrestspokes(int interpolationLevel);

    double moveCrestSpokes(M3DQuadFigure *quadFig, vtkSmartPointer<vtkSRep> srepfig, double *vars, DoubleVec subTVs);
    void shiftCrestSpokeAndSetToPrimitive(vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> shiftCrestSpokes,
                                          M3DQuadFigure *quadFig, double u, double v, int varIndex, double var);
    void updateCrestSpokeInfo(vtkSmartPointer<vtkSRep> srepfig, M3DQuadFigure *quadFig, double * varArraySrep);

    double computeRegularityEntropyOfCrest(M3DQuadFigure *quadFig, DoubleVec subTVs);


    double calculateSpokeLength(vtkSRep::VNLType newSpoke);
    Vector3D calculateUnitSpokeDir(vtkSRep::VNLType newSpoke);


private:    

    int interpolationLevel;
};

#endif // SLIDECRESTSPOKES_H
