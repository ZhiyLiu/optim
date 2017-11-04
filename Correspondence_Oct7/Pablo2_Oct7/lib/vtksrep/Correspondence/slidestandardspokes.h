#ifndef SLIDESTANDARDSPOKES_H
#define SLIDESTANDARDSPOKES_H


#include "toolsfunc.h"
#include "calregularityfeatures.h"
#include "regularityentropy.h"

#include "time.h"



class slidestandardspokes
{
public:
    slidestandardspokes();
    slidestandardspokes(int side, int interpolationLevel);

    double computeRegularityEntropy(M3DQuadFigure* quadFig, DoubleVec subUVs);

    //update the primitive's r, u, x, y, z and its corresponding spoke's deltau, deltav.
    void updateSpokeInfo(M3DQuadFigure *quadFig, double * varArraySrep, int side);

    void setMovedSpokeInfoToPrimitive(M3DQuadFigure *quadFig, double newU, double newV, double primitiveIndexU, double primitiveIndexV);

    //modify m3d file's deltaU, deltaV
    double moveStandardSpokes(M3DQuadFigure* quadFig, double *vars, DoubleVec subUVs);



private:

    int side;
    int interpolationLevel;
};

#endif // SLIDESTANDARDSPOKES_H
