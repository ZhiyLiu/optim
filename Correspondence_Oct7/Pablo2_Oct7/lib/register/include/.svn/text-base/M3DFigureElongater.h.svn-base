#ifndef M3D_FIGURE_ELONGATER_H
#define M3D_FIGURE_ELONGATER_H

#include "M3DQuadFigure.h"

enum M3DElongationMethod
{
    M3D_ELONGATE_ROW,
    M3D_ELONGATE_COLUMN,
    M3D_ELONGATE_LONG_AXIS,
    M3D_ELONGATE_SHORT_AXIS
};

class M3DFigureElongater
{
public:
    M3DFigureElongater(M3DQuadFigure * _figure = NULL);
    ~M3DFigureElongater();

    void initialize(M3DQuadFigure * _figure);

    void setElongationMethod(M3DElongationMethod m);

    void elongate(double value);
	void update();

private:
    void setPrimitiveToPreviousPrediction(M3DQuadFigure * targetFigure,
        int primitiveId, int prevPrimitiveId, double scale);
    void setPrimitiveToAveragePrediction(M3DQuadFigure * targetFigure,
        int primitiveId, int prevPrimitiveId, int neighborPrimitiveId,
        double scale, double elongation);
    void averagePrimitives(M3DPrimitive & avePrim, const M3DPrimitive & prim1,
        const M3DPrimitive & prim2, double weight);

    M3DQuadFigure * figure;
    M3DQuadFigure * referenceFigure;

    M3DElongationMethod elongationMethod;
    bool rowOrientedLinks;
};

#endif

