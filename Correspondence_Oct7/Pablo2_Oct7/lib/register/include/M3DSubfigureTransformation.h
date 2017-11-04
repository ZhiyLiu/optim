#ifndef M3D_SUBFIGURE_TRANSFORMATION_H
#define M3D_SUBFIGURE_TRANSFORMATION_H

#include "M3DObject.h"
#include "M3DFigureTreeNode.h"
#include "M3DQuadFigure.h"

extern const double LOG_BASE;
extern const double INVERSED_LOG_BASE;


class M3DSubfigureTransformation
{
public:
    M3DSubfigureTransformation(M3DObject * _object = NULL,
        M3DFigureTreeNode * _figureTreeNode = NULL);

    ~M3DSubfigureTransformation();

    void initStep();
	void init(M3DObject * _object, M3DFigureTreeNode * _figureTreeNode);
	void reInit();

    void translateStep(double deltaU, double deltaV);
    void rotateStep(double angle);
    void widenStep(double value);
    void elongateStep(double value);
    void hingeStep(double angle);

	void translate(double deltaU, double deltaV);
    void rotate(double angle);
    void widen(double value);
    void elongate(double value);
    void hinge(double angle);
	void updateSubfigure(M3DQuadFigure * oldParrentFigure);

	void getSurfacePointAndNormal(ThallCode::Pointlist_server2 *pList2, 
		Vector3D & point, Vector3D & normal, double u, double v, double t);

//	void translateHingeAtom(double theDeltaU, double theDeltaV, int hingeIndex);

    bool getInterpolatedPrimitive(M3DPrimitive & prim, M3DQuadFigure * figure,
                                  double u, double v, double t);

private:
    void checkForValidLinks();

    void setPrimitiveToPreviousPrediction(M3DQuadFigure * figure,
        int primitiveId, int prevPrimitiveId, double scale);

    void setPrimitiveToAveragePrediction(M3DQuadFigure * figure,
        int primitiveId, int prevPrimitiveId, int neighborPrimitiveId,
        double scale, double elongation);

    void averagePrimitives(M3DPrimitive & avePrim, const M3DPrimitive & prim1,
        const M3DPrimitive & prim2, double weight);

    M3DObject * object;
    M3DFigureTreeNode * figureTreeNode;

    M3DFigure * referenceFigure;
	M3DFigure * referenceFigureOld;

    double lastScaleValue;
    double lastElongateValue;

    // True if links are oriented along a row, false if column orientation
    bool rowOrientedLinks;

    bool invalidLinks;

	double maxParentU, maxParentV;
	double factorT2U, factorT2V;

};

#endif

