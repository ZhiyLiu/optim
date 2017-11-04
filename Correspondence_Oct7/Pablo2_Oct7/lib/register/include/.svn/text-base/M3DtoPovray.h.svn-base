#ifndef M3D_TO_POVRAY_H
#define M3D_TO_POVRAY_H

#include <fstream>
#include "M3DObject.h"
#include "M3DQuadFigure.h"

class M3DtoPovray
{
public:
	M3DtoPovray(M3DObject * _object = NULL);

    void setObject(M3DObject * _object) { object = _object; }
    void output(const char * filename);

private:
	M3DObject * object;
    std::ofstream out;

    void outputHeader();
    void outputFigure(M3DQuadFigure * figure);
	void outputPrimitive(M3DPrimitive * prim);

    void outputCylinder(const Vector3D & a, const Vector3D & b,
                        const char * width, const char * surface);
    void outputSphere(const Vector3D & pos, const char * width,
                      const char * surface);

    void outputVector(const Vector3D & v);

    void outputSurface(M3DQuadFigure * figure);
};

#endif

