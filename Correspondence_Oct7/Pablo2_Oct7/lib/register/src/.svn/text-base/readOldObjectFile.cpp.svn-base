#include <stdarg.h>
#include "M3DObject.h"
#include "M3DQuadFigure.h"
#include "Registry.h"

M3DPrimitive * readPrimitive(Registry &reg, char * regStr, ...)
{
	char newStr[1024];

	va_list val;
	va_start(val, regStr);
	vsprintf(newStr, regStr, val);
	va_end(val);

	strcat(newStr, ".%s");

	double x, y, z, r, theta;
	Quat q;
	bool selected;

	M3DPrimitive * prim;

	x =  reg.getDoubleValue(newStr, 0.0, "x");
	y =  reg.getDoubleValue(newStr, 0.0, "y");
	z =  reg.getDoubleValue(newStr, 0.0, "z");
	r =  reg.getDoubleValue(newStr, 0.0, "r");
	q.setX(reg.getDoubleValue(newStr, 0.0, "qx"));
	q.setY(reg.getDoubleValue(newStr, 0.0, "qy"));
	q.setZ(reg.getDoubleValue(newStr, 0.0, "qz"));
	q.setW(reg.getDoubleValue(newStr, 0.0, "qw"));
	theta = reg.getDoubleValue(newStr, 0.0, "theta") * R_DEGREES_TO_RADIANS;
	selected = reg.getBooleanValue(newStr, 0, "selected");

    if(theta > R_HALF_PI)
    {
        // Set theta to supplementary angle
        theta = R_PI - theta;

        // Rotate the primitive 180 degrees about N
        Quat rotQ;
        rotQ.setAxisAngle(Vector3D(0.0, 1.0, 0.0), R_PI);

        q = q * rotQ;
    }

	prim = new M3DQuadPrimitive(x, y, z, r, q, theta);

	if(selected)
		prim->select();

	return prim;
}

M3DObject * readOldObjectFile(const char * filename)
{
    Registry reg;
    M3DObject * object;
    M3DQuadFigure * fig;
    M3DPrimitive * prim;

    int figcount;
    int slicecount;
    int primcount;
    int i, j, k;

    float color[3];

    double elongation;

	try {
		reg.readFromFile(filename);
	}
	catch (RException excp) {
		excp.print(std::cout);
		return NULL;
	}

    object = new M3DObject;

    figcount = reg.getIntValue("model.figureCount", 0);
	for(i = 0; i < figcount; i++)
	{
        slicecount = reg.getIntValue("model.figure[%d].sliceCount", 0, i);
        primcount = reg.getIntValue("model.figure[%d].slice[0].primitiveCount", 0, i);

        fig = new M3DQuadFigure(slicecount, primcount);

        color[0] = reg.getDoubleValue("model.figure[%d].color.red", DEFAULT_FIGURE_COLOR[0], i);
        color[1] = reg.getDoubleValue("model.figure[%d].color.green", DEFAULT_FIGURE_COLOR[1], i);
        color[2] = reg.getDoubleValue("model.figure[%d].color.blue", DEFAULT_FIGURE_COLOR[2], i);

        fig->setColor(&(color[0]));

		for(j = 0; j < slicecount; j++)
		{
			primcount = reg.getIntValue("model.figure[%d].slice[%d].primitiveCount", 0,i,j);
	        for(k = 0; k < primcount; k++)
			{
				prim = readPrimitive(reg, "model.figure[%d].slice[%d].primitive[%d]", i, j, k);
				if(j == 0 || j == slicecount - 1 || k == 0 || k == primcount - 1)
                {
                    elongation = reg.getDoubleValue("model.figure[%d].slice[%d].primitive[%d].alpha", 1.0, i, j, k);
					// FIXME: This code might break upon execution - rrs
					std::cout << __FILE__ << ":" << __LINE__ << ": code might break." << std::endl;
					prim = new M3DQuadEndPrimitive(prim->getX(), prim->getR(),
						prim->getQ(),prim->getTheta(), elongation);
                }

				fig->setPrimitivePtr(j, k, prim);
			}

			fig->setModified(true);
		}

        object->addFigure(fig);
	}

    return object;
}

