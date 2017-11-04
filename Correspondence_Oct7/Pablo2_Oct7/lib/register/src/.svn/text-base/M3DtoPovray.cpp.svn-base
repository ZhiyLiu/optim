#include "M3DtoPovray.h"
#define D_DIATOM
#define D_XFERLIST
#define D_POINTLIST_SERVER2
#include "Shapedepend.h"

using namespace std;

const char * STD_SURFACE = "pigment{ stdColor }\nfinish{ stdFinish }\nno_shadow";
const char * CONNECTOR_SURFACE = "pigment{ connectorColor }\nfinish{ stdFinish }\nno_shadow";
const char * ATOM_VEC_WIDTH = "atomVecWidth";
const char * ATOM_CENTER_WIDTH = "atomCenterWidth";
const char * CONNECTOR_WIDTH = "connectorWidth";

M3DtoPovray::M3DtoPovray(M3DObject * _object)
{
    object = _object;
}

void M3DtoPovray::output(const char * filename)
{
    int figureCount,
        figureId;
    M3DQuadFigure * figure;

    if(object == NULL)
        return;

    out.open(filename);

    outputHeader();

    figureCount = object->getFigureCount();
    for(figureId = 0; figureId < figureCount; figureId++)
    {
      // FIXME: not really defined for M3DTubeFigure
        figure = dynamic_cast<M3DQuadFigure *>(object->getFigurePtr(figureId));
        if(figure == NULL)
            return;
        outputFigure(figure);
    }

    out.close();
}

void M3DtoPovray::outputHeader()
{
    out << "#declare stdColor = color red .6 green .6 blue .6;" << endl
        << "#declare stdFinish = finish{ambient 0.4 diffuse 0.8 phong 0.74};" << endl
        << "#declare atomVecWidth = 0.0075;" << endl
        << "#declare atomCenterWidth = 0.025;" << endl
        << "#declare connectorWidth = 0.005;" << endl
        << "#declare connectorColor = color red .4 green .4 blue .4;" << endl
        << "global_settings { assumed_gamma 2.2 }" << endl
        << "background{ color rgb<1,1,1> }" << endl
        << "camera {"<< endl
        << "    orthographic" << endl
        << "    right -2.66 * x" << endl
        << "    up 2 * y" << endl
        << "    location  <0.5, 0.5, 1000>" << endl
        << "    look_at   <0.5, 0.5, 0.5>" << endl
        << "}" << endl
        << "light_source { <0, 30, 25> color red 1 green 1 blue 1 }" << endl;
}

void M3DtoPovray::outputFigure(M3DQuadFigure * figure)
{
    M3DPrimitive * prim;
    int i, j,
        primitiveCount;

    if(figure == NULL)
        return;

    primitiveCount = figure->getPrimitiveCount();
    for(i = 0; i < primitiveCount; i++)
    {
        out << "// Start Primitive #" << i << endl;
        outputPrimitive(figure->getPrimitivePtr(i));
        out << "// End Primitive #" << i << endl;
    }

    int numRows, numColumns;
    numRows = figure->getRowCount();
    numColumns = figure->getColumnCount();
    Vector3D a, b;
    for(i = 0; i < numRows; i++)
    {
        for(j = 0; j < numColumns - 1; j++)
        {
            prim = figure->getPrimitivePtr(i, j);
            a = prim->getX();
            prim = figure->getPrimitivePtr(i, j + 1);
            b = prim->getX();

            out << "// Connector (" << i << ", " << j << ") -> (" << i << ", " << j + 1 << ")" << endl;
            outputCylinder(a, b, CONNECTOR_WIDTH, CONNECTOR_SURFACE);
        }
    }

    for(j = 0; j < numColumns; j++)
    {
        for(i = 0; i < numRows - 1; i++)
        {
            prim = figure->getPrimitivePtr(i, j);
            a = prim->getX();
            prim = figure->getPrimitivePtr(i + 1, j);
            b = prim->getX();

            out << "// Connector (" << i << ", " << j << ") -> (" << i + 1 << ", " << j << ")" << endl;
            outputCylinder(a, b, CONNECTOR_WIDTH, CONNECTOR_SURFACE);
        }
    }

    outputSurface(figure);
}

void M3DtoPovray::outputPrimitive(M3DPrimitive * prim)
{
    Vector3D x;
    Vector3D y0, y1;
    Vector3D exB;

    x = prim->getX();
    exB = prim->getExtendedB();
    y0 = prim->getY0();
    y1 = prim->getY1();

    outputCylinder(x, x + y0, ATOM_VEC_WIDTH, STD_SURFACE);
    outputCylinder(x, x + y1, ATOM_VEC_WIDTH, STD_SURFACE);
    outputSphere(x, ATOM_CENTER_WIDTH, STD_SURFACE);

    if(prim->type() == M3D_END_PRIMITIVE)
        outputCylinder(x, x + exB, ATOM_VEC_WIDTH, STD_SURFACE);
}

void M3DtoPovray::outputCylinder(const Vector3D & a, const Vector3D & b, 
                                 const char * width, const char * surface)
{
    out << "cylinder{" << flush;
    outputVector(a);
    out << ", " << flush;
    outputVector(b);
    out << ", " << width << endl
        << surface << endl
        << "}" << endl;
}

void M3DtoPovray::outputSphere(const Vector3D & pos, const char * width,
                               const char * surface)
{
    out << "sphere{" << flush;
    outputVector(pos);
    out << ", " << width << endl << surface << endl << "}" << endl;
}

void M3DtoPovray::outputSurface(M3DQuadFigure * figure)
{
    ThallCode::Pointlist_server2 pList;
    Bpoint * tileset;
    int numtiles;
    int i;

    Xferlist * xferList = convertM3DtoXfer(figure);
    pList.init(xferList);
    pList.ComputeSubdivBoundaryTiles(3);

    pList.subdivtileinfo(&numtiles, &tileset);

    out << "mesh2 {" << endl;

    // Print out vertices
    out << "vertex_vectors {" << endl;
    out << 4*numtiles << endl;
    for(i = 0; i < 4*numtiles; i++)
    {
        outputVector(Vector3D(tileset[i].pnt));
        out << "," << endl;
    }
    out << "}" << endl;

    // Print out normals
    out << "normal_vectors {" << endl;
    out << 4*numtiles << endl;
    for(i = 0; i < 4*numtiles; i++)
    {
        outputVector(Vector3D(tileset[i].norm));
        out << "," << endl;
    }
    out << "}" << endl;

    // Print out faces (as vertex indices)
    out << "face_indices {" << endl;
    out << 2*numtiles << endl;
    for(i = 0; i < 4*numtiles; i+=4)
    {
        outputVector(Vector3D(i+3, i, i+1));
        out << "," << endl;
        outputVector(Vector3D(i+3, i+1, i+2));
        out << "," << endl;
    }
    out << "}" << endl;

    // Print out indices into normals
    out << "normal_indices {" << endl;
    out << 2*numtiles << endl;
    for(i = 0; i < 4*numtiles; i+=4)
    {
        outputVector(Vector3D(i+3, i, i+1));
        out << "," << endl;
        outputVector(Vector3D(i+3, i+1, i+2));
        out << "," << endl;
    }
    out << "}" << endl;

    out << "texture {" << endl
        << "    pigment { color rgbf <1,1,1,0.7> }" << endl
        << "    finish { phong 0.4 phong_size 20 }" << endl
        << "}" << endl
        << "no_shadow" << endl
        << "}" << endl;
}

void M3DtoPovray::outputVector(const Vector3D & v)
{
    out << "<" << v.getX() << ", " << v.getY() << ", " << v.getZ() << ">"
        << flush;
}

