#ifndef M3D_PRIMITIVE_RENDERER_H
#define M3D_PRIMITIVE_RENDERER_H

#include "M3DPrimitive.h"

extern const float DEFAULT_POINT_COLOR[3];
extern const float DEFAULT_SELECTED_POINT_COLOR[3];

class M3DPrimitiveRenderer
{
public:
	M3DPrimitiveRenderer();
	~M3DPrimitiveRenderer();

    void setPrimitive(M3DPrimitive * prim) { primitive = prim; }
    void setSelectionId(int id) { selectionId = id; }
    void setPointScaleFactor(double factor) { pointScaleFactor = factor; }

    void setPointColor(float r, float g, float b)
    {
        pointColor[0] = r;
        pointColor[1] = g;
        pointColor[2] = b;
    }

    void setSelectedPointColor(float r, float g, float b)
    {
        selectedPointColor[0] = r;
        selectedPointColor[1] = g;
        selectedPointColor[2] = b;
    }

    void setTipPointColor(float r, float g, float b)
    {
        tipPointColor[0] = r;
        tipPointColor[1] = g;
        tipPointColor[2] = b;
    }

	void setTipName(char *n) {tipName = n;}
	void setTipT(double t) {tipT = t;}
	void setTipV(double v) {tipV = v;}

	void render();
	void renderLandmark();

protected:
    void drawPoint();

	M3DPrimitive * primitive;

    int selectionId;
    double pointScaleFactor;

    float selectedPointColor[3];
    float pointColor[3];

    bool axesOn;	// Draw BPerpN and N
    bool fancyAxes;	// Use proportional B + eta extension
	int drawBAxes;	// 0 == B vectors off; 1 == All on; 2 == Crest only on
	float lineWidth;

	char *tipName;
	double tipT;	// which spoke the tip is on: -1=Y0 0=b 1=Y1 for quads;
	double tipV;	// which spoke the tip is on: needed for tubes.
    float tipPointColor[3];
};

#endif

