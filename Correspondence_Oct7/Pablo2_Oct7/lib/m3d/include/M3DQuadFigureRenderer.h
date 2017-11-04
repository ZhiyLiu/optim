#ifndef M3D_QUAD_FIGURE_RENDERER_H
#define M3D_QUAD_FIGURE_RENDERER_H

#include "M3DFigureRenderer.h"

class M3DQuadFigureRenderer : public M3DFigureRenderer
{
public:
    M3DQuadFigureRenderer(M3DQuadFigure * fig = NULL);
    ~M3DQuadFigureRenderer();

//	void setFigurePtr(M3DQuadFigure * fig);

    virtual void render();
    virtual void renderLandmarks();
};

#endif

