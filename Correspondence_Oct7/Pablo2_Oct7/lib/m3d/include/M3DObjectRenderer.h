#ifndef M3D_OBJECT_RENDERER_H
#define M3D_OBJECT_RENDERER_H

#include "M3DObject.h"
#include "M3DFigureRenderer.h"

#undef drawString

class M3DObjectRenderer
{

public:

    M3DObjectRenderer(M3DObject * obj = NULL);
    ~M3DObjectRenderer();

    void setObjectPtr(M3DObject * obj);

    void setPrimitiveRenderer(M3DPrimitiveRenderer * pr);
    void setMarkedPrimitiveRenderer(M3DPrimitiveRenderer * mpr)
    {
        markedPrimitiveRenderer = mpr;
    }

    void setMarkedPrimitiveId(int id) { markedPrimitiveId = id; }

    static void setMeshVisibility(bool toggle) { mesh_visibility = toggle; }
	static void setMeshConnectorsType(bool toggle) { mesh_connectors_type = toggle; }
	static void setMeshConnectorsLineWidth(int width) { mesh_connectors_line_width = (float) width; }
	static void setMeshConnectorsColor(float color[3]) {
		mesh_connectors_color[0] = color[0];
		mesh_connectors_color[1] = color[1];
		mesh_connectors_color[2] = color[2];
	}
    static void setAtomVectorVisibility(bool toggle) { atom_vector_visibility = toggle; }
    static void setAtomVisibility(bool toggle) { atom_visibility = toggle; }
	static void setExtraAtomVectors(bool toggle) { extra_atom_vectors = toggle; }
	static void setAtomVectorsType(bool toggle) { atom_vectors_type = toggle; }
	static void setBVectorsType(int toggle) { b_vectors_type = toggle; }
	static void setAtomVectorsLineWidth(float width) { atom_vectors_line_width = width; }
	static void setUseFigureColor(bool toggle) { use_figure_color = toggle; }
	static void setConstraintArrows(bool toggle) { draw_constraint_arrows = toggle; }
    static void setLandmarkVisibility(bool toggle) { landmark_visibility = toggle; }

    void setCurrentScale(double scale);

    void render();
    void renderLandmarks();
	void renderFigureNames();

	static bool meshVisibility() { return mesh_visibility; }
	static bool meshConnectorsType() { return mesh_connectors_type; }
	static float meshConnectorsLineWidth() { return mesh_connectors_line_width; }
	static float * meshConnectorsColor() { return mesh_connectors_color; }
	static bool atomVectorVisibility() { return atom_vector_visibility; }
	static bool atomVisibility() { return atom_visibility; }
	static bool showExtraAtomVectors() { return extra_atom_vectors; }
	static bool atomVectorsType() { return atom_vectors_type; }
	static int atomB_VectorsType() { return b_vectors_type; }
	static float atomVectorsLineWidth() { return atom_vectors_line_width; }
	static bool useFigureColor() { return use_figure_color; }
	static bool constraintArrows() { return draw_constraint_arrows; }

	static bool landmark_visibility;
		// public: defined at object level; used in primitive rendering
	static void drawString(const char * s);
	static int stringDrawingWidth(const char * s);
	
protected:

    M3DObject * object;

    double currentScale;

    int markedPrimitiveId;
    M3DPrimitiveRenderer * markedPrimitiveRenderer;

	static bool mesh_visibility;
	static bool mesh_connectors_type;
	static float mesh_connectors_line_width;
	static float mesh_connectors_color[3];
	static bool atom_vector_visibility;
	static bool atom_visibility;
	static bool extra_atom_vectors;
	static bool atom_vectors_type;
	static int b_vectors_type;
	static float atom_vectors_line_width;
	static bool use_figure_color;
	static bool draw_constraint_arrows;

    std::vector<M3DFigureRenderer *> figureRenderers;

	void renderConstraintArrows();
};


#endif

