/*********************************************************\
 * Mesh.cpp
 * Rohit
 * 2005/9/25
 *
 *  The implementation for the base Mesh class.
 *	Modifications:
 *
\********************************************************************************/

#define D_MESH
#define D_XFERLIST
#define D_PSEUDOSET
#define D_DIATOMGRID
#include "Shapedepend.h"
#include <stdio.h>

using namespace ThallCode;

/*********************************************************\
 * printvals()
\*********************************************************/
void Mesh::printvals(char *message)
{
	if (message != NULL)
		fprintf(stderr, "%s\n");

	fprintf(stderr, "Numverts = %d, Numpolygons = %d\n", nverts, npolys);

	fprintf(stderr, "Vertices are:\n");
	for (int vdex = 0; vdex < nverts; vdex++)
		v_list[vdex].print();

	fprintf(stderr, "Polygons are:\n");
	for (int pdex = 0; pdex < npolys; pdex++) {
		fprintf(stderr, "  poly %d = [", pdex);
		for (int pvdex = 0; pvdex < pi_list[pdex].Length(); pvdex++)
			fprintf(stderr, " %d", pi_list[pdex][pvdex]);
		fprintf(stderr, "]\n");
	}
}

/*********************************************************\
 * glRender()
 * render the tiles of the Mesh object
 * according to input thisstyle (enum in CCSubdivsurf.h)
\*********************************************************/
void Mesh::glRender(RenderStyle thisstyle)
{
	// Save current color
	glPushAttrib(GL_CURRENT_BIT);

	GLenum glstyle;
	switch (thisstyle) {
	case WIREFRAME:
		glstyle = GL_LINE_LOOP;
		break;
	case FILLEDPOLY:
		glstyle = GL_POLYGON;
		break;
	case VPOINTS:
		glstyle = GL_POINTS;
		break;
	case NO_RENDER:
		return;
		break;
	}

	glColor3d(1.0, 0.0, 0.0);
	glEnable(GL_LIGHTING);
	int polydex;
	for (polydex = 0; polydex < npolys; polydex++) {

		glBegin(glstyle);
		for (int vertdex = 0; vertdex < pi_list[polydex].Length() ; vertdex++) {
			Bpoint *thisvert = &v_list[pi_list[polydex][vertdex]];
			glNormal3dv(thisvert->norm);
			glVertex3dv(thisvert->pnt);
		}
		glEnd();
	}

	glDisable(GL_LIGHTING);

	glPopAttrib();
}

void Mesh::CopyXferList(Xferlist *thisList)
{
	//std::cout << "**********************************************" << std::endl;
	fig_grid->CopyXferlist(thisList);
}

void Mesh::CopyDiatomGrid(Diatomgrid *thisGrid)
{
	for(int cnt = 0; cnt < thisGrid->rows*thisGrid->cols; cnt++)
		fig_grid->dlist[cnt] = thisGrid->dlist[cnt];
}

void Mesh::UpdateVertexList()
{
	constructMesh0();
}

void Mesh::UpdateFaceList()
{
	computePolyList();
}

Mesh::~Mesh()
{						
	if (fig_grid != NULL) 
		delete fig_grid;
	if (v_list != NULL) 
		delete []v_list;
	if (pi_list != NULL)
		delete []pi_list;
}

