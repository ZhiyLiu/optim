#ifdef WIN32
# include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include <stdio.h>
#include <iostream>
#include "Vector3D.h"
#include "TileSet.h"
#include "TileSetRenderer.h"



//#define DEBUG


using namespace std;

TileSetRenderer::TileSetRenderer()
{
	nlists = 0;
	clear();
}

TileSetRenderer::~TileSetRenderer()
{
	if (nlists > 0) {
		glDeleteLists(initialDispList, nlists);
#ifdef DEBUG
		cout << "Discarding " << nlists << " display lists\n";
#endif
	}
}

void TileSetRenderer::clear()
{
#ifdef DEBUG
		cout << "TileSetRenderer::clear()\n";
#endif
	if (nlists > 0) {
		glDeleteLists(initialDispList, nlists);
#ifdef DEBUG
		cout << "Discarding " << nlists << " display lists\n";
#endif
		nlists = 0;
	}
    numTiles = 0;
	dflt_color[0] = 0.0f;	// Same as in control_parms_defaults.cpp
	dflt_color[1] = 0.0f;
	dflt_color[2] = 1.0f;
	pDfltColor[0] = dflt_color;
}

void TileSetRenderer::compileFigure(int figureId)
{
	int lo, hi;

	lo = 0;
	for (int i = 0; i < figureId; i++)
		lo += tileCounts[i];
	hi = lo + tileCounts[figureId] - 1;
	compilePart(lo, hi);
}

void TileSetRenderer::compileFigures(int numFigures)
{
	int lo, hi;

	lo = 0;
	for (int i = 0; i < numFigures; i++) {
		hi = lo + tileCounts[i] - 1;
		compilePart(lo, hi);
		lo = hi + 1;
	}
}

// Lo and hi are polygon indexes
void TileSetRenderer::compilePart(int lo, int hi)
{
    int i;

	if (quads) {
		lo *= 4;
		hi = 4*hi + 3;
		glBegin(GL_QUADS);
		for(i = lo; i < hi; i += 4)
		{
			Vector3D v1, v2, norm;

			v1 = coords[i + 1] - coords[i];
			v2 = coords[i + 3] - coords[i];
			norm = v1.cross(v2);
			if(norm.norm() != 0)
				norm.normalize();
			glNormal3d(norm.getX(), norm.getY(), norm.getZ());
			glVertex3d(coords[i].getX(), coords[i].getY(), coords[i].getZ());

			v1 = coords[i] - coords[i + 1];
			v2 = coords[i + 2] - coords[i + 1];
			norm = v2.cross(v1);
			if(norm.norm() != 0)
				norm.normalize();
			glNormal3d(norm.getX(), norm.getY(), norm.getZ());
			glVertex3d(coords[i + 1].getX(), coords[i + 1].getY(), coords[i + 1].getZ());

			v1 = coords[i + 1] - coords[i + 2];
			v2 = coords[i + 3] - coords[i + 2];
			norm = v2.cross(v1);
			if(norm.norm() != 0)
				norm.normalize();
			glNormal3d(norm.getX(), norm.getY(), norm.getZ());
			glVertex3d(coords[i + 2].getX(), coords[i + 2].getY(), coords[i + 2].getZ());

			v1 = coords[i] - coords[i + 3];
			v2 = coords[i + 2] - coords[i + 3];
			norm = v1.cross(v2);
			if(norm.norm() != 0)
				norm.normalize();
			glNormal3d(norm.getX(), norm.getY(), norm.getZ());
			glVertex3d(coords[i + 3].getX(), coords[i + 3].getY(), coords[i + 3].getZ());
		}
	}
    else {
		lo *= 3;
		hi = 3*hi + 2;
		glBegin(GL_TRIANGLES);
		for(i = lo; i < hi; i += 3)
		{
			Vector3D v1, v2, norm;

			v1 = coords[i + 1] - coords[i];
			v2 = coords[i + 2] - coords[i];
			norm = v1.cross(v2);
			if(norm.norm() != 0)
				norm.normalize();
			glNormal3d(norm.getX(), norm.getY(), norm.getZ());
			glVertex3d(coords[i].getX(), coords[i].getY(), coords[i].getZ());

			v1 = coords[i] - coords[i + 1];
			v2 = coords[i + 2] - coords[i + 1];
			norm = v2.cross(v1);
			if(norm.norm() != 0)
				norm.normalize();
			glNormal3d(norm.getX(), norm.getY(), norm.getZ());
			glVertex3d(coords[i + 1].getX(), coords[i + 1].getY(), coords[i + 1].getZ());

			v1 = coords[i] - coords[i + 2];
			v2 = coords[i + 1] - coords[i + 2];
			norm = v1.cross(v2);
			if(norm.norm() != 0)
				norm.normalize();
			glNormal3d(norm.getX(), norm.getY(), norm.getZ());
			glVertex3d(coords[i + 2].getX(), coords[i + 2].getY(), coords[i + 2].getZ());
		}
	}
    glEnd();
}

void TileSetRenderer::render(int index, float alpha, bool useDefault, bool solid)
{
#ifdef DEBUG
//      cout << "TileSetRenderer::render(" << index << ", " << alpha << ", "
//			<< useDefault << ", " << solid << ')' << endl;
#endif

    glEnable(GL_CULL_FACE);
    glEnable(GL_LIGHTING);

    if (! solid) {
		glDepthMask(GL_FALSE);
		glEnable(GL_BLEND);						// Does not seem to matter
		glBlendFunc(GL_SRC_ALPHA, GL_DST_ALPHA);// Does not seem to matter
	}
    else
        glDisable(GL_BLEND);

//	if (colors != NULL)
//		cout << "Tileset color = (" << colors[index][0] << ", " << colors[index][1]
//			<< ", " << colors[index][2] << ")\n";
//	cout << "Tileset default color = (" << 1.0f - dflt_color[0] << ", "
//		<< 1.0f - dflt_color[1] << ", " << 1.0f - dflt_color[2] << ")\n";

    glCullFace(GL_FRONT);
	if (useDefault || index >= ncolors)
		glColor4f(1.0f - dflt_color[0], 1.0f - dflt_color[1], 1.0f - dflt_color[2], alpha);
	else
		glColor4f(1.0f - colors[index][0], 1.0f - colors[index][1],
			1.0f - colors[index][2], alpha);
    glColorMaterial(GL_BACK, GL_AMBIENT_AND_DIFFUSE);
	glCallList(initialDispList + index);

    glCullFace(GL_BACK);
	if (useDefault || index >= ncolors)
		glColor4f(dflt_color[0], dflt_color[1], dflt_color[2], alpha*alpha);
	else
		glColor4f(colors[index][0], colors[index][1], colors[index][2], alpha);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glCallList(initialDispList + index);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    if (! solid) {
		glDepthMask(GL_TRUE);
		glDisable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
    else
        glEnable(GL_BLEND);

	// why's this here? it messes up clip planes! -gst
    //glDisable(GL_CLIP_PLANE0);

    glDisable(GL_CULL_FACE);
    glDisable(GL_LIGHTING);
}

void TileSetRenderer::render(float alpha, bool useDefault, bool solid)
{
#ifdef DEBUG
//      cout << "TileSetRenderer::render()" << endl;
#endif
	if (numTiles <= 0)
		return;

	for (int i = 0; i < nlists; i++)
		render(i, alpha, useDefault, solid);
}

void TileSetRenderer::compileDispList(TileSet * ts, int numFigColors,
	const float ** figureColors)
{
	int numFigures;

#ifdef DEBUG
	cout << "TileSetRenderer::compileDispList(" << hex << "0x"
		<< ts << dec << ")" << endl;
#endif
	if (ts == NULL)
		return;

	numFigures = ts->figureCount();
	numTiles = ts->tileCount();
	coords = ts->getTiles();
	tileCounts = ts->getTileCounts();
	quads = ts->haveQuads();

	if (nlists > 0) {
		glDeleteLists(initialDispList, nlists);
#ifdef DEBUG
		cout << "Discarding " << nlists << " display lists\n";
#endif
		nlists = 0;
	}

	// numFigures is the number of parts in the BYU file
	if (numFigures > 0) {
		// The number of display lists equals the number of colors or figures,
		// which ever is less.  Thus, extra colors are ignored.  On the other
		// hand, if there are more figures than colors, the default color will
		// be used for rendering the extra figures.  This ignores the
		// possibility that some colors may be duplicated.
		if (numFigColors <= 0)
			nlists = 1;
		else if (numFigColors >= numFigures)
			nlists = numFigures;
		else	// numFigColors < numFigures
			nlists = numFigColors + 1;
		// The use of different colors requires multiple display lists
		initialDispList = glGenLists(nlists);
#ifdef DEBUG
		cout << "Creating " << nlists << " display lists\n";
#endif

		colors = figureColors;	// May be NULL
		ncolors = numFigColors;

		for (int n = 0; n < numFigures; n++) {
#ifdef DEBUG
			if (n < nlists)
				glNewList(initialDispList + n, GL_COMPILE_AND_EXECUTE);
#else
			if (n < nlists)
				glNewList(initialDispList + n, GL_COMPILE);
#endif
			compileFigure(n);
#ifdef DEBUG
			GLenum e = glGetError();
			if (e != GL_NO_ERROR) cout << "ERROR: GL err #" << e << endl;
#endif
			if (n < nlists - 1)
				glEndList();
		}
		glEndList();
	}
	if (colors == NULL)
		colors = (const float **) pDfltColor;
}

void TileSetRenderer::setDefaultColor(const float * color)
{
	dflt_color[0] = color[0];
	dflt_color[1] = color[1];
	dflt_color[2] = color[2];
#ifdef DEBUG
	cout << "Default tileset color = (" << color[0] << ", " << color[1]
		<< ", " << color[2] << ")\n";
#endif
}


