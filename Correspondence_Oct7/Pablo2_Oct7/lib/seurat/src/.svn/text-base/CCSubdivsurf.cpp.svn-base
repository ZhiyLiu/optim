
/*
2004/05/16
	Updates - ::Interpolate(...) to support the subdivision of tube m-rep mesh
				extra codes to deal with the new uvt space of the tube
*/
/****************************************************************************/
/*																			*/
/*  	File	:  CCSubdivsurf.cpp											*/
/*																			*/
/*	Description:  class functions for Catmull-Clark subdivision class		*/
/*																			*/
/*		This is meant to replace Subdivsurf.h., ultimately, but for now,	*/
/*		uses Subdivsurf output to Pointlist_server2 to set adjacency info	*/
/*		for quads and vertices for the initial mesh.						*/
/*																			*/
/*	Project :  Seurat														*/
/*																			*/
/*	Author  :  A. Thall														*/
/*																			*/
/*	Date	:  5. July 2002													*/
/*																			*/
/*	Modifications:															*/
/****************************************************************************/
#define D_CCSUBDIVSURF
#define D_XFERLIST
#define D_PSEUDOSET
#include "Shapedepend.h"
#include <stdio.h>
#include <math.h>
#include <map>
#include <float.h>
#include <time.h>
//#include <sys/timeb.h>
#include <assert.h>

#include <stdlib.h>

using namespace ThallCode;

#define QIONG_DEBUG_LINALG

#define QIONG_DEBUG_20030510
#define QIONG_DEBUG_20030425

//#define QIONG_DEBUG_02232003

// This will cause normal-interpolation to be used at level-2 and higher subdivision
const bool NORMINTERP = true;

//#define DEBUG_QIONG_BLENDING

#ifdef DEBUG_QIONG_BLENDING
int dPNum;
int polyList[200];
#endif

// 02052003
// to testify the correctness of the inverse algorithm on arbitrary quad mesh!
// it proves to be well!
//#define DEBUG_QIONG_TESTING_INVERSE_ALGORITHM_ON_QUAD_MESH

// 02152003
// to try to negatify the displacement of the vertex when calculating the IIS mesh
// what i did is to negatify the displacement for the points with valence of 4, so
// this is only to test simple quadmesh, where i inserted a vertex in the middle of
// each quad, unfortunately this DOES NOT work!
//#define DEBUG_QIONG_TESTING_INVERSE_DISPLACEMENT_IISMESH

// 02162003
// to try to mirror the displacement of the points inserted in the middle of each quad
//#define DEBUG_QIONG_TESTING_MIRROR

// 02172003
// work in progress
// to do: to try my own method
// 02172003
#define QIONG_DEBUG_BLENDING_02192003

//#define QIONG_DEBUG_BLENDING

//#define QIONG_DEBUG_BLENDED_SURFACE_NORMALS
#ifdef QIONG_DEBUG_BLENDED_SURFACE_NORMALS
#define QIONG_DEBUG_BLENDED_VERTEX_NORMALS
/*
DbVector3 Cross(DbVector3 v0, DbVector3 v1)
{
	v0.cross(v1);
	return v0;
}
*/
DbVector3 CalculateNormal(CCVertex **thisVert)
{
	DbVector3 normal;
	DbVector3 v0, v1;
	v0=thisVert[0]->limitloc-thisVert[1]->limitloc;
	v1=thisVert[2]->limitloc-thisVert[1]->limitloc;
	v0.selfnormalize(); //v0=v0.normalize();
	v1.selfnormalize(); //v1=v1.normalize();
	//v0=v0.cross(v1);
	return (v0.cross(v1)).normalize();//Cross(v0, v1);
}
//void CCSubdivsurf::CalculateNormalsFromVertices(int meshLvl, CCMesh *thisMesh)
void CalculateNormalsFromVertices(int meshLvl, CCMesh *thisMesh)
{
	int i; //, j;
	//int pDex, vDex;
	//int pdex, vdex;
	int pDex, vdex;
	DbVector3 *faceNormals=new DbVector3[thisMesh->numfaces()], c3tangent;
	CCPolygon *thisPoly;
	CCVertex *thisVert[3];

	for(pDex=0; pDex<thisMesh->numfaces(); pDex++)
	{
		thisPoly=thisMesh->FaceList(pDex);
		thisPoly->Pvertices(0);
		for(i=0; i<3; i++)
			thisVert[i]=thisMesh->VertList(thisPoly->Pvertices[i]);

		faceNormals[pDex]=-CalculateNormal(thisVert);

		//if(meshLvl=0)
		//	level0SurfaceNormals[pDex]=faceNormals[pDex];
		/*
		if(thisPoly->numsides()<4)
			faceNormals[pDex]=-faceNormals[pDex];
		*/
	}

	for(vdex = 0; vdex<thisMesh->numverts(); vdex++) 
	{
		CCVertex *thisVert = thisMesh->VertList(vdex);

		int vertvalence = thisVert->v_valence;

		// for each polygonal neighbor, form the weighted sum of the e_i edge
		//    and the f_i face (possibly computed, for a non-quad) into the
		//    limitloc and limitnorm DbVectors
		//thisVert->limitloc = thisVert->loc * (vertvalence*vertvalence);
		//c2tangent.set(0.0, 0.0, 0.0);
		//c3tangent.set(0.0, 0.0, 0.0);

		thisVert->limitnorm.set(0, 0, 0);
		CCPindex pdex;
		for (pdex = 0; pdex < vertvalence; pdex++) 
		{
			CCPindex thispolydex = thisMesh->VertFaceNeighbors[vdex][pdex];

			//CCPolygon *thispoly = thisMesh->FaceList(thispolydex);
			//int polyvdex = thislevel->VertFaceNeighborsVert[vdex][pdex];
			//DbVector3 e_j, f_j;
			//int pnumsides = thispoly->numsides();

			// The edge vertex is the next vertex counterclockwise from vpolydex
			//CCVertex *e_pdex = thislevel->VertList(thispoly->Pvertices[INTMOD(polyvdex + 1, pnumsides)]);
			//e_j = e_pdex->loc;

			// If pnumsides == 4, 2nd vertex counterclockwise from vpolydex is face vertex
			//    else have to regularize and just compute regloc
			//if (pnumsides == 4) 
			{
			//	CCVertex *f_pdex = thislevel->VertList(thispoly->Pvertices[INTMOD(polyvdex + 2, pnumsides)]);
			//	f_j = f_pdex->loc;
			}
			//else
			//	f_j = thislevel->regularize(thispolydex, polyvdex);

			// Now, add weighted values to limitloc and {c2,c3}tangent
			//thisVert->limitnorm += (e_j * 4.0) + f_j;
			thisVert->limitnorm+=faceNormals[thispolydex];
			//c2tangent += (e_j * Etable[vertvalence][pdex + 1]) + (f_j * Ftable[vertvalence][pdex + 1]);
			//c3tangent += (e_j * Etable[vertvalence][pdex    ]) + (f_j * Ftable[vertvalence][pdex    ]);
		}
		thisVert->limitnorm/=vertvalence;
	}

	delete []faceNormals;
}

//#else
//void CalculateNormalsFromVertices(int meshLvl, CCMesh *thisMesh)
//{
//}
#endif

void CCSubdivsurf::CompareNormalsFromVerticesAndFromSubdivision(int meshLvl)
{
#ifdef QIONG_DEBUG_BLENDED_SURFACE_NORMALS
	CCMesh *thisMesh=submeshes(meshLvl);

	int i; //, j;
	//int pDex, vDex;
	//int pdex, vdex;
	int pDex, vdex;
	DbVector3 *faceNormals=new DbVector3[thisMesh->numfaces()], c3tangent;
	CCPolygon *thisPoly;
	CCVertex *thisVert[3];

	MyFloat *accumulatedMag=new MyFloat[thisMesh->numverts()], magTotal=0, magAverage;

	fprintf(stderr, "The current level is %d\n", meshLvl);

	for(pDex=0; pDex<thisMesh->numfaces(); pDex++)
	{
		thisPoly=thisMesh->FaceList(pDex);
		thisPoly->Pvertices(0);
		for(i=0; i<3; i++)
			thisVert[i]=thisMesh->VertList(thisPoly->Pvertices[i]);

		faceNormals[pDex]=-CalculateNormal(thisVert);

		//if(meshLvl=0)
		//	level0SurfaceNormals[pDex]=faceNormals[pDex];
		/*
		if(thisPoly->numsides()<4)
			faceNormals[pDex]=-faceNormals[pDex];
		*/
	}

	for(vdex = 0; vdex<thisMesh->numverts(); vdex++) 
	{
		CCVertex *trueVert = thisMesh->VertList(vdex);
		CCVertex tmpVert, trueNorm;

		int vertvalence = trueVert->v_valence;

		// for each polygonal neighbor, form the weighted sum of the e_i edge
		//    and the f_i face (possibly computed, for a non-quad) into the
		//    limitloc and limitnorm DbVectors
		//thisVert->limitloc = thisVert->loc * (vertvalence*vertvalence);
		//c2tangent.set(0.0, 0.0, 0.0);
		//c3tangent.set(0.0, 0.0, 0.0);

		tmpVert.limitnorm.set(0, 0, 0);
		CCPindex pdex;
		for (pdex = 0; pdex < vertvalence; pdex++) 
		{
			CCPindex thispolydex = thisMesh->VertFaceNeighbors[vdex][pdex];

			//CCPolygon *thispoly = thisMesh->FaceList(thispolydex);
			//int polyvdex = thislevel->VertFaceNeighborsVert[vdex][pdex];
			//DbVector3 e_j, f_j;
			//int pnumsides = thispoly->numsides();

			// The edge vertex is the next vertex counterclockwise from vpolydex
			//CCVertex *e_pdex = thislevel->VertList(thispoly->Pvertices[INTMOD(polyvdex + 1, pnumsides)]);
			//e_j = e_pdex->loc;

			// If pnumsides == 4, 2nd vertex counterclockwise from vpolydex is face vertex
			//    else have to regularize and just compute regloc
			//if (pnumsides == 4) 
			{
			//	CCVertex *f_pdex = thislevel->VertList(thispoly->Pvertices[INTMOD(polyvdex + 2, pnumsides)]);
			//	f_j = f_pdex->loc;
			}
			//else
			//	f_j = thislevel->regularize(thispolydex, polyvdex);

			// Now, add weighted values to limitloc and {c2,c3}tangent
			//thisVert->limitnorm += (e_j * 4.0) + f_j;
			tmpVert.limitnorm+=faceNormals[thispolydex];
			//c2tangent += (e_j * Etable[vertvalence][pdex + 1]) + (f_j * Ftable[vertvalence][pdex + 1]);
			//c3tangent += (e_j * Etable[vertvalence][pdex    ]) + (f_j * Ftable[vertvalence][pdex    ]);
		}
		//thisVert->limitnorm/=vertvalence;
		tmpVert.limitnorm=tmpVert.limitnorm.normalize();
		trueNorm.limitnorm=trueVert->limitnorm.normalize();
		tmpVert.limitloc=tmpVert.limitnorm-trueNorm.limitnorm;
		accumulatedMag[vdex]=sqrt(tmpVert.limitloc.x()*tmpVert.limitloc.x()+tmpVert.limitloc.y()*tmpVert.limitloc.y()+tmpVert.limitloc.z()*tmpVert.limitloc.z());
		magTotal+=accumulatedMag[vdex];
		//fprintf(stderr, "The current level is %d\n", meshLvl);
	}

	MyFloat SD=0, magMin=999, magMax=-1;
	magAverage=magTotal/thisMesh->numverts();
	for(vdex=0; vdex<thisMesh->numverts(); vdex++)
	{
		if(magMin>accumulatedMag[vdex]) magMin=accumulatedMag[vdex];
		if(magMax<accumulatedMag[vdex]) magMax=accumulatedMag[vdex];
		SD+=(accumulatedMag[vdex]-magAverage)*(accumulatedMag[vdex]-magAverage);
	}
	SD/=thisMesh->numverts()-1;

	fprintf(stderr, "Normal difference: MIN=%f, MAX=%f, AVG=%f, SD=%f\n", magMin, magMax, magAverage, SD);

	delete []accumulatedMag;
	delete []faceNormals;
#endif
}

void CCSubdivsurf::glRenderSurfaceNormals(int nLvl)
{
#ifdef QIONG_DEBUG_BLENDED_SURFACE_NORMALS
	glPushAttrib(GL_CURRENT_BIT);
	glColor3d(1.0, 0, 0);
#ifdef QIONG_DEBUG_BLENDED_VERTEX_NORMALS
	CCMesh *thisMesh=submeshes(nLvl);
	//int pDex, vDex;
	int vDex;
	DbVector3 head, tail;

	glBegin(GL_LINES);

// to display the numerically calculated vertex-normals!
//*
	int i, pDex;
	DbVector3 *faceNormals=new DbVector3[thisMesh->numfaces()], c3tangent;
	CCPolygon *thisPoly;
	CCVertex *thisVert[3];
	for(pDex=0; pDex<thisMesh->numfaces(); pDex++)
	{
		thisPoly=thisMesh->FaceList(pDex);
		thisPoly->Pvertices(0);
		for(i=0; i<3; i++)
			thisVert[i]=thisMesh->VertList(thisPoly->Pvertices[i]);

		faceNormals[pDex]=-CalculateNormal(thisVert);
	}
	for(vDex=0; vDex<thisMesh->numverts(); vDex++)
	{
		CCVertex *thisVert=thisMesh->VertList(vDex);
		DbVector3 norm;
		norm.set(0, 0, 0);
		for(i=0; i<thisVert->v_valence; i++)
		{
			pDex=thisMesh->VertFaceNeighbors[vDex][i];
			norm+=faceNormals[pDex];
		}
		norm/=thisVert->v_valence;
		head=thisVert->limitloc;
		tail=head+0.03*norm.normalize();
		head.glVertex();
		tail.glVertex();
	}
	delete []faceNormals;	//*/

// to display the vertex-normals by the subdivision algorithm,
//	including either the control mesh's normal or the limit normal!
/*
	for(vDex=0; vDex<thisMesh->numverts(); vDex++)
	{
		CCVertex *thisVert=thisMesh->VertList(vDex);
		head=thisVert->limitloc;
		//tail=head+0.03*thisVert->limitnorm.normalize();
		tail=head+0.03*thisVert->limitnorm.normalize();
		head.glVertex();
		tail.glVertex();
	}						//*/
	glEnd();
	//delete []faceNormals;
#else

// to display the numerically calculated surface-normals
	CCMesh *thisMesh=submeshes(nLvl);
	int pDex, vDex;
	DbVector3 center, tail;

	int i; //, j;
	//int pdex, vdex;
	DbVector3 *faceNormals=new DbVector3[thisMesh->numfaces()], c3tangent;
	CCPolygon *thisPoly;
	CCVertex *thisVert[3];

	for(pDex=0; pDex<thisMesh->numfaces(); pDex++)
	{
		thisPoly=thisMesh->FaceList(pDex);
		thisPoly->Pvertices(0);
		for(i=0; i<3; i++)
			thisVert[i]=thisMesh->VertList(thisPoly->Pvertices[i]);

		faceNormals[pDex]=-CalculateNormal(thisVert);
		//if(meshLvl=0)
		//	level0SurfaceNormals[pDex]=faceNormals[pDex];
		/*
		if(thisPoly->numsides()<4)
			faceNormals[pDex]=-faceNormals[pDex];
		*/
	}

	glBegin(GL_LINES);
	for(pDex=0; pDex<thisMesh->numfaces(); pDex++)
	{
		CCPolygon *thisPoly=thisMesh->FaceList(pDex);
		center.set(0, 0, 0);
		for(vDex=0; vDex<thisPoly->numsides(); vDex++)
		{
			CCVertex *thisVert=thisMesh->VertList(thisPoly->Pvertices[vDex]);
			center+=thisVert->limitloc;
		}
		center/=thisPoly->numsides();
		tail=center+0.03*faceNormals[pDex];
		center.glVertex();
		tail.glVertex();
	}
	glEnd();
	delete []faceNormals;
#endif
	glPopAttrib();
#endif
}


/****************************************************************************/
/* printvals() -- Print values in CCVertex
/****************************************************************************/
void CCVertex::printvals(char *message)
{
	if (message != NULL)
		fprintf(stderr, "%s\n", message);

    loc.printvals("loc");
	norm.printvals("norm");
	limitloc.printvals("limitloc");
	limitnorm.printvals("limitnorm");
	true_loc.printvals("trueloc");
	true_norm.printvals("truenorm");

	fprintf(stderr, "rad = %f, uvt = [%f, %f, %f], valence = %d\n",
		    rad, uvt.u, uvt.v, uvt.t, v_valence);
}

/****************************************************************************/
/* CCPolygon::operator=()
/****************************************************************************/
CCPolygon& CCPolygon::operator=(CCPolygon& poly1)
{
	if (Pvertices.Max_Length() != 0)
		Pvertices.Destroy();
	if (Pneighbors.Max_Length() != 0)
		Pneighbors.Destroy();
	if (Pneighboredges.Max_Length() != 0)
		Pneighboredges.Destroy();
	if (childPolygons.Max_Length() != 0)
		childPolygons.Destroy();
	if (newedgeverts.Max_Length() != 0)
		newedgeverts.Destroy();

	Pvertices.Create(poly1.Pvertices.Length());
	Pneighbors.Create(poly1.Pneighbors.Length());
	Pneighboredges.Create(poly1.Pneighboredges.Length());
	childPolygons.Create(poly1.childPolygons.Length());
	newedgeverts.Create(poly1.newedgeverts.Length());

	parentPolygon = poly1.parentPolygon;
	newfacevert = poly1.newfacevert;
	// Now copy data from all lists to their associated in the new Polygon
	CCVindex vdex;

	for (vdex = 0; vdex < Pvertices.Length(); vdex++) {
		Pvertices[vdex] = poly1.Pvertices[vdex];
		Pneighbors[vdex] = poly1.Pneighbors[vdex];
		Pneighboredges[vdex] = poly1.Pneighboredges[vdex];
	}

	for (vdex = 0; vdex < childPolygons.Length(); vdex++)
		childPolygons[vdex] = poly1.childPolygons[vdex];

	for (vdex = 0; vdex < newedgeverts.Length(); vdex++)
		newedgeverts[vdex] = poly1.newedgeverts[vdex];

    return (*this);
}

/********************************************************************************/
/* Compute the vertex that regularizes the given polygon
/*   with <> 4 sides, about the indexed Pvertex
/* See Ch. 4 of Thall Diss for details
/********************************************************************************/
DbVector3 CCMesh::regularize(CCPindex thispolydex, int pvdex)
{
	CCPolygon *thispoly = FaceList(thispolydex);

	int nval = thispoly->numsides();

	DbVector3 regval;
	regval = VertList[thispoly->Pvertices[pvdex]].loc 
		   + VertList[thispoly->Pvertices[INTMOD(pvdex + 1, nval)]].loc
		   + VertList[thispoly->Pvertices[INTMOD(pvdex - 1, nval)]].loc;

	regval *= (4.0 - nval);

	for (int cnt = 2; cnt < nval - 1; cnt++)
		regval += VertList[thispoly->Pvertices[INTMOD(pvdex + cnt, nval)]].loc * 4.0;

	return regval/nval;
}
/********************************************************************************/
/* Compute the r of the vertex that regularizes the given polygon
/*   with <> 4 sides, about the indexed Pvertex
/* See Ch. 4 of Thall Diss for details
/********************************************************************************/
double CCMesh::regularizeR(CCPindex thispolydex, int pvdex)
{
	CCPolygon *thispoly = FaceList(thispolydex);

	int nval = thispoly->numsides();

	double regval;
	regval = VertList[thispoly->Pvertices[pvdex]].rad
		   + VertList[thispoly->Pvertices[INTMOD(pvdex + 1, nval)]].rad
		   + VertList[thispoly->Pvertices[INTMOD(pvdex - 1, nval)]].rad;

	regval *= (4.0 - nval);

	for (int cnt = 2; cnt < nval - 1; cnt++)
		regval += VertList[thispoly->Pvertices[INTMOD(pvdex + cnt, nval)]].rad * 4.0;

	return regval/nval;
}

/****************************************************************************/
/* printvals() -- Print values in CCPolygon
/****************************************************************************/
void CCPolygon::printvals(char *message)
{
	if (message != NULL)
		fprintf(stderr, "%s\n", message);

	int dex;
	fprintf(stderr, "vertex indices are:\n   ");
	for (dex = 0; dex < Pvertices.Length(); dex++)
		fprintf(stderr, "%d ", Pvertices[dex]);

	fprintf(stderr, "\nneighboring polygons are:\n   ");
	for (dex = 0; dex < Pneighbors.Length(); dex++)
        fprintf(stderr, "%d ", Pneighbors[dex]);

	fprintf(stderr, "\nneighboredges are:\n   ");
	for (dex = 0; dex < Pneighboredges.Length(); dex++)
		fprintf(stderr, "%d ", Pneighboredges[dex]);

	fprintf(stderr, "\nparent polygon is: %d\n", parentPolygon);
	fprintf(stderr, "newfacevert is: %d\n", newfacevert);

	fprintf(stderr, "\nchildpolygons are:\n   ");
	for (dex = 0; dex < childPolygons.Length(); dex++)
		fprintf(stderr, "%d ", childPolygons[dex]);

	fprintf(stderr, "\nnewedgeverts are:\n   ");
	for (dex = 0; dex < newedgeverts.Length(); dex++)
		fprintf(stderr, "%d ", newedgeverts[dex]);


}

bool tablesinitFLAG = false;
double CCSubdivsurf::Etable[MAXVALENCE + 1][MAXVALENCE + 1];
double CCSubdivsurf::Ftable[MAXVALENCE + 1][MAXVALENCE + 1];
/********************************************************************************/
/* init_tables() -- initialize limit vertex tables for edge and face weights	
/*    E[n][j] = weight for edge_j for valence n, for j=1 to n
/*    F[n][j] = weight for face_j (counterclockwise from edge i)
/*                for valence n, for j = 1 to n
/* (some wasted space, but don't want to do more indexing arithmetic)
/********************************************************************************/
void CCSubdivsurf::init_tables()
{
	if (tablesinitFLAG)
		return;
	tablesinitFLAG = true;

	for (int n = 3; n < MAXVALENCE + 1; n++) {

		double twopi_over_n = 2*M_PI/n;
		double A_n = 1.0 + cos(twopi_over_n) + cos(M_PI/n)*sqrt(18 + cos(twopi_over_n));

		for (int j = 0; j < MAXVALENCE + 1; j++) {
			Etable[n][j] = A_n * cos(twopi_over_n * j);
			Ftable[n][j] = cos(twopi_over_n * j)  + cos(twopi_over_n * (j + 1));
		}
	}
}

void CCSubdivsurf::glRenderBlended(int nivelo, RenderStyle thisstyle)
{
	int vertdex;
	//int polydex, vertdex;
	// Render the Polygons in submeshes[nivelo]
	//    if submeshes

	// Save current color
	glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_POINT_BIT);

	GLenum glstyle;
	switch (thisstyle) {
	case WIREFRAME:
		glLineWidth(lineDotWidth);
		glstyle = GL_LINE_LOOP;
		break;
	case FILLEDPOLY:
		glstyle = GL_POLYGON;
		break;
	case VPOINTS:
		glPointSize(lineDotWidth);
		//glColor3f(1.0, 0, 0.0);
		glstyle = GL_POINTS;
		break;
	case NO_RENDER:
		return;
		break;
	}

	if(glstyle==GL_POLYGON)	glEnable(GL_LIGHTING);
#ifdef QIONG_DEBUG_02232003
	int badVertexCounter=0;
	for(vertdex=0; vertdex<submeshes[nivelo].numverts(); vertdex++)
	{
		if(nivelo==1&&vRedex[vertdex]==-1)
			continue;
		CCVertex *thisVert = submeshes[nivelo].VertList(vertdex);
		if(thisVert->v_valence==3)
			badVertexCounter++;
	}
	badVertexCounter=badVertexCounter;
	fprintf(stderr, "Level=%d\tBad Vertices: %d\n", nivelo, badVertexCounter);
#endif
	glBegin(glstyle);
	//for(vertdex=submeshes[nivelo].numverts()/2; vertdex<submeshes[nivelo].numverts(); vertdex++)//submeshes[nivelo].numverts(); vertdex++) 
	for(vertdex=0; vertdex<submeshes[nivelo].numverts()/2; vertdex++)//submeshes[nivelo].numverts(); vertdex++) 
	{
		if(nivelo==1&&vRedex[vertdex]==-1)
			continue;
		CCVertex *thisVert = submeshes[nivelo].VertList(vertdex);
		thisVert->loc.glVertex();
	}
	glEnd();

	if(glstyle==GL_POLYGON)	glDisable(GL_LIGHTING);

	glPopAttrib();
}
//glRenderBlended(surfaceLevel, VPOINTS);


void CCSubdivsurf::glRender(int vNum, int fNum, Bpoint2 *vList, PindexList *fList, RenderStyle thisStyle, MyFloat *c)
{
	glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_POINT_BIT);

	GLenum glstyle;
	switch (thisStyle) {
	case WIREFRAME:
		glLineWidth(lineDotWidth);
		glstyle = GL_LINE_LOOP;
		break;
	case FILLEDPOLY:
		glstyle = GL_POLYGON;
		break;
	case VPOINTS:
		glPointSize(lineDotWidth);
		glstyle = GL_POINTS;
		break;
	case NO_RENDER:
		return;
		break;
	}

	glColor3dv(c);
	if(glstyle==GL_POLYGON)	glEnable(GL_LIGHTING);
	for(int pDex=0; pDex<fNum; pDex++) 
	{
		glBegin(glstyle);
		for(int vDex=0; vDex<fList[pDex].Length(); vDex++) 
		{
			glNormal3d(vList[fList[pDex][vDex]].norm[0], vList[fList[pDex][vDex]].norm[1], vList[fList[pDex][vDex]].norm[2]);
			glVertex3d(vList[fList[pDex][vDex]].pnt[0], vList[fList[pDex][vDex]].pnt[1], vList[fList[pDex][vDex]].pnt[2]);

		}
		glEnd();
	}

	if(glstyle==GL_POLYGON)	glDisable(GL_LIGHTING);
	glPopAttrib();
}

void CCSubdivsurf::glRenderSegments(int startPDex, int endPDex, Bpoint2 *pntList, MyFloat *c)
{
	int i;

	glPushAttrib(GL_CURRENT_BIT | GL_POINT_BIT);

	glColor3dv(c);
	glBegin(GL_LINES);
	for(i=startPDex; i<=endPDex; i++)
	{
		glVertex3d(pntList[2*i].pnt[0], pntList[2*i].pnt[1], pntList[2*i].pnt[2]);
		glVertex3d(pntList[2*i+1].pnt[0], pntList[2*i+1].pnt[1], pntList[2*i+1].pnt[2]);
	}
	glEnd();
	glColor3d(1.0, 0.0, 0.0);
	glPointSize(lineDotWidth);
	glBegin(GL_POINTS);
	for(i=startPDex; i<=endPDex; i++)
	{
		glVertex3d(pntList[2*i].pnt[0], pntList[2*i].pnt[1], pntList[2*i].pnt[2]);
		glVertex3d(pntList[2*i+1].pnt[0], pntList[2*i+1].pnt[1], pntList[2*i+1].pnt[2]);
	}
	glEnd();
	glPointSize(1);
	glPopAttrib();
}

void CCSubdivsurf::glRenderPoint(Bpoint2 pnt, MyFloat *c)
{
	glPushAttrib(GL_CURRENT_BIT | GL_POINT_BIT);

	glPointSize(lineDotWidth);
	glColor3dv(c);
	glBegin(GL_POINTS);
	//for(i=startPDex; i<=endPDex; i++)
	//{
		glVertex3d(pnt.pnt[0], pnt.pnt[1], pnt.pnt[2]);
	//}
	glEnd();

	glPopAttrib();
}

void CCSubdivsurf::glRenderCurve(int startPDex, int endPDex, Bpoint2 *pntList, MyFloat *c)
{
	int i;

	// Save current color/line/point information
	glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_POINT_BIT);

	glColor3dv(c);
	glBegin(GL_LINE_STRIP);
	for(i=startPDex; i<=endPDex; i++)
	{
		glVertex3d(pntList[i].pnt[0], pntList[i].pnt[1], pntList[i].pnt[2]);
	}
	glEnd();
	glColor3d(1.0, 0.0, 0.0);
	glPointSize(lineDotWidth);
	glBegin(GL_POINTS);
	for(i=startPDex; i<=endPDex; i++)
	{
		glVertex3d(pntList[i].pnt[0], pntList[i].pnt[1], pntList[i].pnt[2]);
	}
	glEnd();
	glPopAttrib();
}

/****************************************************************************/
/* glRender() render the tiles of submeshes(nivelo) according to input
/*		thisstyle (enum in CCSubdivsurf.h)
/****************************************************************************/
void CCSubdivsurf::glRender(int nivelo, int polyDex, RenderStyle thisstyle, MyFloat *c)
{
	// Render the Polygons in submeshes[nivelo]
	//    if submeshes

	// Save current color/line/point information
	glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_POINT_BIT);

	GLenum glstyle;
	switch (thisstyle) {
	case WIREFRAME:
		glLineWidth(lineDotWidth);
		glstyle = GL_LINE_LOOP;
		break;
	case FILLEDPOLY:
		glstyle = GL_POLYGON;
		break;
	case VPOINTS:
		glPointSize(lineDotWidth);
		glstyle = GL_POINTS;
		break;
	case NO_RENDER:
		return;
		break;
	}

	//fprintf(stderr, "\nNumpolygons = %d\n", submeshes[nivelo].numfaces());

	glColor3dv(c);
	if(glstyle==GL_POLYGON)	glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHTING);
	CCPolygon *thispoly = submeshes[nivelo].FaceList(polyDex);

	glBegin(glstyle);
	for(int vertdex=0; vertdex<thispoly->numsides(); vertdex++) 
	{
		CCVertex *thisVert = submeshes[nivelo].VertList(thispoly->Pvertices[vertdex]);
		thisVert->limitnorm.glNormal();
		thisVert->limitloc.glVertex();
	}
	glEnd();
	if(glstyle==GL_POLYGON)	glDisable(GL_LIGHTING);

	glPopAttrib();
}

/****************************************************************************/
/* glRender() render the tiles of submeshes(nivelo) according to input
/*		thisstyle (enum in CCSubdivsurf.h)
/****************************************************************************/
void CCSubdivsurf::glRender(int nivelo, RenderStyle thisstyle)
{
	// Render the Polygons in submeshes[nivelo]
	//    if submeshes

	// Save current color/line/point information
	glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_POINT_BIT);

	GLenum glstyle;
	switch (thisstyle) {
	case WIREFRAME:
		glLineWidth(lineDotWidth);
		glstyle = GL_LINE_LOOP;
		break;
	case FILLEDPOLY:
		glstyle = GL_POLYGON;
		break;
	case VPOINTS:
		glPointSize(lineDotWidth);
		glstyle = GL_POINTS;
		break;
	case NO_RENDER:
		return;
		break;
	}

	int polydex;

#ifdef QIONG_DEBUG_BLENDED_SURFACE_NORMALS
	glRenderSurfaceNormals(nivelo);
#endif

if(nivelo>=0 && nivelo<=99)
{
	if(glstyle==GL_POLYGON)	glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHTING);
	int polydex;

#ifdef DEBUG_QIONG_BLENDING
	double ext, ext2, ext3;
	if(dPNum>0)
	{
		for(polydex = 0; polydex < dPNum; polydex++) 
		{
			CCPolygon *thispoly = submeshes[nivelo].FaceList(polyList[polydex]);

			glBegin(glstyle);
			for(int vertdex = 0; vertdex < thispoly->numsides(); vertdex++) 
			{
				CCVertex *thisVert = submeshes[nivelo].VertList(thispoly->Pvertices[vertdex]);
				thisVert->limitnorm.glNormal();
				thisVert->limitloc.glVertex();
			}
			glEnd();
		}
	}

	for(polydex=56; polydex<submeshes[nivelo].numfaces(); polydex++) 
	{
		CCPolygon *thispoly = submeshes[nivelo].FaceList(polydex);
		CCVertex *thisVert;

		ext=999;
			for(int vertdex=0; vertdex<thispoly->numsides(); vertdex++) 
			{
				thisVert = submeshes[nivelo].VertList(thispoly->Pvertices[vertdex]);
				if(thisVert->uvt.u<ext)
					ext=thisVert->uvt.u;
			}
			if(ext>=1.5)
			{
				glBegin(glstyle);
				for(vertdex=0; vertdex<thispoly->numsides(); vertdex++) 
				{
					thisVert = submeshes[nivelo].VertList(thispoly->Pvertices[vertdex]);
					thisVert->limitnorm.glNormal();
					thisVert->limitloc.glVertex();
				}
				glEnd();
			}
	}
#else
	for(polydex = 0; polydex < submeshes[nivelo].numfaces(); polydex++) 
	{
		CCPolygon *thispoly = submeshes[nivelo].FaceList(polydex);
		glBegin(glstyle);
		for (int vertdex = 0; vertdex < thispoly->numsides(); vertdex++) 
		{
			CCVertex *thisVert = submeshes[nivelo].VertList(thispoly->Pvertices[vertdex]);

			thisVert->limitnorm.glNormal();
			thisVert->limitloc.glVertex();
		}
		glEnd();
	}
#endif

	if(glstyle==GL_POLYGON)	glDisable(GL_LIGHTING);
	/*
	for (polydex = 0; polydex < submeshes[nivelo].numfaces(); polydex++) {

		CCPolygon *thispoly = submeshes[nivelo].FaceList(polydex);

		glBegin(glstyle);
		for (int vertdex = 0; vertdex < thispoly->numsides(); vertdex++) {
			CCVertex *thisVert = submeshes[nivelo].VertList(thispoly->Pvertices[vertdex]);
			thisVert->limitnorm.glNormal();
			thisVert->limitloc.glVertex();
		}
		glEnd();
	}
	if(glstyle==GL_POLYGON)	glDisable(GL_LIGHTING);
	*/

	/*
	glColor3d(0.0, 1.0, 1.0);
	for (polydex = 0; polydex < submeshes[nivelo].numfaces(); polydex++) {

		CCPolygon *thispoly = submeshes[nivelo].FaceList(polydex);

		glBegin(GL_LINE_LOOP);
		for (int vertdex = 0; vertdex < thispoly->numsides(); vertdex++) {
			CCVertex *thisVert = submeshes[nivelo].VertList(thispoly->Pvertices[vertdex]);
			thisVert->norm.glNormal();
			thisVert->loc.glVertex();
		}
		glEnd();
	}
	glDisable(GL_LIGHTING);
	*/

}
else
{
	glstyle = GL_LINE_LOOP;
	if(nivelo==-1)
	{
#if (defined(DEBUG_QIONG_TESTING_INVERSE_ALGORITHM_ON_QUAD_MESH) || defined(QIONG_DEBUG_BLENDING_02192003))
		nivelo=1;
#else
		nivelo++;
#endif

	fprintf(stderr, "\nNumpolygons = %d\n", submeshes[nivelo].numfaces());
	glColor3d(0.0, 1.0, 0.0);
	for(polydex = 0; polydex < submeshes[nivelo].numfaces(); polydex++) 
	{
		CCPolygon *thispoly = submeshes[nivelo].FaceList(polydex);
		glBegin(glstyle);
		for (int vertdex = 0; vertdex < thispoly->numsides(); vertdex++) 
		{
			CCVertex *thisVert = submeshes[nivelo].VertList(thispoly->Pvertices[vertdex]);

			thisVert->limitnorm.glNormal();
			thisVert->limitloc.glVertex();
			//thisVert->loc.glVertex();

		}
		glEnd();
	}
	}
}
	glPopAttrib();
}

void CCSubdivsurf::glRenderBlendedObject(int nivelo, RenderStyle thisstyle)
{
	// Render the Polygons in submeshes[nivelo]
	//    if submeshes

	// Save current color/line/point information
	glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_POINT_BIT);

	GLenum glstyle;
	switch (thisstyle) {
	case WIREFRAME:
		glLineWidth(lineDotWidth);
		glstyle = GL_LINE_LOOP;
		break;
	case FILLEDPOLY:
		glstyle = GL_POLYGON;
		break;
	case VPOINTS:
		glPointSize(lineDotWidth);
		glstyle = GL_POINTS;
		break;
	case NO_RENDER:
		return;
		break;
	}

	int polydex;

#ifdef QIONG_DEBUG_BLENDED_SURFACE_NORMALS
	glRenderSurfaceNormals(nivelo);
#endif

if(nivelo>=0 && nivelo<=99)
{
	if(glstyle==GL_POLYGON)	glEnable(GL_LIGHTING);
	int polydex;

	for(polydex = 0; polydex < submeshes[nivelo].numfaces(); polydex++) 
	{
		CCPolygon *thispoly = submeshes[nivelo].FaceList(polydex);
		glBegin(glstyle);
		for (int vertdex = 0; vertdex < thispoly->numsides(); vertdex++) 
		{
			CCVertex *thisVert = submeshes[nivelo].VertList(thispoly->Pvertices[vertdex]);

			thisVert->norm.normalize().glNormal();
			thisVert->limitloc.glVertex();
		}
		glEnd();
	}

	if(glstyle==GL_POLYGON)	glDisable(GL_LIGHTING);
}
else
{
	glstyle = GL_LINE_LOOP;
	if(nivelo==-1)
	{
#if (defined(DEBUG_QIONG_TESTING_INVERSE_ALGORITHM_ON_QUAD_MESH) || defined(QIONG_DEBUG_BLENDING_02192003))
		nivelo=1;
#else
		nivelo++;
#endif

	fprintf(stderr, "\nNumpolygons = %d\n", submeshes[nivelo].numfaces());
	glColor3d(0.0, 1.0, 0.0);
	for(polydex = 0; polydex < submeshes[nivelo].numfaces(); polydex++) 
	{
		CCPolygon *thispoly = submeshes[nivelo].FaceList(polydex);
		glBegin(glstyle);
		for (int vertdex = 0; vertdex < thispoly->numsides(); vertdex++) 
		{
			CCVertex *thisVert = submeshes[nivelo].VertList(thispoly->Pvertices[vertdex]);

			thisVert->limitnorm.glNormal();
			thisVert->limitloc.glVertex();
		}
		glEnd();
	}

	}
}
	glPopAttrib();
}

/********************************************************************************/
/* initdummy() -- create parent mesh for tritorus or pyramid with 4-sided base	*/
/********************************************************************************/
void CCSubdivsurf::initdummy(TestMeshShape thisshape)
{
	switch(thisshape)
	{
	case TRITORUS:
		inittritorus();
		break;
	case FIGURE8:
		initfigure8();
		break;
	case PYRAMID:
		initpyramid();
		break;
	case CUSTOMIZED:
		initcustomized();
		break;
	case FROMFILE:
		initfromfile();
		break;
	default:
		break;
	}
	loadlevelzero();
}


CCMesh *CCSubdivsurf::SubMeshPointer(int subLvl)
{
	return submeshes(subLvl);
}

void CCSubdivsurf::LoadLevelZeroValue()
{
	//int tempVValence;
	int numverts1=fathermesh.numverts();
	CCMesh *thisMesh=submeshes(0);
	CCVindex vdex;
	for (vdex = 0; vdex < numverts1; vdex++)
	{
		//tempVValence=thisMesh->VertList[vdex].v_valence;
		thisMesh->VertList[vdex] = fathermesh.VertList[vdex];
	}
	computeIISmesh(5);
	computelimitmesh(0);
}
/********************************************************************************/
/* loadlevelzero() -- load values of fathermesh into submeshes[0]				*/
/********************************************************************************/
void CCSubdivsurf::loadlevelzero()
{
	CCMesh *thisMesh = submeshes(0);

	int numfaces1, numverts1, nedges1;
	numfaces1 = fathermesh.numfaces();
	numverts1 = fathermesh.numverts();
	nedges1 = fathermesh.numedges();

	thisMesh->VertList.Create(numverts1);
	thisMesh->FaceList.Create(numfaces1);
	thisMesh->nedges = nedges1;

	CCVindex vdex;
	//for (vdex = 0; vdex < numverts1; vdex++)
	//	thisMesh->VertList[vdex] = fathermesh.VertList[vdex];

	CCPindex fdex;
	for (fdex = 0; fdex < numfaces1; fdex++) {

		CCPolygon *fathersface = fathermesh.FaceList(fdex);

		// Copy the Polygons
		thisMesh->FaceList[fdex] = fathermesh.FaceList[fdex];
	}

	// Now that VertexList and FaceList are done, do VertFaceNeighbors and FaceFaceNeighbors
	//   though we may not need these both --- For now, skip FaceFaceNeighbors
	thisMesh->VertFaceNeighbors.Create(numverts1);
	thisMesh->VertFaceNeighborsVert.Create(numverts1);

	for (vdex = 0; vdex < numverts1; vdex++) {

		int numvneighbors = fathermesh.VertFaceNeighbors[vdex].Length();
		thisMesh->VertFaceNeighbors[vdex].Create(numvneighbors);
		thisMesh->VertFaceNeighborsVert[vdex].Create(numvneighbors);

		for (fdex = 0; fdex < numvneighbors; fdex++) {
			thisMesh->VertFaceNeighbors[vdex][fdex] = fathermesh.VertFaceNeighbors[vdex][fdex];
			thisMesh->VertFaceNeighborsVert[vdex][fdex] = fathermesh.VertFaceNeighborsVert[vdex][fdex];
		}
	}

	// Now, compute interpolating mesh, then limitmesh
	/*
	glPushAttrib(GL_CURRENT_BIT);
	glColor3d(1.0, 1.0, 0.0);
	for (int polydex = 0; polydex < submeshes[0].numfaces(); polydex++) {

		CCPolygon *thispoly = submeshes[0].FaceList(polydex);

		glBegin(GL_LINE_LOOP);
		for (int vertdex = 0; vertdex < thispoly->numsides(); vertdex++) {
			CCVertex *thisVert = submeshes[0].VertList(thispoly->Pvertices[vertdex]);
			thisVert->norm.glNormal();
			thisVert->loc.glVertex();
		}
		glEnd();
	}
	glPopAttrib();
	*/

	/*
	computeIISmesh(5);
	computelimitmesh(0);
	*/
	LoadLevelZeroValue();
}

/*
//******************************************************************************
// SplitandAverage() -- split nivelo to give the next finer level (nivelo + 1)	
//******************************************************************************
bool CCSubdivsurf::splitandaverage(int nivelo)
{
	if (nivelo + 1 >= NUMLEVELS) {
		fprintf(stderr, "ERROR: CCSubdivsurf::splitandaverage() -- Number of levels too high for split.\n");
		return false;
	}

	CCMesh *oldlevel = submeshes(nivelo);
	CCMesh *newlevel = submeshes(nivelo + 1);

	// The new mesh has the even vertices + the new face vertices + the new edge vertices
	int newnumverts = oldlevel->numverts() + oldlevel->numfaces() + oldlevel->numedges();

	// Each new edge-vertex generates 2 new faces
	//   (As proof, each new edge from the edge-vertex to the 2 face-vertices forms the
	//    lead edge (in counterclockwise order) of one of the new faces)
	int newnumfaces = 2*oldlevel->numedges();

	newlevel->VertList.Create(newnumverts);
	newlevel->FaceList.Create(newnumfaces);
	newlevel->VertFaceNeighbors.Create(newnumverts);
	newlevel->VertFaceNeighborsVert.Create(newnumverts);

	int numoldfaces = oldlevel->numfaces();
	int numoldverts = oldlevel->numverts();

	// Copy even vertices into new vertex list
	// NO, begin by zeroing out the even vertices...use parent's even vertices
	//    (which have the same indices) for computation, and use the new ones
	//    as accumulators for computing the new even vertices.
	// Formula is neweven = oldweight*v_old + neighborweight*(sum[4*newedgeneighbors - newfaceneighbors])
	//    Derivation of this is in Thall diss, App.A.
	int evenvindex;
	for (evenvindex = 0; evenvindex < numoldverts; evenvindex++) {
		newlevel->VertList[evenvindex] = oldlevel->VertList[evenvindex];
		newlevel->VertList[evenvindex].is_oddvert = false;

		newlevel->VertList[evenvindex].loc = NULLVEC;
		newlevel->VertList[evenvindex].norm = NULLVEC;
	}

	DbVector3 newloc;
	DbVector3 newnorm;
	Medcoord newuvt;
	double newradius;

	int polydex;
	// Create new face vertices at centroids of their polygons
	for (polydex = 0; polydex < numoldfaces; polydex++) {

		CCPolygon *thispoly = oldlevel->FaceList(polydex);

//		fprintf(stderr, "old polygon #%d is:\n", polydex);
//		thispoly->printvals();
//		fprintf(stderr, "\n");

		newloc.set(0.0, 0.0, 0.0);
		newnorm.set(0.0, 0.0, 0.0);
		newuvt.setvals(0.0, 0.0, 0.0);
		newradius = 0.0;

		int v_index;
		for (v_index = 0; v_index < thispoly->numsides(); v_index++) {

			CCVertex *thisVertex = oldlevel->VertList(thispoly->Pvertices[v_index]);

			// Sum the node values
			newloc += thisVertex->loc;
			newnorm += thisVertex->norm;
			newradius += thisVertex->rad;
			newuvt += thisVertex->uvt;
		}

		newloc /= thispoly->numsides();
		newnorm /= thispoly->numsides();
		newradius /= thispoly->numsides();
		newuvt /= thispoly->numsides();

		newlevel->VertList[numoldverts + polydex].setvals(newloc, newnorm, newradius,
			                                              newuvt, thispoly->numsides(), true);
		thispoly->newfacevert = numoldverts + polydex;
		thispoly->newedgeverts.Create(thispoly->numsides()); // go ahead and initialize newedgevert list

		// Now, subtract values of the new face vertex into loc and norm of surrounding even vertices
		for (v_index = 0; v_index < thispoly->numsides(); v_index++) {

			CCVertex *thisVertex = newlevel->VertList(thispoly->Pvertices[v_index]);

			thisVertex->loc -= newloc;
			thisVertex->norm -= newnorm;
			// NOTE:  if you make this an addition instead of a subtraction,
			//   it screws everything.  :)  Makes some nice sea anemonies, though
		}
	}

	// For each edge of polygon, create new edge vertex at averaged position of
	//    of old edge endpoints and new face vertices at centroids
	// For r and uvt interpolation, just average old endpoints
	CCVindex currentnewvert = numoldverts + numoldfaces;

	for (polydex = 0; polydex < numoldfaces; polydex++) {

		CCPolygon *thispoly = oldlevel->FaceList(polydex);

		// Now, need those pointers from polygon to new vertices at next level for
		//   face and for edges
		for (int edgedex = 0; edgedex < thispoly->numsides(); edgedex++) {

			CCPindex neighbordex = thispoly->Pneighbors[edgedex];

			CCPolygon *neighborpoly = oldlevel->FaceList(neighbordex);

			// If the across-edge neighbor polygon has a lower index, its edges
			//   have already been computed, so just store pointer to new edgevertex in parent Polygon
			if (polydex > neighbordex) {
				thispoly->newedgeverts[edgedex] = neighborpoly->newedgeverts[thispoly->Pneighboredges[edgedex]];
			}
			// else create a new vertex by averaging the new face and old edge-endpoints
			else {
				CCVertex *facevert1 = newlevel->VertList(thispoly->newfacevert);
				CCVertex *facevert2 = newlevel->VertList(neighborpoly->newfacevert);
				CCVertex *edgevert1 = oldlevel->VertList(thispoly->Pvertices[edgedex]);
				CCVertex *edgevert2 = oldlevel->VertList(thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]);


				newloc = edgevert1->loc + edgevert2->loc + facevert1->loc + facevert2->loc ;
				newnorm = edgevert1->norm + edgevert2->norm + facevert1->norm + facevert2->norm;
				newuvt = edgevert1->uvt + edgevert2->uvt;
				newradius = edgevert1->rad + edgevert2->rad;

				newloc /= 4;
				newnorm /= 4;
				newradius /= 2;
				newuvt /= 2;

				newlevel->VertList[currentnewvert].setvals(newloc, newnorm, newradius, newuvt, 4, true);
				thispoly->newedgeverts[edgedex] = currentnewvert++;

				// Now, add 4x the loc and norm values to each of the two edgeverts
				//   This happens just once, when this vertex is created
				CCVertex *thisVertex1 = newlevel->VertList(thispoly->Pvertices[edgedex]);
				CCVertex *thisVertex2 = newlevel->VertList(thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]);

				thisVertex1->loc += newloc*4.0;
				thisVertex1->norm += newnorm*4.0;
				thisVertex2->loc += newloc*4.0;
				thisVertex2->norm += newnorm*4.0;
			}
		}
	}

	// Now, form new polygons and compute all connectivity
	//    Do this in two stages:
	//    1) create new polygons and assign their
	//            numbers to parent polygons' childpolygon lists
	//    2) Use this info to set neighbor and neighboredge lists
	CCPindex newpolydex = 0;
	for (polydex = 0; polydex < numoldfaces; polydex++) {

		CCPolygon *thispoly = oldlevel->FaceList(polydex);
		int nsides = thispoly->numsides();

		thispoly->childPolygons.Create(nsides);

		for (int vdex = 0; vdex < nsides; vdex++) 
			thispoly->childPolygons[vdex] = newpolydex++;
	}

	newpolydex = 0;
	int totalsides = 0;
	for (polydex = 0; polydex < numoldfaces; polydex++) {

		CCPolygon *thispoly = oldlevel->FaceList(polydex);
		int nsides = thispoly->numsides();

		// count number of edges for total
		totalsides += nsides*4;

		int firstnewpolydex = newpolydex;
		int lastnewpolydex = newpolydex + nsides - 1;

		for (evenvindex = 0; evenvindex < nsides; evenvindex++) {

			CCPolygon *newpoly = newlevel->FaceList(newpolydex);
			newpoly->Pvertices.Create(4);

			newpoly->Pvertices[0] = thispoly->newfacevert;
			newpoly->Pvertices[1] = thispoly->newedgeverts[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pvertices[2] = thispoly->Pvertices[evenvindex];
			newpoly->Pvertices[3] = thispoly->newedgeverts[evenvindex];

			// Now do polygon neighbors

			int facenum;
			newpoly->Pneighbors.Create(4);
			// The edge number is clockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex - 1 < 0 ? nsides - 1 : evenvindex - 1);
			newpoly->Pneighbors[0] = thispoly->childPolygons[facenum];

			// The edge number of the neighbor polygon gives the polygon-index
			facenum = thispoly->Pneighboredges[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pneighbors[1] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[facenum];

			// The edge number + 1 MOD neighbors-numsides of the neighboredge gives the polygon-index
			facenum = INTMOD(thispoly->Pneighboredges[evenvindex] + 1,
				             oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].numsides());
			newpoly->Pneighbors[2] = oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].childPolygons[facenum];

			// The edge number is counterclockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex + 1 > nsides - 1 ? 0 : evenvindex + 1);
			newpoly->Pneighbors[3] = thispoly->childPolygons[facenum];

			// -- if all polygons start at the new face vertex, then the edge#'s are all 1
			//    or 2 for neighboring polygons.  (2 facing 1, 1 facing 2)
			newpoly->Pneighboredges.Create(4);
			newpoly->Pneighboredges[0] = 3;
			newpoly->Pneighboredges[1] = 2;
			newpoly->Pneighboredges[2] = 1;
			newpoly->Pneighboredges[3] = 0;

			newpoly->parentPolygon = polydex;
			newpoly->newfacevert = NULLVERTDEX;
//			newpoly->printvals();
//			fprintf(stderr, "\n");
			newpolydex++;
		}
	}

	newlevel->nedges = totalsides/2;

	// Now, average the even vertices according to their connectivity
	// ----no, this can be done on the fly.  Zero all nivelo+1 even vertices, and
	// add 4x the new edge values and 1x each new face vertex.
	//    Once here, weight the computed sum and add the weighted parent vertex weight
	// Use the weights:
	//    1/(n^2) for the face-edge sum, and (n-3)/n for the old vertex
	// (See App.A of Thall Diss for the derivation of this.)
	for (evenvindex = 0; evenvindex < numoldverts; evenvindex++) 
	{
		double thisval = newlevel->VertList[evenvindex].v_valence;

		double sumweight, oldweight;

		sumweight = 1.0/(thisval*thisval);
		oldweight = (thisval - 3.0)/thisval;

		newlevel->VertList[evenvindex].loc *= sumweight;
		newlevel->VertList[evenvindex].norm *= sumweight;

		newlevel->VertList[evenvindex].loc += oldlevel->VertList[evenvindex].loc * oldweight;
		newlevel->VertList[evenvindex].norm += oldlevel->VertList[evenvindex].norm * oldweight;
	}

	// Now, compute all the connectivity relationships for the CCmesh, as necessary
	// Can probably do this on the fly, too.
	// For each polygon in old mesh:
	//    new children are neighborhood of newfacevert
	//    two new children and their 2-edge and 1-edge neighbors are neighborhood of new edge
	//    new child and its 1-edge neighbor and ITS 1-edge neighbors are
	//       neighborhood of new even vertex
	//    Vertex 2 of each new polygon is old even vertex
	// If VertFaceNeighbors[vdex] has been initialized, then already computed!
	for (polydex = 0; polydex < numoldfaces; polydex++) {
		CCPolygon *oldpoly = oldlevel->FaceList(polydex);

		CCVindex newvertdex = oldpoly->newfacevert;

		// newvalence is number of edges of old polygon, which will be used
		//    all over the place
		int newvalence = newlevel->VertList[newvertdex].v_valence;

		// Compute neighborhood of new face vertex
		newlevel->VertFaceNeighbors[newvertdex].Create(newvalence);
		newlevel->VertFaceNeighborsVert[newvertdex].Create(newvalence);

		CCPindex newpolydex;
		for (newpolydex = 0; newpolydex < newvalence; newpolydex++) {
			newlevel->VertFaceNeighbors[newvertdex][newpolydex] = oldpoly->childPolygons[newpolydex];
			newlevel->VertFaceNeighborsVert[newvertdex][newpolydex] = 0;
		}

		// Compute neighborhood of each new edge vertex if VertFaceNeighbors[] not already initialized
		for (int newoddedge = 0; newoddedge < newvalence; newoddedge++) {
			newvertdex = oldpoly->newedgeverts[newoddedge];
			if (newlevel->VertFaceNeighbors[newvertdex].Length() == 0) {
				newlevel->VertFaceNeighbors[newvertdex].Create(4);
				newlevel->VertFaceNeighborsVert[newvertdex].Create(4);

				CCPindex pdexA, pdexB, pdexC, pdexD;
				pdexA = oldpoly->childPolygons[INTMOD(newoddedge + 1, newvalence)];
				pdexB = oldpoly->childPolygons[newoddedge];
				pdexC = newlevel->FaceList[pdexB].Pneighbors[2];
				pdexD = newlevel->FaceList[pdexA].Pneighbors[1];

				newlevel->VertFaceNeighbors[newvertdex][0] = pdexA;
				newlevel->VertFaceNeighbors[newvertdex][1] = pdexB;
				newlevel->VertFaceNeighbors[newvertdex][2] = pdexC;
				newlevel->VertFaceNeighbors[newvertdex][3] = pdexD;
				newlevel->VertFaceNeighborsVert[newvertdex][0] = 1;
				newlevel->VertFaceNeighborsVert[newvertdex][1] = 3;
				newlevel->VertFaceNeighborsVert[newvertdex][2] = 1;
				newlevel->VertFaceNeighborsVert[newvertdex][3] = 3;
			}
		}

		// Compute neighborhood of each new even vertex if VertFaceNeighbors[] not already initialized
		for (int childpolydex = 0; childpolydex < newvalence; childpolydex++) {

			CCPindex newface = oldpoly->childPolygons[childpolydex];
			newvertdex = newlevel->FaceList[newface].Pvertices[2];

			if (newlevel->VertFaceNeighbors[newvertdex].Length() == 0) {

				int vertvalence = newlevel->VertList[newvertdex].v_valence;
				newlevel->VertFaceNeighbors[newvertdex].Create(vertvalence);
				newlevel->VertFaceNeighborsVert[newvertdex].Create(vertvalence);

				CCPindex pdexE = newface;

				newlevel->VertFaceNeighbors[newvertdex][0] = pdexE;
				newlevel->VertFaceNeighborsVert[newvertdex][0] = 2;

				for (int facedex = 1; facedex < vertvalence; facedex++) {

					pdexE = newlevel->FaceList[pdexE].Pneighbors[1];
					newlevel->VertFaceNeighbors[newvertdex][facedex] = pdexE;
					newlevel->VertFaceNeighborsVert[newvertdex][facedex] = 2;
				}
			}
		}
	}

	// Lastly, compute the limit mesh
	computelimitmesh(nivelo + 1);

	// This will cause normal-interpolation to be used at level-2 and higher subdivision
	if(nivelo==1 && NORMINTERP) 
		perturbleveltwo();

	return true;
}
*/

//******************************************************************************
// SplitandAverage() -- split nivelo to give the next finer level (nivelo + 1)	
//******************************************************************************
bool CCSubdivsurf::splitandaverage(int nivelo)
{
	if (nivelo + 1 >= NUMLEVELS) 
	{
		fprintf(stderr, "ERROR: CCSubdivsurf::splitandaverage() -- Number of levels too high for split.\n");
		return false;
	}
#ifdef QIONG_DEBUG_SEURAT_OPTIMIZATION
if(nivelo==2) //if(nivelo==NUMLEVELS-2)
{
    struct _timeb thistime;
	int i;
	int totalRuns=50;
	double ltime1, ltime2;
for(int ii=0; ii<10; ii++)
{
    _ftime(&thistime);    
    ltime1=thistime.time + thistime.millitm * 1e-3;
	for(i=0; i<totalRuns; i++)
	{
		Subdivide(0);
		Subdivide(1);
		//Subdivide(2);
	}
    _ftime(&thistime);    
    ltime2=thistime.time + thistime.millitm * 1e-3;
    fprintf(stderr, "Time for running subdividing for %d times is \t%f secs.\n", totalRuns, ltime2-ltime1);

    _ftime(&thistime);    
    ltime1=thistime.time + thistime.millitm * 1e-3;
	for(i=0; i<totalRuns; i++)
	{
		Interpolate(0);
		Interpolate(1);
		//Interpolate(2);
	}
    _ftime(&thistime);    
    ltime2=thistime.time + thistime.millitm * 1e-3;
    fprintf(stderr, "Time for running interpolation for %d times is \t%f secs.\n", totalRuns, ltime2-ltime1);
}
	return true;
}
#endif	

	Subdivide(nivelo);
	Interpolate(nivelo);
	return true;
}

void CCSubdivsurf::Subdivide(int nivelo)
{

	CCMesh *oldlevel = submeshes(nivelo);
	CCMesh *newlevel = submeshes(nivelo + 1);

	// The new mesh has the even vertices + the new face vertices + the new edge vertices
	int newnumverts = oldlevel->numverts() + oldlevel->numfaces() + oldlevel->numedges();

	// Each new edge-vertex generates 2 new faces
	//   (As proof, each new edge from the edge-vertex to the 2 face-vertices forms the
	//    lead edge (in counterclockwise order) of one of the new faces)
	int newnumfaces = 2*oldlevel->numedges();

	newlevel->VertList.Create(newnumverts);
	newlevel->FaceList.Create(newnumfaces);
	newlevel->VertFaceNeighbors.Create(newnumverts);
	newlevel->VertFaceNeighborsVert.Create(newnumverts);

	int numoldfaces = oldlevel->numfaces();
	int numoldverts = oldlevel->numverts();

	// Copy even vertices into new vertex list
	// NO, begin by zeroing out the even vertices...use parent's even vertices
	//    (which have the same indices) for computation, and use the new ones
	//    as accumulators for computing the new even vertices.
	// Formula is neweven = oldweight*v_old + neighborweight*(sum[4*newedgeneighbors - newfaceneighbors])
	//    Derivation of this is in Thall diss, App.A.
	int evenvindex;
	for (evenvindex = 0; evenvindex < numoldverts; evenvindex++) 
	{
		newlevel->VertList[evenvindex] = oldlevel->VertList[evenvindex];
		newlevel->VertList[evenvindex].is_oddvert = false;
	}

	DbVector3 newloc;
	DbVector3 newnorm;
	Medcoord newuvt;
	//double newradius;

	int polydex;

	// Create new face vertices at centroids of their polygons
	for(polydex = 0; polydex < numoldfaces; polydex++) 
	{

		CCPolygon *thispoly = oldlevel->FaceList(polydex);

		newlevel->VertList[numoldverts + polydex].setvalenceandoddflag(thispoly->numsides(), true);

		thispoly->newfacevert = numoldverts + polydex;
		thispoly->newedgeverts.Create(thispoly->numsides()); // go ahead and initialize newedgevert list

		// Now, subtract values of the new face vertex into loc and norm of surrounding even vertices
	}

	// For each edge of polygon, create new edge vertex at averaged position of
	//    of old edge endpoints and new face vertices at centroids
	// For r and uvt interpolation, just average old endpoints
	CCVindex currentnewvert = numoldverts + numoldfaces;
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		CCPolygon *thispoly = oldlevel->FaceList(polydex);

		// Now, need those pointers from polygon to new vertices at next level for
		//   face and for edges
		for (int edgedex = 0; edgedex < thispoly->numsides(); edgedex++) 
		{
			CCPindex neighbordex = thispoly->Pneighbors[edgedex];

			CCPolygon *neighborpoly = oldlevel->FaceList(neighbordex);

			// If the across-edge neighbor polygon has a lower index, its edges
			//   have already been computed, so just store pointer to new edgevertex in parent Polygon
			if (polydex > neighbordex) 
			{
				thispoly->newedgeverts[edgedex] = neighborpoly->newedgeverts[thispoly->Pneighboredges[edgedex]];
			}

			// else create a new vertex by averaging the new face and old edge-endpoints
			else 
			{

				CCVertex *facevert1 = newlevel->VertList(thispoly->newfacevert);
				CCVertex *facevert2 = newlevel->VertList(neighborpoly->newfacevert);
				CCVertex *edgevert1 = oldlevel->VertList(thispoly->Pvertices[edgedex]);
				CCVertex *edgevert2 = oldlevel->VertList(thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]);
				/*
				newloc = edgevert1->loc + edgevert2->loc + facevert1->loc + facevert2->loc ;
				newnorm = edgevert1->norm + edgevert2->norm + facevert1->norm + facevert2->norm;
				newuvt = edgevert1->uvt + edgevert2->uvt;
				newradius = edgevert1->rad + edgevert2->rad;
				newloc /= 4;
				newnorm /= 4;
				newradius /= 2;
				newuvt /= 2;
				newlevel->VertList[currentnewvert].setvals(newloc, newnorm, newradius, newuvt, 4, true);
				*/
				newlevel->VertList[currentnewvert].setvalenceandoddflag(4, true);
				thispoly->newedgeverts[edgedex] = currentnewvert++;
				// Now, add 4x the loc and norm values to each of the two edgeverts
				//   This happens just once, when this vertex is created
				CCVertex *thisVertex1 = newlevel->VertList(thispoly->Pvertices[edgedex]);
				CCVertex *thisVertex2 = newlevel->VertList(thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]);
				//thisVertex1->loc += newloc*4.0;
				//thisVertex1->norm += newnorm*4.0;
				//thisVertex2->loc += newloc*4.0;
				//thisVertex2->norm += newnorm*4.0;

			}
		}
	}

	// Now, form new polygons and compute all connectivity
	//    Do this in two stages:
	//    1) create new polygons and assign their
	//            numbers to parent polygons' childpolygon lists
	//    2) Use this info to set neighbor and neighboredge lists
	CCPindex newpolydex = 0;
	for(polydex = 0; polydex < numoldfaces; polydex++) 
	{
		CCPolygon *thispoly = oldlevel->FaceList(polydex);
		int nsides = thispoly->numsides();

		thispoly->childPolygons.Create(nsides);

		for (int vdex = 0; vdex < nsides; vdex++) 
			thispoly->childPolygons[vdex] = newpolydex++;
	}

	newpolydex = 0;
	int totalsides = 0;
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		CCPolygon *thispoly = oldlevel->FaceList(polydex);
		int nsides = thispoly->numsides();

		// count number of edges for total
		totalsides += nsides*4;

		int firstnewpolydex = newpolydex;
		int lastnewpolydex = newpolydex + nsides - 1;

		for (evenvindex = 0; evenvindex < nsides; evenvindex++) {

			CCPolygon *newpoly = newlevel->FaceList(newpolydex);
			newpoly->Pvertices.Create(4);

			newpoly->Pvertices[0] = thispoly->newfacevert;
			newpoly->Pvertices[1] = thispoly->newedgeverts[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pvertices[2] = thispoly->Pvertices[evenvindex];
			newpoly->Pvertices[3] = thispoly->newedgeverts[evenvindex];

			// Now do polygon neighbors

			int facenum;
			newpoly->Pneighbors.Create(4);
			// The edge number is clockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex - 1 < 0 ? nsides - 1 : evenvindex - 1);
			newpoly->Pneighbors[0] = thispoly->childPolygons[facenum];

			// The edge number of the neighbor polygon gives the polygon-index
			facenum = thispoly->Pneighboredges[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pneighbors[1] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[facenum];

			// The edge number + 1 MOD neighbors-numsides of the neighboredge gives the polygon-index
			facenum = INTMOD(thispoly->Pneighboredges[evenvindex] + 1,
				             oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].numsides());
			newpoly->Pneighbors[2] = oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].childPolygons[facenum];

			// The edge number is counterclockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex + 1 > nsides - 1 ? 0 : evenvindex + 1);
			newpoly->Pneighbors[3] = thispoly->childPolygons[facenum];

			// -- if all polygons start at the new face vertex, then the edge#'s are all 1
			//    or 2 for neighboring polygons.  (2 facing 1, 1 facing 2)
			newpoly->Pneighboredges.Create(4);
			newpoly->Pneighboredges[0] = 3;
			newpoly->Pneighboredges[1] = 2;
			newpoly->Pneighboredges[2] = 1;
			newpoly->Pneighboredges[3] = 0;

			newpoly->parentPolygon = polydex;
			newpoly->newfacevert = NULLVERTDEX;
//			newpoly->printvals();
//			fprintf(stderr, "\n");
			newpolydex++;
		}
	}

	newlevel->nedges = totalsides/2;

	/*
	// Now, average the even vertices according to their connectivity
	// ----no, this can be done on the fly.  Zero all nivelo+1 even vertices, and
	// add 4x the new edge values and 1x each new face vertex.
	//    Once here, weight the computed sum and add the weighted parent vertex weight
	// Use the weights:
	//    1/(n^2) for the face-edge sum, and (n-3)/n for the old vertex
	// (See App.A of Thall Diss for the derivation of this.)
	*/

	// Now, compute all the connectivity relationships for the CCmesh, as necessary
	// Can probably do this on the fly, too.
	// For each polygon in old mesh:
	//    new children are neighborhood of newfacevert
	//    two new children and their 2-edge and 1-edge neighbors are neighborhood of new edge
	//    new child and its 1-edge neighbor and ITS 1-edge neighbors are
	//       neighborhood of new even vertex
	//    Vertex 2 of each new polygon is old even vertex
	// If VertFaceNeighbors[vdex] has been initialized, then already computed!
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		CCPolygon *oldpoly = oldlevel->FaceList(polydex);

		CCVindex newvertdex = oldpoly->newfacevert;

		// newvalence is number of edges of old polygon, which will be used
		//    all over the place
		int newvalence = newlevel->VertList[newvertdex].v_valence;

		// Compute neighborhood of new face vertex
		newlevel->VertFaceNeighbors[newvertdex].Create(newvalence);
		newlevel->VertFaceNeighborsVert[newvertdex].Create(newvalence);

		CCPindex newpolydex;
		for (newpolydex = 0; newpolydex < newvalence; newpolydex++) {
			newlevel->VertFaceNeighbors[newvertdex][newpolydex] = oldpoly->childPolygons[newpolydex];
			newlevel->VertFaceNeighborsVert[newvertdex][newpolydex] = 0;
		}

		// Compute neighborhood of each new edge vertex if VertFaceNeighbors[] not already initialized
		for (int newoddedge = 0; newoddedge < newvalence; newoddedge++) {
			newvertdex = oldpoly->newedgeverts[newoddedge];
			if (newlevel->VertFaceNeighbors[newvertdex].Length() == 0) 
			{
				newlevel->VertFaceNeighbors[newvertdex].Create(4);
				newlevel->VertFaceNeighborsVert[newvertdex].Create(4);

				CCPindex pdexA, pdexB, pdexC, pdexD;
				pdexA = oldpoly->childPolygons[INTMOD(newoddedge + 1, newvalence)];
				pdexB = oldpoly->childPolygons[newoddedge];
				pdexC = newlevel->FaceList[pdexB].Pneighbors[2];
				pdexD = newlevel->FaceList[pdexA].Pneighbors[1];

				newlevel->VertFaceNeighbors[newvertdex][0] = pdexA;
				newlevel->VertFaceNeighbors[newvertdex][1] = pdexB;
				newlevel->VertFaceNeighbors[newvertdex][2] = pdexC;
				newlevel->VertFaceNeighbors[newvertdex][3] = pdexD;
				newlevel->VertFaceNeighborsVert[newvertdex][0] = 1;
				newlevel->VertFaceNeighborsVert[newvertdex][1] = 3;
				newlevel->VertFaceNeighborsVert[newvertdex][2] = 1;
				newlevel->VertFaceNeighborsVert[newvertdex][3] = 3;
			}
		}

		// Compute neighborhood of each new even vertex if VertFaceNeighbors[] not already initialized
		for (int childpolydex = 0; childpolydex < newvalence; childpolydex++) 
		{
			CCPindex newface = oldpoly->childPolygons[childpolydex];
			newvertdex = newlevel->FaceList[newface].Pvertices[2];

			if (newlevel->VertFaceNeighbors[newvertdex].Length() == 0) {

				int vertvalence = newlevel->VertList[newvertdex].v_valence;
				newlevel->VertFaceNeighbors[newvertdex].Create(vertvalence);
				newlevel->VertFaceNeighborsVert[newvertdex].Create(vertvalence);

				CCPindex pdexE = newface;

				newlevel->VertFaceNeighbors[newvertdex][0] = pdexE;
				newlevel->VertFaceNeighborsVert[newvertdex][0] = 2;

				for (int facedex = 1; facedex < vertvalence; facedex++) {

					pdexE = newlevel->FaceList[pdexE].Pneighbors[1];
					newlevel->VertFaceNeighbors[newvertdex][facedex] = pdexE;
					newlevel->VertFaceNeighborsVert[newvertdex][facedex] = 2;
				}
			}
		}
	}
}


//	Updates - ::Interpolate(...) to support the subdivision of tube m-rep mesh
//				extra codes to deal with the new uvt space of the tube
void CCSubdivsurf::Interpolate(int nivelo)
{
	//std::cout << "Interpolate" << std::endl;
	CCMesh *oldlevel = submeshes(nivelo);
	CCMesh *newlevel = submeshes(nivelo + 1);

	// The new mesh has the even vertices + the new face vertices + the new edge vertices
	int newnumverts = oldlevel->numverts() + oldlevel->numfaces() + oldlevel->numedges();

	// Each new edge-vertex generates 2 new faces
	//   (As proof, each new edge from the edge-vertex to the 2 face-vertices forms the
	//    lead edge (in counterclockwise order) of one of the new faces)
	int newnumfaces = 2*oldlevel->numedges();

	int numoldfaces = oldlevel->numfaces();
	int numoldverts = oldlevel->numverts();

	// Copy even vertices into new vertex list
	// NO, begin by zeroing out the even vertices...use parent's even vertices
	//    (which have the same indices) for computation, and use the new ones
	//    as accumulators for computing the new even vertices.
	// Formula is neweven = oldweight*v_old + neighborweight*(sum[4*newedgeneighbors - newfaceneighbors])
	//    Derivation of this is in Thall diss, App.A.
	int evenvindex;
	for (evenvindex = 0; evenvindex < numoldverts; evenvindex++) 
	{
		newlevel->VertList[evenvindex].loc = NULLVEC;
		newlevel->VertList[evenvindex].norm = NULLVEC;

		newlevel->VertList[evenvindex].rad=0;
	}

	DbVector3 newloc;
	DbVector3 newnorm;
	Medcoord newuvt;
	double newradius;

	int polydex;

	// Create new face vertices at centroids of their polygons
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{

		CCPolygon *thispoly = oldlevel->FaceList(polydex);

		newloc.set(0.0, 0.0, 0.0);
		newnorm.set(0.0, 0.0, 0.0);
		newuvt.setvals(0.0, 0.0, 0.0);
		newradius = 0.0;
		int v_index;
		for (v_index = 0; v_index < thispoly->numsides(); v_index++) {

			CCVertex *thisVertex = oldlevel->VertList(thispoly->Pvertices[v_index]);

			// Sum the node values
			newloc += thisVertex->loc;
			newnorm += thisVertex->norm;
			newradius += thisVertex->rad;
		}
		newloc /= thispoly->numsides();
		newnorm /= thispoly->numsides();
		newradius /= thispoly->numsides();

		Medcoord* uvt	= new Medcoord[thispoly->numsides()];
		for (v_index = 0; v_index < thispoly->numsides(); v_index++) {
			uvt[v_index] = oldlevel->VertList(thispoly->Pvertices[v_index])->uvt;
		}
		newuvt.average(thispoly->numsides(), uvt,
			(isTubeMesh) ? Medcoord::TubeMedCoord : Medcoord::SlabMedCoord);
		delete[] uvt;

		newlevel->VertList[numoldverts + polydex].setvals(newloc, newnorm, newradius,
			                                              newuvt, thispoly->numsides(), true);

		// Now, subtract values of the new face vertex into loc and norm of surrounding even vertices
		for (v_index = 0; v_index < thispoly->numsides(); v_index++) 
		{
			CCVertex *thisVertex = newlevel->VertList(thispoly->Pvertices[v_index]);
			thisVertex->loc -= newloc;
			thisVertex->norm -= newnorm;
			thisVertex->rad -= newradius;
			// NOTE:  if you make this an addition instead of a subtraction,
			//   it screws everything.  :)  Makes some nice sea anemonies, though
		}
	}

	// For each edge of polygon, create new edge vertex at averaged position of
	//    of old edge endpoints and new face vertices at centroids
	// For r and uvt interpolation, just average old endpoints
	CCVindex currentnewvert = numoldverts + numoldfaces;
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		CCPolygon *thispoly = oldlevel->FaceList(polydex);

		// Now, need those pointers from polygon to new vertices at next level for
		//   face and for edges
		for (int edgedex = 0; edgedex < thispoly->numsides(); edgedex++) 
		{
			CCPindex neighbordex = thispoly->Pneighbors[edgedex];

			CCPolygon *neighborpoly = oldlevel->FaceList(neighbordex);

			// If the across-edge neighbor polygon has a lower index, its edges
			//   have already been computed, so just store pointer to new edge vertex in parent Polygon
			if (polydex > neighbordex) 
			{
				thispoly->newedgeverts[edgedex] = neighborpoly->newedgeverts[thispoly->Pneighboredges[edgedex]];
			}
			// else create a new vertex by averaging the new face and old edge-endpoints
			else 
			{
				CCVertex *facevert1 = newlevel->VertList(thispoly->newfacevert);
				CCVertex *facevert2 = newlevel->VertList(neighborpoly->newfacevert);
				CCVertex *edgevert1 = oldlevel->VertList(thispoly->Pvertices[edgedex]);
				CCVertex *edgevert2 = oldlevel->VertList(thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]);

				newloc = edgevert1->loc + edgevert2->loc + facevert1->loc + facevert2->loc ;
				newnorm = edgevert1->norm + edgevert2->norm + facevert1->norm + facevert2->norm;
				newradius = edgevert1->rad + edgevert2->rad + facevert1->rad + facevert2->rad;

				newloc /= 4;
				newnorm /= 4;
				newradius /= 4;

				Medcoord uvt[2];
				uvt[0]	= edgevert1->uvt;
				uvt[1]	= edgevert2->uvt;
				newuvt.average(2, uvt, (isTubeMesh) ? Medcoord::TubeMedCoord : Medcoord::SlabMedCoord);
/*
				if(!isTubeMesh) {
					// for quad m-rep mesh
					newuvt = edgevert1->uvt + edgevert2->uvt;
					newuvt /= 2;
				}
				else {
					// for tube m-rep mesh
					newuvt.u = edgevert1->uvt.u + edgevert2->uvt.u;
					newuvt.t = edgevert1->uvt.t + edgevert2->uvt.t;
					newuvt.u /= 2;
					newuvt.t /= 2;

					double minV, maxV;
					if(edgevert1->uvt.v<edgevert2->uvt.v) {
						minV = edgevert1->uvt.v;
						maxV = edgevert2->uvt.v;
					}
					else {
						minV = edgevert2->uvt.v;
						maxV = edgevert1->uvt.v;
					}
					if(minV >= 0.0) {
						// no end vertices were involved.
						if((1.0-maxV) + minV < maxV - minV) {
							// map minV to 1.0+minV
							minV += 1.0;
						}
						newuvt.v = minV+maxV;
						if(newuvt.v >= 1.0) {
							newuvt.v -= 1.0;
						}
					}
					else {
						newuvt.v = maxV;
					}
				}
*/
				newlevel->VertList[currentnewvert].setvals(newloc, newnorm, newradius, newuvt, 4, true);
				thispoly->newedgeverts[edgedex] = currentnewvert++;

				// Now, add 4x the loc and norm values to each of the two edgeverts
				//   This happens just once, when this vertex is created
				CCVertex *thisVertex1 = newlevel->VertList(thispoly->Pvertices[edgedex]);
				CCVertex *thisVertex2 = newlevel->VertList(thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]);

				thisVertex1->loc += newloc*4.0;
				thisVertex1->norm += newnorm*4.0;
				thisVertex2->loc += newloc*4.0;
				thisVertex2->norm += newnorm*4.0;

				thisVertex1->rad+=newradius*4.0;
				thisVertex2->rad+=newradius*4.0;
			}
		}
	}

	// Now, average the even vertices according to their connectivity
	// ----no, this can be done on the fly.  Zero all nivelo+1 even vertices, and
	// add 4x the new edge values and 1x each new face vertex.
	//    Once here, weight the computed sum and add the weighted parent vertex weight
	// Use the weights:
	//    1/(n^2) for the face-edge sum, and (n-3)/n for the old vertex
	// (See App.A of Thall Diss for the derivation of this.)
	for (evenvindex = 0; evenvindex < numoldverts; evenvindex++) 
	{
		double thisval = newlevel->VertList[evenvindex].v_valence;

		double sumweight, oldweight;

		sumweight = 1.0/(thisval*thisval);
		oldweight = (thisval - 3.0)/thisval;

		newlevel->VertList[evenvindex].loc *= sumweight;
		newlevel->VertList[evenvindex].norm *= sumweight;

		newlevel->VertList[evenvindex].loc += oldlevel->VertList[evenvindex].loc * oldweight;
		newlevel->VertList[evenvindex].norm += oldlevel->VertList[evenvindex].norm * oldweight;

		newlevel->VertList[evenvindex].rad *= sumweight;
		newlevel->VertList[evenvindex].rad += oldlevel->VertList[evenvindex].rad * oldweight;
	}

	// Lastly, compute the limit mesh
	computelimitmesh(nivelo + 1);
	// This will cause normal-interpolation to be used at level-2 and higher subdivision
	if(nivelo==1 && NORMINTERP) 
	{
		perturbleveltwo();
		computelimitmesh(nivelo+1);
	}
	//return true;
}


/********************************************************************************/
/* computelimitmesh() -- Compute the limit mesh for the given submesh level
/*							--returns false if level not computed
/********************************************************************************/
bool CCSubdivsurf::computelimitmesh(int nivelo)
{
	CCMesh *thislevel = submeshes(nivelo);

	int nverts = thislevel->numverts();
	if (nverts == 0)
		return false;

	DbVector3 c2tangent, c3tangent;
	double e_r, f_r;
	int vertvalence;
	// For each vertex in thislevel, send it to its limit pposition
	//    and compute the limit normal as well
	for (CCVindex vdex = 0; vdex < nverts; vdex++) 
	{

		CCVertex *thisVert = thislevel->VertList(vdex);

		vertvalence = thisVert->v_valence;

		if( vertvalence > MAXVALENCE ) {
			printf("CCSubdivsurf encountered a vertex with valency %d, while the tables"
				"are designed only for a maximum of %d. Increase MAXVALENCE in CCSubdivsurf.h"
				"and recompile.\n", vertvalence, MAXVALENCE );
			assert(false);
		}

		// for each polygonal neighbor, form the weighted sum of the e_i edge
		//    and the f_i face (possibly computed, for a non-quad) into the
		//    limitloc and limitnorm DbVectors
		thisVert->limitloc = thisVert->loc * (vertvalence*vertvalence);
		thisVert->limitRad = thisVert->rad * (vertvalence*vertvalence);
		c2tangent.set(0.0, 0.0, 0.0);
		c3tangent.set(0.0, 0.0, 0.0);

		CCPindex pdex;
		for (pdex = 0; pdex < vertvalence; pdex++) 
		{

			CCPindex thispolydex = thislevel->VertFaceNeighbors[vdex][pdex];
			CCPolygon *thispoly = thislevel->FaceList(thispolydex);

			int polyvdex = thislevel->VertFaceNeighborsVert[vdex][pdex];

			DbVector3 e_j, f_j;
			int pnumsides = thispoly->numsides();

			// The edge vertex is the next vertex counterclockwise from vpolydex
			CCVertex *e_pdex = thislevel->VertList(thispoly->Pvertices[INTMOD(polyvdex + 1, pnumsides)]);
			e_j = e_pdex->loc;
			e_r = e_pdex->rad;

			// If pnumsides == 4, 2nd vertex counterclockwise from vpolydex is face vertex
			//    else have to regularize and just compute regloc
			if (pnumsides == 4) 
			{
				CCVertex *f_pdex = thislevel->VertList(thispoly->Pvertices[INTMOD(polyvdex + 2, pnumsides)]);
				f_j = f_pdex->loc;
				f_r = f_pdex->rad;
			}
			else
			{
				f_j = thislevel->regularize(thispolydex, polyvdex);
				f_r = thislevel->regularizeR(thispolydex, polyvdex);
			}

			// Now, add weighted values to limitloc and {c2,c3}tangent
			thisVert->limitloc += (e_j * 4.0) + f_j;
			thisVert->limitRad += (e_r * 4.0) + f_r;

			c2tangent += (e_j * Etable[vertvalence][pdex + 1]) + (f_j * Ftable[vertvalence][pdex + 1]);
			c3tangent += (e_j * Etable[vertvalence][pdex    ]) + (f_j * Ftable[vertvalence][pdex    ]);
		}

		thisVert->limitloc /= (double) vertvalence*(vertvalence + 5);
		thisVert->limitRad /= (double) (vertvalence*(vertvalence + 5));
		thisVert->limitnorm = c2tangent.cross(c3tangent);
		thisVert->limitnorm.selfnormalize();
	}

	return true;
}

/********************************************************************************/
/* computeIISmesh() -- Iteratively compute a submeshes(0) that interpolates the
/*	 fathermesh in the limit.  numiter gives the number of iterations.
/* --returns false if level 0 not set not computed
/********************************************************************************/
bool CCSubdivsurf::computeIISmesh(int numiter)
{
	CCMesh *IISmesh = submeshes(0);
	CCMesh *dadmesh = &fathermesh;

	int nverts = IISmesh->numverts();
	if (nverts == 0)
		return false;

	double e_r, f_r;

	// For each iteration, perturb vertex submeshes[0].VertList[j] from their
	//    initialized values (same as values in fathermesh)
	for (int icnt = 0; icnt < numiter; icnt++) {
		int vertvalence;

		// For each vertex in thislevel, send it to its inverse limit position
		// --This can't be done in-place, so use limitloc as a temporary variable
		// (should experiment---perhaps CAN do in-place)
		// REWRITE:  can do in-place, and thus, get Gauss-Seidel which should
		//    converge faster than Jacobi, and is faster per iteration
		for (CCVindex vdex = 0; vdex < nverts; vdex++) {

			CCVertex *thisVert = IISmesh->VertList(vdex);
			CCVertex *dadvert = dadmesh->VertList(vdex);
			vertvalence = thisVert->v_valence;
			DbVector3 iterant;
			double iterantR;

			// Formula:  as described in ThallDiss Ch. 4
			//    v_j^[i+1] = 1/2 * [ v_j^[i]  + (5+valence)/valence * v_limit_j
			//                                 - 1/(valence^2) * SUM[4*e_k + f_k]
	        //						]
			thisVert->limitloc = thisVert->loc;
			thisVert->limitloc += dadvert->loc * (vertvalence + 5.0)/((double) vertvalence);
			thisVert->limitRad=thisVert->rad;
			thisVert->limitRad += dadvert->rad * (vertvalence + 5.0)/((double) vertvalence);

			iterant.set(0.0, 0.0, 0.0);
			iterantR=0;

			// for each polygonal neighbor, form the weighted sum of the e_i edge
			//    and the f_i face (possibly computed, for a non-quad) into the
			//    limitloc and limitnorm DbVectors
			CCPindex pdex;
			for (pdex = 0; pdex < vertvalence; pdex++) {

				CCPindex thispolydex = IISmesh->VertFaceNeighbors[vdex][pdex];
				CCPolygon *thispoly = IISmesh->FaceList(thispolydex);

				int polyvdex = IISmesh->VertFaceNeighborsVert[vdex][pdex];

				DbVector3 e_j, f_j;
				int pnumsides = thispoly->numsides();

				// The edge vertex is the next vertex counterclockwise from polyvdex
				CCVertex *e_pdex = IISmesh->VertList(thispoly->Pvertices[INTMOD(polyvdex + 1, pnumsides)]);
				e_j = e_pdex->loc;
				e_r = e_pdex->rad;

				// If pnumsides == 4, 2nd vertex counterclockwise from vpolydex is face vertex
				//    else have to regularize and just compute regloc
				if (pnumsides == 4) {
					CCVertex *f_pdex = IISmesh->VertList(thispoly->Pvertices[INTMOD(polyvdex + 2, pnumsides)]);
					f_j = f_pdex->loc;
					f_r = f_pdex->rad;
				}
				else
				{
					f_j = IISmesh->regularize(thispolydex, polyvdex);
					f_r = IISmesh->regularizeR(thispolydex, polyvdex);
				}

				// Now, add weighted values to limitloc and {c2,c3}tangent
				iterant += (e_j * 4.0) + f_j;
				iterantR += (e_r*4.0)+f_r;
			}

			iterant /= (double) vertvalence*vertvalence;
			thisVert->limitloc -= iterant;
			thisVert->limitloc /= 2.0;
			iterantR/=(double)(vertvalence*vertvalence);
			thisVert->limitRad-=iterantR;
			thisVert->limitRad/=2.0;
		}
		// After each complete iteration through all vertices,
		//    copy values back from limitloc to loc
		for (CCVindex vdex2 = 0; vdex2 < nverts; vdex2++) {
			CCVertex *thisVert = IISmesh->VertList(vdex2);

			thisVert->loc = thisVert->limitloc;
			thisVert->rad = thisVert->limitRad;
		}
	}
	return true;	
}

bool CCSubdivsurf::WriteCMesh(int lvl, int vNumBase)
{
	FILE *fp;
	char filename[20];
	int vNums, fNums;

	if(lvl<0||lvl>NUMLEVELS-1) 
		return false;

	CCMesh *thisMesh=SubMeshPointer(lvl);

	vNums=thisMesh->VertList.Length();
	fNums=thisMesh->FaceList.Length();

	sprintf(filename, "%s%d%d%s", "blending", vNums, fNums, ".dat");

	fp=fopen(filename, "wb");

	fprintf(fp, "%d\n", vNums);
	for(int vdex=0; vdex<vNums; vdex++)
	{
		CCVertex *thisVert=thisMesh->VertList(vdex);
		fprintf(fp, "%f, %f, %f\n", thisVert->limitloc.x(), thisVert->limitloc.y(), thisVert->limitloc.z());
		fprintf(fp, "%f, %f, %f\n", thisVert->limitnorm.x(), thisVert->limitnorm.y(), thisVert->limitnorm.z());
		fprintf(fp, "%f, %f, %f, %f\n", thisVert->rad, thisVert->uvt.u, thisVert->uvt.v, thisVert->uvt.t);
	}

	fprintf(fp, "%d\n", fNums);
	for(int fdex=0; fdex<fNums; fdex++)
	{
		CCPolygon *thisPoly=thisMesh->FaceList(fdex);
		fprintf(fp, "%d", thisPoly->numsides());
		for(int fvdex=0; fvdex<thisPoly->numsides(); fvdex++)
			fprintf(fp, ", %d", thisPoly->Pvertices[fvdex]+vNumBase);
		fprintf(fp, "\n");
	}
	fclose(fp);
	return true;
}

/********************************************************************************/
/* init() -- Pass in a Bpoint list of initial boundary involutes and an array
/*		of Pindexlists giving the polygons associated with the pntlist vertices.
/*		This will load fathermesh with the vertices and polygons, and compute
/*		all connectivity information for the initial CCMesh objects.
/********************************************************************************/

typedef struct {
	int facenum, edgenum;
} edge_element;

typedef std::map<int, edge_element> EdgeSet;

void CCSubdivsurf::InitTopology(int numverts, int numfaces, Bpoint *pntlist, PindexList *faces)
{
	int vdex, pdex;

	// Create a map for the edges
	EdgeSet Eset;
	// Create another for the vertices, key is vertex number
	//    edgenum gives order of vertex in face
	EdgeSet Vset;

	// Load data into FaceLists---compute VertList valences at same time
	//  First, allocate space for all faces
	for (pdex = 0; pdex < numfaces; pdex++) {

		int nsides = faces[pdex].Length();

		// Create FaceList fields
		fathermesh.FaceList[pdex].Pvertices.Create(nsides);
		fathermesh.FaceList[pdex].Pneighbors.Create(nsides);
		fathermesh.FaceList[pdex].Pneighboredges.Create(nsides);
		fathermesh.FaceList[pdex].parentPolygon = NULLPOLYDEX;
		fathermesh.FaceList[pdex].newfacevert = NULLVERTDEX;
	}
	
	for(vdex=0; vdex<numverts; vdex++)
		fathermesh.VertList[vdex].v_valence = 0;

	//  Next, load data
	for (pdex = 0; pdex < numfaces; pdex++) {

		int nsides = faces[pdex].Length();

		for (int sidecnt = 0; sidecnt < nsides; sidecnt++) {

			// Set CCVindex of FaceList to index in faces[]
			CCVindex thisVertdex = faces[pdex][sidecnt];
			fathermesh.FaceList[pdex].Pvertices[sidecnt] = thisVertdex;

			// Increment valence counter of VertList[thisVertdex]
			fathermesh.VertList[thisVertdex].v_valence += 1;

			// if vertex hasn't been added to Vset, add polygon # and vindex in poly
			//    to Vset.  Just use one.  FASTER if rewrite PSEUDOSET for this.
			//    Can make O(n) instead of O(n logn)
			if (Vset.find(thisVertdex) == Vset.end()) {
				edge_element v_el;
				v_el.facenum = pdex;
				v_el.edgenum = sidecnt;
				Vset[thisVertdex] = v_el;
			}

			// FIGURE HOW TO DO THIS:  add this edge to set along with
			//    name of Face and location in Face
			//    If edge already in set, include Face and location in Face
			//    of this polygon as well.
			// WHEN ALL DONE, go through all edges and add Pneighbor and
			//    Pneighborfields to each edge in set

			// OKAY:  for edge vb-ve, where vb and ve are vertex indices
			//    let v1 = vb > ve ? vb : ve
			//        v2 = vb > ve ? ve : vb
			// Let key = v1 * numverts + v2
			// To retrieve, v1 = key/numverts (with int divide)
			//              v2 = key % numverts
			CCVindex nextvertdex = faces[pdex][(sidecnt + 1) % nsides];
			int v1 = thisVertdex > nextvertdex ? thisVertdex : nextvertdex;
			int v2 = thisVertdex > nextvertdex ? nextvertdex : thisVertdex;
			int thiskey = v1*numverts + v2;

			// If key already present, copy values to this FaceList and
			//    and those of key in Eset.
			// Else add to Eset
			if (Eset.find(thiskey) != Eset.end()) {
				CCPindex neighface = Eset[thiskey].facenum;
				int neighedgenum = Eset[thiskey].edgenum;
				fathermesh.FaceList[pdex].Pneighbors[sidecnt] = neighface;
				fathermesh.FaceList[pdex].Pneighboredges[sidecnt] = neighedgenum;

				fathermesh.FaceList[neighface].Pneighbors[neighedgenum] = pdex;
				fathermesh.FaceList[neighface].Pneighboredges[neighedgenum] = sidecnt;
			}
			else {
				edge_element e_el;
				e_el.facenum = pdex;
				e_el.edgenum = sidecnt;
				Eset[thiskey] = e_el;
			}
		}
	}

	// Set numedges from distinct edges in Eset
	fathermesh.nedges = Eset.size();

	// Now, compute VertFaceNeighbors and VertFaceNeighborsVert
	fathermesh.VertFaceNeighbors.Create(numverts);
	fathermesh.VertFaceNeighborsVert.Create(numverts);

	for (vdex = 0; vdex < numverts; vdex++) {
		int thisvalence = fathermesh.VertList[vdex].v_valence;

		fathermesh.VertFaceNeighbors[vdex].Create(thisvalence);
		fathermesh.VertFaceNeighborsVert[vdex].Create(thisvalence);

		CCPindex firstpolydex = Vset[vdex].facenum;
		int firstvertpos = Vset[vdex].edgenum;

		fathermesh.VertFaceNeighbors[vdex][0] = firstpolydex;
		fathermesh.VertFaceNeighborsVert[vdex][0] = firstvertpos;

		// Add neighbor polygon of (firstvertpos - 1 % numedge) edge
		//    as next VertFaceNeighbors[] and compute position of
		//    vertex in it as next VertFaceNeighborsVert[]
		CCPindex nextpolydex;
		int nextvertpos;
		for (int pcnt = 1; pcnt < thisvalence; pcnt++) {

			int npedges = fathermesh.FaceList[firstpolydex].numsides();

			// neighbor polygon is along (firstvertpos - 1 % numedges) edge
			// position of vertex in neighbor is simply the number of the edge in neighbor
			nextpolydex = fathermesh.FaceList[firstpolydex].Pneighbors[INTMOD(firstvertpos - 1, npedges)];
			nextvertpos = fathermesh.FaceList[firstpolydex].Pneighboredges[INTMOD(firstvertpos - 1, npedges)];
			fathermesh.VertFaceNeighbors[vdex][pcnt] = nextpolydex;
			fathermesh.VertFaceNeighborsVert[vdex][pcnt] = nextvertpos;
			firstpolydex = nextpolydex;
			firstvertpos = nextvertpos;
		}
	}
}
void CCSubdivsurf::InitValue(int numverts, Bpoint *pntlist)
{
	// Load data into VertLists
	for(int vdex = 0; vdex < numverts; vdex++) {
		fathermesh.VertList[vdex].loc.set(pntlist[vdex].pnt);
		fathermesh.VertList[vdex].norm.set(pntlist[vdex].norm);
		fathermesh.VertList[vdex].rad = pntlist[vdex].rad;
		fathermesh.VertList[vdex].uvt.setvals(pntlist[vdex].u, pntlist[vdex].v, pntlist[vdex].t);
 		//fathermesh.VertList[vdex].v_valence = 0;		// NEED TO COMPUTE VALENCE, using as counter
		fathermesh.VertList[vdex].is_oddvert = false;
	}
}
void CCSubdivsurf::init(int numverts, int numfaces, Bpoint *pntlist, PindexList *faces)
{
	fathermesh.VertList.Create(numverts);
	fathermesh.FaceList.Create(numfaces);

	InitValue(numverts, pntlist);
	InitTopology(numverts, numfaces, pntlist, faces);
}

inline double AngleBetweenVectors(DbVector3 v1, DbVector3 v2)
{
	double dotP = v1.dot(v2);
	double c = acos(dotP);
#ifdef _MSC_VER
	if (_isnan(c)) {
		// For an unknown reason, this actually happens on Intel-based PC's running MS Visual Studio
		//printf("Error: Subdivsurf.cpp/AngleBetweenVectors, v1.dot(v2) is out of range [-1, 1]!\n");
		if(dotP > 1)
			return 0;
		else
			return M_PI;
	}
	else
#endif
		return c;
}


// What should the tolerance be on theta?  Chose 4 degrees, just arbitrarily.
//    Less than this, and no action taken.
//    to Interoplate normals in lvl2 to those in fathermesh
void CCSubdivsurf::perturbleveltwo()//create_normalinterpolating_mesh()
{
    //int nrows, ncols, nsrows, nscols;

	CCMesh *thisMesh = submeshes(2);
	// For involute vertices
	//    (1) find their limit positions and limit normals
	//    (2) compute a rotation about the limit position (as a quaternion)
	//           bringing it to location of true_nodenorm
	//    (3) rotate the involute vertex about the true_nodenorm by given rotation
	//    (4) rotate all neighbors of the involute vertex by this same amount

	CCVertex *thisVert, *adjVert;// = submeshes[nivelo].VertList(thispoly->Pvertices[vertdex]);
	int vIndex, vAdjIndex;
    //int trow, tcol, meshdex;
	//DbVector3 pval, pvalplusn;
	DbVector3 true_loc, true_norm;
	DbVector3 limit_loc, limit_norm;
	double theta;

	DbVector3 rotaxis;
	DbVector3 posvec;
	Quat rot_true2limitnorm;
	//subgridnode limitnode;

	// Save current color
	//glPushAttrib(GL_CURRENT_BIT);

	for(vIndex=0; vIndex<fathermesh.numverts(); vIndex++)
	{
		thisVert=thisMesh->VertList(vIndex);

		//oldlevel->VertList(thispoly->Pvertices[v_index]);
		limit_loc=thisVert->limitloc;
		limit_norm=thisVert->limitnorm;
		true_norm=fathermesh.VertList[vIndex].norm;
		limit_norm.selfnormalize();
		true_norm.selfnormalize();

		theta=AngleBetweenVectors(limit_norm, true_norm);
		//theta = acos(limit_norm.dot(true_norm));

		// the angle diffference between the actural/real normals is big enough
		if (theta > theta_tol) 
		{
			//Vertexneighbors *nearvert = &nearverts->topneighbors[meshdex];
			// we need to rotate the vertex and its neighbor vertices to approximate the real normal

			rotaxis = limit_norm.cross(true_norm);

			//rot_true2limitnorm.rot(rotaxis, theta);
			rot_true2limitnorm.setAxisAngle(Vector3D(rotaxis.X(), rotaxis.Y(), rotaxis.Z()), theta);

			posvec = thisVert->loc - limit_loc;

			//rot_true2limitnorm.rotate_vec(posvec);
			Vector3D posvecTmp(posvec.X(), posvec.Y(), posvec.Z()); 
			rot_true2limitnorm.rotateVector(posvecTmp);
			posvec.set(posvecTmp.getX(), posvecTmp.getY(), posvecTmp.getZ());

			thisVert->loc = posvec + limit_loc;

			for(int i=0; i<thisVert->v_valence; i++)
			{
				// Find out the adjacent neighbor vertices one by one
				//   and rotote them into along with the same angle
				CCPindex thisPolyIndex = thisMesh->VertFaceNeighbors[vIndex][i];
				CCPolygon *thisPoly = thisMesh->FaceList(thisPolyIndex);
				int polyVIndex=thisMesh->VertFaceNeighborsVert[vIndex][i];
				int pnumsides=thisPoly->numsides();
				// The adjacent vertex is the next vertex counterclockwise from vPolyIndex
				//CCVertex *e_pdex=thisMesh->VertList(thisPoly->Pvertices[INTMOD(polyVIndex+1, pnumsides)]);

				vAdjIndex=thisPoly->Pvertices[INTMOD(polyVIndex+1, pnumsides)];//thisMesh->VertFaceNeighborsVert[vIndex];

				adjVert=thisMesh->VertList(vAdjIndex);
				posvec = adjVert->loc-limit_loc;

				//rot_true2limitnorm.rotate_vec(posvec);
				Vector3D posvecTmp(posvec.X(), posvec.Y(), posvec.Z()); 
				rot_true2limitnorm.rotateVector(posvecTmp);
				posvec.set(posvecTmp.getX(), posvecTmp.getY(), posvecTmp.getZ());

				adjVert->loc=posvec+limit_loc;
			}
			/*
		for (pdex = 0; pdex < vertvalence; pdex++) {

			CCPindex thispolydex = thislevel->VertFaceNeighbors[vdex][pdex];
			CCPolygon *thispoly = thislevel->FaceList(thispolydex);

			int polyvdex = thislevel->VertFaceNeighborsVert[vdex][pdex];

			DbVector3 e_j, f_j;
			int pnumsides = thispoly->numsides();

			// The edge vertex is the next vertex counterclockwise from vpolydex
			CCVertex *e_pdex = thislevel->VertList(thispoly->Pvertices[INTMOD(polyvdex + 1, pnumsides)]);
			e_j = e_pdex->loc;

			// If pnumsides == 4, 2nd vertex counterclockwise from vpolydex is face vertex
			//    else have to regularize and just compute regloc
			if (pnumsides == 4) {
				CCVertex *f_pdex = thislevel->VertList(thispoly->Pvertices[INTMOD(polyvdex + 2, pnumsides)]);
				f_j = f_pdex->loc;
			}
			else
				f_j = thislevel->regularize(thispolydex, polyvdex);

			// Now, add weighted values to limitloc and {c2,c3}tangent
			thisVert->limitloc += (e_j * 4.0) + f_j;

			c2tangent += (e_j * Etable[vertvalence][pdex + 1]) + (f_j * Ftable[vertvalence][pdex + 1]);
			c3tangent += (e_j * Etable[vertvalence][pdex    ]) + (f_j * Ftable[vertvalence][pdex    ]);
		}
		*/


			/*
			// Draw perturbed neighborhood if DRAWNORMINTERP set
			if (DRAWNORMINTERP) {
				glColor3d(1.0, 0.0, 1.0);
				glBegin(GL_POINTS);
				submesh->topgrid[meshdex].nodeloc.glVertex();
				glEnd();
				glBegin(GL_LINE_LOOP);
				for (i = 0; i < nearvert->degree(); i++) {

					switch ((*nearvert)[i].whichmesh) {
					case TOP:
						submesh->topgrid[(*nearvert)[i].meshindex].nodeloc.glVertex();
						break;
					case BOTTOM:
						submesh->bottomgrid[(*nearvert)[i].meshindex].nodeloc.glVertex();
						break;
					case SIDE:
						submesh->sidegrid[(*nearvert)[i].meshindex].nodeloc.glVertex();
						break;
					}
				}
				glEnd();
			}
			*/
		}
	}
	// Save current color
	//glPopAttrib();
}




















/********************************************************************************/
/* initpyramid() -- create a Bpoint list and an array of PindexLists with the	*/
/*		the vertex and polygon information, then call init() to create the		*/
/*		fathermesh.																*/
/********************************************************************************/
void CCSubdivsurf::initpyramid()
{
	Bpoint vList[5];
	PindexList pList[5];

	vList[0].setvals(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.0, 1.0);
	vList[1].setvals(1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0);
	vList[2].setvals(0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.0, 1.0, 0.0);
	vList[3].setvals(-1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.5, -1.0, 0.0, 0.0);
	vList[4].setvals(0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.5, 0.0, -1.0, 0.0);

	pList[0].Create(4);
	pList[0][0] = 1;
	pList[0][1] = 4;
	pList[0][2] = 3;
	pList[0][3] = 2;

	pList[1].Create(3);
	pList[1][0] = 1;
	pList[1][1] = 0;
	pList[1][2] = 4;

	pList[2].Create(3);
	pList[2][0] = 4;
	pList[2][1] = 0;
	pList[2][2] = 3;

	pList[3].Create(3);
	pList[3][0] = 3;
	pList[3][1] = 0;
	pList[3][2] = 2;

	pList[4].Create(3);
	pList[4][0] = 2;
	pList[4][1] = 0;
	pList[4][2] = 1;

	init(5, 5, vList, pList);
}

/********************************************************************************/
/* inittritorus() -- create a Bpoint list and an array of PindexLists with the	*/
/*		the vertex and polygon information, then call init() to create the		*/
/*		fathermesh.																*/
/********************************************************************************/
void CCSubdivsurf::inittritorus()
{
	Bpoint vList[12];
	PindexList pList[12];

	// Vertex 0 is apex, Vertex 1-4 are cardinal points of axes in counterclockwise
	//    order -- NOTE: (uvt are wrong---these need to be moved to polygons, maybe)
	double sqrt3 = sqrt(3.0);
	double radval = sqrt(2.0)/2.0;		

	DbVector3 thisloc[12];
	DbVector3 thisnorm[12];

	thisloc[0].set(2.0, 0.0, 0.5);
	thisloc[1].set(-1.0, sqrt3, 0.5);
	thisloc[2].set(-1.0, -sqrt3, 0.5);
	thisloc[3].set(1.0, 0.0, 0.5);
	thisloc[4].set(-0.5, sqrt3/2.0, 0.5);
	thisloc[5].set(-0.5, -sqrt3/2.0, 0.5);
	thisloc[6].set(2.0, 0.0, 0.0);
	thisloc[7].set(-1.0, sqrt3, 0.0);
	thisloc[8].set(-1.0, -sqrt3, 0.0);
	thisloc[9].set(1.0, 0.0, 0.0);
	thisloc[10].set(-0.5, sqrt3/2.0, 0.0);
	thisloc[11].set(-0.5, -sqrt3/2.0, 0.0);

	// make normals vectors from each vertex to centroid of its corner
	DbVector3 centroidA = thisloc[0] + thisloc[3]
		                + thisloc[6] + thisloc[9];
	centroidA /= 4.0;
	DbVector3 centroidB = thisloc[1] + thisloc[4]
		                + thisloc[7] + thisloc[10];
	centroidB /= 4.0;
	DbVector3 centroidC = thisloc[2] + thisloc[5]
		                + thisloc[8] + thisloc[11];
	centroidC /= 4.0;

	thisnorm[0] = thisloc[0] - centroidA;
	thisnorm[1] = thisloc[1] - centroidB;
	thisnorm[2] = thisloc[2] - centroidC;
	thisnorm[3] = thisloc[3] - centroidA;
	thisnorm[4] = thisloc[4] - centroidB;
	thisnorm[5] = thisloc[5] - centroidC;
	thisnorm[6] = thisloc[6] - centroidA;
	thisnorm[7] = thisloc[7] - centroidB;
	thisnorm[8] = thisloc[8] - centroidC;
	thisnorm[9] = thisloc[9] - centroidA;
	thisnorm[10] = thisloc[10] - centroidB;
	thisnorm[11] = thisloc[11] - centroidC;

	vList[0].setvals(2.0, 0.0, 0.5,			thisnorm[0].x(), thisnorm[0].y(), thisnorm[0].z(), radval, 0.0, 0.0, 1.0);
	vList[1].setvals(-1.0, sqrt3, 0.5,		thisnorm[1].x(), thisnorm[1].y(), thisnorm[1].z(), radval, 1.0, 0.0, 1.0);
	vList[2].setvals(-1.0, -sqrt3, 0.5,		thisnorm[2].x(), thisnorm[2].y(), thisnorm[2].z(), radval, 2.0, 0.0, 1.0);
	vList[3].setvals(1.0, 0.0, 0.5,			thisnorm[3].x(), thisnorm[3].y(), thisnorm[3].z(), radval, 0.0, 1.0, 1.0);
	vList[4].setvals(-0.5, sqrt3/2.0, 0.5,	thisnorm[4].x(), thisnorm[4].y(), thisnorm[4].z(), radval, 1.0, 1.0, 1.0);
	vList[5].setvals(-0.5, -sqrt3/2.0, 0.5, thisnorm[5].x(), thisnorm[5].y(), thisnorm[5].z(), radval, 2.0, 1.0, 1.0);
	vList[6].setvals(2.0, 0.0, 0.0,			thisnorm[6].x(), thisnorm[6].y(), thisnorm[6].z(), radval, 0.0, 0.0, -1.0);
	vList[7].setvals(-1.0, sqrt3, 0.0,		thisnorm[7].x(), thisnorm[7].y(), thisnorm[7].z(), radval, 1.0, 0.0, -1.0);
	vList[8].setvals(-1.0, -sqrt3, 0.0,		thisnorm[8].x(), thisnorm[8].y(), thisnorm[8].z(), radval, 2.0, 0.0, -1.0);
	vList[9].setvals(1.0, 0.0, 0.0,			thisnorm[9].x(), thisnorm[9].y(), thisnorm[9].z(), radval, 0.0, 1.0, -1.0);
	vList[10].setvals(-0.5, sqrt3/2.0, 0.0, thisnorm[10].x(), thisnorm[10].y(), thisnorm[10].z(), radval, 1.0, 1.0, -1.0);
	vList[11].setvals(-0.5, -sqrt3/2.0, 0.0,thisnorm[11].x(), thisnorm[11].y(), thisnorm[11].z(), radval, 2.0, 1.0, -1.0);


	// Face 0-2 are the top in counterclockwise order,
	//    then 3-5 are the adjacent interior, then 6-8 the adjacent bottom,
	//    then 9-11 are the adjacent exterior
	pList[0].Create(4);
	pList[0][0] = 1;
	pList[0][1] = 4;
	pList[0][2] = 3;
	pList[0][3] = 2;

	pList[0].Create(4);
	pList[0][0] = 0;
	pList[0][1] = 1;
	pList[0][2] = 4;
	pList[0][3] = 3;

	pList[1].Create(4);
	pList[1][0] = 1;
	pList[1][1] = 2;
	pList[1][2] = 5;
	pList[1][3] = 4;

	pList[2].Create(4);
	pList[2][0] = 2;
	pList[2][1] = 0;
	pList[2][2] = 3;
	pList[2][3] = 5;

	pList[3].Create(4);
	pList[3][0] = 3;
	pList[3][1] = 4;
	pList[3][2] = 10;
	pList[3][3] = 9;

	pList[4].Create(4);
	pList[4][0] = 4;
	pList[4][1] = 5;
	pList[4][2] = 11;
	pList[4][3] = 10;


	pList[5].Create(4);
	pList[5][0] = 5;
	pList[5][1] = 3;
	pList[5][2] = 9;
	pList[5][3] = 11;

	pList[6].Create(4);
	pList[6][0] = 9;
	pList[6][1] = 10;
	pList[6][2] = 7;
	pList[6][3] = 6;

	pList[7].Create(4);
	pList[7][0] = 10;
	pList[7][1] = 11;
	pList[7][2] = 8;
	pList[7][3] = 7;

	pList[8].Create(4);
	pList[8][0] = 11;
	pList[8][1] = 9;
	pList[8][2] = 6;
	pList[8][3] = 8;

	pList[9].Create(4);
	pList[9][0] = 6;
	pList[9][1] = 7;
	pList[9][2] = 1;
	pList[9][3] = 0;

	pList[10].Create(4);
	pList[10][0] = 7;
	pList[10][1] = 8;
	pList[10][2] = 2;
	pList[10][3] = 1;

	pList[11].Create(4);
	pList[11][0] = 8;
	pList[11][1] = 6;
	pList[11][2] = 0;
	pList[11][3] = 2;

	init(12, 12, vList, pList);
}

/********************************************************************************/
/* initfigure8() -- create a Bpoint list and an array of PindexLists with the	*/
/*		the vertex and polygon information, then call init() to create the		*/
/*		fathermesh.																*/
/********************************************************************************/
void CCSubdivsurf::initfigure8()
{
	Bpoint vList[28];
	PindexList pList[28];

	vList[0].setvals(-2.0,   3.0, 0.5, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[1].setvals(-2.0,   0.0, 0.5, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);
	vList[2].setvals(-2.0,  -3.0, 0.5, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0);
	vList[3].setvals(2.0,   -3.0, 0.5, -1.0, 0.0, 0.0, 0.25, -1.0, 0.0, 0.0);
	vList[4].setvals(2.0,    0.0, 0.5, 0.0, -1.0, 0.0, 0.25, 0.0, -1.0, 0.0);
	vList[5].setvals(2.0,    3.0, 0.5, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[6].setvals(-1.0,   2.0, 0.5, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);
	vList[7].setvals(-1.0,   1.0, 0.5, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0);
	vList[8].setvals(1.0,    1.0, 0.5, -1.0, 0.0, 0.0, 0.25, -1.0, 0.0, 0.0);
	vList[9].setvals(1.0,    2.0, 0.5, 0.0, -1.0, 0.0, 0.25, 0.0, -1.0, 0.0);
	vList[10].setvals(-1.0, -1.0, 0.5, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[11].setvals(-1.0, -2.0, 0.5, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);
	vList[12].setvals(1.0,  -2.0, 0.5, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0);
	vList[13].setvals(1.0,  -1.0, 0.5, -1.0, 0.0, 0.0, 0.25, -1.0, 0.0, 0.0);
	vList[14].setvals(-2.0,  3.0, 0.0, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[15].setvals(-2.0,  0.0, 0.0, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);
	vList[16].setvals(-2.0, -3.0, 0.0, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0);
	vList[17].setvals(2.0,  -3.0, 0.0, -1.0, 0.0, 0.0, 0.25, -1.0, 0.0, 0.0);
	vList[18].setvals(2.0,   0.0, 0.0, 0.0, -1.0, 0.0, 0.25, 0.0, -1.0, 0.0);
	vList[19].setvals(2.0,   3.0, 0.0, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[20].setvals(-1.0,  2.0, 0.0, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);
	vList[21].setvals(-1.0,  1.0, 0.0, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0);
	vList[22].setvals(1.0,   1.0, 0.0, -1.0, 0.0, 0.0, 0.25, -1.0, 0.0, 0.0);
	vList[23].setvals(1.0,   2.0, 0.0, 0.0, -1.0, 0.0, 0.25, 0.0, -1.0, 0.0);
	vList[24].setvals(-1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[25].setvals(-1.0, -2.0, 0.0, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);
	vList[26].setvals(1.0,  -2.0, 0.0, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0);
	vList[27].setvals(1.0,  -1.0, 0.0, -1.0, 0.0, 0.0, 0.25, -1.0, 0.0, 0.0);

	pList[0].Create(6);
	pList[0][0] = 1;
	pList[0][1] = 10;
	pList[0][2] = 13;
	pList[0][3] = 4;
	pList[0][4] = 8;
	pList[0][5] = 7;


	pList[1].Create(6);
	pList[1][0] = 15;
	pList[1][1] = 21;
	pList[1][2] = 22;
	pList[1][3] = 18;
	pList[1][4] = 27;
	pList[1][5] = 24;

	pList[2].Create(4);
	pList[2][0] = 0;
	pList[2][1] = 14;
	pList[2][2] = 15;
	pList[2][3] = 1;

	pList[3].Create(4);
	pList[3][0] = 1;
	pList[3][1] = 15;
	pList[3][2] = 16;
	pList[3][3] = 2;

	pList[4].Create(4);
	pList[4][0] = 2;
	pList[4][1] = 16;
	pList[4][2] = 17;
	pList[4][3] = 3;

	pList[5].Create(4);
	pList[5][0] = 3;
	pList[5][1] = 17;
	pList[5][2] = 18;
	pList[5][3] = 4;

	pList[6].Create(4);
	pList[6][0] = 4;
	pList[6][1] = 18;
	pList[6][2] = 19;
	pList[6][3] = 5;

	pList[7].Create(4);
	pList[7][0] = 5;
	pList[7][1] = 19;
	pList[7][2] = 14;
	pList[7][3] = 0;

	pList[8].Create(4);
	pList[8][0] = 6;
	pList[8][1] = 7;
	pList[8][2] = 21;
	pList[8][3] = 20;

	pList[9].Create(4);
	pList[9][0] = 7;
	pList[9][1] = 8;
	pList[9][2] = 22;
	pList[9][3] = 21;

	pList[10].Create(4);
	pList[10][0] = 8;
	pList[10][1] = 9;
	pList[10][2] = 23;
	pList[10][3] = 22;

	pList[11].Create(4);
	pList[11][0] = 9;
	pList[11][1] = 6;
	pList[11][2] = 20;
	pList[11][3] = 23;

	pList[12].Create(4);
	pList[12][0] = 10;
	pList[12][1] = 11;
	pList[12][2] = 25;
	pList[12][3] = 24;

	pList[13].Create(4);
	pList[13][0] = 11;
	pList[13][1] = 12;
	pList[13][2] = 26;
	pList[13][3] = 25;

	pList[14].Create(4);
	pList[14][0] = 12;
	pList[14][1] = 13;
	pList[14][2] = 27;
	pList[14][3] = 26;

	pList[15].Create(4);
	pList[15][0] = 13;
	pList[15][1] = 10;
	pList[15][2] = 24;
	pList[15][3] = 27;

	pList[16].Create(4);
	pList[16][0] = 0;
	pList[16][1] = 6;
	pList[16][2] = 9;
	pList[16][3] = 5;

	pList[17].Create(4);
	pList[17][0] = 1;
	pList[17][1] = 7;
	pList[17][2] = 6;
	pList[17][3] = 0;

	pList[18].Create(4);
	pList[18][0] = 2;
	pList[18][1] = 11;
	pList[18][2] = 10;
	pList[18][3] = 1;

	pList[19].Create(4);
	pList[19][0] = 2;
	pList[19][1] = 3;
	pList[19][2] = 12;
	pList[19][3] = 11;

	pList[20].Create(4);
	pList[20][0] = 3;
	pList[20][1] = 4;
	pList[20][2] = 13;
	pList[20][3] = 12;

	pList[21].Create(4);
	pList[21][0] = 4;
	pList[21][1] = 5;
	pList[21][2] = 9;
	pList[21][3] = 8;

	pList[22].Create(4);
	pList[22][0] = 14;
	pList[22][1] = 19;
	pList[22][2] = 23;
	pList[22][3] = 20;

	pList[23].Create(4);
	pList[23][0] = 15;
	pList[23][1] = 14;
	pList[23][2] = 20;
	pList[23][3] = 21;

	pList[24].Create(4);
	pList[24][0] = 16;
	pList[24][1] = 15;
	pList[24][2] = 24;
	pList[24][3] = 25;

	pList[25].Create(4);
	pList[25][0] = 17;
	pList[25][1] = 16;
	pList[25][2] = 25;
	pList[25][3] = 26;

	pList[26].Create(4);
	pList[26][0] = 18;
	pList[26][1] = 17;
	pList[26][2] = 26;
	pList[26][3] = 27;

	pList[27].Create(4);
	pList[27][0] = 19;
	pList[27][1] = 18;
	pList[27][2] = 22;
	pList[27][3] = 23;

	init(28, 28, vList, pList);
}

double CCSubdivsurf::MyRand(double low, double high)
{
	double range=high-low;
	double each=1.0*RAND_MAX/range;
	return (rand()/each+low);
}

/********************************************************************************/
/* initfigure8() -- create a Bpoint list and an array of PindexLists with the	*/
/*		the vertex and polygon information, then call init() to create the		*/
/*		fathermesh.																*/
/********************************************************************************/
void CCSubdivsurf::initcustomized()
{
	//srand((unsigned)time(NULL));
	Bpoint vList[36];
	PindexList pList[34];
	int vNum=36, pNum=34;

	vList[0].setvals(6,  6.0, -1.1, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[1].setvals(-6.0,   6, -1.2, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);
	vList[2].setvals(-6.0,  -6.0, -1.2, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0);
	vList[3].setvals(6.0,   -6.0, -1.2, -1.0, 0.0, 0.0, 0.25, -1.0, 0.0, 0.0);

	vList[32].setvals(4,   4.0, -0.08, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[33].setvals(-4.0,   4, -0.08, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);
	vList[34].setvals(-4.0,  -4.0, -0.08, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0);
	vList[35].setvals(4.0,   -4.0, -0.08, -1.0, 0.0, 0.0, 0.25, -1.0, 0.0, 0.0);

	vList[4].setvals(3.0,    3.0, -0.05, 0.0, -1.0, 0.0, 0.25, 0.0, -1.0, 0.0);
	vList[5].setvals(-3.0,    3.0, -0.05, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[6].setvals(-3.0,   -3.0, -0.05, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);
	vList[7].setvals(3.0,   -3.0, -0.05, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0);

	vList[24].setvals(2.0,    2.0, -0.05, 0.0, -1.0, 0.0, 0.25, 0.0, -1.0, 0.0);
	vList[25].setvals(-2.0,    2.0, -0.05, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[26].setvals(-2.0,   -2.0, -0.05, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);
	vList[27].setvals(2.0,   -2.0, -0.05, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0);

	vList[8].setvals(1.5,    1.5, 0, -1.0, 0.0, 0.0, 0.25, -1.0, 0.0, 0.0);
	vList[9].setvals(-1.5,    1.5, 0, 0.0, -1.0, 0.0, 0.25, 0.0, -1.0, 0.0);
	vList[10].setvals(-1.5, -1.5, 0, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[11].setvals(1.5, -1.5, 0, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);

	vList[20].setvals(1.15,    1.15, 0.5, -1.0, 0.0, 0.0, 0.25, -1.0, 0.0, 0.0);
	vList[21].setvals(-1.15,    1.15, 0.5, 0.0, -1.0, 0.0, 0.25, 0.0, -1.0, 0.0);
	vList[22].setvals(-1.15, -1.15, 0.5, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[23].setvals(1.15, -1.15, 0.5, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);

	vList[12].setvals(1.0,  1.0, 1.0, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0);
	vList[13].setvals(-1.0,  1.0, 1.0, -1.0, 0.0, 0.0, 0.25, -1.0, 0.0, 0.0);
	vList[14].setvals(-1.0,  -1.0, 1.0, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[15].setvals(1.0,  -1.0, 1.0, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);

	vList[28].setvals(1.0,  1.0, 1.5, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0);
	vList[29].setvals(-1.0,  1.0, 1.5, -1.0, 0.0, 0.0, 0.25, -1.0, 0.0, 0.0);
	vList[30].setvals(-1.0,  -1.0, 1.5, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);
	vList[31].setvals(1.0,  -1.0, 1.5, 1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.0);

	vList[16].setvals(1.0, 1.0, 2.0, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0);
	vList[17].setvals(-1.0,  1.0, 2.0, -1.0, 0.0, 0.0, 0.25, -1.0, 0.0, 0.0);
	vList[18].setvals(-1.0,   -1.0, 2.0, 0.0, -1.0, 0.0, 0.25, 0.0, -1.0, 0.0);
	vList[19].setvals(1.0,   -1.0, 2.0, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 1.0);

	int i;
	double xx, yy, zz;
	for(i=0; i<20; i++)
	{
		xx=vList[i].pnt[0];
		if(xx<0) xx=-xx;
		xx=xx*0.05;
		yy=vList[i].pnt[1];
		if(yy<0) yy=-yy;
		yy=yy*0.05;
		zz=vList[i].pnt[2];
		if(zz<0) zz=-zz;
		zz=zz*0.05;
		vList[i].pnt[0]+=MyRand(-xx, xx);
		vList[i].pnt[1]+=MyRand(-yy, yy);
		vList[i].pnt[2]+=MyRand(-zz, zz);
	}

	for(i=0; i<pNum; i++)
		pList[i].Create(4);

	pList[0][0] = 0;
	pList[0][1] = 3;
	pList[0][2] = 2;
	pList[0][3] = 1;

/////////////////////////////////
	pList[1][0] = 0;
	pList[1][1] = 32;
	pList[1][2] = 35;
	pList[1][3] = 3;

	pList[2][0] = 1;
	pList[2][1] = 33;
	pList[2][2] = 32;
	pList[2][3] = 0;

	pList[3][0] = 2;
	pList[3][1] = 34;
	pList[3][2] = 33;
	pList[3][3] = 1;

	pList[4][0] = 3;
	pList[4][1] = 35;
	pList[4][2] = 34;
	pList[4][3] = 2;
/////////////////////////////////
	pList[30][0] = 32;
	pList[30][1] = 4;
	pList[30][2] = 7;
	pList[30][3] = 35;

	pList[31][0] = 33;
	pList[31][1] = 5;
	pList[31][2] = 4;
	pList[31][3] = 32;

	pList[32][0] = 34;
	pList[32][1] = 6;
	pList[32][2] = 5;
	pList[32][3] = 33;

	pList[33][0] = 35;
	pList[33][1] = 7;
	pList[33][2] = 6;
	pList[33][3] = 34;
/////////////////////////////////
	pList[5][0] = 4;
	pList[5][1] = 24;
	pList[5][2] = 27;
	pList[5][3] = 7;

	pList[6][0] = 5;
	pList[6][1] = 25;
	pList[6][2] = 24;
	pList[6][3] = 4;

	pList[7][0] = 6;
	pList[7][1] = 26;
	pList[7][2] = 25;
	pList[7][3] = 5;

	pList[8][0] = 7;
	pList[8][1] = 27;
	pList[8][2] = 26;
	pList[8][3] = 6;
/////////////////////////////////
	pList[22][0] = 24;
	pList[22][1] = 8;
	pList[22][2] = 11;
	pList[22][3] = 27;

	pList[23][0] = 25;
	pList[23][1] = 9;
	pList[23][2] = 8;
	pList[23][3] = 24;

	pList[24][0] = 26;
	pList[24][1] = 10;
	pList[24][2] = 9;
	pList[24][3] = 25;

	pList[25][0] = 27;
	pList[25][1] = 11;
	pList[25][2] = 10;
	pList[25][3] = 26;
/////////////////////////////////
	pList[9][0] = 8;
	pList[9][1] = 20;
	pList[9][2] = 23;
	pList[9][3] = 11;

	pList[10][0] = 9;
	pList[10][1] = 21;
	pList[10][2] = 20;
	pList[10][3] = 8;

	pList[11][0] = 10;
	pList[11][1] = 22;
	pList[11][2] = 21;
	pList[11][3] = 9;

	pList[12][0] = 11;
	pList[12][1] = 23;
	pList[12][2] = 22;
	pList[12][3] = 10;
/////////////////////////////////
	pList[18][0] = 20;
	pList[18][1] = 12;
	pList[18][2] = 15;
	pList[18][3] = 23;

	pList[19][0] = 21;
	pList[19][1] = 13;
	pList[19][2] = 12;
	pList[19][3] = 20;

	pList[20][0] = 22;
	pList[20][1] = 14;
	pList[20][2] = 13;
	pList[20][3] = 21;

	pList[21][0] = 23;
	pList[21][1] = 15;
	pList[21][2] = 14;
	pList[21][3] = 22;
//////////////////////////////////
	pList[13][0] = 12;
	pList[13][1] = 28;
	pList[13][2] = 31;
	pList[13][3] = 15;

	pList[14][0] = 13;
	pList[14][1] = 29;
	pList[14][2] = 28;
	pList[14][3] = 12;

	pList[15][0] = 14;
	pList[15][1] = 30;
	pList[15][2] = 29;
	pList[15][3] = 13;

	//pList[16].Create(4);
	pList[16][0] = 15;
	pList[16][1] = 31;
	pList[16][2] = 30;
	pList[16][3] = 14;
//////////////////////////////////
	pList[26][0] = 28;
	pList[26][1] = 16;
	pList[26][2] = 19;
	pList[26][3] = 31;

	pList[27][0] = 29;
	pList[27][1] = 17;
	pList[27][2] = 16;
	pList[27][3] = 28;

	pList[28][0] = 30;
	pList[28][1] = 18;
	pList[28][2] = 17;
	pList[28][3] = 29;

	//pList[16].Create(4);
	pList[29][0] = 31;
	pList[29][1] = 19;
	pList[29][2] = 18;
	pList[29][3] = 30;
/////////////////////////////////
	pList[17][0] = 16;
	pList[17][1] = 17;
	pList[17][2] = 18;
	pList[17][3] = 19;

	init(vNum, pNum, vList, pList);
}

/********************************************************************************/
/* initfigure8() -- create a Bpoint list and an array of PindexLists with the	*/
/*		the vertex and polygon information, then call init() to create the		*/
/*		fathermesh.																*/
/********************************************************************************/
void CCSubdivsurf::initfromfile()
{
	FILE *fp;
	char filename[20];
	int vNum, fNum;
	double vec[10];

	Bpoint *vList;
	PindexList *pList;

	sprintf(filename, "%s", "test.dat");
	if((fp=fopen(filename, "rb"))==NULL)
		return;

	fscanf(fp, "%d\n", &vNum);
	vList=new Bpoint[vNum];
	for(int vdex=0; vdex<vNum; vdex++)
	{
		fscanf(fp, "%lf, %lf, %lf\n", &vec[0], &vec[1], &vec[2]);
		fscanf(fp, "%lf, %lf, %lf\n", &vec[3], &vec[4], &vec[5]);
		fscanf(fp, "%lf, %lf, %lf, %lf\n", &vec[6], &vec[7], &vec[8], &vec[9]);
		vList[vdex].setvals(vec[0], vec[1], vec[2], vec[3], vec[4], vec[5], vec[6], vec[7], vec[8], vec[9]);

		/*
		CCVertex *thisVert=new CCVertex();
		fscanf(fp, "%f, %f, %f\n", &vec[0], &vec[1], &vec[2]);
		thisVert->loc.set(vec);		//, thisVert->loc.y(), thisVert->loc.z());
		fscanf(fp, "%f, %f, %f\n", &vec[0], &vec[1], &vec[2]);
		thisVert->norm.set(vec);
		//fscanf(fp, "%f, %f, %f\n", thisVert->norm.x(), thisVert->norm.y(), thisVert->norm.z());
		fscanf(fp, "%f, %f, %f, %f\n", &vec[0], &vec[1], &vec[2], &vec[3]);

		//thisVert->rad, thisVert->uvt.u, thisVert->uvt.v, thisVert->uvt.t);
		*/
	}

	int sNum;
	fscanf(fp, "%d\n", &fNum);
	pList=new PindexList[fNum];
	for(int fdex=0; fdex<fNum; fdex++)
	{
		//CCPolygon *thisPoly=thisMesh->FaceList(fdex);
		fscanf(fp, "%d", &sNum);
		pList[fdex].Create(sNum);
		for(int fvdex=0; fvdex<sNum; fvdex++)
			fscanf(fp, ", %d", &pList[fdex][fvdex]);
		fscanf(fp, "\n");
	}

#ifdef DEBUG_QIONG_BLENDING
	fscanf(fp, "%d", &dPNum);
	if(dPNum>0)
	{
		for(int fdex=0; fdex<dPNum; fdex++)
		{
			fscanf(fp, ", %d", &polyList[fdex]);
		}
	}
#endif

	fclose(fp);

	init(vNum, fNum, vList, pList);
	//return true;
}

/********************************************************************************/
/* returns pointer to  list of vertices, each stored as a Bpoint.				*/
/*           where npoints gives # of bpoints on boundary.						*/
/*   NOTE:  Bptlistptr is created here, but control (and deletion of) passed to */
/*    calling routine															*/
/* NOTE:  to use, first call GLrendersubdiv(subdivlevel, NO_RENDER)				*/
/********************************************************************************/
void CCSubdivsurf::SubdivBoundaryInfo(int subLvl, int *npoints, Bpoint **bptlistptr)
{
	CCMesh *thisMesh=submeshes(subLvl);
	int numPoints=thisMesh->numverts();

	*npoints = numPoints;

	// deal with NULL pointers here instead of in the Ps2
	if(*bptlistptr==NULL)
		*bptlistptr = new Bpoint[numPoints];

	int vnum = 0;
	for(vnum=0; vnum<numPoints; vnum++)
	{
		CCVertex *thisVert=thisMesh->VertList(vnum);
		(*bptlistptr)[vnum].pnt[0] = thisVert->limitloc.x();
		(*bptlistptr)[vnum].pnt[1] = thisVert->limitloc.y();
		(*bptlistptr)[vnum].pnt[2] = thisVert->limitloc.z();
		(*bptlistptr)[vnum].norm[0] = thisVert->limitnorm.x();
		(*bptlistptr)[vnum].norm[1] = thisVert->limitnorm.y();
		(*bptlistptr)[vnum].norm[2] = thisVert->limitnorm.z();
		(*bptlistptr)[vnum].rad = thisVert->limitRad;
		(*bptlistptr)[vnum].u = thisVert->uvt.u;
		(*bptlistptr)[vnum].v = thisVert->uvt.v;
		(*bptlistptr)[vnum].t = thisVert->uvt.t;
	}
}

/********************************************************************************/
/* create neighborlist associated with vertices produced by mod_subdiv()
/*    for current activesubmesh.  Contains a list of nverts Tileneighbors arrays,
/*    each giving index values associated with neighbors of indexing vertexnumber.
/*    Tileneighbors (int[8]) typedefed above.
/*    If fewer than 8 neighbors, last elements in list will be == -1.
/* Vertices are in counterclockwise order starting with an edge vertex
/*    (See def of class Tileneighbor above)
/* Returns false if no vertex-list exists at that subdivlevel.
/********************************************************************************/
bool CCSubdivsurf::GetVertexNeighbors(int subLvl, int *nverts, Tileneighbors **neighborlist)
{
	//int numrows, numcols, numsiderows, numsidecols, lengthedge;
	int vcount, i;

	CCMesh *thisMesh=submeshes(subLvl);

		// number of points is top & bottom grids + sidegrid minus top and bottom edge,
		//   which were already listed in top & bottom grids
		int numpoints = thisMesh->numverts(); //2*numrows*numcols + (numsiderows*numsidecols - 2*lengthedge);

		*nverts = numpoints;

		// deal with NULL pointers here instead of in the Ps2
		if(*neighborlist==NULL)
			*neighborlist = new Tileneighbors[numpoints];

		CCVertex *thisVert;
		for(vcount = 0; vcount<numpoints; vcount++)
		{
			thisVert=thisMesh->VertList(vcount);
			for(i=0; i<thisVert->v_valence; i++)
			{
				// Find out the adjacent neighbor vertices one by one
				//   and rotote them into along with the same angle
				CCPindex thisPolyIndex = thisMesh->VertFaceNeighbors[vcount][i];
				CCPolygon *thisPoly=thisMesh->FaceList(thisPolyIndex);
				int polyVIndex=thisMesh->VertFaceNeighborsVert[vcount][i];
				int pnumsides=thisPoly->numsides();

				// The adjacent vertex is the next vertex counterclockwise from vPolyIndex
				(*neighborlist)[vcount][i]=thisPoly->Pvertices[INTMOD(polyVIndex+1, pnumsides)];
			}
			for(i=thisVert->v_valence; i<MAX_TILE_NEIGHBORS; i++)
				(*neighborlist)[vcount][i]=-1;

		}

	return true;
}

// return the new neighboring vertex list for the area calculation in Mask
bool CCSubdivsurf::GetNeighboringVertices(int subLvl, NeighborVertices **nVertices)
{
	int vcount, i;

	CCMesh *thisMesh=submeshes(subLvl);

	// number of points
	int numpoints = thisMesh->numverts(); //2*numrows*numcols + (numsiderows*numsidecols - 2*lengthedge);

	// deal with NULL pointers here instead of in the Ps2
	if(*nVertices==NULL)
		*nVertices = new NeighborVertices[numpoints];

	CCVertex *thisVert;
	CCVertex *thatVert;
	CCVertex avgVert;
	// there are 2*v_valence vertices for the area calculations for each vertex
	for(vcount = 0; vcount<numpoints; vcount++)
	{
		thisVert=thisMesh->VertList(vcount);
		(*nVertices)[vcount].nNumber=thisVert->v_valence;
		for(i=0; i<thisVert->v_valence; i++)
		{
			// Find out the adjacent neighbor vertices one by one
			//   and rotote them into along with the same angle
			CCPindex thisPolyIndex = thisMesh->VertFaceNeighbors[vcount][i];
			avgVert.limitloc.set(0, 0, 0);
			CCPolygon *thisPoly=thisMesh->FaceList(thisPolyIndex);
			int pnumsides=thisPoly->numsides();

			// average of the #pnumsides vertices of a neighboring face
			for(int j=0; j<pnumsides; j++) {
				avgVert.limitloc+=thisMesh->VertList(thisPoly->Pvertices[j])->limitloc;
			}
			avgVert.limitloc /= pnumsides;
			(*nVertices)[vcount][2*i+1].pnt[0]=avgVert.limitloc.X();
			(*nVertices)[vcount][2*i+1].pnt[1]=avgVert.limitloc.Y();
			(*nVertices)[vcount][2*i+1].pnt[2]=avgVert.limitloc.Z();

			// the neighboring vertex
			int polyVIndex=thisMesh->VertFaceNeighborsVert[vcount][i];
			// The adjacent vertex is the next vertex counterclockwise from vPolyIndex
			thatVert=thisMesh->VertList(thisPoly->Pvertices[INTMOD(polyVIndex+1, pnumsides)]);
			// average of the current and the neiboring vertices
			avgVert.limitloc=(thisVert->limitloc+thatVert->limitloc)/2;
			(*nVertices)[vcount][2*i].pnt[0]=avgVert.limitloc.X();
			(*nVertices)[vcount][2*i].pnt[1]=avgVert.limitloc.Y();
			(*nVertices)[vcount][2*i].pnt[2]=avgVert.limitloc.Z();
		}
	}

	return true;
}

// New function to return the neighboring-polygon-indices of each tile
void CCSubdivsurf::GetTileNeighbors(int subLvl, int *ntiles, Tileneighbors **neighborlist)
{
	int tcount, i;

	CCMesh *thisMesh=submeshes(subLvl);
	int numTiles = thisMesh->numfaces(); //2*numrows*numcols + (numsiderows*numsidecols - 2*lengthedge);
	*ntiles = numTiles;

	// deal with NULL pointers here instead of in the Ps2
	if(*neighborlist==NULL)
		*neighborlist = new Tileneighbors[numTiles];

	CCPolygon *thisPoly;
	for(tcount = 0; tcount<numTiles; tcount++)
	{
		thisPoly=thisMesh->FaceList(tcount);
		for(i = 0; i < thisPoly->numsides(); i++)
		{
			(*neighborlist)[tcount][i]=thisPoly->Pneighbors[i];
		}
		for(i = thisPoly->numsides(); i < MAX_TILE_NEIGHBORS; i++)
			(*neighborlist)[tcount][i]=-1;
	}
}

// New function to return the vertex-indices of each tile
bool CCSubdivsurf::GetTileIndices(int subLevel, int **subdivtileindices)
{
	CCMesh *thisMesh=submeshes(subLevel);
	int numTiles = thisMesh->numfaces(); //2*numrows*numcols + (numsiderows*numsidecols - 2*lengthedge);

	// deal with NULL pointers here instead of in the Ps2
	if(*subdivtileindices==NULL)
		*subdivtileindices = new int[4*numTiles];

	CCPolygon *thisPoly;
	for(int tcount = 0; tcount<numTiles; tcount++)
	{
		thisPoly=thisMesh->FaceList(tcount);

		// The assumption is, each tile in the higher level is a quad due to the subdivision
		if(thisPoly->numsides()!=4)
		{
			delete []subdivtileindices;
			subdivtileindices=NULL;
			return false;
		}

		for(int i=0; i<4; i++) //for(int i=0; i<thisPoly->numsides(); i++)
			(*subdivtileindices)[4*tcount+i]=thisPoly->Pvertices[i];
	}
	return true;
}

/********************************************************************************/
/* returns pointer to list of quad-tiles, each stored as 4 sequential Bpoints.  */
/*     where numtiles gives # of such 4-Bpoint tiles							*/
/*   NOTE:  Bptlistptr is created here, but control (and deletion of) passed to */
/*    calling routine															*/
/* NOTE:  to use, first call GLrendersubdiv(subdivlevel, NO_RENDER)				*/
/********************************************************************************/
void CCSubdivsurf::SubdivTileInfo(int subdivlevel, int *ntiles, Bpoint **bptlistptr)
//								 int *num_ucoltiles, int *num_vrowtiles, int *num_t_tiles)
{
	CCMesh *thisMesh=submeshes(subdivlevel);

	int numtiles = thisMesh->numfaces();

	*ntiles=numtiles;

	// deal with NULL pointers here instead of in the Ps2
	if(*bptlistptr==NULL)
		*bptlistptr=new Bpoint[numtiles*4];

	int vNum=0, pDex, vDex;
	CCPolygon *thisPoly;
	CCVertex *thisVert;
	for(pDex=0; pDex<numtiles; pDex++) 
	{
		thisPoly=thisMesh->FaceList(pDex);
		for(int i=0; i<4; i++)
		{
			vDex=thisPoly->Pvertices[i];
			thisVert=thisMesh->VertList(vDex);
			(*bptlistptr)[vNum].pnt[0]=thisVert->limitloc.x();
			(*bptlistptr)[vNum].pnt[1]=thisVert->limitloc.y();
			(*bptlistptr)[vNum].pnt[2]=thisVert->limitloc.z();
			(*bptlistptr)[vNum].norm[0]=thisVert->limitnorm.x();
			(*bptlistptr)[vNum].norm[1]=thisVert->limitnorm.y();
			(*bptlistptr)[vNum].norm[2]=thisVert->limitnorm.z();
			(*bptlistptr)[vNum].rad=thisVert->limitRad;
			(*bptlistptr)[vNum].u=thisVert->uvt.u;
			(*bptlistptr)[vNum].v=thisVert->uvt.v;
			(*bptlistptr)[vNum].t=thisVert->uvt.t;
			vNum++;
		}
	}

}

void CCSubdivsurf::loadlevelzeroBlending()
{
	CCMesh *thisMesh = submeshes(0);

	int numfaces1, numverts1, nedges1;
	numfaces1 = fathermesh.numfaces();
	numverts1 = fathermesh.numverts();
	nedges1 = fathermesh.numedges();

	thisMesh->VertList.Create(numverts1);
	thisMesh->FaceList.Create(numfaces1);
	thisMesh->nedges = nedges1;

	CCVindex vdex;
	for (vdex = 0; vdex < numverts1; vdex++)
		thisMesh->VertList[vdex] = fathermesh.VertList[vdex];

	CCPindex fdex;
	for (fdex = 0; fdex < numfaces1; fdex++) {

		CCPolygon *fathersface = fathermesh.FaceList(fdex);

		// Copy the Polygons
		thisMesh->FaceList[fdex] = fathermesh.FaceList[fdex];
	}

	// Now that VertexList and FaceList are done, do VertFaceNeighbors and FaceFaceNeighbors
	//   though we may not need these both --- For now, skip FaceFaceNeighbors
	thisMesh->VertFaceNeighbors.Create(numverts1);
	thisMesh->VertFaceNeighborsVert.Create(numverts1);

	for (vdex = 0; vdex < numverts1; vdex++) {

		int numvneighbors = fathermesh.VertFaceNeighbors[vdex].Length();
		thisMesh->VertFaceNeighbors[vdex].Create(numvneighbors);
		thisMesh->VertFaceNeighborsVert[vdex].Create(numvneighbors);

		for (fdex = 0; fdex < numvneighbors; fdex++) {
			thisMesh->VertFaceNeighbors[vdex][fdex] = fathermesh.VertFaceNeighbors[vdex][fdex];
			thisMesh->VertFaceNeighborsVert[vdex][fdex] = fathermesh.VertFaceNeighborsVert[vdex][fdex];
		}
	}

#if !(defined(DEBUG_QIONG_TESTING_INVERSE_ALGORITHM_ON_QUAD_MESH) || defined(QIONG_DEBUG_BLENDING_02192003))
	computeIISmesh(0, 12);
	computelimitmesh(0);
#else
	for (vdex = 0; vdex < numverts1; vdex++)
	{
		CCVertex *thisVert;
		thisVert=thisMesh->VertList(vdex);
		thisVert->limitloc=thisVert->loc;
		thisVert->limitnorm=thisVert->norm;
	}
#endif
}


/********************************************************************************/
/* loadlevelzero() -- load values of fathermesh into submeshes[0]				*/
/********************************************************************************/
void CCSubdivsurf::loadlevelzeroPureSmoothing()
{
	CCMesh *thisMesh = submeshes(0);

	int numfaces1, numverts1, nedges1;
	numfaces1 = fathermesh.numfaces();
	numverts1 = fathermesh.numverts();
	nedges1 = fathermesh.numedges();

	thisMesh->VertList.Create(numverts1);
	thisMesh->FaceList.Create(numfaces1);
	thisMesh->nedges = nedges1;

	CCVindex vdex;
	for (vdex = 0; vdex < numverts1; vdex++)
		thisMesh->VertList[vdex] = fathermesh.VertList[vdex];

	CCPindex fdex;
	for (fdex = 0; fdex < numfaces1; fdex++) {

		CCPolygon *fathersface = fathermesh.FaceList(fdex);

		// Copy the Polygons
		thisMesh->FaceList[fdex] = fathermesh.FaceList[fdex];
	}

	// Now that VertexList and FaceList are done, do VertFaceNeighbors and FaceFaceNeighbors
	//   though we may not need these both --- For now, skip FaceFaceNeighbors
	thisMesh->VertFaceNeighbors.Create(numverts1);
	thisMesh->VertFaceNeighborsVert.Create(numverts1);

	for (vdex = 0; vdex < numverts1; vdex++) {

		int numvneighbors = fathermesh.VertFaceNeighbors[vdex].Length();
		thisMesh->VertFaceNeighbors[vdex].Create(numvneighbors);
		thisMesh->VertFaceNeighborsVert[vdex].Create(numvneighbors);

		for (fdex = 0; fdex < numvneighbors; fdex++) {
			thisMesh->VertFaceNeighbors[vdex][fdex] = fathermesh.VertFaceNeighbors[vdex][fdex];
			thisMesh->VertFaceNeighborsVert[vdex][fdex] = fathermesh.VertFaceNeighborsVert[vdex][fdex];
		}
	}
}

/********************************************************************************/
/* SplitandAverage() -- split nivelo to give the next finer level (nivelo + 1)	*/
/********************************************************************************/
bool CCSubdivsurf::splitandaverageBlending(int nivelo, int *vTag, int *fTag)
{
	int i;
#if (defined(DEBUG_QIONG_TESTING_INVERSE_ALGORITHM_ON_QUAD_MESH) || defined(QIONG_DEBUG_BLENDING_02192003))
if (nivelo==1)
{
	computeIISmeshBlending(1, 12);
}
#endif
if(nivelo==0)
{
	if (nivelo + 1 >= NUMLEVELS) 
	{
		fprintf(stderr, "ERROR: CCSubdivsurf::splitandaverage() -- Number of levels too high for split.\n");
		return false;
	}

	CCMesh *oldlevel = submeshes(nivelo);
	CCMesh *newlevel = submeshes(nivelo + 1);

	// The new mesh has the even vertices + the new face vertices + the new edge vertices
	int newnumverts = oldlevel->numverts() + oldlevel->numfaces() + oldlevel->numedges();

	int interFaceCounter=0;
	for(i=0; i<oldlevel->numfaces(); i++)
		if(fTag[i]==2) 
			interFaceCounter++;

	// Each new edge-vertex generates 2 new faces
	//   (As proof, each new edge from the edge-vertex to the 2 face-vertices forms the
	//    lead edge (in counterclockwise order) of one of the new faces)
	int newnumfaces = 2*oldlevel->numedges()-2*interFaceCounter;
	//int newnumfaces = 2*oldlevel->numedges();

	// work to do
	newlevel->VertList.Create(newnumverts);
	newlevel->FaceList.Create(newnumfaces);
	newlevel->VertFaceNeighbors.Create(newnumverts);
	newlevel->VertFaceNeighborsVert.Create(newnumverts);

	int numoldfaces = oldlevel->numfaces();
	int numoldverts = oldlevel->numverts();

#if defined(QIONG_DEBUG_BLENDING_02192003)
	//int vRedex[1000];			// a re-index of all the vertices in the new mesh
	//int vFlag[1000];
	vFaceReduction=0;		// the # of face vertices in the new mesh that should be removed
	vEdgeReduction=0;		// the # of edge vertices in the new mesh that should be removed
#endif

	// Some constantly used variables are defined here!
	CCPolygon *thispoly;
	CCVertex *thisVertex;
	int polydex;

	// Copy even vertices into new vertex list
	// NO, begin by zeroing out the even vertices...use parent's even vertices
	//    (which have the same indices) for computation, and use the new ones
	//    as accumulators for computing the new even vertices.
	// Formula is neweven = oldweight*v_old + neighborweight*(sum[4*newedgeneighbors - newfaceneighbors])
	//    Derivation of this is in Thall diss, App.A.
	int evenvindex;
	for (evenvindex = 0; evenvindex < numoldverts; evenvindex++)
	{
		newlevel->VertList[evenvindex] = oldlevel->VertList[evenvindex];
		newlevel->VertList[evenvindex].is_oddvert = false;

		newlevel->VertList[evenvindex].loc = NULLVEC;
		newlevel->VertList[evenvindex].norm = NULLVEC;

	// so all the even vertices should be KEPT INTACT
#if defined(QIONG_DEBUG_BLENDING_02192003)
		vRedex[evenvindex]=evenvindex;
		// changed changed changed
		// 0 means main figure and 1 means subfigure!
		vFlag[evenvindex]=vTag[evenvindex]>=2?0:1;
#endif
	}

#ifdef QIONG_DEBUG_20030425
	for(evenvindex=0; evenvindex<numoldverts; evenvindex++)
	{
		thisVertex=newlevel->VertList(evenvindex);
	}
#endif

	DbVector3 newloc;
	DbVector3 newnorm;
	Medcoord newuvt;
	double newradius;

#if defined(QIONG_DEBUG_BLENDING_02192003)
#endif

#ifdef QIONG_DEBUG_20030425
	bool badPoly;
	float whereIsThePoly;					// 1 on top, -1 on bottom, otherwise on crest!
#endif

	// Create new face vertices at centroids of their polygons
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		thispoly = oldlevel->FaceList(polydex);

		newloc.set(0.0, 0.0, 0.0);
		newnorm.set(0.0, 0.0, 0.0);
		newuvt.setvals(0.0, 0.0, 0.0);
		newradius = 0.0;

#ifdef QIONG_DEBUG_20030425
	CCVertex *thisOldVertex = oldlevel->VertList(thispoly->Pvertices[0]);
	whereIsThePoly=thisOldVertex->uvt.t;
	badPoly=false;
#endif

		int v_index;
		for (v_index = 0; v_index < thispoly->numsides(); v_index++) 
		{

			//CCVertex *thisVertex = oldlevel->VertList(thispoly->Pvertices[v_index]);
			thisVertex = oldlevel->VertList(thispoly->Pvertices[v_index]);
#ifdef QIONG_DEBUG_20030425
			if(thisVertex->uvt.t!=whereIsThePoly)
			{
				if(vTag[thispoly->Pvertices[v_index]]>=2 && vTag[thispoly->Pvertices[0]]>=2 &&
					thisVertex->uvt.u!=thisOldVertex->uvt.u && thisVertex->uvt.v!=thisOldVertex->uvt.v)
					badPoly=true;
			}
			if(badPoly)
				printf("thisVert->uvt.t=%5.4f\t", thisVertex->uvt.t);
#endif

			// Sum the node values
			newloc += thisVertex->loc;
			newnorm += thisVertex->norm;
			newradius += thisVertex->rad;
			newuvt += thisVertex->uvt;
		}
#ifdef QIONG_DEBUG_20030425
		if(badPoly)
			printf("\n");
#endif

		newloc /= thispoly->numsides();
		newnorm /= thispoly->numsides();
		newradius /= thispoly->numsides();
		newuvt /= thispoly->numsides();

		newlevel->VertList[numoldverts + polydex].setvals(newloc, newnorm, newradius,
			                                              newuvt, thispoly->numsides(), true);

		// all the face vertex corresponding to the INTER face are tagged
#if defined(QIONG_DEBUG_BLENDING_02192003)
		if(fTag[polydex]==2)
		{
			vRedex[numoldverts + polydex]=-1;
			vFlag[numoldverts + polydex]=-1;
			vFaceReduction++;
			//continue;
		}
		else
		{
			vRedex[numoldverts + polydex]=numoldverts+polydex-vFaceReduction;
			vFlag[numoldverts + polydex]=fTag[polydex]>=3?0:1;
		}
		thispoly->newfacevert = numoldverts + polydex;
		thispoly->newedgeverts.Create(thispoly->numsides()); // go ahead and initialize newedgevert list
#else
		thispoly->newfacevert = numoldverts + polydex;
		thispoly->newedgeverts.Create(thispoly->numsides()); // go ahead and initialize newedgevert list
#endif

		// Now, subtract values of the new face vertex into loc and norm of surrounding even vertices
		for (v_index = 0; v_index < thispoly->numsides(); v_index++) 
		{
			// CCVertex *thisVertex = newlevel->VertList(thispoly->Pvertices[v_index]);
			thisVertex = newlevel->VertList(thispoly->Pvertices[v_index]);
			thisVertex->loc -= newloc;
			thisVertex->norm -= newnorm;
			// NOTE:  if you make this an addition instead of a subtraction,
			//   it screws everything.  :)  Makes some nice sea anemonies, though
		}

	}
#ifdef QIONG_DEBUG_20030425
	for(polydex=numoldverts; polydex<numoldverts+numoldfaces; polydex++)
	{
		thisVertex=newlevel->VertList(polydex);
	}
#endif

	// need to re-index all the edge vertices in the middle!

	// For each edge of polygon, create new edge vertex at averaged position of
	//    of old edge endpoints and new face vertices at centroids
	// For r and uvt interpolation, just average old endpoints
	CCVindex currentnewvert = numoldverts + numoldfaces;

	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{

		//CCPolygon *thispoly = oldlevel->FaceList(polydex);
		thispoly = oldlevel->FaceList(polydex);

		// Now, need those pointers from polygon to new vertices at next level for
		//   face and for edges

		CCPolygon *neighborpoly;// = oldlevel->FaceList(neighbordex);

		for (int edgedex = 0; edgedex < thispoly->numsides(); edgedex++) 
		{
			CCPindex neighbordex = thispoly->Pneighbors[edgedex];

			//CCPolygon *neighborpoly = oldlevel->FaceList(neighbordex);
			neighborpoly = oldlevel->FaceList(neighbordex);

			// If the across-edge neighbor polygon has a lower index, its edges
			//   have already been computed, so just store pointer to new edgevertex in parent Polygon
			if (polydex > neighbordex) 
			{
				thispoly->newedgeverts[edgedex] = neighborpoly->newedgeverts[thispoly->Pneighboredges[edgedex]];
			}
			// else create a new vertex by averaging the new face and old edge-endpoints
			else 
			{
				CCVertex *facevert1 = newlevel->VertList(thispoly->newfacevert);
				CCVertex *facevert2 = newlevel->VertList(neighborpoly->newfacevert);
				CCVertex *edgevert1 = oldlevel->VertList(thispoly->Pvertices[edgedex]);
				CCVertex *edgevert2 = oldlevel->VertList(thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]);

				newloc = edgevert1->loc + edgevert2->loc + facevert1->loc + facevert2->loc ;
				newnorm = edgevert1->norm + edgevert2->norm + facevert1->norm + facevert2->norm;
				newuvt = edgevert1->uvt + edgevert2->uvt;
				newradius = edgevert1->rad + edgevert2->rad;

				newloc /= 4;
				newnorm /= 4;
				newradius /= 2;
				newuvt /= 2;

				newlevel->VertList[currentnewvert].setvals(newloc, newnorm, newradius, newuvt, 4, true);
				thispoly->newedgeverts[edgedex] = currentnewvert++;

				// all the edge vertex corresponding to the INTER face are tagged
#if defined(QIONG_DEBUG_BLENDING_02192003)
			if(vTag[thispoly->Pvertices[edgedex]]==1&&vTag[thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]]==2 ||
			   vTag[thispoly->Pvertices[edgedex]]==2&&vTag[thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]]==1)
			{
				vRedex[currentnewvert-1]=-1;
				vFlag[currentnewvert-1]=-1;
				vEdgeReduction++;
				//continue;
			}
			else
			{
				vRedex[currentnewvert-1]=currentnewvert-1-vFaceReduction-vEdgeReduction;
				//changed changed changed
				//vFlag[currentnewvert-1]=fTag[polydex];
				vFlag[currentnewvert-1]=vTag[thispoly->Pvertices[edgedex]]>=2?0:1;
			}

//		thispoly->newfacevert = vRedex[numoldverts + polydex];
//#else
//		thispoly->newfacevert = numoldverts + polydex;
#endif

				// Now, add 4x the loc and norm values to each of the two edgeverts
				//   This happens just once, when this vertex is created
				CCVertex *thisVertex1 = newlevel->VertList(thispoly->Pvertices[edgedex]);
				CCVertex *thisVertex2 = newlevel->VertList(thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]);

				thisVertex1->loc += newloc*4.0;
				thisVertex1->norm += newnorm*4.0;
				thisVertex2->loc += newloc*4.0;
				thisVertex2->norm += newnorm*4.0;
			}
		}
	}


	// need to prevent from splitting INTER triangles into 4 sub-quads, instead, we only split it to one quad
#if defined(QIONG_DEBUG_BLENDING_02192003)
	int fReduction=0;		// the # of INTER triangles that should be tagged		
#endif
	// Now, form new polygons and compute all connectivity
	//    Do this in two stages:
	//    1) create new polygons and assign their
	//            numbers to parent polygons' childpolygon lists
	//    2) Use this info to set neighbor and neighboredge lists
	CCPindex newpolydex = 0;

	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		thispoly = oldlevel->FaceList(polydex);

		int nsides = thispoly->numsides();

#if defined(QIONG_DEBUG_BLENDING_02192003)
		if(fTag[polydex]==2)
		{
			thispoly->childPolygons.Create(1);
			thispoly->childPolygons[0] = newpolydex++;
			fReduction++;
		}
		else
		{
			thispoly->childPolygons.Create(nsides);
			for (int vdex = 0; vdex < nsides; vdex++) 
				thispoly->childPolygons[vdex] = newpolydex++;
		}
#else
		thispoly->childPolygons.Create(nsides);

		for (int vdex = 0; vdex < nsides; vdex++) 
			thispoly->childPolygons[vdex] = newpolydex++;
#endif
	}

	CCPolygon * newpoly;

	newpolydex = 0;
	int totalsides = 0;
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		thispoly = oldlevel->FaceList(polydex);

		int nsides = thispoly->numsides();
#if defined(QIONG_DEBUG_BLENDING_02192003)
{
		if(fTag[polydex]==2)
		{
		int sharedTag;
		totalsides += 4;
		//int firstnewpolydex = newpolydex;
		//int lastnewpolydex = newpolydex + 1 - 1;

		for (evenvindex = 0; evenvindex < 1; evenvindex++) 
		{
			int edgedex;

			newpoly = newlevel->FaceList(newpolydex);

			newpoly->Pvertices.Create(4);

			for(edgedex=0; edgedex<3; edgedex++)
				if(vTag[thispoly->Pvertices[edgedex]]==vTag[thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]])
				{
					sharedTag=vTag[thispoly->Pvertices[edgedex]];
					break;
				}

			edgedex=INTMOD(edgedex+2, thispoly->numsides());

			/*
			newpoly->Pvertices[0] = thispoly->newfacevert;
			newpoly->Pvertices[1] = thispoly->newedgeverts[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pvertices[2] = thispoly->Pvertices[evenvindex];
			newpoly->Pvertices[3] = thispoly->newedgeverts[evenvindex];
			*/
			newpoly->Pvertices[0]=thispoly->Pvertices[edgedex];
			newpoly->Pvertices[1]=thispoly->Pvertices[INTMOD(edgedex+1, thispoly->numsides())];
			// work to do
			// need to figure which edge vertex it is!
			newpoly->Pvertices[2]=thispoly->newedgeverts[INTMOD(edgedex+1, thispoly->numsides())];
			//newpoly->Pvertices[2]=thispoly->Pvertices[edgedex];
			newpoly->Pvertices[3]=thispoly->Pvertices[INTMOD(edgedex+2, thispoly->numsides())];

//work to do!
{
			// Now do polygon neighbors
			int facenum;
			newpoly->Pneighbors.Create(4);

			// new 0, old 1
			// The edge number + 1 MOD neighbors-numsides of the neighboredge gives the polygon-index
			//facenum = INTMOD(thispoly->Pneighboredges[evenvindex] + 1,
			//	             oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].numsides());
			//newpoly->Pneighbors[0] = oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].childPolygons[facenum];
			newpoly->Pneighbors[0] = oldlevel->FaceList[thispoly->Pneighbors[edgedex]].childPolygons[0];

			// new 1
			// The edge number of the neighbor polygon gives the polygon-index
			facenum = INTMOD(thispoly->Pneighboredges[INTMOD(edgedex+1, nsides)]+1,
				             oldlevel->FaceList[thispoly->Pneighbors[INTMOD(edgedex+1, nsides)]].numsides());
			//facenum = thispoly->Pneighboredges[INTMOD(evenvindex + 1, nsides)];
			//newpoly->Pneighbors[3] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[facenum];
			newpoly->Pneighbors[1]=oldlevel->FaceList[thispoly->Pneighbors[INTMOD(edgedex+1, nsides)]].childPolygons[facenum];
				//oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[0];

			// new 2
			// The edge number of the neighbor polygon gives the polygon-index
			facenum = thispoly->Pneighboredges[INTMOD(edgedex + 1, nsides)];
			//newpoly->Pneighbors[3] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[facenum];
			newpoly->Pneighbors[2]=oldlevel->FaceList[thispoly->Pneighbors[INTMOD(edgedex+1, nsides)]].childPolygons[facenum];
				//oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[0];

			// new 3, old 2
			// The edge number of the neighbor polygon gives the polygon-index
			//facenum = thispoly->Pneighboredges[INTMOD(evenvindex - 1, nsides)];
			//newpoly->Pneighbors[3] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[facenum];
			newpoly->Pneighbors[3] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(edgedex - 1, nsides)]].childPolygons[0];

			/*
			// The edge number is clockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex - 1 < 0 ? nsides - 1 : evenvindex - 1);
			newpoly->Pneighbors[0] = thispoly->childPolygons[facenum];

			// The edge number of the neighbor polygon gives the polygon-index
			facenum = thispoly->Pneighboredges[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pneighbors[1] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[facenum];

			// The edge number + 1 MOD neighbors-numsides of the neighboredge gives the polygon-index
			facenum = INTMOD(thispoly->Pneighboredges[evenvindex] + 1,
				             oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].numsides());
			newpoly->Pneighbors[2] = oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].childPolygons[facenum];

			// The edge number is counterclockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex + 1 > nsides - 1 ? 0 : evenvindex + 1);
			newpoly->Pneighbors[3] = thispoly->childPolygons[facenum];
			*/

			// -- if all polygons start at the new face vertex, then the edge#'s are all 1
			//    or 2 for neighboring polygons.  (2 facing 1, 1 facing 2)

			// work to do
			CCPolygon *nPoly;
			int startTag, endTag, nsidesNeighbor, i;
			newpoly->Pneighboredges.Create(4);

			// work to do right now
			nPoly=oldlevel->FaceList(thispoly->Pneighbors[edgedex]);
			nsidesNeighbor=nPoly->numsides();
			startTag=vTag[thispoly->Pvertices[INTMOD(edgedex+1, nsides)]];
			endTag=vTag[thispoly->Pvertices[edgedex]];
			for(i=0; i<nsidesNeighbor; i++)
			{
				if(vTag[nPoly->Pvertices[i]]==startTag						   &&
				   vTag[nPoly->Pvertices[INTMOD(i+1, nsidesNeighbor)]]==endTag   )
					break;
			}
			if(i>=nsidesNeighbor)
			{
				int i;
				fprintf(stderr, "\nevenvindex");
				for(i=0; i<nsidesNeighbor; i++)
				{
					fprintf(stderr, "%d \t", vTag[thispoly->Pvertices[i]]);
					/*
					if(vTag[nPoly->Pvertices[i]]==startTag						   &&
						vTag[nPoly->Pvertices[INTMOD(i+1, nsidesNeighbor)]]==endTag   )
						break;
					*/
				}
				fprintf(stderr, "\n");
				for(i=0; i<nsidesNeighbor; i++)
				{
					fprintf(stderr, "%d \t", vTag[nPoly->Pvertices[i]]);
					/*
					if(vTag[nPoly->Pvertices[i]]==startTag						   &&
						vTag[nPoly->Pvertices[INTMOD(i+1, nsidesNeighbor)]]==endTag   )
						break;
					*/
				}
				exit(0);
			}
			int sharedTag2;
			int edgedex2;

			for(edgedex2=0; edgedex2<3; edgedex2++)
				if(vTag[nPoly->Pvertices[edgedex2]]==vTag[nPoly->Pvertices[INTMOD(edgedex2+1, nPoly->numsides())]])
				{
					sharedTag2=vTag[nPoly->Pvertices[edgedex2]];
					break;
				}
			if(sharedTag==sharedTag2)
				newpoly->Pneighboredges[0] = 3;
			else
				newpoly->Pneighboredges[0] = 0;

			newpoly->Pneighboredges[1] = 1;
			newpoly->Pneighboredges[2] = 2;

			// work to do right now
			nPoly=oldlevel->FaceList(thispoly->Pneighbors[INTMOD(edgedex+2, nsides)]);
			nsidesNeighbor=nPoly->numsides();
			startTag=vTag[thispoly->Pvertices[INTMOD(edgedex+2, nsides)]];
			endTag=vTag[thispoly->Pvertices[edgedex]];
			for(i=0; i<nsidesNeighbor; i++)
			{
				if(vTag[nPoly->Pvertices[i]]==startTag					    && 
				   vTag[nPoly->Pvertices[INTMOD(i+1, nsidesNeighbor)]]==endTag   )
					break;
			}
			if(i>=nsidesNeighbor)
			{
				for(int i=0; i<nsidesNeighbor; i++)
				{
					fprintf(stderr, "%d \t", vTag[nPoly->Pvertices[i]]);
					/*
					if(vTag[nPoly->Pvertices[i]]==startTag						   &&
						vTag[nPoly->Pvertices[INTMOD(i+1, nsidesNeighbor)]]==endTag   )
						break;
					*/
				}
				exit(0);
			}
			for(edgedex2 = 0; edgedex2 < 3; edgedex2++)
				if(vTag[nPoly->Pvertices[edgedex2]]==vTag[nPoly->Pvertices[INTMOD(edgedex2+1, nPoly->numsides())]])
				{
					sharedTag2=vTag[nPoly->Pvertices[edgedex2]];
					break;
				}
			if(sharedTag==sharedTag2)
				newpoly->Pneighboredges[3] = 0;
			else
				newpoly->Pneighboredges[3] = 3;
			//newpoly->Pneighboredges[3] = 0;
}

			newpoly->parentPolygon = polydex;
			newpoly->newfacevert = NULLVERTDEX;

			newpolydex++;
		}
		}

		else
		{
		totalsides += nsides*4;

		for (evenvindex = 0; evenvindex < nsides; evenvindex++) 
		{
			newpoly = newlevel->FaceList(newpolydex);

			newpoly->Pvertices.Create(4);

			newpoly->Pvertices[0] = thispoly->newfacevert;
			newpoly->Pvertices[1] = thispoly->newedgeverts[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pvertices[2] = thispoly->Pvertices[evenvindex];
			newpoly->Pvertices[3] = thispoly->newedgeverts[evenvindex];

			// Now do polygon neighbors
			int facenum;
			newpoly->Pneighbors.Create(4);

// work to do
			if(vTag[thispoly->Pvertices[evenvindex]]==1&&vTag[thispoly->Pvertices[INTMOD(evenvindex+1, nsides)]]==1 ||
			   vTag[thispoly->Pvertices[evenvindex]]==2&&vTag[thispoly->Pvertices[INTMOD(evenvindex+1, nsides)]]==2   )
			{
			// The edge number is clockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex - 1 < 0 ? nsides - 1 : evenvindex - 1);
			newpoly->Pneighbors[0] = thispoly->childPolygons[facenum];

			// The edge number of the neighbor polygon gives the polygon-index
			facenum = thispoly->Pneighboredges[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pneighbors[1] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[facenum];

			// The edge number + 1 MOD neighbors-numsides of the neighboredge gives the polygon-index
			//facenum = INTMOD(thispoly->Pneighboredges[evenvindex] + 1,
			//	             oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].numsides());
			newpoly->Pneighbors[2] = oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].childPolygons[0];

			// The edge number is counterclockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex + 1 > nsides - 1 ? 0 : evenvindex + 1);
			newpoly->Pneighbors[3] = thispoly->childPolygons[facenum];

			// -- if all polygons start at the new face vertex, then the edge#'s are all 1
			//    or 2 for neighboring polygons.  (2 facing 1, 1 facing 2)
			newpoly->Pneighboredges.Create(4);
			newpoly->Pneighboredges[0] = 3;
			newpoly->Pneighboredges[1] = 2;
			newpoly->Pneighboredges[2] = 2;
			newpoly->Pneighboredges[3] = 0;
			}
			else
		{
			// bug fixed 02202003
			if(vTag[thispoly->Pvertices[evenvindex]]==1&&vTag[thispoly->Pvertices[INTMOD(evenvindex-1, nsides)]]==1 || 
			   vTag[thispoly->Pvertices[evenvindex]]==2&&vTag[thispoly->Pvertices[INTMOD(evenvindex-1, nsides)]]==2   )
			{
			// The edge number is clockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex - 1 < 0 ? nsides - 1 : evenvindex - 1);
			newpoly->Pneighbors[0] = thispoly->childPolygons[facenum];

			// The edge number of the neighbor polygon gives the polygon-index
			//facenum = thispoly->Pneighboredges[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pneighbors[1] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[0];

			// The edge number + 1 MOD neighbors-numsides of the neighboredge gives the polygon-index
			facenum = INTMOD(thispoly->Pneighboredges[evenvindex] + 1,
				             oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].numsides());
			newpoly->Pneighbors[2] = oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].childPolygons[facenum];

			// The edge number is counterclockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex + 1 > nsides - 1 ? 0 : evenvindex + 1);
			newpoly->Pneighbors[3] = thispoly->childPolygons[facenum];

			// -- if all polygons start at the new face vertex, then the edge#'s are all 1
			//    or 2 for neighboring polygons.  (2 facing 1, 1 facing 2)
			newpoly->Pneighboredges.Create(4);
			newpoly->Pneighboredges[0] = 3;
			newpoly->Pneighboredges[1] = 1;
			newpoly->Pneighboredges[2] = 1;
			newpoly->Pneighboredges[3] = 0;
			}
			else
			{
			// The edge number is clockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex - 1 < 0 ? nsides - 1 : evenvindex - 1);
			newpoly->Pneighbors[0] = thispoly->childPolygons[facenum];

			// The edge number of the neighbor polygon gives the polygon-index
			facenum = thispoly->Pneighboredges[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pneighbors[1] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[facenum];

			// The edge number + 1 MOD neighbors-numsides of the neighboredge gives the polygon-index
			facenum = INTMOD(thispoly->Pneighboredges[evenvindex] + 1,
				             oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].numsides());
			newpoly->Pneighbors[2] = oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].childPolygons[facenum];

			// The edge number is counterclockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex + 1 > nsides - 1 ? 0 : evenvindex + 1);
			newpoly->Pneighbors[3] = thispoly->childPolygons[facenum];

			// -- if all polygons start at the new face vertex, then the edge#'s are all 1
			//    or 2 for neighboring polygons.  (2 facing 1, 1 facing 2)
			newpoly->Pneighboredges.Create(4);
			newpoly->Pneighboredges[0] = 3;
			newpoly->Pneighboredges[1] = 2;
			newpoly->Pneighboredges[2] = 1;
			newpoly->Pneighboredges[3] = 0;
			}
		}
			newpoly->parentPolygon = polydex;
			newpoly->newfacevert = NULLVERTDEX;
//			newpoly->printvals();
//			fprintf(stderr, "\n");
			newpolydex++;
		}
		}
}
#else
		// count number of edges for total
		totalsides += nsides*4;

		for (evenvindex = 0; evenvindex < nsides; evenvindex++) 
		{
			newpoly = newlevel->FaceList(newpolydex);

			newpoly->Pvertices.Create(4);

			newpoly->Pvertices[0] = thispoly->newfacevert;
			newpoly->Pvertices[1] = thispoly->newedgeverts[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pvertices[2] = thispoly->Pvertices[evenvindex];
			newpoly->Pvertices[3] = thispoly->newedgeverts[evenvindex];

			// Now do polygon neighbors

			int facenum;
			newpoly->Pneighbors.Create(4);
			// The edge number is clockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex - 1 < 0 ? nsides - 1 : evenvindex - 1);
			newpoly->Pneighbors[0] = thispoly->childPolygons[facenum];

			// The edge number of the neighbor polygon gives the polygon-index
			facenum = thispoly->Pneighboredges[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pneighbors[1] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[facenum];

			// The edge number + 1 MOD neighbors-numsides of the neighboredge gives the polygon-index
			facenum = INTMOD(thispoly->Pneighboredges[evenvindex] + 1,
				             oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].numsides());
			newpoly->Pneighbors[2] = oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].childPolygons[facenum];

			// The edge number is counterclockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex + 1 > nsides - 1 ? 0 : evenvindex + 1);
			newpoly->Pneighbors[3] = thispoly->childPolygons[facenum];

			// -- if all polygons start at the new face vertex, then the edge#'s are all 1
			//    or 2 for neighboring polygons.  (2 facing 1, 1 facing 2)
			newpoly->Pneighboredges.Create(4);
			newpoly->Pneighboredges[0] = 3;
			newpoly->Pneighboredges[1] = 2;
			newpoly->Pneighboredges[2] = 1;
			newpoly->Pneighboredges[3] = 0;

			newpoly->parentPolygon = polydex;
			newpoly->newfacevert = NULLVERTDEX;
//			newpoly->printvals();
//			fprintf(stderr, "\n");
			newpolydex++;
		}
#endif
	}

	newlevel->nedges = totalsides/2;

	// Now, average the even vertices according to their connectivity
	// ----no, this can be done on the fly.  Zero all nivelo+1 even vertices, and
	// add 4x the new edge values and 1x each new face vertex.
	//    Once here, weight the computed sum and add the weighted parent vertex weight
	// Use the weights:
	//    1/(n^2) for the face-edge sum, and (n-3)/n for the old vertex
	// (See App.A of Thall Diss for the derivation of this.)
	for (evenvindex = 0; evenvindex < numoldverts; evenvindex++)
	{
		double thisval = newlevel->VertList[evenvindex].v_valence;

		double sumweight, oldweight;

		sumweight = 1.0/(thisval*thisval);
		oldweight = (thisval - 3.0)/thisval;

		newlevel->VertList[evenvindex].loc *= sumweight;
		newlevel->VertList[evenvindex].norm *= sumweight;

		newlevel->VertList[evenvindex].loc += oldlevel->VertList[evenvindex].loc * oldweight;
		newlevel->VertList[evenvindex].norm += oldlevel->VertList[evenvindex].norm * oldweight;
	}

	// Now, compute all the connectivity relationships for the CCmesh, as necessary
	// Can probably do this on the fly, too.
	// For each polygon in old mesh:
	//    new children are neighborhood of newfacevert
	//    two new children and their 2-edge and 1-edge neighbors are neighborhood of new edge
	//    new child and its 1-edge neighbor and ITS 1-edge neighbors are
	//       neighborhood of new even vertex
	//    Vertex 2 of each new polygon is old even vertex
	// If VertFaceNeighbors[vdex] has been initialized, then already computed!

	CCPolygon *oldpoly;

	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		oldpoly = oldlevel->FaceList(polydex);

		CCVindex newvertdex = oldpoly->newfacevert;

		// newvalence is number of edges of old polygon, which will be used
		//    all over the place
		int newvalence = newlevel->VertList[newvertdex].v_valence;

		// Compute neighborhood of new face vertex
		/*
		newlevel->VertFaceNeighbors[newvertdex].Create(newvalence);
		newlevel->VertFaceNeighborsVert[newvertdex].Create(newvalence);
		CCPindex newpolydex;
		for (newpolydex = 0; newpolydex < newvalence; newpolydex++) 
		{
			newlevel->VertFaceNeighbors[newvertdex][newpolydex] = oldpoly->childPolygons[newpolydex];
			newlevel->VertFaceNeighborsVert[newvertdex][newpolydex] = 0;
		}
		*/
		newlevel->VertFaceNeighbors[newvertdex].Create(newvalence);
		newlevel->VertFaceNeighborsVert[newvertdex].Create(newvalence);
		CCPindex newpolydex;
		// those face vertices in the INTER triangles should be ignored (eventually removed)
		if(fTag[polydex]==2)
			continue;

		if(fTag[polydex]!=2)
		{
			for (newpolydex = 0; newpolydex < newvalence; newpolydex++) 
			{
				newlevel->VertFaceNeighbors[newvertdex][newpolydex] = oldpoly->childPolygons[newpolydex];
				newlevel->VertFaceNeighborsVert[newvertdex][newpolydex] = 0;
			}
		}
		else
		{
			fprintf(stderr, "error at line 1636 in CCSubdivsurf.cpp!\n");
		}

		// Compute neighborhood of each new edge vertex if VertFaceNeighbors[] not already initialized
		for (int newoddedge = 0; newoddedge < newvalence; newoddedge++) 
		{
			newvertdex = oldpoly->newedgeverts[newoddedge];
			if (newlevel->VertFaceNeighbors[newvertdex].Length() == 0) 
			{
				if(vRedex[newvertdex]!=-1)
				{	
					if(vTag[oldpoly->Pvertices[newoddedge]]==1 && vTag[oldpoly->Pvertices[INTMOD(newoddedge+1, newvalence)]]==1 ||
					   vTag[oldpoly->Pvertices[newoddedge]]==2 && vTag[oldpoly->Pvertices[INTMOD(newoddedge+1, newvalence)]]==2   )
					{
						newlevel->VertList[newvertdex].v_valence=3;
						newlevel->VertFaceNeighbors[newvertdex].Create(3);
						newlevel->VertFaceNeighborsVert[newvertdex].Create(3);

						CCPindex pdexA, pdexB, pdexC; //, pdexD;
						pdexA = oldpoly->childPolygons[INTMOD(newoddedge + 1, newvalence)];
						pdexB = oldpoly->childPolygons[newoddedge];
						pdexC = newlevel->FaceList[pdexB].Pneighbors[2];
						//pdexD = newlevel->FaceList[pdexA].Pneighbors[1];

						newlevel->VertFaceNeighbors[newvertdex][0] = pdexA;
						newlevel->VertFaceNeighbors[newvertdex][1] = pdexB;
						newlevel->VertFaceNeighbors[newvertdex][2] = pdexC;
						//newlevel->VertFaceNeighbors[newvertdex][3] = pdexD;
						newlevel->VertFaceNeighborsVert[newvertdex][0] = 1;
						newlevel->VertFaceNeighborsVert[newvertdex][1] = 3;
						newlevel->VertFaceNeighborsVert[newvertdex][2] = 2;
						//newlevel->VertFaceNeighborsVert[newvertdex][3] = 3;
					}
					else
					{
						newlevel->VertFaceNeighbors[newvertdex].Create(4);
						newlevel->VertFaceNeighborsVert[newvertdex].Create(4);

						CCPindex pdexA, pdexB, pdexC, pdexD;
						pdexA = oldpoly->childPolygons[INTMOD(newoddedge + 1, newvalence)];
						pdexB = oldpoly->childPolygons[newoddedge];
						pdexC = newlevel->FaceList[pdexB].Pneighbors[2];
						pdexD = newlevel->FaceList[pdexA].Pneighbors[1];

						newlevel->VertFaceNeighbors[newvertdex][0] = pdexA;
						newlevel->VertFaceNeighbors[newvertdex][1] = pdexB;
						newlevel->VertFaceNeighbors[newvertdex][2] = pdexC;
						newlevel->VertFaceNeighbors[newvertdex][3] = pdexD;
						newlevel->VertFaceNeighborsVert[newvertdex][0] = 1;
						newlevel->VertFaceNeighborsVert[newvertdex][1] = 3;
						newlevel->VertFaceNeighborsVert[newvertdex][2] = 1;
						newlevel->VertFaceNeighborsVert[newvertdex][3] = 3;
					}
				}
				else
				{
					fprintf(stderr, "error at line 1691 in CCSubdivsurf.cpp!\n");
				}
			}
		}

		//*
		// Compute neighborhood of each new even vertex if VertFaceNeighbors[] not already initialized
		int oldCur, oldCurEdge;
		int oldNeighbor, oldNeighborEdge;

		if(fTag[polydex]!=2)
		{
			for(int childpolydex = 0; childpolydex < newvalence; childpolydex++) 
			{
				CCPindex newface = oldpoly->childPolygons[childpolydex];
				newvertdex = newlevel->FaceList[newface].Pvertices[2];
				if(newlevel->VertFaceNeighbors[newvertdex].Length() == 0) 
				{
					if(vTag[newvertdex]==1 || vTag[newvertdex]==2)
					{
						oldCur=polydex;
						oldCurEdge=INTMOD(childpolydex-1, newvalence);
						oldNeighbor=oldpoly->Pneighbors[oldCurEdge];
						oldNeighborEdge=oldpoly->Pneighboredges[oldCurEdge];
						thispoly=oldpoly;

						int vertvalence = newlevel->VertList[newvertdex].v_valence;
						newlevel->VertFaceNeighbors[newvertdex].Create(vertvalence);
						newlevel->VertFaceNeighborsVert[newvertdex].Create(vertvalence);

						CCPindex pdexE = newface;
						newlevel->VertFaceNeighbors[newvertdex][0] = pdexE;
						newlevel->VertFaceNeighborsVert[newvertdex][0] = 2;

						// work to do right now
						for(int facedex = 1; facedex < vertvalence; facedex++) 
						{
							// case 1
							if(fTag[oldCur]!=2&&fTag[oldNeighbor]!=2)
							{
								pdexE = newlevel->FaceList[pdexE].Pneighbors[1];
								newlevel->VertFaceNeighbors[newvertdex][facedex] = pdexE;
								newlevel->VertFaceNeighborsVert[newvertdex][facedex] = 2;
							}
							else
								// case 2
								if(fTag[oldCur]!=2&&fTag[oldNeighbor]==2)
								{
#ifdef QIONG_DEBUG_BLENDING
	fprintf(stderr, "\nThe Pvertices are: ");
	for(int edgedex=0; edgedex<newlevel->FaceList[pdexE].numsides(); edgedex++)
	{
		fprintf(stderr, "%d\t", newlevel->FaceList[pdexE].Pvertices[edgedex]);
	}
	fprintf(stderr, "\nThe Pneighbors are: ");
	for(edgedex=0; edgedex<newlevel->FaceList[pdexE].numsides(); edgedex++)
	{
		fprintf(stderr, "%d\t", newlevel->FaceList[pdexE].Pneighbors[edgedex]);
	}
#endif
									pdexE = newlevel->FaceList[pdexE].Pneighbors[1];
									newlevel->VertFaceNeighbors[newvertdex][facedex] = pdexE;
									newlevel->VertFaceNeighborsVert[newvertdex][facedex] = 1;
								}
								else
									// case 3
									if(fTag[oldCur]==2&&fTag[oldNeighbor]!=2)
									{
										pdexE = newlevel->FaceList[pdexE].Pneighbors[2];
										newlevel->VertFaceNeighbors[newvertdex][facedex] = pdexE;
										newlevel->VertFaceNeighborsVert[newvertdex][facedex] = 2;
									}
									else
										// case 4
										if(fTag[oldCur]==2&&fTag[oldNeighbor]==2)
										{
											CCPolygon *newPoly=oldlevel->FaceList(oldNeighbor);
											int sharedTag3, sharedTag4;
											for(int edgedex3=0; edgedex3<3; edgedex3++)
												if(vTag[thispoly->Pvertices[edgedex3]]==vTag[thispoly->Pvertices[INTMOD(edgedex3+1, thispoly->numsides())]])
												{
													sharedTag3=vTag[thispoly->Pvertices[edgedex3]];
													break;
												}
											for(int edgedex4=0; edgedex4<3; edgedex4++)
												if(vTag[newPoly->Pvertices[edgedex4]]==vTag[newPoly->Pvertices[INTMOD(edgedex4+1, newPoly->numsides())]])
												{
													sharedTag4=vTag[newPoly->Pvertices[edgedex4]];
													break;
												}
											// case 4.1
											if(sharedTag3==vTag[newvertdex] && sharedTag4==vTag[newvertdex])
											{
												pdexE = newlevel->FaceList[pdexE].Pneighbors[0];
												newlevel->VertFaceNeighbors[newvertdex][facedex] = pdexE;
												newlevel->VertFaceNeighborsVert[newvertdex][facedex] = 3;
											}
											else
												// case 4.2
												if(sharedTag3==vTag[newvertdex] && sharedTag4!=vTag[newvertdex])
												{
													pdexE = newlevel->FaceList[pdexE].Pneighbors[0];
													newlevel->VertFaceNeighbors[newvertdex][facedex] = pdexE;
													newlevel->VertFaceNeighborsVert[newvertdex][facedex] = 0;
												}
												else
													// case 4.3
													if(sharedTag3!=vTag[newvertdex] && sharedTag4==vTag[newvertdex])
													{
														pdexE = newlevel->FaceList[pdexE].Pneighbors[3];
														newlevel->VertFaceNeighbors[newvertdex][facedex] = pdexE;
														newlevel->VertFaceNeighborsVert[newvertdex][facedex] = 3;
													}
													else
														// case 4.4
														if(sharedTag3!=vTag[newvertdex] && sharedTag4!=vTag[newvertdex])
														{
															pdexE = newlevel->FaceList[pdexE].Pneighbors[3];
															newlevel->VertFaceNeighbors[newvertdex][facedex] = pdexE;
															newlevel->VertFaceNeighborsVert[newvertdex][facedex] = 0;
														}
										}

							//thispoly=
							oldCur=oldNeighbor;
							//oldCurEdge=thispoly->Pneighboredges[oldCurEdge];
							thispoly=oldlevel->FaceList(oldNeighbor);
							oldCurEdge=INTMOD(oldCurEdge-1, thispoly->numsides());
							//oldNeighborEdge=thispoly;
							//thispoly=oldlevel->FaceList(oldCurPoly);
							oldNeighbor=thispoly->Pneighbors[oldCurEdge];
							oldNeighborEdge=oldpoly->Pneighboredges[oldCurEdge];
#ifdef QIONG_DEBUG_BLENDING
	for(int edgedex=0; edgedex<thispoly->numsides(); edgedex++)
	{
		fprintf(stderr, "\n%d\t", thispoly->Pvertices[edgedex]);
	}
#endif
						}
					}
					else
					{
						int vertvalence = newlevel->VertList[newvertdex].v_valence;
						newlevel->VertFaceNeighbors[newvertdex].Create(vertvalence);
						newlevel->VertFaceNeighborsVert[newvertdex].Create(vertvalence);

						CCPindex pdexE = newface;
						newlevel->VertFaceNeighbors[newvertdex][0] = pdexE;
						newlevel->VertFaceNeighborsVert[newvertdex][0] = 2;

						for(int facedex = 1; facedex < vertvalence; facedex++) 
						{
							pdexE = newlevel->FaceList[pdexE].Pneighbors[1];
							newlevel->VertFaceNeighbors[newvertdex][facedex] = pdexE;
							newlevel->VertFaceNeighborsVert[newvertdex][facedex] = 2;
						}
					}
				}
			}
		/*
		for (int childpolydex = 0; childpolydex < newvalence; childpolydex++) 
		{
			CCPindex newface = oldpoly->childPolygons[childpolydex];
			newvertdex = newlevel->FaceList[newface].Pvertices[2];
			if (newlevel->VertFaceNeighbors[newvertdex].Length() == 0) 
			{
				int vertvalence = newlevel->VertList[newvertdex].v_valence;
				newlevel->VertFaceNeighbors[newvertdex].Create(vertvalence);
				newlevel->VertFaceNeighborsVert[newvertdex].Create(vertvalence);

				CCPindex pdexE = newface;

				newlevel->VertFaceNeighbors[newvertdex][0] = pdexE;
				newlevel->VertFaceNeighborsVert[newvertdex][0] = 2;

				for (int facedex = 1; facedex < vertvalence; facedex++) {

					pdexE = newlevel->FaceList[pdexE].Pneighbors[1];
					newlevel->VertFaceNeighbors[newvertdex][facedex] = pdexE;
					newlevel->VertFaceNeighborsVert[newvertdex][facedex] = 2;
				}
			}
		}
		*/
		}
		// */
	}

#ifdef QIONG_DEBUG_20030425
	for(polydex=numoldverts+92; polydex<numoldverts+101; polydex++)
	{
		thisVertex=newlevel->VertList(polydex);
		//printf("The vertex %d is: loc=%f %f %f, norm=%f %f %f\n", polydex, thisVertex->loc.x(), thisVertex->loc.y(), \
		//	thisVertex->loc.z(), thisVertex->norm.x(), thisVertex->norm.y(), thisVertex->norm.z());

		//printf("The vertex %d is: loc=%f %f %f\n", polydex, thisVertex->loc.x(), thisVertex->loc.y(), thisVertex->loc.z());
	}
#endif

	// This will cause normal-interpolation to be used at level-2 and higher subdivision
	if(nivelo==1 && NORMINTERP) 
		perturbleveltwo();

	return true;
}
else
	if(nivelo==1)
	{
		return splitandaverageBlending_1(1);	
	}
	else
		return splitandaverage(nivelo);
}


bool CCSubdivsurf::splitandaverageBlending_1(int nivelo)
{
	if (nivelo + 1 >= NUMLEVELS) 
	{
		fprintf(stderr, "ERROR: CCSubdivsurf::splitandaverage() -- Number of levels too high for split.\n");
		return false;
	}

	if(nivelo>1)
		return splitandaverage(nivelo);

	//computeIISmeshBlending(1, 1);

	CCMesh *oldlevel = submeshes(nivelo);
	CCMesh *newlevel = submeshes(nivelo + 1);

	// The new mesh has the even vertices + the new face vertices + the new edge vertices
	int newnumverts = (oldlevel->numverts()-vFaceReduction-vEdgeReduction) + oldlevel->numfaces() + oldlevel->numedges();

	// Each new edge-vertex generates 2 new faces
	//   (As proof, each new edge from the edge-vertex to the 2 face-vertices forms the
	//    lead edge (in counterclockwise order) of one of the new faces)
	int newnumfaces = 2*oldlevel->numedges();

	newlevel->VertList.Create(newnumverts);
	newlevel->FaceList.Create(newnumfaces);
	newlevel->VertFaceNeighbors.Create(newnumverts);
	newlevel->VertFaceNeighborsVert.Create(newnumverts);

	int numoldfaces = oldlevel->numfaces();
	int numoldverts = oldlevel->numverts()-vFaceReduction-vEdgeReduction;

#if defined(QIONG_DEBUG_BLENDING_02192003)
	//int vRedex[1000];			// a re-index of all the vertices in the new mesh
	//int vFlag[1000];
	//vFaceReduction=0;		// the # of face vertices in the new mesh that should be removed
	//vEdgeReduction=0;		// the # of edge vertices in the new mesh that should be removed
#endif

	// Copy even vertices into new vertex list
	// NO, begin by zeroing out the even vertices...use parent's even vertices
	//    (which have the same indices) for computation, and use the new ones
	//    as accumulators for computing the new even vertices.
	// Formula is neweven = oldweight*v_old + neighborweight*(sum[4*newedgeneighbors - newfaceneighbors])
	//    Derivation of this is in Thall diss, App.A.

	// changed changed changed
	//for(i=0; i<numoldverts+vFaceReduction+vEdgeReduction; i++)
	//{
	//}
	int evenvindex;

#ifdef QIONG_DEBUG_02232003
	fprintf(stderr, "\n");
	for (evenvindex = 0; evenvindex < numoldverts+vFaceReduction+vEdgeReduction; evenvindex++)
		fprintf(stderr, "%d: vrdx=%d\t", evenvindex, vRedex[evenvindex]);
#endif
	// COME BACK HERE
	for (evenvindex = 0; evenvindex < numoldverts+vFaceReduction+vEdgeReduction; evenvindex++)
	{
		while(vRedex[evenvindex]==-1)
			evenvindex++;

		// changed changed changed
		newlevel->VertList[vRedex[evenvindex]] = oldlevel->VertList[evenvindex];
		newlevel->VertList[vRedex[evenvindex]].is_oddvert = false;
		newlevel->VertList[vRedex[evenvindex]].loc = NULLVEC;
		newlevel->VertList[vRedex[evenvindex]].norm = NULLVEC;

	// so all the even vertices should be KEPT INTACT
#if defined(QIONG_DEBUG_BLENDING_02192003)
		//vRedex[evenvindex]=evenvindex;
		//vFlag[evenvindex]=vTag[evenvindex];
#endif
	}

	DbVector3 newloc;
	DbVector3 newnorm;
	Medcoord newuvt;
	double newradius;

#if defined(QIONG_DEBUG_BLENDING_02192003)
#endif

	int polydex;
	// Create new face vertices at centroids of their polygons
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		CCPolygon *thispoly = oldlevel->FaceList(polydex);

//		fprintf(stderr, "old polygon #%d is:\n", polydex);
//		thispoly->printvals();
//		fprintf(stderr, "\n");

		newloc.set(0.0, 0.0, 0.0);
		newnorm.set(0.0, 0.0, 0.0);
		newuvt.setvals(0.0, 0.0, 0.0);
		newradius = 0.0;

		int v_index;
		for (v_index = 0; v_index < thispoly->numsides(); v_index++) 
		{
			CCVertex *thisVertex = oldlevel->VertList(thispoly->Pvertices[v_index]);

			// Sum the node values
			newloc += thisVertex->loc;
			newnorm += thisVertex->norm;
			newradius += thisVertex->rad;
			newuvt += thisVertex->uvt;
		}

		newloc /= thispoly->numsides();
		newnorm /= thispoly->numsides();
		newradius /= thispoly->numsides();
		newuvt /= thispoly->numsides();

		newlevel->VertList[numoldverts + polydex].setvals(newloc, newnorm, newradius,
			                                              newuvt, thispoly->numsides(), true);

		// all the face vertex corresponding to the INTER face are tagged
#if defined(QIONG_DEBUG_BLENDING_02192003)
		//if(fTag[polydex]==1)
		//{
		//	vRedex[numoldverts + polydex]=-1;
		//	vFaceReduction++;
			//continue;
		//}
		//else
		//{
		//	vRedex[numoldverts + polydex]=numoldverts + polydex;
		//	vFlag[numoldverts + polydex]=fTag[polydex];
		//}
//		thispoly->newfacevert = vRedex[numoldverts + polydex];
//#else
//		thispoly->newfacevert = numoldverts + polydex;
#endif

		thispoly->newfacevert = numoldverts + polydex;
		thispoly->newedgeverts.Create(thispoly->numsides()); // go ahead and initialize newedgevert list

		// Now, subtract values of the new face vertex into loc and norm of surrounding even vertices
		for (v_index = 0; v_index < thispoly->numsides(); v_index++) 
		{
			// changed changed changed
			CCVertex *thisVertex = newlevel->VertList(vRedex[thispoly->Pvertices[v_index]]);
			thisVertex->loc -= newloc;
			thisVertex->norm -= newnorm;
			// NOTE:  if you make this an addition instead of a subtraction,
			//   it screws everything.  :)  Makes some nice sea anemonies, though
		}
	}

	// For each edge of polygon, create new edge vertex at averaged position of
	//    of old edge endpoints and new face vertices at centroids
	// For r and uvt interpolation, just average old endpoints
	CCVindex currentnewvert = numoldverts + numoldfaces;

	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{

		CCPolygon *thispoly = oldlevel->FaceList(polydex);

		// Now, need those pointers from polygon to new vertices at next level for
		//   face and for edges
		for (int edgedex = 0; edgedex < thispoly->numsides(); edgedex++) 
		{
			CCPindex neighbordex = thispoly->Pneighbors[edgedex];

			CCPolygon *neighborpoly = oldlevel->FaceList(neighbordex);

			// If the across-edge neighbor polygon has a lower index, its edges
			//   have already been computed, so just store pointer to new edgevertex in parent Polygon
			if (polydex > neighbordex) 
			{
				thispoly->newedgeverts[edgedex] = neighborpoly->newedgeverts[thispoly->Pneighboredges[edgedex]];
			}
			// else create a new vertex by averaging the new face and old edge-endpoints
			else 
			{
				CCVertex *facevert1 = newlevel->VertList(thispoly->newfacevert);
				CCVertex *facevert2 = newlevel->VertList(neighborpoly->newfacevert);
				CCVertex *edgevert1 = oldlevel->VertList(thispoly->Pvertices[edgedex]);
				CCVertex *edgevert2 = oldlevel->VertList(thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]);

				newloc = edgevert1->loc + edgevert2->loc + facevert1->loc + facevert2->loc ;
				newnorm = edgevert1->norm + edgevert2->norm + facevert1->norm + facevert2->norm;
				newuvt = edgevert1->uvt + edgevert2->uvt;
				newradius = edgevert1->rad + edgevert2->rad;

				newloc /= 4;
				newnorm /= 4;
				newradius /= 2;
				newuvt /= 2;

				newlevel->VertList[currentnewvert].setvals(newloc, newnorm, newradius, newuvt, 4, true);
				thispoly->newedgeverts[edgedex] = currentnewvert++;

				// all the edge vertex corresponding to the INTER face are tagged

				// changed changed changed
				// Now, add 4x the loc and norm values to each of the two edgeverts
				//   This happens just once, when this vertex is created
				CCVertex *thisVertex1 = newlevel->VertList(vRedex[thispoly->Pvertices[edgedex]]);
				CCVertex *thisVertex2 = newlevel->VertList(vRedex[thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]]);

				thisVertex1->loc += newloc*4.0;
				thisVertex1->norm += newnorm*4.0;
				thisVertex2->loc += newloc*4.0;
				thisVertex2->norm += newnorm*4.0;
			}
		}
	}

	// need to prevent from splitting INTER triangles into 4 sub-quads, instead, we only split it to one quad
#if defined(QIONG_DEBUG_BLENDING_02192003)
	//int fReduction=0;		// the # of INTER triangles that should be tagged		
#endif
	// Now, form new polygons and compute all connectivity
	//    Do this in two stages:
	//    1) create new polygons and assign their
	//            numbers to parent polygons' childpolygon lists
	//    2) Use this info to set neighbor and neighboredge lists
	CCPindex newpolydex = 0;
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		CCPolygon *thispoly = oldlevel->FaceList(polydex);
		int nsides = thispoly->numsides();

#if defined(QIONG_DEBUG_BLENDING_02192003)
		//if(fTag[polydex]==1)
		//{
		//	thispoly->childPolygons.Create(1);
		//	thispoly->childPolygons[0] = newpolydex++;
		//	fReduction++;
		//}
		//else
		{
			thispoly->childPolygons.Create(nsides);
			for (int vdex = 0; vdex < nsides; vdex++) 
				thispoly->childPolygons[vdex] = newpolydex++;
		}
#else
		thispoly->childPolygons.Create(nsides);
		for (int vdex = 0; vdex < nsides; vdex++) 
			thispoly->childPolygons[vdex] = newpolydex++;
#endif
	}

	newpolydex = 0;
	int totalsides = 0;
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		CCPolygon *thispoly = oldlevel->FaceList(polydex);
		int nsides = thispoly->numsides();
#if defined(QIONG_DEBUG_BLENDING_02192003)
{
		// right now
		// work to do
		totalsides += nsides*4;
		//int firstnewpolydex = newpolydex;
		//int lastnewpolydex = newpolydex + nsides - 1;

		for (evenvindex = 0; evenvindex < nsides; evenvindex++) 
		{
			CCPolygon *newpoly = newlevel->FaceList(newpolydex);
			newpoly->Pvertices.Create(4);

			// changed changed changed
			// back here
			newpoly->Pvertices[0] = thispoly->newfacevert;
			newpoly->Pvertices[1] = thispoly->newedgeverts[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pvertices[2] = vRedex[thispoly->Pvertices[evenvindex]];
			newpoly->Pvertices[3] = thispoly->newedgeverts[evenvindex];

			// Now do polygon neighbors
			int facenum;
			newpoly->Pneighbors.Create(4);
// work to do
			{
			// The edge number is clockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex - 1 < 0 ? nsides - 1 : evenvindex - 1);
			newpoly->Pneighbors[0] = thispoly->childPolygons[facenum];

			// The edge number of the neighbor polygon gives the polygon-index
			facenum = thispoly->Pneighboredges[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pneighbors[1] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[facenum];

			// The edge number + 1 MOD neighbors-numsides of the neighboredge gives the polygon-index
			facenum = INTMOD(thispoly->Pneighboredges[evenvindex] + 1,
				             oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].numsides());
			newpoly->Pneighbors[2] = oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].childPolygons[facenum];

			// The edge number is counterclockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex + 1 > nsides - 1 ? 0 : evenvindex + 1);
			newpoly->Pneighbors[3] = thispoly->childPolygons[facenum];

			// -- if all polygons start at the new face vertex, then the edge#'s are all 1
			//    or 2 for neighboring polygons.  (2 facing 1, 1 facing 2)
			newpoly->Pneighboredges.Create(4);
			newpoly->Pneighboredges[0] = 3;
			newpoly->Pneighboredges[1] = 2;
			newpoly->Pneighboredges[2] = 1;
			newpoly->Pneighboredges[3] = 0;
			}
			newpoly->parentPolygon = polydex;
			newpoly->newfacevert = NULLVERTDEX;
//			newpoly->printvals();
//			fprintf(stderr, "\n");
			newpolydex++;
		}
}
#else
		// count number of edges for total
		totalsides += nsides*4;

		//int firstnewpolydex = newpolydex;
		//int lastnewpolydex = newpolydex + nsides - 1;

		for (evenvindex = 0; evenvindex < nsides; evenvindex++) 
		{
			CCPolygon *newpoly = newlevel->FaceList(newpolydex);
			newpoly->Pvertices.Create(4);

			newpoly->Pvertices[0] = thispoly->newfacevert;
			newpoly->Pvertices[1] = thispoly->newedgeverts[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pvertices[2] = thispoly->Pvertices[evenvindex];
			newpoly->Pvertices[3] = thispoly->newedgeverts[evenvindex];

			// Now do polygon neighbors

			int facenum;
			newpoly->Pneighbors.Create(4);
			// The edge number is clockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex - 1 < 0 ? nsides - 1 : evenvindex - 1);
			newpoly->Pneighbors[0] = thispoly->childPolygons[facenum];

			// The edge number of the neighbor polygon gives the polygon-index
			facenum = thispoly->Pneighboredges[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pneighbors[1] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[facenum];

			// The edge number + 1 MOD neighbors-numsides of the neighboredge gives the polygon-index
			facenum = INTMOD(thispoly->Pneighboredges[evenvindex] + 1,
				             oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].numsides());
			newpoly->Pneighbors[2] = oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].childPolygons[facenum];

			// The edge number is counterclockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex + 1 > nsides - 1 ? 0 : evenvindex + 1);
			newpoly->Pneighbors[3] = thispoly->childPolygons[facenum];

			// -- if all polygons start at the new face vertex, then the edge#'s are all 1
			//    or 2 for neighboring polygons.  (2 facing 1, 1 facing 2)
			newpoly->Pneighboredges.Create(4);
			newpoly->Pneighboredges[0] = 3;
			newpoly->Pneighboredges[1] = 2;
			newpoly->Pneighboredges[2] = 1;
			newpoly->Pneighboredges[3] = 0;

			newpoly->parentPolygon = polydex;
			newpoly->newfacevert = NULLVERTDEX;
//			newpoly->printvals();
//			fprintf(stderr, "\n");
			newpolydex++;
		}
#endif
	}

	newlevel->nedges = totalsides/2;

	// Now, average the even vertices according to their connectivity
	// ----no, this can be done on the fly.  Zero all nivelo+1 even vertices, and
	// add 4x the new edge values and 1x each new face vertex.
	//    Once here, weight the computed sum and add the weighted parent vertex weight
	// Use the weights:
	//    1/(n^2) for the face-edge sum, and (n-3)/n for the old vertex
	// (See App.A of Thall Diss for the derivation of this.)

	// changed changed changed
	for (evenvindex = 0; evenvindex < numoldverts+vFaceReduction+vEdgeReduction; evenvindex++)
	{
		while(vRedex[evenvindex]==-1)
			evenvindex++;
	//or (evenvindex = 0; evenvindex < numoldverts; evenvindex++)
	//
		double thisval = newlevel->VertList[vRedex[evenvindex]].v_valence;
		double sumweight, oldweight;

		sumweight = 1.0/(thisval*thisval);
		oldweight = (thisval - 3.0)/thisval;

		newlevel->VertList[vRedex[evenvindex]].loc *= sumweight;
		newlevel->VertList[vRedex[evenvindex]].norm *= sumweight;
		newlevel->VertList[vRedex[evenvindex]].loc += oldlevel->VertList[evenvindex].loc * oldweight;
		newlevel->VertList[vRedex[evenvindex]].norm += oldlevel->VertList[evenvindex].norm * oldweight;
	}

	// Now, compute all the connectivity relationships for the CCmesh, as necessary
	// Can probably do this on the fly, too.
	// For each polygon in old mesh:
	//    new children are neighborhood of newfacevert
	//    two new children and their 2-edge and 1-edge neighbors are neighborhood of new edge
	//    new child and its 1-edge neighbor and ITS 1-edge neighbors are
	//       neighborhood of new even vertex
	//    Vertex 2 of each new polygon is old even vertex
	// If VertFaceNeighbors[vdex] has been initialized, then already computed!
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		CCPolygon *oldpoly = oldlevel->FaceList(polydex);

		// // changed changed changed
		CCVindex newvertdex = oldpoly->newfacevert;
		// newvalence is number of edges of old polygon, which will be used
		//    all over the place
		int newvalence = newlevel->VertList[newvertdex].v_valence;

		// Compute neighborhood of new face vertex
#ifdef QIONG_DEBUG_02232003
		//fprintf(stderr, "newvertdex=%d, newvalence=%d\n", newvertdex, newvalence);
#endif
		newlevel->VertFaceNeighbors[newvertdex].Create(newvalence);
		newlevel->VertFaceNeighborsVert[newvertdex].Create(newvalence);

		CCPindex newpolydex;
		for (newpolydex = 0; newpolydex < newvalence; newpolydex++) {
			newlevel->VertFaceNeighbors[newvertdex][newpolydex] = oldpoly->childPolygons[newpolydex];
			newlevel->VertFaceNeighborsVert[newvertdex][newpolydex] = 0;
		}

		// Compute neighborhood of each new edge vertex if VertFaceNeighbors[] not already initialized
		for (int newoddedge = 0; newoddedge < newvalence; newoddedge++) 
		{
			newvertdex = oldpoly->newedgeverts[newoddedge];
			if (newlevel->VertFaceNeighbors[newvertdex].Length() == 0) 
			{
#ifdef QIONG_DEBUG_02232003
		//fprintf(stderr, "newvertdex=%d, newvalence=%d\n", newvertdex, 4);
#endif
				newlevel->VertFaceNeighbors[newvertdex].Create(4);
				newlevel->VertFaceNeighborsVert[newvertdex].Create(4);

				CCPindex pdexA, pdexB, pdexC, pdexD;
				pdexA = oldpoly->childPolygons[INTMOD(newoddedge + 1, newvalence)];
				pdexB = oldpoly->childPolygons[newoddedge];
				pdexC = newlevel->FaceList[pdexB].Pneighbors[2];
				pdexD = newlevel->FaceList[pdexA].Pneighbors[1];

				newlevel->VertFaceNeighbors[newvertdex][0] = pdexA;
				newlevel->VertFaceNeighbors[newvertdex][1] = pdexB;
				newlevel->VertFaceNeighbors[newvertdex][2] = pdexC;
				newlevel->VertFaceNeighbors[newvertdex][3] = pdexD;
				newlevel->VertFaceNeighborsVert[newvertdex][0] = 1;
				newlevel->VertFaceNeighborsVert[newvertdex][1] = 3;
				newlevel->VertFaceNeighborsVert[newvertdex][2] = 1;
				newlevel->VertFaceNeighborsVert[newvertdex][3] = 3;
			}
		}

		// Compute neighborhood of each new even vertex if VertFaceNeighbors[] not already initialized
		for (int childpolydex = 0; childpolydex < newvalence; childpolydex++) 
		{
			CCPindex newface = oldpoly->childPolygons[childpolydex];
			newvertdex = newlevel->FaceList[newface].Pvertices[2];

			if (newlevel->VertFaceNeighbors[newvertdex].Length() == 0) 
			{
				int vertvalence = newlevel->VertList[newvertdex].v_valence;
#ifdef QIONG_DEBUG_02232003
		//fprintf(stderr, "newvertdex=%d, newvalence=%d\n", newvertdex, vertvalence);
#endif
				newlevel->VertFaceNeighbors[newvertdex].Create(vertvalence);
				newlevel->VertFaceNeighborsVert[newvertdex].Create(vertvalence);

				CCPindex pdexE = newface;

				newlevel->VertFaceNeighbors[newvertdex][0] = pdexE;
				newlevel->VertFaceNeighborsVert[newvertdex][0] = 2;

				for (int facedex = 1; facedex < vertvalence; facedex++) 
				{
					pdexE = newlevel->FaceList[pdexE].Pneighbors[1];
					newlevel->VertFaceNeighbors[newvertdex][facedex] = pdexE;
					newlevel->VertFaceNeighborsVert[newvertdex][facedex] = 2;
				}
			}
		}
	}

	// Lastly, compute the limit mesh
	computelimitmesh(nivelo + 1);

#define QIONG_DEBUG_CALCULATING_NORMALS_1
#ifdef QIONG_DEBUG_CALCULATING_NORMALS_1
	for(int vDex=0; vDex<newlevel->numverts(); vDex++)
	{
		CCPindex thisPDex;
		int neigVDex;
		CCPolygon *thisP;
		CCVertex *neigV;
		DbVector3 newNorm, tempNorm, normalizedNormal;
		DbVector3 newLoc1, newLoc2;

		newNorm.set(0, 0, 0);
		//tempLoc.set(0, 0, 0);

		CCVertex *thisV=newlevel->VertList(vDex);
		//double projLength;
		for(int nNeigh=0; nNeigh<thisV->v_valence; nNeigh++)
		{
			thisPDex=newlevel->VertFaceNeighbors[vDex][nNeigh];
			thisP=newlevel->FaceList(thisPDex);
			/*
			if(vDex==314)
			{
				for(int iii=0; iii<thisP->numsides(); iii++)
				{
					fprintf(stderr, "%d\t", newlevel->VertFaceNeighborsVert[vDex][iii]);
				}
				fprintf(stderr, "\n");
			}
			*/
			neigVDex=thisP->Pvertices[INTMOD(newlevel->VertFaceNeighborsVert[vDex][nNeigh]-1, thisP->numsides())];
			neigV=newlevel->VertList(neigVDex);
			newLoc1=neigV->limitloc;

			thisPDex=newlevel->VertFaceNeighbors[vDex][(nNeigh+1)%thisV->v_valence];
			thisP=newlevel->FaceList(thisPDex);
			/*
			if(vDex==314)
			{
				for(int iii=0; iii<thisP->numsides(); iii++)
				{
					fprintf(stderr, "%d\t", newlevel->VertFaceNeighborsVert[vDex][iii]);
				}
				fprintf(stderr, "\n");
			}
			*/
			neigVDex=thisP->Pvertices[INTMOD(newlevel->VertFaceNeighborsVert[vDex][(nNeigh+1)%thisV->v_valence]-1, thisP->numsides())];
			neigV=newlevel->VertList(neigVDex);
			newLoc2=neigV->limitloc;

			newLoc1=newLoc1-thisV->limitloc;
			newLoc2=newLoc2-thisV->limitloc;
			tempNorm=newLoc1.cross(newLoc2);
			newNorm+=tempNorm;
		}
		//if(vDex==313)
		//	int trap=true;
		//newNorm=newNorm/(thisV->v_valence);
		thisV->limitnorm=newNorm.normalize();
		/*
		newLoc=newLoc-thisV->loc;
		normalizedNormal=thisV->limitnorm.normalize();
		projLength=newLoc.dot(normalizedNormal);
		tempLoc=projLength*normalizedNormal;
		newLoc=newLoc-tempLoc;
		newLoc=thisV->loc+newLoc;
		//thisV->loc.set();
		tempVerts[vDex].loc.set(newLoc.x(), newLoc.y(), newLoc.z());
		*/
	}
#endif

//#define QIONG_DEBUG_SMOOTHING_BLENDING_1
#ifdef QIONG_DEBUG_SMOOTHING_BLENDING_1
	CCVertex *tempVerts=new CCVertex[newlevel->numverts()];
for(int iteration=0; iteration<10; iteration++)
{
	for(vDex=0; vDex<newlevel->numverts(); vDex++)
	//for(int vDex=0; vDex<newlevel->numverts(); vDex++)
	{
		CCPindex thisPDex;
		int neigVDex;
		CCPolygon *thisP;
		CCVertex *neigV;
		DbVector3 newLoc, tempLoc, normalizedNormal;
		newLoc.set(0, 0, 0);
		//tempLoc.set(0, 0, 0);

		CCVertex *thisV=newlevel->VertList(vDex);
		double projLength;
		for(int nNeigh=0; nNeigh<thisV->v_valence; nNeigh++)
		{
			thisPDex=newlevel->VertFaceNeighbors[vDex][nNeigh];
			thisP=newlevel->FaceList(thisPDex);
			/*
			if(vDex==314)
			{
				for(int iii=0; iii<thisP->numsides(); iii++)
				{
					fprintf(stderr, "%d\t", newlevel->VertFaceNeighborsVert[vDex][iii]);
				}
				fprintf(stderr, "\n");
			}
			*/
			neigVDex=thisP->Pvertices[INTMOD(newlevel->VertFaceNeighborsVert[vDex][nNeigh]-1, thisP->numsides())];
			neigV=newlevel->VertList(neigVDex);
			newLoc=newLoc+neigV->loc;
		}
		if(vDex==313)
			int trap=true;
		newLoc=newLoc/(thisV->v_valence);
		newLoc=newLoc-thisV->loc;
		normalizedNormal=thisV->limitnorm.normalize();
		projLength=newLoc.dot(normalizedNormal);
		tempLoc=projLength*normalizedNormal;
		newLoc=newLoc-tempLoc;
		newLoc=thisV->loc+newLoc;
		//thisV->loc.set();
		tempVerts[vDex].loc.set(newLoc.x(), newLoc.y(), newLoc.z());
	}
	for(vDex=0; vDex<newlevel->numverts(); vDex++)
	{
		CCVertex *thisV=newlevel->VertList(vDex);
		thisV->loc.set(tempVerts[vDex].loc.x(), tempVerts[vDex].loc.y(), tempVerts[vDex].loc.z());
	}
}
	delete [] tempVerts;
	// Lastly, compute the limit mesh
	computelimitmesh(nivelo + 1);
#endif

	// This will cause normal-interpolation to be used at level-2 and higher subdivision
	if(nivelo==1 && NORMINTERP) 
		perturbleveltwo();

	return true;
}

/********************************************************************************/
/* SplitandAverage() -- split nivelo to give the next finer level (nivelo + 1)	*/
/********************************************************************************/
bool CCSubdivsurf::splitandaverageBlending_2(int nivelo)
{
	if (nivelo + 1 >= NUMLEVELS) {
		fprintf(stderr, "ERROR: CCSubdivsurf::splitandaverage() -- Number of levels too high for split.\n");
		return false;
	}

	CCMesh *oldlevel = submeshes(nivelo);
	CCMesh *newlevel = submeshes(nivelo + 1);

	// work to do
	// The new mesh has the even vertices + the new face vertices + the new edge vertices
	int newnumverts = oldlevel->numverts() + oldlevel->numfaces() + oldlevel->numedges();

	// work to do
	// Each new edge-vertex generates 2 new faces
	//   (As proof, each new edge from the edge-vertex to the 2 face-vertices forms the
	//    lead edge (in counterclockwise order) of one of the new faces)
	int newnumfaces = 2*oldlevel->numedges();

	// work to do
	newlevel->VertList.Create(newnumverts);
	newlevel->FaceList.Create(newnumfaces);
	newlevel->VertFaceNeighbors.Create(newnumverts);
	newlevel->VertFaceNeighborsVert.Create(newnumverts);

	int numoldfaces = oldlevel->numfaces();
	int numoldverts = oldlevel->numverts();

	// Copy even vertices into new vertex list
	// NO, begin by zeroing out the even vertices...use parent's even vertices
	//    (which have the same indices) for computation, and use the new ones
	//    as accumulators for computing the new even vertices.
	// Formula is neweven = oldweight*v_old + neighborweight*(sum[4*newedgeneighbors - newfaceneighbors])
	//    Derivation of this is in Thall diss, App.A.
	int evenvindex;
	for (evenvindex = 0; evenvindex < numoldverts; evenvindex++) {
		newlevel->VertList[evenvindex] = oldlevel->VertList[evenvindex];
		newlevel->VertList[evenvindex].is_oddvert = false;

		newlevel->VertList[evenvindex].loc = NULLVEC;
		newlevel->VertList[evenvindex].norm = NULLVEC;

	// so all the even vertices should be KEPT INTACT
	}

	DbVector3 newloc;
	DbVector3 newnorm;
	Medcoord newuvt;
	double newradius;

	int polydex;
	// Create new face vertices at centroids of their polygons
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		CCPolygon *thispoly = oldlevel->FaceList(polydex);

//		fprintf(stderr, "old polygon #%d is:\n", polydex);
//		thispoly->printvals();
//		fprintf(stderr, "\n");

		newloc.set(0.0, 0.0, 0.0);
		newnorm.set(0.0, 0.0, 0.0);
		newuvt.setvals(0.0, 0.0, 0.0);
		newradius = 0.0;

		int v_index;
		for (v_index = 0; v_index < thispoly->numsides(); v_index++) {

			CCVertex *thisVertex = oldlevel->VertList(thispoly->Pvertices[v_index]);

			// Sum the node values
			newloc += thisVertex->loc;
			newnorm += thisVertex->norm;
			newradius += thisVertex->rad;
			newuvt += thisVertex->uvt;
		}

		newloc /= thispoly->numsides();
		newnorm /= thispoly->numsides();
		newradius /= thispoly->numsides();
		newuvt /= thispoly->numsides();
		newlevel->VertList[numoldverts + polydex].setvals(newloc, newnorm, newradius,
			                                              newuvt, thispoly->numsides(), true);
		// Now, subtract values of the new face vertex into loc and norm of surrounding even vertices
		for (v_index = 0; v_index < thispoly->numsides(); v_index++) {

			CCVertex *thisVertex = newlevel->VertList(thispoly->Pvertices[v_index]);

			thisVertex->loc -= newloc;
			thisVertex->norm -= newnorm;
			// NOTE:  if you make this an addition instead of a subtraction,
			//   it screws everything.  :)  Makes some nice sea anemonies, though
		}

		// all the face vertex corresponding to the INTER face are tagged
		thispoly->newfacevert = numoldverts + polydex;
		thispoly->newedgeverts.Create(thispoly->numsides()); // go ahead and initialize newedgevert list
	}

	// 02192003
	// need to re-index all the edge vertices in the middle!

	// For each edge of polygon, create new edge vertex at averaged position of
	//    of old edge endpoints and new face vertices at centroids
	// For r and uvt interpolation, just average old endpoints
	CCVindex currentnewvert = numoldverts + numoldfaces;

	for (polydex = 0; polydex < numoldfaces; polydex++) {

		CCPolygon *thispoly = oldlevel->FaceList(polydex);

		// Now, need those pointers from polygon to new vertices at next level for
		//   face and for edges
		for (int edgedex = 0; edgedex < thispoly->numsides(); edgedex++) 
		{
			CCPindex neighbordex = thispoly->Pneighbors[edgedex];

			CCPolygon *neighborpoly = oldlevel->FaceList(neighbordex);

			// If the across-edge neighbor polygon has a lower index, its edges
			//   have already been computed, so just store pointer to new edgevertex in parent Polygon
			if (polydex > neighbordex) 
			{
				thispoly->newedgeverts[edgedex] = neighborpoly->newedgeverts[thispoly->Pneighboredges[edgedex]];
			}
			// else create a new vertex by averaging the new face and old edge-endpoints
			else 
			{
				CCVertex *facevert1 = newlevel->VertList(thispoly->newfacevert);
				CCVertex *facevert2 = newlevel->VertList(neighborpoly->newfacevert);
				CCVertex *edgevert1 = oldlevel->VertList(thispoly->Pvertices[edgedex]);
				CCVertex *edgevert2 = oldlevel->VertList(thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]);

				newloc = edgevert1->loc + edgevert2->loc + facevert1->loc + facevert2->loc ;
				newnorm = edgevert1->norm + edgevert2->norm + facevert1->norm + facevert2->norm;
				newuvt = edgevert1->uvt + edgevert2->uvt;
				newradius = edgevert1->rad + edgevert2->rad;

				newloc /= 4;
				newnorm /= 4;
				newradius /= 2;
				newuvt /= 2;

				newlevel->VertList[currentnewvert].setvals(newloc, newnorm, newradius, newuvt, 4, true);
				thispoly->newedgeverts[edgedex] = currentnewvert++;

				// all the edge vertex corresponding to the INTER face are tagged

				// Now, add 4x the loc and norm values to each of the two edgeverts
				//   This happens just once, when this vertex is created
				CCVertex *thisVertex1 = newlevel->VertList(thispoly->Pvertices[edgedex]);
				CCVertex *thisVertex2 = newlevel->VertList(thispoly->Pvertices[INTMOD(edgedex + 1, thispoly->numsides())]);

				thisVertex1->loc += newloc*4.0;
				thisVertex1->norm += newnorm*4.0;
				thisVertex2->loc += newloc*4.0;
				thisVertex2->norm += newnorm*4.0;
			}
		}
	}

	// Now, form new polygons and compute all connectivity
	//    Do this in two stages:
	//    1) create new polygons and assign their
	//            numbers to parent polygons' childpolygon lists
	//    2) Use this info to set neighbor and neighboredge lists
	CCPindex newpolydex = 0;
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		CCPolygon *thispoly = oldlevel->FaceList(polydex);
		int nsides = thispoly->numsides();

		thispoly->childPolygons.Create(nsides);

		for (int vdex = 0; vdex < nsides; vdex++) 
			thispoly->childPolygons[vdex] = newpolydex++;
	}

	newpolydex = 0;
	int totalsides = 0;
	for (polydex = 0; polydex < numoldfaces; polydex++) 
	{
		CCPolygon *thispoly = oldlevel->FaceList(polydex);
		int nsides = thispoly->numsides();

		// count number of edges for total
		totalsides += nsides*4;

		//int firstnewpolydex = newpolydex;
		//int lastnewpolydex = newpolydex + nsides - 1;

		for (evenvindex = 0; evenvindex < nsides; evenvindex++) {

			CCPolygon *newpoly = newlevel->FaceList(newpolydex);
			newpoly->Pvertices.Create(4);

			newpoly->Pvertices[0] = thispoly->newfacevert;
			newpoly->Pvertices[1] = thispoly->newedgeverts[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pvertices[2] = thispoly->Pvertices[evenvindex];
			newpoly->Pvertices[3] = thispoly->newedgeverts[evenvindex];

			// Now do polygon neighbors

			int facenum;
			newpoly->Pneighbors.Create(4);
			// The edge number is clockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex - 1 < 0 ? nsides - 1 : evenvindex - 1);
			newpoly->Pneighbors[0] = thispoly->childPolygons[facenum];

			// The edge number of the neighbor polygon gives the polygon-index
			facenum = thispoly->Pneighboredges[INTMOD(evenvindex - 1, nsides)];
			newpoly->Pneighbors[1] = oldlevel->FaceList[thispoly->Pneighbors[INTMOD(evenvindex - 1, nsides)]].childPolygons[facenum];

			// The edge number + 1 MOD neighbors-numsides of the neighboredge gives the polygon-index
			facenum = INTMOD(thispoly->Pneighboredges[evenvindex] + 1,
				             oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].numsides());
			newpoly->Pneighbors[2] = oldlevel->FaceList[thispoly->Pneighbors[evenvindex]].childPolygons[facenum];

			// The edge number is counterclockwise neighbor in new neighborhood about face vertex
			facenum = (evenvindex + 1 > nsides - 1 ? 0 : evenvindex + 1);
			newpoly->Pneighbors[3] = thispoly->childPolygons[facenum];

			// -- if all polygons start at the new face vertex, then the edge#'s are all 1
			//    or 2 for neighboring polygons.  (2 facing 1, 1 facing 2)
			newpoly->Pneighboredges.Create(4);
			newpoly->Pneighboredges[0] = 3;
			newpoly->Pneighboredges[1] = 2;
			newpoly->Pneighboredges[2] = 1;
			newpoly->Pneighboredges[3] = 0;

			newpoly->parentPolygon = polydex;
			newpoly->newfacevert = NULLVERTDEX;
//			newpoly->printvals();
//			fprintf(stderr, "\n");
			newpolydex++;
		}
	}

	newlevel->nedges = totalsides/2;

	// Now, average the even vertices according to their connectivity
	// ----no, this can be done on the fly.  Zero all nivelo+1 even vertices, and
	// add 4x the new edge values and 1x each new face vertex.
	//    Once here, weight the computed sum and add the weighted parent vertex weight
	// Use the weights:
	//    1/(n^2) for the face-edge sum, and (n-3)/n for the old vertex
	// (See App.A of Thall Diss for the derivation of this.)
	for (evenvindex = 0; evenvindex < numoldverts; evenvindex++) {
		double thisval = newlevel->VertList[evenvindex].v_valence;

		double sumweight, oldweight;

		sumweight = 1.0/(thisval*thisval);
		oldweight = (thisval - 3.0)/thisval;

		newlevel->VertList[evenvindex].loc *= sumweight;
		newlevel->VertList[evenvindex].norm *= sumweight;

		newlevel->VertList[evenvindex].loc += oldlevel->VertList[evenvindex].loc * oldweight;
		newlevel->VertList[evenvindex].norm += oldlevel->VertList[evenvindex].norm * oldweight;
	}

	// Now, compute all the connectivity relationships for the CCmesh, as necessary
	// Can probably do this on the fly, too.
	// For each polygon in old mesh:
	//    new children are neighborhood of newfacevert
	//    two new children and their 2-edge and 1-edge neighbors are neighborhood of new edge
	//    new child and its 1-edge neighbor and ITS 1-edge neighbors are
	//       neighborhood of new even vertex
	//    Vertex 2 of each new polygon is old even vertex
	// If VertFaceNeighbors[vdex] has been initialized, then already computed!
	for (polydex = 0; polydex < numoldfaces; polydex++) {
		CCPolygon *oldpoly = oldlevel->FaceList(polydex);

		CCVindex newvertdex = oldpoly->newfacevert;

		// newvalence is number of edges of old polygon, which will be used
		//    all over the place
		int newvalence = newlevel->VertList[newvertdex].v_valence;

		// Compute neighborhood of new face vertex
		newlevel->VertFaceNeighbors[newvertdex].Create(newvalence);
		newlevel->VertFaceNeighborsVert[newvertdex].Create(newvalence);

		CCPindex newpolydex;
		for (newpolydex = 0; newpolydex < newvalence; newpolydex++) {
			newlevel->VertFaceNeighbors[newvertdex][newpolydex] = oldpoly->childPolygons[newpolydex];
			newlevel->VertFaceNeighborsVert[newvertdex][newpolydex] = 0;
		}

		// Compute neighborhood of each new edge vertex if VertFaceNeighbors[] not already initialized
		for (int newoddedge = 0; newoddedge < newvalence; newoddedge++) {
			newvertdex = oldpoly->newedgeverts[newoddedge];
			if (newlevel->VertFaceNeighbors[newvertdex].Length() == 0) {
				newlevel->VertFaceNeighbors[newvertdex].Create(4);
				newlevel->VertFaceNeighborsVert[newvertdex].Create(4);

				CCPindex pdexA, pdexB, pdexC, pdexD;
				pdexA = oldpoly->childPolygons[INTMOD(newoddedge + 1, newvalence)];
				pdexB = oldpoly->childPolygons[newoddedge];
				pdexC = newlevel->FaceList[pdexB].Pneighbors[2];
				pdexD = newlevel->FaceList[pdexA].Pneighbors[1];

				newlevel->VertFaceNeighbors[newvertdex][0] = pdexA;
				newlevel->VertFaceNeighbors[newvertdex][1] = pdexB;
				newlevel->VertFaceNeighbors[newvertdex][2] = pdexC;
				newlevel->VertFaceNeighbors[newvertdex][3] = pdexD;
				newlevel->VertFaceNeighborsVert[newvertdex][0] = 1;
				newlevel->VertFaceNeighborsVert[newvertdex][1] = 3;
				newlevel->VertFaceNeighborsVert[newvertdex][2] = 1;
				newlevel->VertFaceNeighborsVert[newvertdex][3] = 3;
			}
		}

		// Compute neighborhood of each new even vertex if VertFaceNeighbors[] not already initialized
		for (int childpolydex = 0; childpolydex < newvalence; childpolydex++) 
		{

			CCPindex newface = oldpoly->childPolygons[childpolydex];
			newvertdex = newlevel->FaceList[newface].Pvertices[2];

			if (newlevel->VertFaceNeighbors[newvertdex].Length() == 0) {

				int vertvalence = newlevel->VertList[newvertdex].v_valence;
				newlevel->VertFaceNeighbors[newvertdex].Create(vertvalence);
				newlevel->VertFaceNeighborsVert[newvertdex].Create(vertvalence);

				CCPindex pdexE = newface;

				newlevel->VertFaceNeighbors[newvertdex][0] = pdexE;
				newlevel->VertFaceNeighborsVert[newvertdex][0] = 2;

				for (int facedex = 1; facedex < vertvalence; facedex++) {

					pdexE = newlevel->FaceList[pdexE].Pneighbors[1];
					newlevel->VertFaceNeighbors[newvertdex][facedex] = pdexE;
					newlevel->VertFaceNeighborsVert[newvertdex][facedex] = 2;
				}
			}
		}
	}

	// Lastly, compute the limit mesh
	computelimitmesh(nivelo + 1);

#define QIONG_DEBUG_CALCULATING_NORMALS_2
#ifdef QIONG_DEBUG_CALCULATING_NORMALS_2
	for(int vDex=0; vDex<newlevel->numverts(); vDex++)
	//for(int vDex=0; vDex<newlevel->numverts(); vDex++)
	{
		CCPindex thisPDex;
		int neigVDex;
		CCPolygon *thisP;
		CCVertex *neigV;
		DbVector3 newNorm, tempNorm, normalizedNormal;
		DbVector3 newLoc1, newLoc2;
		newNorm.set(0, 0, 0);
		//tempLoc.set(0, 0, 0);

		CCVertex *thisV=newlevel->VertList(vDex);
		//double projLength;
		for(int nNeigh=0; nNeigh<thisV->v_valence; nNeigh++)
		{
			thisPDex=newlevel->VertFaceNeighbors[vDex][nNeigh];
			thisP=newlevel->FaceList(thisPDex);
			/*
			if(vDex==314)
			{
				for(int iii=0; iii<thisP->numsides(); iii++)
				{
					fprintf(stderr, "%d\t", newlevel->VertFaceNeighborsVert[vDex][iii]);
				}
				fprintf(stderr, "\n");
			}
			*/
			neigVDex=thisP->Pvertices[INTMOD(newlevel->VertFaceNeighborsVert[vDex][nNeigh]-1, thisP->numsides())];
			neigV=newlevel->VertList(neigVDex);
			newLoc1=neigV->limitloc;

			thisPDex=newlevel->VertFaceNeighbors[vDex][(nNeigh+1)%thisV->v_valence];
			thisP=newlevel->FaceList(thisPDex);
			/*
			if(vDex==314)
			{
				for(int iii=0; iii<thisP->numsides(); iii++)
				{
					fprintf(stderr, "%d\t", newlevel->VertFaceNeighborsVert[vDex][iii]);
				}
				fprintf(stderr, "\n");
			}
			*/
			neigVDex=thisP->Pvertices[INTMOD(newlevel->VertFaceNeighborsVert[vDex][(nNeigh+1)%thisV->v_valence]-1, thisP->numsides())];
			neigV=newlevel->VertList(neigVDex);
			newLoc2=neigV->limitloc;

			newLoc1=newLoc1-thisV->limitloc;
			newLoc2=newLoc2-thisV->limitloc;
			tempNorm=newLoc1.cross(newLoc2);
			newNorm+=tempNorm;
		}
		//if(vDex==313)
		//	int trap=true;
		//newNorm=newNorm/(thisV->v_valence);
		thisV->limitnorm=newNorm.normalize();
		/*
		newLoc=newLoc-thisV->loc;
		normalizedNormal=thisV->limitnorm.normalize();
		projLength=newLoc.dot(normalizedNormal);
		tempLoc=projLength*normalizedNormal;
		newLoc=newLoc-tempLoc;
		newLoc=thisV->loc+newLoc;
		//thisV->loc.set();
		tempVerts[vDex].loc.set(newLoc.x(), newLoc.y(), newLoc.z());
		*/
	}
#endif

//#define QIONG_DEBUG_SMOOTHING_BLENDING_2
#ifdef QIONG_DEBUG_SMOOTHING_BLENDING_2
	//int vDex;
	CCVertex *tempVerts=new CCVertex[newlevel->numverts()];
for(int iteration=0; iteration<2; iteration++)
{
	for(vDex=0; vDex<newlevel->numverts(); vDex++)
	//for(int vDex=0; vDex<newlevel->numverts(); vDex++)
	{
		CCPindex thisPDex;
		int neigVDex;
		CCPolygon *thisP;
		CCVertex *neigV;
		DbVector3 avgLoc, tempLoc, normalizedNormal;
		avgLoc.set(0, 0, 0);
		//tempLoc.set(0, 0, 0);

		CCVertex *thisV=newlevel->VertList(vDex);
		double projLength;
		for(int nNeigh=0; nNeigh<thisV->v_valence; nNeigh++)
		{
			thisPDex=newlevel->VertFaceNeighbors[vDex][nNeigh];
			thisP=newlevel->FaceList(thisPDex);
			/*
			if(vDex==314)
			{
				for(int iii=0; iii<thisP->numsides(); iii++)
				{
					fprintf(stderr, "%d\t", newlevel->VertFaceNeighborsVert[vDex][iii]);
				}
				fprintf(stderr, "\n");
			}
			*/
			neigVDex=thisP->Pvertices[INTMOD(newlevel->VertFaceNeighborsVert[vDex][nNeigh]+1, thisP->numsides())];
			neigV=newlevel->VertList(neigVDex);
			avgLoc=avgLoc+neigV->loc;
		}
		if(vDex==313)
			int trap=true;
		avgLoc=avgLoc/(thisV->v_valence);
		avgLoc=avgLoc-thisV->loc;
		normalizedNormal=thisV->limitnorm.normalize();
		projLength=avgLoc.dot(normalizedNormal);
		tempLoc=projLength*normalizedNormal;
		avgLoc=avgLoc-tempLoc;
		avgLoc=thisV->loc+avgLoc;
		//thisV->loc.set();
		tempVerts[vDex].loc.set(avgLoc.x(), avgLoc.y(), avgLoc.z());
	}
	for(vDex=0; vDex<newlevel->numverts(); vDex++)
	{
		CCVertex *thisV=newlevel->VertList(vDex);
		thisV->loc.set(tempVerts[vDex].loc.x(), tempVerts[vDex].loc.y(), tempVerts[vDex].loc.z());
	}
}
	delete [] tempVerts;
	// Lastly, compute the limit mesh
	computelimitmesh(nivelo + 1);
#endif

	// This will cause normal-interpolation to be used at level-2 and higher subdivision
	if(nivelo==1 && NORMINTERP) 
		perturbleveltwo();
	return true;
}







bool CCSubdivsurf::computelimitmeshBlending(int nivelo)
{
	CCMesh *thislevel = submeshes(nivelo);

	int nverts = thislevel->numverts();
	if (nverts == 0)
		return false;

	DbVector3 c2tangent, c3tangent;
	int vertvalence;
	// For each vertex in thislevel, send it to its limit pposition
	//    and compute the limit normal as well
	for (CCVindex vdex = 0; vdex < nverts; vdex++) 
	{
		if(vRedex[vdex]==-1)
			continue;

		CCVertex *thisVert = thislevel->VertList(vdex);
		vertvalence = thisVert->v_valence;

		// for each polygonal neighbor, form the weighted sum of the e_i edge
		//    and the f_i face (possibly computed, for a non-quad) into the
		//    limitloc and limitnorm DbVectors
		thisVert->limitloc = thisVert->loc * (vertvalence*vertvalence);
		c2tangent.set(0.0, 0.0, 0.0);
		c3tangent.set(0.0, 0.0, 0.0);

#ifdef QIONG_DEBUG_02232003
			if(vdex==100)
			{
				int i;
				fprintf(stderr, "\n");
				for(i=0; i<vertvalence; i++)
					fprintf(stderr, "thislevel->VertFaceNeighbors[vdex][i]=%d, thislevel->FaceList[thispolydex]=%d\n", thislevel->VertFaceNeighbors[vdex][i]), thislevel->FaceList[thislevel->VertFaceNeighbors[vdex][i]];
			}
#endif

		CCPindex pdex;
		for (pdex = 0; pdex < vertvalence; pdex++) 
		{
			CCPindex thispolydex = thislevel->VertFaceNeighbors[vdex][pdex];
			CCPolygon *thispoly = thislevel->FaceList(thispolydex);

			int polyvdex = thislevel->VertFaceNeighborsVert[vdex][pdex];

			DbVector3 e_j, f_j;
			int pnumsides = thispoly->numsides();

			// The edge vertex is the next vertex counterclockwise from vpolydex
			CCVertex *e_pdex = thislevel->VertList(thispoly->Pvertices[INTMOD(polyvdex + 1, pnumsides)]);
			e_j = e_pdex->loc;

			// If pnumsides == 4, 2nd vertex counterclockwise from vpolydex is face vertex
			//    else have to regularize and just compute regloc
			if (pnumsides == 4) {
				CCVertex *f_pdex = thislevel->VertList(thispoly->Pvertices[INTMOD(polyvdex + 2, pnumsides)]);
				f_j = f_pdex->loc;
			}
			else
				f_j = thislevel->regularize(thispolydex, polyvdex);

			// Now, add weighted values to limitloc and {c2,c3}tangent
			thisVert->limitloc += (e_j * 4.0) + f_j;

			c2tangent += (e_j * Etable[vertvalence][pdex + 1]) + (f_j * Ftable[vertvalence][pdex + 1]);
			c3tangent += (e_j * Etable[vertvalence][pdex    ]) + (f_j * Ftable[vertvalence][pdex    ]);
		}

		thisVert->limitloc /= (double) vertvalence*(vertvalence + 5);
		thisVert->limitnorm = c2tangent.cross(c3tangent);
	}
	return true;
}



// 02052003
bool CCSubdivsurf::computeIISmesh2(int numlvl, int numiter)
{
	if(numlvl==0)
	{
		return computeIISmesh();
	}
	else
	{
	CCVertex *fatherVerts;
	CCMesh *IISmesh = submeshes(numlvl);
	//CCMesh *dadmesh = &fathermesh;

	int nverts = IISmesh->numverts();
	if (nverts == 0)
		return false;

	CCVindex vdex;
	fatherVerts=new CCVertex[nverts];
	Medcoord a(0, 0, 0);
	for(vdex=0; vdex<nverts; vdex++)
	{
		CCVertex *thisVert = IISmesh->VertList(vdex);
		fatherVerts[vdex].setvals(thisVert->limitloc, thisVert->limitnorm, 0, a, 0, 0);
	}

	// For each iteration, perturb vertex submeshes[0].VertList[j] from their
	//    initialized values (same as values in fathermesh)
	for (int icnt = 0; icnt < numiter; icnt++) 
	{
		int vertvalence;

		// For each vertex in thislevel, send it to its inverse limit position
		// --This can't be done in-place, so use limitloc as a temporary variable
		// (should experiment---perhaps CAN do in-place)
		// REWRITE:  can do in-place, and thus, get Gauss-Seidel which should
		//    converge faster than Jacobi, and is faster per iteration
		for(vdex = 0; vdex < nverts; vdex++) {

			CCVertex *thisVert = IISmesh->VertList(vdex);
			//CCVertex *dadvert = dadmesh->VertList(vdex);
			vertvalence = thisVert->v_valence;
			DbVector3 iterant;

			// Formula:  as described in ThallDiss Ch. 4
			//    v_j^[i+1] = 1/2 * [ v_j^[i]  + (5+valence)/valence * v_limit_j
			//                                 - 1/(valence^2) * SUM[4*e_k + f_k]
	        //						]
			thisVert->limitloc = thisVert->loc;

			thisVert->limitloc += fatherVerts[vdex].loc * (vertvalence + 5.0)/((double) vertvalence);
			//thisVert->limitloc += dadvert->loc * (vertvalence + 5.0)/((double) vertvalence);
			iterant.set(0.0, 0.0, 0.0);

			// for each polygonal neighbor, form the weighted sum of the e_i edge
			//    and the f_i face (possibly computed, for a non-quad) into the
			//    limitloc and limitnorm DbVectors
			CCPindex pdex;
			for (pdex = 0; pdex < vertvalence; pdex++) {

				CCPindex thispolydex = IISmesh->VertFaceNeighbors[vdex][pdex];
				CCPolygon *thispoly = IISmesh->FaceList(thispolydex);

				int polyvdex = IISmesh->VertFaceNeighborsVert[vdex][pdex];

				DbVector3 e_j, f_j;
				int pnumsides = thispoly->numsides();

				// The edge vertex is the next vertex counterclockwise from polyvdex
				CCVertex *e_pdex = IISmesh->VertList(thispoly->Pvertices[INTMOD(polyvdex + 1, pnumsides)]);
				e_j = e_pdex->loc;

				// If pnumsides == 4, 2nd vertex counterclockwise from vpolydex is face vertex
				//    else have to regularize and just compute regloc
				if (pnumsides == 4) {
					CCVertex *f_pdex = IISmesh->VertList(thispoly->Pvertices[INTMOD(polyvdex + 2, pnumsides)]);
					f_j = f_pdex->loc;
				}
				else
					f_j = IISmesh->regularize(thispolydex, polyvdex);

				// Now, add weighted values to limitloc and {c2,c3}tangent
				iterant += (e_j * 4.0) + f_j;
			}

			iterant /= (double) vertvalence*vertvalence;
			thisVert->limitloc -= iterant;

			thisVert->limitloc /= 2.0;
		}
		// After each complete iteration through all vertices,
		//    copy values back from limitloc to loc
		for (CCVindex vdex2 = 0; vdex2 < nverts; vdex2++) {
			CCVertex *thisVert = IISmesh->VertList(vdex2);

			thisVert->loc = thisVert->limitloc;
		}
	}
	delete []fatherVerts;
	return true;	
	}
}


// 02052003
bool CCSubdivsurf::computeIISmesh(int numlvl, int numiter)
{
	if(numlvl==0)
	{
		return computeIISmesh();
	}
	else
	{
	CCMesh *IISmesh = submeshes(numlvl);
	CCMesh *dadmesh = &fathermesh;

	int nverts = IISmesh->numverts();
	if (nverts == 0)
		return false;

	// For each iteration, perturb vertex submeshes[0].VertList[j] from their
	//    initialized values (same as values in fathermesh)
	for (int icnt = 0; icnt < numiter; icnt++) {
		int vertvalence;

		// For each vertex in thislevel, send it to its inverse limit position
		// --This can't be done in-place, so use limitloc as a temporary variable
		// (should experiment---perhaps CAN do in-place)
		// REWRITE:  can do in-place, and thus, get Gauss-Seidel which should
		//    converge faster than Jacobi, and is faster per iteration
		for (CCVindex vdex = 0; vdex < nverts; vdex++) {

			CCVertex *thisVert = IISmesh->VertList(vdex);
			CCVertex *dadvert = dadmesh->VertList(vdex);
			vertvalence = thisVert->v_valence;
			DbVector3 iterant;

			// Formula:  as described in ThallDiss Ch. 4
			//    v_j^[i+1] = 1/2 * [ v_j^[i]  + (5+valence)/valence * v_limit_j
			//                                 - 1/(valence^2) * SUM[4*e_k + f_k]
	        //						]
			thisVert->limitloc = thisVert->loc;

			thisVert->limitloc += dadvert->loc * (vertvalence + 5.0)/((double) vertvalence);
			iterant.set(0.0, 0.0, 0.0);

			// for each polygonal neighbor, form the weighted sum of the e_i edge
			//    and the f_i face (possibly computed, for a non-quad) into the
			//    limitloc and limitnorm DbVectors
			CCPindex pdex;
			for (pdex = 0; pdex < vertvalence; pdex++) {

				CCPindex thispolydex = IISmesh->VertFaceNeighbors[vdex][pdex];
				CCPolygon *thispoly = IISmesh->FaceList(thispolydex);

				int polyvdex = IISmesh->VertFaceNeighborsVert[vdex][pdex];

				DbVector3 e_j, f_j;
				int pnumsides = thispoly->numsides();

				// The edge vertex is the next vertex counterclockwise from polyvdex
				CCVertex *e_pdex = IISmesh->VertList(thispoly->Pvertices[INTMOD(polyvdex + 1, pnumsides)]);
				e_j = e_pdex->loc;

				// If pnumsides == 4, 2nd vertex counterclockwise from vpolydex is face vertex
				//    else have to regularize and just compute regloc
				if (pnumsides == 4) {
					CCVertex *f_pdex = IISmesh->VertList(thispoly->Pvertices[INTMOD(polyvdex + 2, pnumsides)]);
					f_j = f_pdex->loc;
				}
				else
					f_j = IISmesh->regularize(thispolydex, polyvdex);

				// Now, add weighted values to limitloc and {c2,c3}tangent
				iterant += (e_j * 4.0) + f_j;
			}

			iterant /= (double) vertvalence*vertvalence;
			thisVert->limitloc -= iterant;
			thisVert->limitloc /= 2.0;
		}
		// After each complete iteration through all vertices,
		//    copy values back from limitloc to loc
		for (CCVindex vdex2 = 0; vdex2 < nverts; vdex2++) {
			CCVertex *thisVert = IISmesh->VertList(vdex2);

			thisVert->loc = thisVert->limitloc;
		}
	}
	return true;	
	}
}

/*
bool CCSubdivsurf::computeFatherMeshBlending(Pointlist_server2 **pList)
{
	int vertdex;
	CCMesh *thisMesh=smeshes(1);
	CCVertex *thisVertex;
	Bpoint bp;
	for(vertdex=0; vertdex<thisMesh->numverts(); vertdex++)
	{
		while(vRedex[vertdex]==-1)
			vertdex++;
		thisVertex=thisMesh->VertList(vertdex);
		// the vertex is on the sub-figure
		if(vFlag[vertdex]<=1)
		{
			pList[1]->subdivBPosition(&bp, thisVertex->uvt.u, thisVertex->uvt.v, thisVertex->uvt.t);
			thisVertex->loc.set(bp.pnt[0], bp.pnt[1], bp.pnt[2]);
			thisVertex->norm.set(bp.norm[0], bp.norm[1], bp.norm[2]);
			thisVertex->rad=bp.rad;
		}
		// the vertex is on the main-figure
		else
		{
			pList[0]->subdivBPosition(&bp, thisVertex->uvt.u, thisVertex->uvt.v, thisVertex->uvt.t);
			thisVertex->loc.set(bp.pnt[0], bp.pnt[1], bp.pnt[2]);
			thisVertex->norm.set(bp.norm[0], bp.norm[1], bp.norm[2]);
			thisVertex->rad=bp.rad;
		}
	}
	return true;
}
*/

// 02052003
bool CCSubdivsurf::computeIISmeshBlending(int numlvl, int numiter)
{
	if(numlvl==0)
	{
		return computeIISmesh();
	}
	else
	{
	CCMesh *IISmesh = submeshes(numlvl);
	//CCMesh *dadmesh = &fathermesh;

	int nverts = IISmesh->numverts();
	if (nverts == 0)
		return false;

	int i;
	CCVertex *dadVerts=new CCVertex[nverts-vFaceReduction-vEdgeReduction];
	for(i=0; i<nverts; i++)
	{
		while(vRedex[i]==-1) i++;
		dadVerts[vRedex[i]].loc=IISmesh->VertList[i].loc;
	}

	// For each iteration, perturb vertex submeshes[0].VertList[j] from their
	//    initialized values (same as values in fathermesh)
	for (int icnt = 0; icnt < numiter; icnt++) {
		int vertvalence;

		// For each vertex in thislevel, send it to its inverse limit position
		// --This can't be done in-place, so use limitloc as a temporary variable
		// (should experiment---perhaps CAN do in-place)
		// REWRITE:  can do in-place, and thus, get Gauss-Seidel which should
		//    converge faster than Jacobi, and is faster per iteration
		for (CCVindex vdex = 0; vdex < nverts; vdex++) 
		{
			if(vRedex[vdex]==-1)
				continue;

			CCVertex *thisVert = IISmesh->VertList(vdex);
			//CCVertex *dadvert = dadmesh->VertList(vdex);
			vertvalence = thisVert->v_valence;
			DbVector3 iterant;

			// Formula:  as described in ThallDiss Ch. 4
			//    v_j^[i+1] = 1/2 * [ v_j^[i]  + (5+valence)/valence * v_limit_j
			//                                 - 1/(valence^2) * SUM[4*e_k + f_k]
	        //						]
			thisVert->limitloc = thisVert->loc;

			//thisVert->limitloc += dadvert->loc * (vertvalence + 5.0)/((double) vertvalence);
			thisVert->limitloc += dadVerts[vRedex[vdex]].loc * (vertvalence + 5.0)/((double) vertvalence);
			iterant.set(0.0, 0.0, 0.0);

			// for each polygonal neighbor, form the weighted sum of the e_i edge
			//    and the f_i face (possibly computed, for a non-quad) into the
			//    limitloc and limitnorm DbVectors
			CCPindex pdex;
			for (pdex = 0; pdex < vertvalence; pdex++) {

				CCPindex thispolydex = IISmesh->VertFaceNeighbors[vdex][pdex];
				CCPolygon *thispoly = IISmesh->FaceList(thispolydex);

				int polyvdex = IISmesh->VertFaceNeighborsVert[vdex][pdex];

				DbVector3 e_j, f_j;
				int pnumsides = thispoly->numsides();

				// The edge vertex is the next vertex counterclockwise from polyvdex
				CCVertex *e_pdex = IISmesh->VertList(thispoly->Pvertices[INTMOD(polyvdex + 1, pnumsides)]);
				e_j = e_pdex->loc;

				// If pnumsides == 4, 2nd vertex counterclockwise from vpolydex is face vertex
				//    else have to regularize and just compute regloc
				if (pnumsides == 4) {
					CCVertex *f_pdex = IISmesh->VertList(thispoly->Pvertices[INTMOD(polyvdex + 2, pnumsides)]);
					f_j = f_pdex->loc;
				}
				else
				{
					// error!
					fprintf(stderr, "In computeIISMeshBlending, ahouldn't be any shapes but quads! \n");
					f_j = IISmesh->regularize(thispolydex, polyvdex);
				}

				// Now, add weighted values to limitloc and {c2,c3}tangent
				iterant += (e_j * 4.0) + f_j;
			}

			iterant /= (double) vertvalence*vertvalence;
			thisVert->limitloc -= iterant;
			thisVert->limitloc /= 2.0;
		}
		// After each complete iteration through all vertices,
		//    copy values back from limitloc to loc
		for (CCVindex vdex2 = 0; vdex2 < nverts; vdex2++) 
		{
			if(vRedex[vdex2]==-1) 
			{
				CCVertex *thisVert = IISmesh->VertList(vdex2);
				thisVert->loc.set(-99999, -99999, -99999);
				continue;
			}

			CCVertex *thisVert = IISmesh->VertList(vdex2);
			thisVert->loc = thisVert->limitloc;
		}
	}
	delete []dadVerts;
	return true;
	}
}








void CCSubdivsurf::initBlending(int numverts, int numfaces, Bpoint2 *pntlist, PindexList *faces)
{
	fathermesh.VertList.Create(numverts);
	fathermesh.FaceList.Create(numfaces);

	int vdex, pdex;

	// Load data into VertLists
	for (vdex = 0; vdex < numverts; vdex++) {
		fathermesh.VertList[vdex].loc.set(pntlist[vdex].pnt);
		fathermesh.VertList[vdex].norm.set(pntlist[vdex].norm);
		fathermesh.VertList[vdex].rad = pntlist[vdex].rad;
		fathermesh.VertList[vdex].uvt.setvals(pntlist[vdex].u, pntlist[vdex].v, pntlist[vdex].t);
 		fathermesh.VertList[vdex].v_valence = 0;		// NEED TO COMPUTE VALENCE, using as counter
		fathermesh.VertList[vdex].is_oddvert = false;

		//fathermesh.VertList[vdex].rad_2 = pntlist[vdex].rad_2;
		//fathermesh.VertList[vdex].uvt_2.setvals(pntlist[vdex].u_2, pntlist[vdex].v_2, pntlist[vdex].t_2);
	}

	// Create a map for the edges
	EdgeSet Eset;
	// Create another for the vertices, key is vertex number
	//    edgenum gives order of vertex in face
	EdgeSet Vset;

	// Load data into FaceLists---compute VertList valences at same time
	//  First, allocate space for all faces
#ifdef QIONG_DEBUG_02232003
	int triangleCounter=0;
	for (pdex = 0; pdex < numfaces; pdex++) 
		if(faces[pdex].Length()==3)
			triangleCounter++;
	fprintf(stderr, "numverts=%d\t numfaces=%d\t triangleCounter=%d\n", numverts, numfaces, triangleCounter);
#endif
	for (pdex = 0; pdex < numfaces; pdex++) 
	{
		int nsides = faces[pdex].Length();
		// Create FaceList fields
		fathermesh.FaceList[pdex].Pvertices.Create(nsides);
		fathermesh.FaceList[pdex].Pneighbors.Create(nsides);
		fathermesh.FaceList[pdex].Pneighboredges.Create(nsides);
		fathermesh.FaceList[pdex].parentPolygon = NULLPOLYDEX;
		fathermesh.FaceList[pdex].newfacevert = NULLVERTDEX;
	}

	//  Next, load data
	for (pdex = 0; pdex < numfaces; pdex++) {

		int nsides = faces[pdex].Length();

		for (int sidecnt = 0; sidecnt < nsides; sidecnt++) {

			// Set CCVindex of FaceList to index in faces[]
			CCVindex thisVertdex = faces[pdex][sidecnt];

#ifdef QIONG_DEBUG_20030510
	if(thisVertdex==59)
		bool trap=true;
#endif

			fathermesh.FaceList[pdex].Pvertices[sidecnt] = thisVertdex;

			// Increment valence counter of VertList[thisVertdex]
			fathermesh.VertList[thisVertdex].v_valence += 1;

			// if vertex hasn't been added to Vset, add polygon # and vindex in poly
			//    to Vset.  Just use one.  FASTER if rewrite PSEUDOSET for this.
			//    Can make O(n) instead of O(n logn)
			if (Vset.find(thisVertdex) == Vset.end()) {
				edge_element v_el;
				v_el.facenum = pdex;
				v_el.edgenum = sidecnt;
				Vset[thisVertdex] = v_el;
			}

			// FIGURE HOW TO DO THIS:  add this edge to set along with
			//    name of Face and location in Face
			//    If edge already in set, include Face and location in Face
			//    of this polygon as well.
			// WHEN ALL DONE, go through all edges and add Pneighbor and
			//    Pneighborfields to each edge in set

			// OKAY:  for edge vb-ve, where vb and ve are vertex indices
			//    let v1 = vb > ve ? vb : ve
			//        v2 = vb > ve ? ve : vb
			// Let key = v1 * numverts + v2
			// To retrieve, v1 = key/numverts (with int divide)
			//              v2 = key % numverts
			CCVindex nextvertdex = faces[pdex][(sidecnt + 1) % nsides];
			int v1 = thisVertdex > nextvertdex ? thisVertdex : nextvertdex;
			int v2 = thisVertdex > nextvertdex ? nextvertdex : thisVertdex;
			int thiskey = v1*numverts + v2;

			// If key already present, copy values to this FaceList and
			//    and those of key in Eset.
			// Else add to Eset
			if (Eset.find(thiskey) != Eset.end()) {
				CCPindex neighface = Eset[thiskey].facenum;
				int neighedgenum = Eset[thiskey].edgenum;
				fathermesh.FaceList[pdex].Pneighbors[sidecnt] = neighface;
				fathermesh.FaceList[pdex].Pneighboredges[sidecnt] = neighedgenum;

				fathermesh.FaceList[neighface].Pneighbors[neighedgenum] = pdex;
				fathermesh.FaceList[neighface].Pneighboredges[neighedgenum] = sidecnt;
			}
			else {
				edge_element e_el;
				e_el.facenum = pdex;
				e_el.edgenum = sidecnt;
				Eset[thiskey] = e_el;
			}
		}
	}

#ifdef QIONG_DEBUG_20030510
	//for(vdex=0; vdex<numverts; vdex++)
	//{
	//	fprintf(stderr, "%d=f%d,e%d\t", vdex, Vset[vdex].facenum, Vset[vdex].edgenum);
	//}
	//fprintf(stderr, "\n");
#endif

	// Set numedges from distinct edges in Eset
	fathermesh.nedges = Eset.size();

	// Now, compute VertFaceNeighbors and VertFaceNeighborsVert
	fathermesh.VertFaceNeighbors.Create(numverts);
	fathermesh.VertFaceNeighborsVert.Create(numverts);

	for (vdex = 0; vdex < numverts; vdex++) {
		int thisvalence = fathermesh.VertList[vdex].v_valence;

		fathermesh.VertFaceNeighbors[vdex].Create(thisvalence);
		fathermesh.VertFaceNeighborsVert[vdex].Create(thisvalence);

		CCPindex firstpolydex = Vset[vdex].facenum;
		int firstvertpos = Vset[vdex].edgenum;

		fathermesh.VertFaceNeighbors[vdex][0] = firstpolydex;
		fathermesh.VertFaceNeighborsVert[vdex][0] = firstvertpos;

		// Add neighbor polygon of (firstvertpos - 1 % numedge) edge
		//    as next VertFaceNeighbors[] and compute position of
		//    vertex in it as next VertFaceNeighborsVert[]
		CCPindex nextpolydex;
		int nextvertpos;
		for (int pcnt = 1; pcnt < thisvalence; pcnt++) 
		{
			int npedges = fathermesh.FaceList[firstpolydex].numsides();

			// neighbor polygon is along (firstvertpos - 1 % numedges) edge
			// position of vertex in neighbor is simply the number of the edge in neighbor
			nextpolydex = fathermesh.FaceList[firstpolydex].Pneighbors[INTMOD(firstvertpos - 1, npedges)];
			nextvertpos = fathermesh.FaceList[firstpolydex].Pneighboredges[INTMOD(firstvertpos - 1, npedges)];
			fathermesh.VertFaceNeighbors[vdex][pcnt] = nextpolydex;
			fathermesh.VertFaceNeighborsVert[vdex][pcnt] = nextvertpos;
			firstpolydex = nextpolydex;
			firstvertpos = nextvertpos;
		}
	}

	/*
	glPushAttrib(GL_CURRENT_BIT);
	glColor3d(1.0, 1.0, 0.0);
	for (int polydex = 0; polydex < fathermesh.numfaces(); polydex++) {

		CCPolygon *thispoly = fathermesh.FaceList(polydex);

		glBegin(GL_LINE_LOOP);
		for (int vertdex = 0; vertdex < thispoly->numsides(); vertdex++) {
			CCVertex *thisVert = fathermesh.VertList(thispoly->Pvertices[vertdex]);
			thisVert->norm.glNormal();
			thisVert->loc.glVertex();
		}
		glEnd();
	}
	glPopAttrib();
	*/
}

bool CCSubdivsurf::computeLvl0VertAtomList(int numRows, int numCols, int * list)
{
    int row, col, vertDex;
    CCVertex * thisVert;

    CCMesh * thisMesh = submeshes(0);
    int vertNum = thisMesh->numverts();

    if (list == NULL)
        list = new int[vertNum];

    for (vertDex = 0; vertDex<vertNum; vertDex++) {
        thisVert = thisMesh->VertList(vertDex);

        col = (int) thisVert->uvt.u;
        row = (int) thisVert->uvt.v;

        if (row >= numRows || col >= numCols) {
//			printf("Returning false\n");
			return false;
        }

      //  list[vertDex] = col * numRows + row;
		list[vertDex] = row * numCols + col;
	//	printf("%d: %lf %lf, %d\n", vertDex, thisVert->uvt.u, thisVert->uvt.v, list[vertDex]);


    }
    return true;
}


