#include <iostream>
#include <vector>
#include "Image3D.h"


#ifdef AE2_BUILD

#include "Anastruct.h"


// Class Anastruct3D functions


Anastruct3D::Anastruct3D() {
    contours = NULL;
    name = NULL;
    contour_count = 0;
}

Anastruct3D::~Anastruct3D() {
    for (int c = 0; c < contour_count; c++) {
        delete [] contours[c].modelX;
        delete [] contours[c].modelY;
        // contours[c].x and contours[c].y are not owned by this class
    }
    delete [] contours;
}


// Class AnastructList functions


AnastructList::~AnastructList() {
    for (int a = 0; a < ana.size(); a++)
        delete ana[a];
    ana.clear();
}

bool AnastructList::convertAnastructCoords(Image3D * im)
{
    Vector3D coord;

    if (im == NULL)
        return false;

    for (int a = 0; a < ana.size(); a++) {
        if (ana[a]->contours->modelX != NULL)
            return true;

        for (int c = 0; c < ana[a]->contour_count; c++) {
            ana[a]->contours[c].modelX = new double[ana[a]->contours[c].pointCount];
            ana[a]->contours[c].modelY = new double[ana[a]->contours[c].pointCount];

            coord = Vector3D(ana[a]->contours[c].min[0], ana[a]->contours[c].min[1], ana[a]->contours[c].min[2]);
            im->worldToModelCoordinates(coord);
            ana[a]->contours[c].min[0] = coord.getX();
            ana[a]->contours[c].min[1] = coord.getY();
            ana[a]->contours[c].min[2] = coord.getZ();

            coord = Vector3D(ana[a]->contours[c].max[0], ana[a]->contours[c].max[1], ana[a]->contours[c].max[2]);
            im->worldToModelCoordinates(coord);
            ana[a]->contours[c].max[0] = coord.getX();
            ana[a]->contours[c].max[1] = coord.getY();
            ana[a]->contours[c].max[2] = coord.getZ();

            coord = Vector3D(ana[a]->contours[c].x[0], ana[a]->contours[c].y[0], ana[a]->contours[c].z);
            im->worldToModelCoordinates(coord);
            ana[a]->contours[c].modelX[0] = coord.getX();
            ana[a]->contours[c].modelY[0] = coord.getY();
            ana[a]->contours[c].modelZ = coord.getZ();

            for (int i = 1; i < ana[a]->contours[c].pointCount; i++) {
                coord = Vector3D(ana[a]->contours[c].x[i], ana[a]->contours[c].y[i], ana[a]->contours[c].z);
                im->worldToModelCoordinates(coord);
                ana[a]->contours[c].modelX[i] = coord.getX();
                ana[a]->contours[c].modelY[i] = coord.getY();
            }

            // Convert anastruct colors to GL colors
            switch (ana[a]->color) {
	        case 0: ana[a]->glColor[0] = 1.0f;
                        ana[a]->glColor[1] = 0.0f;
                        ana[a]->glColor[2] = 0.0f;
                        break;
	        case 1: ana[a]->glColor[0] = 0.0f;
                        ana[a]->glColor[1] = 1.0f;
                        ana[a]->glColor[2] = 0.0f;
                        break;
	        case 2: ana[a]->glColor[0] = 0.0f;
                        ana[a]->glColor[1] = 0.0f;
                        ana[a]->glColor[2] = 1.0f;
                        break;
	        case 3: ana[a]->glColor[0] = 1.0f;
                        ana[a]->glColor[1] = 1.0f;
                        ana[a]->glColor[2] = 0.0f;
                        break;
	        case 4: ana[a]->glColor[0] = 1.0f;
                        ana[a]->glColor[1] = 0.0f;
                        ana[a]->glColor[2] = 1.0f;
                        break;
	        case 5: ana[a]->glColor[0] = 0.0f;
                        ana[a]->glColor[1] = 1.0f;
                        ana[a]->glColor[2] = 1.0f;
                        break;
	        case 6: ana[a]->glColor[0] = 1.0f;
                        ana[a]->glColor[1] = 1.0f;
                        ana[a]->glColor[2] = 1.0f;
                        break;
	        case 7: ana[a]->glColor[0] = 0.0f;
                        ana[a]->glColor[1] = 0.0f;
                        ana[a]->glColor[2] = 0.0f;
                        break;
            }
        }
    }
    return true;
}

void AnastructList::addAnastruct(Anastruct3D * a)
{
    ana.push_back(a);
}

void AnastructList::deleteAnastruct(int index)
{
    if (index >= ana.size() || index < 0)
	return;
    if (index == ana.size() - 1)
        ana.pop_back();
    else {
        Anastruct3D * a = ana[index];
        ana[index] = ana[ana.size() - 1];
        ana[ana.size() - 1] = a;
        ana.pop_back();
    }
}

Anastruct3D * AnastructList::getAnastruct(int index) const
{
    if (index >= ana.size() || index < 0)
	return NULL;
    return ana[index];
}

const char * AnastructList::getName(int index) const
{
    if (index >= ana.size() || index < 0)
	return NULL;
    return ana[index]->name;
}

int AnastructList::contourCount(int index) const
{
    if (index >= ana.size() || index < 0)
	return 0;
    return ana[index]->contour_count;
}

int AnastructList::pointCount(int index, int contourId) const
{
    if (index >= ana.size() || index < 0)
	return 0;
    return ana[index]->pointCount(contourId);
}

int AnastructList::getColor(int index) const
{
    if (index >= ana.size() || index < 0)
	return 0;
    return ana[index]->getColor();
}


#endif

