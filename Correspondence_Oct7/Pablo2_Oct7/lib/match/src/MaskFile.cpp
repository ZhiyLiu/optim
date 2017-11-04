#include <iostream>
#include <fstream>

#include "MaskFile.h"


using namespace std;

bool MaskFile::read(const char * filename, Match * match)
{
    M3DObject * object;
    Mask * mask;
    int figureCount,
        testFigureCount,
        size;
    int i, j, k;

    double cutoff;
    double sigma;
    bool includeInterior;
    MaskElement * newElement;
    std::vector<MaskElement *> elements;

    std::vector<referencePoint> ** constrainedPointLists;

    ifstream inFile;


    if (filename == NULL)
        return false;

    object = match->getReferenceObject();
    if (object == NULL)
    {
        cout << "MaskFile::read() : Match not initialized." << endl;
        return false;
    }

    inFile.open(filename, ios::in | ios::binary);
    if (! inFile) {
		cout << "Could not open mask file\n";
        return false;
	}

    figureCount = object->getFigureCount();
    inFile.read((char *) (&testFigureCount), sizeof(int));

    if (figureCount != testFigureCount)
    {
        cout << "MaskFile::read() : Number of figures does not match." << endl;
        return false;
    }

    constrainedPointLists = new std::vector<referencePoint> * [figureCount];
    for (i = 0; i < figureCount; i++)
        constrainedPointLists[i] = NULL;

    for (i = 0; i < figureCount; i++)
    {
        mask = match->getMask(i);
        if (mask == NULL)
            continue;

        inFile.read((char *) (&cutoff), sizeof(double));
        inFile.read((char *) (&sigma), sizeof(double));
        inFile.read((char *) (&includeInterior), sizeof(bool));
        inFile.read((char *) (&size), sizeof(int));

        elements.erase(elements.begin(), elements.end());
        for (j = 0; j < size; j++)
        {
            newElement = new MaskElement;

            inFile.read((char *) (newElement), sizeof(MaskElement));
            elements.push_back(newElement);
        }

        mask->setVals(i, cutoff, sigma, includeInterior, elements);

        constrainedPointLists[i] = new std::vector<referencePoint> [figureCount];

        for (j = 0; j < figureCount; j++)
        {
            inFile.read((char *) (&size), sizeof(int));

            for (k = 0; k < size; k++)
            {
                referencePoint rp;
                inFile.read((char *) (&rp), sizeof(referencePoint));

                constrainedPointLists[i][j].push_back(rp);
            }
        }
    }

    match->setConstrainedPointLists(constrainedPointLists);

    inFile.close();
	return true;
}

void MaskFile::write(const char * filename, Match * match)
{
    M3DObject * object;
    M3DFigure * figure;
    Mask * mask;
    MaskElement * element;
    int figureCount,
        size;
    int i, j, k;
	std::ofstream outfile;

    if(filename == NULL)
        return;

    object = match->getReferenceObject();
    if(object == NULL)
    {
        cout << "MaskFile::write() : Match not initialized." << endl;
        return;
    }

    outfile.open(filename, ios::out | ios::binary);

    figureCount = object->getFigureCount();
    outfile.write((char *) (&figureCount), sizeof(int));

    std::vector<referencePoint> ** constrainedPointLists = match->getConstrainedPointLists();

    for(i = 0; i < figureCount; i++)
    {
        mask = match->getMask(i);
        if(mask == NULL)
            continue;
        size = mask->getSize();

        double cutoff = mask->getCutoff();
        double sigma = mask->getSigma();
        bool interiorIsIncluded = mask->interiorIsIncluded();

        outfile.write((char *) (&cutoff), sizeof(double));
        outfile.write((char *) (&sigma), sizeof(double));
        outfile.write((char *) (&interiorIsIncluded), sizeof(bool));
        outfile.write((char *) (&size), sizeof(int));

        for(j = 0; j < size; j++)
        {
            element = mask->getElement(j);
            if(element == NULL)
                continue;

            outfile.write((char *) (element), sizeof(MaskElement));
        }

        // Output interfigural constraints
        figure = object->getFigurePtr(i);
        InterfiguralConstraints & governors = figure->inverseConstraints();

        if(constrainedPointLists[i] == NULL)
        {
            cout << "UH-OH" << endl;
            continue;
        }

        for(j = 0; j < figureCount; j++)
        {
            size = constrainedPointLists[i][j].size();

            outfile.write((char *) (&size), sizeof(int));

            for(k = 0; k < size; k++)
            {
                referencePoint * rp = &(constrainedPointLists[i][j][k]);

                outfile.write((char *) (rp), sizeof(referencePoint));
            }
        }
    }

    outfile.close();
}

