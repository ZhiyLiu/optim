#include "SimpleMaskFile.h"
#include <iostream>
#include <stdio.h>
#include <fstream>


using namespace std;


// Reads the text version of a smsk file.
void SimpleMaskFile::read(const char * filename)
{
	ifstream fin;

	fin.open(filename);

	if (fin == NULL) {
		cout << "Error SimpleMaskFile: could not open " << filename << ".\n";
		return;
	}

	fin >> numFigures;

	//Now initialize some member variables accordingly.
	values = new double * [numFigures];
	numPoints = new int [numFigures];
	numSamples = new int [numFigures];
	figIds = new int [numFigures];

	figNames = new char * [numFigures];

	int numels, surflevel, rowcount, colcount, i, j, k, vind;

	double tempd;
	char *c = new char[10];

	//Now for each figure, read in the data.
	for (i = 0; i < numFigures; i++)
	{
		figNames[i] = new char[256];

		// Current id and name.
		fin >> figIds[i];
		fin.getline(figNames[i], 256);

		// numpoints, numelements, numsamples, surflevel, rowcount, colcount.
		fin >> numPoints[i] >> numels >> numSamples[i] >> surflevel >> rowcount >> colcount;
		char ignoredChar;
		fin >> ignoredChar;
		// The following doesn't link -
		//fin.ignore();

		// The rest of the data is binary, the numels data doubles per point, with
		// numSamples[i] intensity samples that matter for the moment.
		values[i] = new double [numSamples[i]*numPoints[i]];

		// First, read all the data, then write the relevant part to values.
		for (j = 0; j < numPoints[i]; j ++)
		{
			vind = j*numSamples[i];
			for (k = 0; k < numSamples[i]; k ++, vind ++)
			{
				fin >> tempd;
				values[i][vind] = tempd;
			}
			fin >> tempd;
		}
	}
}

// Read in the information from a smsk file that is binary.
// I'm having trouble alternating between reading text and 
// binary data.  Thus the above read.
void SimpleMaskFile::readBinary(const char *filename)
{
	FILE *fin = NULL;
	fin = fopen(filename, "rb");

	if (fin == NULL)
	{
		cout << "Error SimpleMaskFile: could not open " << filename << ".\n";
		return;
	}

	fscanf(fin, "%d", &numFigures);

	// Now initialize some member variables accordingly.
	values = new double * [numFigures];
	numPoints = new int [numFigures];
	numSamples = new int [numFigures];
	figIds = new int [numFigures];
	taus = new double [2*numFigures];
	figNames = new char * [numFigures];

	int numels, surflevel, i, j, k, dind, vind;

	double * data;

	// Now for each figure, read in the data.
	for (i = 0; i < numFigures; i ++)
	{
		figNames[i] = new char[256];

		//current id and name.
		fscanf(fin, "%d %s;", &figIds[i], figNames[i]);

		// numpoints, numelements, numsamples, surflevel, rowcount, colcount.
		fscanf(fin, "%d %d %d %d %lf %lf", &numPoints[i], &numels, &numSamples[i],
			&surflevel, &taus[2*i], &taus[2*i+1]);
		if ((char) fgetc(fin) != ';')
			cout << "Error SimpleMaskFile: misread.\n";

		//The rest of the data is binary, the numels data doubles per point, with
		// numSamples[i] intensity samples that matter for the moment.
		values[i] = new double [numSamples[i]*numPoints[i]];
		data = new double [numels*numPoints[i]];

		// First, read all the data, then write the relevant part to values.
		fread(&(data[0]), sizeof(double), numels*numPoints[i], fin);

		// Now set values
		for (j = 0; j < numPoints[i]; j ++)
		{
			dind = j*numels;
			vind = j*numSamples[i];

			for (k = 0; k < numSamples[i]; k ++, dind ++, vind ++)
				values[i][vind] = data[dind];
		}

		// That's it for this figure...on to the next.
		delete [] data;
	}

}

// Function for getting an element of values.
bool SimpleMaskFile::getValue(int figureId, int point_ind, int sample, double & value)
{
	int ind = point_ind*numSamples[figureId] + sample;

	if (point_ind >= numPoints[figureId]) {
		cout << "SimpleMaskFile: wrong number of points\n";
		return false;
	}

	if (figureId < 0 || figureId > numFigures)
	{
		cout << "Error SimpleMaskFile getValue: figureId " << figureId 
			 << " or index " << ind << " not valid.\n";
		return false;
	}

	value = values[figureId][ind];
	return true;
}

double SimpleMaskFile::getTau(int figureId, int which)
{
	if (figureId < 0 || figureId > numFigures)
	{
		cout << "Error SimpleMaskFile getTau: figureId " << figureId 
			 << " not valid.\n";
	}

	if (which == 0)
		return taus[2*figureId];
	else
		return taus[2*figureId + 1];
}

int SimpleMaskFile::getNumSamples(int figureId)
{
	if (figureId < 0 || figureId > numFigures)
	{
		cout << "Error SimpleMaskFile getNumSamples: figureId " << figureId 
			 << " not valid.\n";
	}

	return numSamples[figureId];
}

