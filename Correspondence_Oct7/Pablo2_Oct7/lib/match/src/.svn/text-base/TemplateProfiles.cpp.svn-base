#include <math.h>
#include <iostream>
#include "M3DObject.h"
#include "Image3D.h"
#include "TemplateProfiles.h"
#include "DistanceMap3D.h"


using namespace std;

TemplateProfiles::TemplateProfiles()
{
	numTemplateTypes = 0;
	ProfileTypes = NULL;
}

TemplateProfiles::~TemplateProfiles()
{
	delete [] ProfileTypes;
}

void TemplateProfiles::loadTPdata(M3DObject *referenceObject, int _figureId, float _cutoff,
		const char * TPfilename, int surfaceLevel)
{
	FILE *filterfile;
	double tempd;
	int curtype, numtypes, numpoints, numsamples, i,j;

	if (TPfilename != NULL) {

		filterfile = fopen(TPfilename, "r");

		if (filterfile == NULL) {
			cout << "TemplateProfiles.cpp: template file not found.\n";
			TPfilename = NULL;
		}
		else {
			// Read in the info.
			// First the header, describing how many template types, the length
			// of the templates. 

			fscanf(filterfile, "%i %i %i ", &numpoints, &numtypes, &numsamples);

			//Now read in the templates.
			dimension = numsamples;
			numTemplateTypes = numtypes;
			numProfiles = numpoints;

			ProfileTypes = new int[numProfiles];

			fscanf(filterfile, "%lf %lf %lf ", &CenterIntensity, &ProfileMatchVariance, &MeanMatchVariance);

			//Now read the mean offset and std.
			for (i = 0; i < numpoints; i ++)
			{
				fscanf(filterfile, "%lf", &tempd);
				MeanOffset[i] = tempd;
			}

			for (i = 0; i < numpoints; i ++)
			{
				fscanf(filterfile, "%lf", &tempd);
				MeanStd[i] = tempd;
			}

			//Now read the rmses.
			for (i = 0; i < numpoints; i ++)
			{
				fscanf(filterfile, "%lf", &tempd);
				Rmses[i] = tempd;
			}


			//Now read the templates.
			for (i = 0; i < numtypes; i ++)
			{
				for (j = 0; j < numsamples; j ++)
				{
					fscanf(filterfile, "%lf", &tempd);
					Templates[i][j] = tempd;
				}
			}

			//Now for each point, read in the type.

			for (i = 0; i < numpoints; i ++) 
			{
				fscanf(filterfile, "%i", &curtype);
				ProfileTypes[i] = curtype;
			}


			fclose(filterfile);
		}
	}
}

void TemplateProfiles::NormalizeTemplates()
//Mean zero and magnitude 1 all the templates.
{
	//int i,j;
	//compute the rmses. 
	//double *mags = new double[numTemplateTypes];
}

double TemplateProfiles::getTemplateVal(int i, int j)
{
	if (i >= numProfiles)
		return Templates[ProfileTypes[numProfiles - 1]][j];
	else
		return Templates[ProfileTypes[i]][j];
}


