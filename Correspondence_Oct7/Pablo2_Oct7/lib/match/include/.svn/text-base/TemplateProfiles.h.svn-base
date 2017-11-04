#ifndef TEMPLATEPROFILES_H
#define TEMPLATEPROFILES_H



#define MAX_PROFILE_DIM 20
#define MAX_NUM_TEMPLATES 2562
#define MAX_NUM_POINTS 10000


class M3DObject;

struct MTElement //This stuff is for the future.  
{
	//double x[MAX_PROFILE_DIM],y[MAX_PROFILE_DIM],z[MAX_PROFILE_DIM];
	//double u,v,t, tau[MAX_PROFILE_DIM];
	double area, radius, volume, normconst;
	bool inBlendRegion;

	double mixture[MAX_NUM_TEMPLATES];
};

class TemplateProfiles
//Currently (10/03) this class simply contains an array of templates and
//and an array of template index per boundary point.
{
public:

	TemplateProfiles();
	~TemplateProfiles();

	void loadTPdata(M3DObject *referenceObject, int _figureId, float _cutoff,
		const char * TPfilename, int surfaceLevel);

	void NormalizeTemplates();

	int getNumTemplateTypes() { return numTemplateTypes; }
	int getDimension() {return dimension; }
	int getSize() {return MTElementsSize;}
	double getTotalVolume() { return totalVolume; }
	MTElement * getMTElement(int i) {return &MTElements[i];}

	//int GetTemplateType(int i) {return ProfileTypes[i];}
	double getTemplateVal(int i, int j);
	double getTempVal(int i, int j) { return Templates[i][j];}
	double getRms(int i) {return Rmses[i];}
	double getMeanOffset(int i) {return MeanOffset[i];}
	double getMeanStd(int i) {return MeanStd[i];}

	int getProfileType(int i) {return ProfileTypes[i];}

	double getPMatchVar() {return ProfileMatchVariance;}
	double getMMatchVar() {return MeanMatchVariance;}
	double getCenterIntensity() {return CenterIntensity;}


private:
	MTElement * MTElements;

	double Templates[MAX_NUM_TEMPLATES][MAX_PROFILE_DIM];
	double Rmses[MAX_NUM_POINTS];
	double MeanOffset[MAX_NUM_POINTS];
	double MeanStd[MAX_NUM_POINTS];
	int *ProfileTypes;

	int MTElementsSize;

	int figureId;

	int numProfiles;
	int numTemplateTypes;
	int dimension;
	double totalVolume, cutoff;

	double ProfileMatchVariance, MeanMatchVariance, CenterIntensity;
};


#endif

