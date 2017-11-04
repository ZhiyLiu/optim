#ifndef M3DFIGURESTATS_H
#define M3DFIGURESTATS_H

#include <vector>

/*Class M3DFigureStats

  This class is used to store statstics information needed to compute 
  a trained image match.  The idea is that a pointer to this class is
  in M3DFigure, and its information is gathered from the m3d file 
  during read (see M3DObjectFile.readQuadFigure, which declares a
  figure to return to the m3d object, and then reads from the registry
  to fill the figure).

  The kind of training information important to computing a trained 
  image match is quite dynamic in the research.  Check such functions
  as Match::ComputeMainFigureMatch to see what the most up to date 
  method is...
*/


class M3DFigureStats
{

public:
	M3DFigureStats() { }
	~M3DFigureStats();

	//setting member structures.
	void setTemplate(int template_ind, int sample, double val); 
	void setTemplateType(int point_ind, int template_ind);
	void setRmses(int point_ind, double val);
	void setMeanOffset(int point_ind, double val);
	void setMeanStd(int point_ind, double val);

	void setnumPoints(int val) { numPoints = val; }
	void setnumTemplates(int val) { numTemplates = val; }
	void setdimension(int dims) { dimension = dims; }
	void setProfileMatchVar(double variance) { profileMatchVariance = variance; }
	void setMeanMatchVar(double variance) { meanMatchVariance = variance; }

	//accessing member structure values.
	void getTemplate(int point_ind, double v[]);
	double getTemplateVal(int point_ind, int sample);
	double getTemplateValbyType(int temp_type, int sample);
	int getTemplateType(int point_ind) { return template_types[point_ind]; }
	double getRms(int point_ind) { return rmses[point_ind]; }
	double getMeanOffset(int point_ind) { return mean_offset[point_ind]; }
	double getMeanStd(int point_ind) { return mean_std[point_ind]; }

	int getnumPoints() { return numPoints; }
	int getnumTemplates() { return numTemplates; }
	int getdimension() { return dimension; }
	double getProfMatchVar() { return profileMatchVariance; }
	double getMeanMatchVar() { return meanMatchVariance; }

protected:

	std::vector<double> templates;	//the possible template types at a point.
	std::vector<int> template_types; //indexes to type per point.
	std::vector<double> rmses;		//the trained intensity normalization per
										//per point.

	std::vector<double> mean_offset; //the trained intensity per point.
	std::vector<double> mean_std;	//the intensity variance per point.


	int numPoints;
	int numTemplates;
	int dimension;	//number of samples per template (or profile).

	double profileMatchVariance, meanMatchVariance;

private:

	// These members are not implemented
	M3DFigureStats(M3DFigureStats & fs);
	M3DFigureStats & operator=(const M3DFigureStats & fs);

};

#endif

