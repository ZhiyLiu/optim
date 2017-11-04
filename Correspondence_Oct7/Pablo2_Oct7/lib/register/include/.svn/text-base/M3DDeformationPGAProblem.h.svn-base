#ifndef M3D_DEFORMATION_PGA_PROBLEM_H
#define M3D_DEFORMATION_PGA_PROBLEM_H



class M3DPGAPrimitiveStats;

class M3DDeformationPGAProblem : public M3DDeformationProblem  
{

public:

	M3DDeformationPGAProblem(Match * _match, M3DObject * _referenceObject,
		M3DObject * _candidateObject, M3DPGAPrimitiveStats * _primPga,
		int _figId, int _atomId, int _applicationID);

	M3DDeformationPGAProblem(M3DPGAPrimitiveStats * pgaStats) {
		primPga = pgaStats;
	};

	virtual ~M3DDeformationPGAProblem() { };
    void setPGA(M3DPGAPrimitiveStats * primPgaStats);	
	void applyVector(M3DPrimitive & prim, const Vector & x, int order);

	double evaluate(const Vector & x);
	M3DPGAPrimitiveStats * primPga;
};


#endif

