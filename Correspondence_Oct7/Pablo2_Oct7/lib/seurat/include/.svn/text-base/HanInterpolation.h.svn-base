#ifndef HANINTERPOLATION_H
#define HANINTERPOLATION_H

#include "M3DQuadFigure.h"

class HanInterpolation{

public:
	HanInterpolation();
	~HanInterpolation();
	

	void SetM3DQuadFigure(M3DQuadFigure* quadfig){
		m_M3DQuadFigure = quadfig;
	}

	void SetInterpolationLevel(int level){
		m_InterpolationLevel = level;
	}

	/* 
	*	Contains the hub possition of the spoke
	*/
	vector< vector< vector< Vector3D > > > GetInterpolatedPoints(){
		return m_InterpolatedPoints;
	}

	/* 
	*	Contains the U*r direction of the spoke, see GetInterpolatedPoints for the
	*	possition at the medial sheet
	*/
	vector< vector< vector< Vector3D > > > GetInterpolatedSpokes(){
		return m_InterpolatedSpokes;
	}

	

	/*
	*	\fn void Update();
	*	\brief Brings the interpolation up to date
	*		   Set the interpolation level and the quadfigure to interpolate
	*	
	*/
	void Update();

	/*void SetSubMeshes(MyList<CCMesh> submeshes){
		m_Submeshes = submeshes;
	}*/

private:

	M3DQuadFigure* m_M3DQuadFigure;

	int m_InterpolationLevel;

	//to access->  m_InterpolatedSpokes[side][posx][posy]       side [0-1] 
	vector< vector< vector< Vector3D > > > m_InterpolatedSpokes;
	vector< vector< vector< Vector3D > > > m_InterpolatedPoints;

	/*
	*	\brief Contains all the calculation for rSrad Matrices for each spoke and
	*		   rSe if is an edge spoke.
	*		   To access the element m_rSMatrix[i][j][0] -> side0 rSrad
	*		   To access the element m_rSMatrix[i][j][1] -> side1 rSrad
	*		   To access the element m_rSMatrix[i][j][2] -> edge  rSe if edge spoke
	*	
	*/
	vector< vector< vector< Matrix > > > m_rSradMatrix;

	vector< vector< vector< Vector3D > > > m_Derivatives;

	/*
	*   \fn Vector3D* getHermiteMatrix(Vector3D p0, Vector3D p1, Vector3D p2, Vector3D p3, Vector3D n0, Vector3D n1, Vector3D n2, Vector3D n3)
	*	\brief Calculates the corresponding hermiteMatrix from the 4 control points and the 4 given normals
	*	\param Vector3D p0-p1 the control points
	*	\param Vector3D n0-n1 the normals of the control points
	*/
	void GetHermiteMatrix(Vector3D p0, Vector3D p1, Vector3D p2, Vector3D p3, vector< Vector3D > dp0, vector< Vector3D > dp1, vector< Vector3D > dp2, vector< Vector3D > dp3, Vector3D n0, Vector3D n1, Vector3D n2, Vector3D n3, Vector3D H[4][4]);

	/*
	* \fn void GetInterpolatedPoint(double u, double v, Vector3D H[4][4], Vector3D p)
	* \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
	* \param double u scalar for the weight functions
	* \param double v second scalar for the weight functions
	* \param Vector3D H[4][4] cubic hermite matrix
	* \param Vector3D p result of the interpolation
	* \pre u and v should be between [0, 1]
	*/
	void GetInterpolatedPoint(double u, double v, Vector3D H[4][4], Vector3D& p);

	/*
	*	\fn Get the interpolated points
	*/
	void GetInterpolatedEdgePoints(Matrix rSe0, Matrix rSe1, Matrix rSe2, Matrix rSe3);
	

	void GetRSEMatrix(M3DPrimitive* prim0, M3DPrimitive* prim1, Vector3D du, Matrix& rSe);

	/*
	*	\fn GetRSEMatrix
	*	\brief Computes the rSe matrix for the given values
	*	\param Vector3D x0, the position of the primitive at (0,0)
	*	\param Vector3D x1, the next position of the primitive to the edge
	*	\param Vector3D U0	the end spoke0
	*	\param Vector3D N0  the normal to the medial sheet at p0
	*	\param Vector3D U1	the end spoke1
	*	\param Vector3D N1  the normal to the medial sheet at p1
	*	\param double r0 the radius of the end spoke0 
	*	\param double r1 the radius of the end spoke1
	*	\param double rn0 the radius of the normal spoke to the medial sheet at spoke0
	*	\param double rn1 the radius of the normal spoke to the medial sheet at spoke1
	*	
	*	\param Vector3D N	normal to the medial sheet
	*	\param Matrix& rSe	matrix to store the calculated values
	*/
	void GetRSEMatrix(Vector3D x0, Vector3D x1, Vector3D U0, Vector3D N0, Vector3D U1, Vector3D N1, double r0, double rn0, double r1, double rn1, Matrix& rSe);

	void GetRSEMatrix(Vector3D du, Vector3D dndu, Vector3D dUdu, Vector3D U0, double r0, Matrix& rSe);

	void GetdEdu(Vector3D x0, Vector3D x1, Vector3D n, Vector3D U0, Matrix rSe, Matrix& dEdun);

	/**
	*	\fn Fits a cuadratic curve between the two given vectors
	*	\pre	The two vectors share the same starting point
	*	\param  Vector3D r0 starting vector to interpolate, note that the derivative at this point is 0
	*	\param  Vector3D rmax end vector, the derivative at this point comes from sRad interpolation
	*	\post	A list of points according to the interpolation level is returned
	*	\return vector<Vector3D> interpolated points
	*/
	vector<Matrix> FitEdgeCurve(Vector3D r0, Vector3D rmax, double drmax);


	//void GetRSRadMatrix(M3DPrimitive* prim0, M3DPrimitive* prim1, M3DPrimitive* prim2, Matrix& rSrad, bool side = false);
	//void GetRSRadMatrix(Vector3D x0, Vector3D x1, Vector3D x2, double r0, double r1, double r2, Vector3D U0, Vector3D U1, Vector3D U2, Matrix& rSrad);
	void GetRSRadMatrix(Vector3D pu1, Vector3D pu2, Vector3D dUdu1, Vector3D dUdu2, Vector3D U0, double r0, Matrix& rSrad);
	void GetdSdu(Vector3D x0, Vector3D x1, Vector3D x2, double stepu, double stepv, Vector3D U0, Matrix rSrad, Matrix& dSdu);

	/**
	*	\fn double GetLogLambda(double lam)
	*	\brief evaluates if the value lam < 1 and returns a value
	*	\post if the value is greater than 1, returns a large negative value
	*/
	double GetLogLambda(double lam);


	/*
	*	\brief Computes the discrete derivatives of the point at i, j
	*	\post The values du1 and du2 have the values of the derivatives in the u1 and u2 direction
	*	the optional parameters prim1 and prim2 contain the primitives used in the quad
	*/
	void GetPDerivatives(int i, int j, Vector3D& du1, Vector3D& du2, M3DPrimitive** outprim1, M3DPrimitive** outprim2);

	/*
	*/
	void GetEdgeDerivatives(int i, int j, Vector3D& du, Vector3D& dndu, Vector3D& dUdu);

	void CalculateSradSeAndDerivatives();
	void InterpolatePoints();
	void InterpolateMedial(int i, int j, M3DPrimitive* prim0, Vector3D Hm[4][4]);
	void InterpolateCrest(int i,int j, int inext, int jnext, Vector3D Hm[4][4]);
	void ClearVectors();

};

#endif
