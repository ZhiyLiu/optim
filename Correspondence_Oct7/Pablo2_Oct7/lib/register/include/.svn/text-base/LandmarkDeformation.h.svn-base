#ifndef __LANDMARK_DEFORMATION_H__
#define __LANDMARK_DEFORMATION_H__


#include "matrix.h"


class LandmarkDeformation
{
	public:

		LandmarkDeformation(Matrix & m1, Matrix & m2, int rows);
		void printResults();

		double getBeta() { return beta; }
		Matrix getGamma() { return gamma; }
		Vector getOffset() { return final_offset; }

	private:

		int n;
		double beta;

		Matrix C_Matrix;

		Matrix PreCentered_X1;
		Matrix PreCentered_X2;
		Matrix X1;
		Matrix X2;
		Matrix PreCentered_X1MinusX1;
		Matrix PreCentered_X2MinusX2;
		Vector X2_offset;
		Vector X1_offset;
		Vector final_offset;

		Matrix gamma;
		double X1_Euclidean_Norm;
		double X2_Euclidean_Norm;

		void compute_C_matrix();
		void run();
};


#endif

