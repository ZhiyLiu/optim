/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
//   File----landmark.cpp
//   Author----John Glotzer
//   Date------Nov. 3, 1998
//   Function-----Does a procrustes fit to two sets of landmarks. The landmarks
//     are assumed to be in 3D. The cardinality of both sets must be the same.
//     Homology is a necessary precondition.
//     Uses the LaPack numerical library along with a wrapper class built around
//     the library. The wrapper classes are Matrix, and Vector and they
//     are written by Paul Yuskevich. The method is taken from Dryden and Mardia's
//     book and shape analysis and theory, and is called by them 
//
//				Ordinary (Procrustes) Sum of Squares

/*
		The method works as follows. From each of our n sets of points we form two
		n x 3 matrices. We also compute the centering matrix. The centering
		matrix is an n x n matrix.  The difference of the centered matrix and 
		the initial matrix is the optimum translation. We operate on "centered"
		sets of points or matrices.

		The assumption is that we want the matrix that transforms X1 into X2.
		Put another way 

				((beta * X1) * gamma) + Offset 

			should be as close to X2 as possible. 

		Next, the Matrix X2-transpose * X1 is calculated and this matrix is divided
		by the product of the Euclidean norms of X1 and X2.

					i.e.     X2T * X1
							__________

							||X1||*||X2||  (all matrices centered)

		The Euclidean norm of a matrix is just the square root of the trace of the product
		of a matrix with its transpose. In other words it is the square root of the sum
		of the magnitudes of the vectors. 


		Once this matrix is calculated, we run SVD on it to produce matrices
		V*Lambda*Utrans.

		The optimal rotation matrix, gamma, is just U*Vtrans. (note we take trans of U and V)

		The optimal scale (to be applied to X1), is given by 


						trace (X2trans * X1 * gamma)
						____________________________
						trace (X1trans * X1)

		As mentioned, we already have the optimal translation. 

Basically, according to Dryden and Mardia, we can get a close form solution if

	- we operate in 2 dimensions and have N configuration matrices
	- we operate in m dimensions and have 2 configuration matrices

  If we operate in m dimensions and have N configruation matrices, we require
  a numerical solution. 
*/


#include <iostream>
#include <stdlib.h>
#include "LandmarkDeformation.h"
#include "matrix.h"


using namespace std;
using namespace pauly;


//#define DEBUG

LandmarkDeformation::LandmarkDeformation(Matrix & m1, Matrix & m2, int rows)
    : n(rows), X1(n,3), X2(n,3)
{
	compute_C_matrix();

	PreCentered_X1 = m1;
#ifdef DEBUG
	cout << "This is the Precentered X1 matrix." << endl;
	PreCentered_X1.print();
#endif

	X1 = C_Matrix * m1;
#ifdef DEBUG
	cout << "And this is the Centered X1." << endl;
	X1.print();
#endif

	PreCentered_X2 = m2;
#ifdef DEBUG
	cout << "This is the PreCentered X2 matrix." << endl;
	PreCentered_X2.print();
#endif

	X2 = C_Matrix * m2;
#ifdef DEBUG
	cout << "And this is the Centered X2. " << endl;
	X2.print();
#endif

	PreCentered_X1MinusX1 = PreCentered_X1 - X1;
    PreCentered_X2MinusX2 = PreCentered_X2 - X2;
	X1_Euclidean_Norm = X1.twoNorm();
	X2_Euclidean_Norm = X2.twoNorm();

	run();
}

void LandmarkDeformation::run()
{
	Matrix X2T, X1T;
	Matrix X2TX1;
	Matrix SVDMatrix;

	X2T = X2.t();
	X1T = X1.t();
	X2TX1 = X2T * X1;

#ifdef DEBUG
	cout << "The norm is " << X2_Euclidean_Norm << endl;
#endif

	SVDMatrix = X2TX1 / (X1_Euclidean_Norm * X2_Euclidean_Norm);

	Matrix V; 
	Vector Sigma; 
	Matrix Utrans;
	Matrix U;
	Matrix Vtrans;

	SVDMatrix.factorSVD(V,Utrans,Sigma,0);
	U = Utrans.t();
	Vtrans = V.t();

	(gamma) = (U) * (Vtrans);

#ifdef DEBUG
	cout << "The value of gamma is " << endl;
	gamma.print();
#endif

	Matrix numerator = (X2TX1 * gamma);
	Matrix denom = (X1T * X1);

	beta = numerator.trace()/ denom.trace();
#ifdef DEBUG
	cout << "The value of beta is " << beta << endl;
#endif

	X1_offset = PreCentered_X1MinusX1.getRow(0);
	X2_offset = PreCentered_X2MinusX2.getRow(0);
	final_offset = X2_offset - X1_offset;
}

void LandmarkDeformation::printResults()
{
	cout << "The rotation matrix is " << endl;
	gamma.print();
	cout << "The value of beta is " << beta << endl;
	cout << "The value of PrecenteredX1 -  X1 is " << endl;
	PreCentered_X1MinusX1.print();
	cout << "The value of PreCentered X2 -  X2 is " << endl;
	PreCentered_X2MinusX2.print();
	cout << "The value of the final translation vector is " << endl;
	final_offset.print();
}

void LandmarkDeformation::compute_C_matrix()
{
	Matrix Identity(n,n), OneOverN(n,n);
	Identity = 1.0;
//	OneOverN.setAll(1.0 / n); 
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			OneOverN(i, j) = 1.0 / n;

	C_Matrix = Identity - OneOverN;

#ifdef DEBUG
	cout << "Now print the result of the Center matrix." << endl;
	C_Matrix.print();
#endif
}


