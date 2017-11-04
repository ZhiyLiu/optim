/********************************************************************************/
/*																				*/
/*  	File	:  Zerofinder.h													*/
/*																				*/
/*	Description:  declaration for class to find zeros of T(u, v) in the range	*/
/*					[0, 1]x[0, 1] for T as defined below.  This used to find	*/
/*					nearpoint (according to Phong-derived normal-vector) to		*/
/*					a triangle.													*/
/*																				*/
/*	Project :  Seurat															*/
/*																				*/
/*	Dependencies:  Bpoints and DbVectors only, so need Xferlist.h and LinAlg.h	*/
/*																				*/
/*	Author  :  A. Thall															*/
/*																				*/
/*	Date	:  6. August 2001													*/
/*	Modifications:																*/
/********************************************************************************/
class Zerofinder
{
	DbVector3 A, B, C, nA, nB, nC;
	DbVector3 xbar;

	DbVector3 pbar_fun(double u, double v);
	DbVector3 pbar_fun_u(double u, double v);
	DbVector3 pbar_fun_v(double u, double v);

	DbVector3 nvec_fun(double u, double v);
	DbVector3 nvec_fun_u(double u, double v);
	DbVector3 nvec_fun_v(double u, double v);


	void ComputeTandDerivs_V0(double u);
	void ComputeTandDerivs_V1(double u);
	void ComputeTandDerivs_U0(double v);

	DbVector3 pbar_funV0(double u, double v);
	DbVector3 nvec_funV0(double u, double v);
	DbVector3 pbar_funV1(double u, double v);
	DbVector3 nvec_funV1(double u, double v);
	DbVector3 pbar_funU0(double u, double v);
	DbVector3 nvec_funU0(double u, double v);

	double T_V0(double u);
	double dT_V0(double u);
	double T_V1(double u);
	double dT_V1(double u);
	double T_U0(double v);
	double dT_U0(double v);

	double dbrent(double ax, double bx, double cx, double tol,
		          double (Zerofinder::*f)(double), double (Zerofinder::*df)(double), double *xmin);

	// Keep track of recursive depth in findroot4(), so can't get in an infinite regress
	//    (if the test fails, will simply go to one of the slower, non-recursive)
	int rd;

public:
	Zerofinder() { rd = 0;}

	double T_fun(double u, double v);
	double T_fun_u(double u, double v);
	double T_fun_v(double u, double v);

	// same as above, but assumes pbar and nvec already set
	double T_fun(const DbVector3& pbar, const DbVector3& nvec);

	// This is called by the (non-member)function funfunfunFAST() in Zerofinder.cpp,
	//    which is passed to the optimization code in Conjgrad.cpp .
	double Tfast, T_ufast, T_vfast;
	void ComputeTandDerivs(double u, double v);

	// find u and v coordinates
	bool findroot(double &u, double &v, DbVector3& p_bar, DbVector3& n_vec,
		const DbVector3& x_bar, const DbVector3& a,
		const DbVector3& b, const DbVector3& c, const DbVector3& n_A,
		const DbVector3& n_B, const DbVector3& n_C);

	// find u and v coordinates -- adaptively localized brute force search
	bool findroot2(double &u, double &v, DbVector3& p_bar, DbVector3& n_vec,
		const DbVector3& x_bar, const DbVector3& a,
		const DbVector3& b, const DbVector3& c, const DbVector3& n_A,
		const DbVector3& n_B, const DbVector3& n_C);

	// find u and v coordinates	-- most brutal force on 100x100 grid
	bool findroot3(double &u, double &v, DbVector3& p_bar, DbVector3& n_vec,
		const DbVector3& x_bar, const DbVector3& a,
		const DbVector3& b, const DbVector3& c, const DbVector3& n_A,
		const DbVector3& n_B, const DbVector3& n_C);

	// find u and v coordinates	-- subdivide using Loop splitting (no averaging)
	//    and recursive legality check of subtriangles
	bool findroot4(Bpoint *bvals, const DbVector3& x_bar,
				   DbVector3& a, DbVector3& b, DbVector3& c,
				   DbVector3& n_A, DbVector3& n_B, DbVector3& n_C,
				   DbVector3& uvt_A, DbVector3& uvt_B, DbVector3& uvt_C,
				   double r_A, double r_B, double r_C);

	bool findroot_boundary(double &u, double &v, DbVector3& p_bar, DbVector3& n_vec,
		const DbVector3& x_bar, 
		const DbVector3& a, const DbVector3& b, const DbVector3& c,
		const DbVector3& n_A, const DbVector3& n_B, const DbVector3& n_C);

	bool drawpointPandtval(double u, double v,
						   const DbVector3& x_bar,
						   const DbVector3& a, const DbVector3& b, const DbVector3& c,
						   const DbVector3& n_A, const DbVector3& n_B, const DbVector3& n_C);

};

