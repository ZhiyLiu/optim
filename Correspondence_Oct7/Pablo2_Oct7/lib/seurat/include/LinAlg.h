/********************************************************************************/
/*																				*/
/*  	File	:  LinAlg.H														*/
/*																				*/
/*	Description:  declarations for Shapemonger 2D and 3D double matrix			*/
/*		and vector routines.  (a rudimentary class, necessary because			*/
/*		smooch at present handles only float vectors and matrices.)				*/
/*																				*/
/*		contains:	class DbVector2												*/
/*				class DbVector3													*/
/*				class SqDbMatrix2												*/
/*	Project :  Seurat															*/
/*																				*/
/*	Author  :  A. Thall															*/
/*																				*/
/*	Date	:  14. January 1997													*/
/*	Modifications:	THIS CODE SUBSTITUTES FROM NagaGamma, and eliminates 2D		*/
/*				30. May 99 -- added operator /= to DbVector2&3					*/
/*				21. May 2000 -- added unary-() and binary*() for dot product	*/
/*				 1. Aug -- changed void returns on operators to DbVector3&		*/
/*				18. Sept -- readded 2D code needed for Surfmath routines		*/
/*					added back SqDblMatrix2 as well								*/
/*				25. Aug 2001 -- added == operator to DbVector3					*/
/*								added += and unary - to DbVector2				*/
/*								made operator= return (*this) rather than void.	*/
/*				18. March 02 -- added DbVector3::outsideface(), testing self	*/
/*						 against a triangle ABC in counterclockwise order.		*/
/*				4. May 02 -- added .x(), .y(), and .z() const functions, to		*/
/*						allow access when a const variable.						*/
/*				16. May 02 -- added set(double x[]) and put(double x[])			*/
/*						functions to DbVector3, as well as copy constructor.	*/
/*				19. June 02 -- added DbVector3::glVertex() and glNormal() for	*/
/*						output to OpenGL pipeline								*/			
/*																				*/
/********************************************************************************/


class DbVector2
{

    public:

        DbVector2() { ; }
        DbVector2(double xval, double yval) { vec[0] = xval; vec[1] = yval; }
        ~DbVector2() { ; }

        double operator() (int A) const { return vec[A]; }
        inline DbVector2& operator=(const DbVector2& dvec2);
		void operator*=(double scaleval) {
			vec[0] *= scaleval; vec[1] *= scaleval;
		}

		void operator +=(const DbVector2& dvec) {
			vec[0] += dvec.vec[0];
			vec[1] += dvec.vec[1];
		}

		DbVector2& set(double xval, double yval) {
			vec[0] = xval; vec[1] = yval; return (*this);
		}

		double& X() { return vec[0]; }
		double& Y() { return vec[1]; }

		double x() const { return vec[0]; }
		double y() const { return vec[0]; }

		inline DbVector2 operator*(double scaleval) const;
		inline DbVector2 operator/(double scaleval) const;

		DbVector2 operator+(const DbVector2& dvec) {
			return DbVector2(vec[0] + dvec.vec[0], vec[1] + dvec.vec[1]);
		}
		DbVector2 operator-(const DbVector2& dvec) {
			return DbVector2(vec[0] - dvec.vec[0], vec[1] - dvec.vec[1]);
		}

		DbVector2 operator-() {
			return DbVector2(-vec[0], -vec[1]);
		}

		inline DbVector2 normalize() const;

		inline bool selfnormalize();

		inline double length() const;
		inline double dot(const DbVector2& dvec2) const;

    private:

        double vec[2];
    	void error(bool errorNumber) const;
};

inline DbVector2& DbVector2::operator=(const DbVector2& dvec2)
{
    vec[0] = dvec2.vec[0];
    vec[1] = dvec2.vec[1];

	return (*this);
}

inline double DbVector2::length() const
{
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1]);
}

inline double DbVector2::dot(const DbVector2& dvec2) const
{
    return vec[0]*dvec2.vec[0] + vec[1]*dvec2.vec[1];
}

inline DbVector2 DbVector2::operator*(double scaleval) const
{
    return DbVector2(vec[0]*scaleval, vec[1]*scaleval);
}

inline DbVector2 DbVector2::operator/(double scaleval) const
{
    return DbVector2(vec[0]/scaleval, vec[1]/scaleval);
}

inline DbVector2 DbVector2::normalize() const
{
    double dblen = length();
    if (dblen != 0)
        return DbVector2(vec[0]/dblen, vec[1]/dblen); 
	else
		error(false);
    return (*this);
}

inline bool DbVector2::selfnormalize()
{
    double dblen = length();
    if (dblen != 0) 
	{
        vec[0] /= dblen;
        vec[1] /= dblen;
		return true;
    }
	else
	{
		error(true);
		return false;
	}
}    

class DbVector3
{

	public:

		DbVector3() { ; }
		DbVector3(double xval, double yval, double zval) {
			vec[0] = xval; vec[1] = yval; vec[2] = zval;
		}
		DbVector3(const DbVector3& dvec) {
			vec[0] = dvec.vec[0]; vec[1] = dvec.vec[1]; vec[2] = dvec.vec[2];
		}
		~DbVector3() { ; }

		double operator()(int A) const { return vec[A]; }

		bool operator==(const DbVector3& dvec) {
			return (vec[0] == dvec.vec[0] && vec[1] == dvec.vec[1] && vec[2] == dvec.vec[2]);
		}

		inline DbVector3& operator=(const DbVector3& dvec);
		void operator*=(double scaleval) {
			vec[0] *= scaleval; vec[1] *= scaleval; vec[2] *= scaleval;
		}
		void operator/=(double scaleval) {
			vec[0] /= scaleval; vec[1] /= scaleval; vec[2] /= scaleval;
		}

		void operator+=(const DbVector3& dvec) {
			vec[0] += dvec.vec[0]; vec[1] += dvec.vec[1]; vec[2] += dvec.vec[2];
		}

		void operator-=(const DbVector3& dvec) {
			vec[0] -= dvec.vec[0]; vec[1] -= dvec.vec[1]; vec[2] -= dvec.vec[2];
		}

		DbVector3& set(double xval, double yval, double zval) {
			vec[0] = xval; vec[1] = yval; vec[2] = zval; return (*this);
			return (*this);
		}

		DbVector3& set(const double x[]) {
			vec[0] = x[0]; vec[1] = x[1]; vec[2] = x[2];
			return (*this);
		}

		DbVector3 put(double x[]) const {
			x[0] = vec[0]; x[1] = vec[1]; x[2] = vec[2];
			return (*this);
		}

		double& X() { return vec[0]; }
		double& Y() { return vec[1]; }
		double& Z() { return vec[2]; }

		double x() const { return vec[0]; }
		double y() const { return vec[1]; }
		double z() const { return vec[2]; }

		DbVector3 operator*(double scaleval) const {
			return DbVector3(vec[0]*scaleval, vec[1]*scaleval, vec[2]*scaleval);
		}
		inline DbVector3 operator/(double scaleval) const;

		DbVector3 operator+(const DbVector3& dvec) const {
			return DbVector3(vec[0] + dvec.vec[0], vec[1] + dvec.vec[1],
							 vec[2] + dvec.vec[2]);
		}
		DbVector3 operator-(const DbVector3& dvec) const {
			return DbVector3(vec[0] - dvec.vec[0], vec[1] - dvec.vec[1],
							 vec[2] - dvec.vec[2]);
		}

		double operator*(const DbVector3& dvec) const{
			return vec[0]*dvec.vec[0] + vec[1]*dvec.vec[1] + vec[2]*dvec.vec[2];
		}

		DbVector3 operator-() const {
			return DbVector3(-vec[0], -vec[1], -vec[2]);
		}

		inline DbVector3 normalize() const;
		inline bool selfnormalize();

		inline double length() const;
		inline double dot(const DbVector3& dvec) const;
		inline DbVector3 cross(const DbVector3& dvec) const;

		// return 1 if self is on positive side of triangle formed by ABC in counterclockwise
		//   order.  return -1 if on negative side, and return 0 if within e-8 of zero
		inline int outsideface(const DbVector3& A, const DbVector3& B, const DbVector3& C) const;

	    void printvals(char *message = NULL);

		// This is part of a revision I need to make to get better coordination
		//   with standard math operations.  Left multiply by scalar.
		friend inline DbVector3 operator*(const double&, const DbVector3&);

		// call glVertex() on point
		void glVertex() { glVertex3d(vec[0], vec[1], vec[2]); }
		void glNormal() { glNormal3d(vec[0], vec[1], vec[2]); }

    private:

		double vec[3];
    	void error(bool errorNumber) const;
};

inline DbVector3& DbVector3::operator=(const DbVector3& dvec)
{
    vec[0] = dvec.vec[0];
    vec[1] = dvec.vec[1];
    vec[2] = dvec.vec[2];

	return (*this);
}

// See above
inline DbVector3 operator*(const double& scalar, const DbVector3& v)
{
	DbVector3 vect;

	vect = v*scalar;

//	vect.X() = scalar * v.X();
//	vect.Y() = scalar * v.Y();
//	vect.Z() = scalar * v.Z();

	return vect;
}

inline double DbVector3::length() const
{
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

inline double DbVector3::dot(const DbVector3& dvec) const
{
    return vec[0]*dvec.vec[0] + vec[1]*dvec.vec[1] + vec[2]*dvec.vec[2];
}

inline DbVector3 DbVector3::operator/(double scaleval) const
{
    return DbVector3(vec[0]/scaleval, vec[1]/scaleval, vec[2]/scaleval);
}

inline DbVector3 DbVector3::normalize() const
{
    double dblen = length();
    if (dblen != 0)
        return DbVector3(vec[0]/dblen, vec[1]/dblen, vec[2]/dblen); 
	else
		error(false);
    return (*this);
}

inline bool DbVector3::selfnormalize()
{
    double dblen = length();
    if (dblen != 0)
	{
        vec[0] /= dblen;
        vec[1] /= dblen;
        vec[2] /= dblen;
		return true;
    }
	else
	{
		error(true);
		return false;
	}
}

inline DbVector3 DbVector3::cross(const DbVector3& dvec) const
{
    double xval, yval, zval;
    xval = vec[1]*dvec.vec[2] - vec[2]*dvec.vec[1];
    yval = vec[2]*dvec.vec[0] - vec[0]*dvec.vec[2];
    zval = vec[0]*dvec.vec[1] - vec[1]*dvec.vec[0];

    return DbVector3(xval, yval, zval);
}

// return 1 if self is on positive side of triangle formed by ABC in counterclockwise	
//   order.  return -1 if on negative side, and return 0 if within e-10 of zero
inline int DbVector3::outsideface(const DbVector3& A, const DbVector3& B, const DbVector3& C) const
{
	const double ZEROTOL = 1.0e-10;

	DbVector3 AB = (B - A);//.normalize();
	DbVector3 AC = (C - A);//.normalize();
	DbVector3 ABCperp = AB.cross(AC);

	DbVector3 Aself = (*this) - A;

	double testval = Aself.dot(ABCperp);
	if (testval > ZEROTOL)
		return 1;
	else if (testval < -ZEROTOL)
		return -1;
	else
		return 0;
}



class SqDblMatrix2
{
    double mat[2][2];
public:
    SqDblMatrix2() { ; }
    SqDblMatrix2(double a00, double a01,
                 double a10, double a11) {
        mat[0][0] = a00; mat[0][1] = a01; mat[1][0] = a10; mat[1][1] = a11;
    }
    ~SqDblMatrix2() { ; }

    double operator()(int row, int col) const { return mat[row][col]; }

    SqDblMatrix2& set(double a00, double a01, double a10, double a11) {
        mat[0][0] = a00; mat[0][1] = a01; mat[1][0] = a10; mat[1][1] = a11;
        return (*this);
    }

    inline void operator=(const SqDblMatrix2& dmat);
    inline double determinant() const;
    inline double trace() const;
};

inline void SqDblMatrix2::operator=(const SqDblMatrix2& dmat)
{
    mat[0][0] = dmat.mat[0][0];
    mat[0][1] = dmat.mat[0][1];
    mat[1][0] = dmat.mat[1][0];
    mat[1][1] = dmat.mat[1][1];
}

inline double SqDblMatrix2::determinant() const
{
    return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
}

inline double SqDblMatrix2::trace() const
{
    return mat[0][0] + mat[1][1];
}




class DbMatrix3
{
    double mat[3][3];
public:
    DbMatrix3() { ; }
    DbMatrix3(double a00, double a01, double a02,
                 double a10, double a11, double a12,
		 double a20, double a21, double a22) {
        mat[0][0] = a00; mat[0][1] = a01; mat[0][2] = a02;
	mat[1][0] = a10; mat[1][1] = a11; mat[1][2] = a12;
	mat[2][0] = a20; mat[2][1] = a21; mat[2][2] = a22;
    }
    ~DbMatrix3() { ; }

    double& operator()(int row, int col) { return mat[row][col]; }

    DbMatrix3& set(double a00, double a01, double a02,
		      double a10, double a11, double a12,
		      double a20, double a21, double a22) {
        mat[0][0] = a00; mat[0][1] = a01; mat[0][2] = a02;
	mat[1][0] = a10; mat[1][1] = a11; mat[1][2] = a12;
	mat[2][0] = a20; mat[2][1] = a21; mat[2][2] = a22;
        return (*this);
    }

    inline void operator=(const DbMatrix3& dmat);
    inline double determinant() const;
    inline double trace() const;
};

inline void DbMatrix3::operator=(const DbMatrix3& dmat)
{
    mat[0][0] = dmat.mat[0][0];
    mat[0][1] = dmat.mat[0][1];
    mat[0][2] = dmat.mat[0][2];
    mat[1][0] = dmat.mat[1][0];
    mat[1][1] = dmat.mat[1][1];
    mat[1][2] = dmat.mat[1][2];
    mat[2][0] = dmat.mat[2][0];
    mat[2][1] = dmat.mat[2][1];
    mat[2][2] = dmat.mat[2][2];
}

inline double DbMatrix3::determinant() const
{
    double mydet;

    mydet = mat[0][0]*mat[1][1]*mat[2][2] 
          + mat[0][1]*mat[1][2]*mat[2][0]
          + mat[0][2]*mat[1][0]*mat[2][1]
          - mat[0][2]*mat[1][1]*mat[2][0]
          - mat[0][1]*mat[1][0]*mat[2][2]
          - mat[0][0]*mat[1][2]*mat[2][1];

    return mydet;
}

inline double DbMatrix3::trace() const
{
    return mat[0][0] + mat[1][1] + mat[2][2];
}


