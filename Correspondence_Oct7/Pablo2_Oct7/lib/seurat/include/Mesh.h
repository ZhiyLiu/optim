/************************************************************\
 * Mesh.h
 *  2005/9/25: Rohit
 *
 *	A class meant to be used as a base class for all meshes.
 *  Current specializations include -
 *  . QuadMesh
 *  . TubeMesh
 *
 *  This is an abstract base class. The specialiazations
 * have to define the set of functions:
 *  . constructMesh0() - constructs the zero-level vertex list
 *  . computePolyList() - constructs the polygons list
 *  . isValidMedcoord() - returns true/false based on whether
 *   the co-ordinates are valid or not.
 *
 * Further possible improvements:
 *  Encapsulate the process of interpolation/averaging of
 * <u,v,t> coords within this framework.
 *
 *  A better design would implement this whole thing as a
 * pipeline rather than providing base/derived class support
 * to allow for different sub-division methods
 *
\************************************************************/

class Mesh
{
  protected:
	int npolys;
	int nverts;
	Bpoint *v_list;
	PindexList *pi_list;
	Diatomgrid *fig_grid;

	// Functions which work on the entire structure.
	virtual void constructMesh0() = 0;
	virtual void computePolyList() = 0;

  public:
    Mesh() : v_list(NULL), pi_list(NULL), fig_grid(NULL),
	  npolys(0), nverts(0) {}
    virtual ~Mesh();

	int numverts() { return nverts; }
	int numpolys() { return npolys; }

	Bpoint *vertlist() { return v_list; }
	PindexList *polylist() { return pi_list; }

	void CopyDiatomGrid(Diatomgrid *thisGrid);
	void CopyXferList(Xferlist *thisList);
	void UpdateVertexList();
	void UpdateFaceList();

	// Some helpful debugging functions
	// (euphemism for left-overs)
	void printvals(char *message = NULL);
	void glRender(RenderStyle thisstyle);

	// Functions which work on co-ordinates
	virtual bool isValidMedcoord( const Medcoord& uvt) = 0;

	/**
	 * Takes a medial coordinate and says whether it lies on the
	 * ends (whatever the definition of end is) or not.
	 * @param	uvt	a medial co-ordinate
	 * @return	true if on end/false otherwise
	 */
	virtual bool isEndMedcoord( const Medcoord& uvt ) = 0;

	/**
	 * Takes n medial coordinates and computes their average.
	 * @param	n	number of medial coordinates.
	 * @param	uvt	an array of medial coordinates.
	 * @return	the average medial coordinate.
	 */
	virtual Medcoord average( int n, const Medcoord* uvt ) = 0;
};
