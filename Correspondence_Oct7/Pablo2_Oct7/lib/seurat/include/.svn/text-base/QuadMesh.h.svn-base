/************************************************************\
 * QuadMesh.h
 * A. Thall
 * 18. July 2002
 *
 * A quadmesh class as interface from Pointlist_server2 to
 * CCSubdivsurf. It takes a nxm Diatomgrid or Xferlist and
 * outputs vertexlist and tile-index list.  Note:  col > 1
 * and row > 1 is mandatory here.
 *
 * Modifications:
 *  2005/9/25: Rohit - Moved most of the work done in
 *  QuadMesh into Mesh to make a generic structure for tubes.
 * 
 *
\************************************************************/

class QuadMesh : public Mesh
{
  protected:
	int numdivs;	// number of samples along crest
					// between v1 and v2, not including
					// v1 and v2 (typically numdivs == 1,
					// just the bvector)
	int numrows;
	int numcols;
	int numsiderows;
	int numsidecols;

	/* The number of vertex-samples per atom
	 * (not including endcaps)
	 */
	// VC++ 6.0 does not let us define the default value here
	static const int NUM_SIDE_ROW;
	//static const int NUM_SIDE_ROW	= 1;

	// Functions which work on the entire structure.
	virtual void constructMesh0();
	virtual void computePolyList();

  public:
	// Constructors
	QuadMesh(Xferlist *thislist);
	QuadMesh(Diatomgrid *thisgrid);

	// Functions which work on co-ordinates
	virtual bool isValidMedcoord( const Medcoord& uvt);
	/**
	 * Takes a medial coordinate and says whether it lies on the
	 * ends or not.
	 * @param	uvt	a medial co-ordinate
	 * @return	true if on end/false otherwise
	 */
	virtual bool isEndMedcoord( const Medcoord& uvt );
	/**
	 * Takes n medial coordinates and computes their average.
	 * @param	n	number of medial coordinates.
	 * @param	uvt	an array of medial coordinates.
	 * @return	the average medial coordinate.
	 */
	virtual Medcoord average( int n, const Medcoord* uvt );
};
