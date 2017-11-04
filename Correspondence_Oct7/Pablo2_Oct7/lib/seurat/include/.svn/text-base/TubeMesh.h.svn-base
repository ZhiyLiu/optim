/************************************************************\
 * TubeMesh.h
 * 2005/9/25: Rohit
 *
 * A tubemesh class as interface from Pointlist_server2 to
 * CCSubdivsurf. It takes a 1 x n Diatomgrid or Xferlist and
 * outputs vertexlist and tile-index list.  Note:  col=1
 * and row > 1 is mandatory here.
 *
 * Modifications:
 *
\************************************************************/

class TubeMesh : public Mesh
{
  protected:
	int numdivs;	// number of circular samples
	int TubeLength;

	// Functions which work on the entire structure.
	virtual void constructMesh0();
	virtual void computePolyList();

  public:
	// Constructors
	TubeMesh(Xferlist *thislist);
	TubeMesh(Diatomgrid *thisgrid);

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
