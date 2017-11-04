/********************************************************************************/
/*                                                                              */
/*      File    :  Plist_subdivcomp.h											*/
/*                                                                              */
/*      Description:  class definition for new NPC_server, to compare local		*/
/*					regions near a modified Diatom.  Initialize with base		*/
/*					Pointlist_server2 pList, then make queries by passing in	*/
/*					a modified pList and the (u, v) coordinates of the modified	*/
/*					Diatom.														*/
/*																				*/
/*      Project :  Seurat                                                       */
/*                                                                              */
/*      Author  :  A. Thall, Shawn Liu                                          */
/*                                                                              */
/*      Date    :  16. October 2001												*/
/*                                                                              */
/*      Modifications:                                                          */
/********************************************************************************/

/********************************************************************************/
/* Have 9 cases:							
/*   u0v0, u0vn, u0vmax, unvmax, umaxvmax, umaxvn, umaxv0, unv0, unvn, where
/*   {unvn} is a central atom
/*   {u0v0, u0vmax, umaxvmax, umaxv0} are corner atoms
/*   {u0vn, unvmax, umaxvn, unv0} are edge atoms
/********************************************************************************/
typedef enum { u0v0, u0vn, u0vmax, unvmax, umaxvmax, umaxvn, umaxv0, unv0, unvn } atomposition;

class Plist_subdivcomp
{
	Pointlist_server2 *pListbase, *pListmod;

	int tsamp;
	int max_u, max_v;

	double compare_faceatom(int uval, int vval);
	double compare_corneratom(atomposition apos, int uval, int vval);
	double compare_edgeatom(atomposition apos, int uval, int vval);

public:
	// inline constructor and destructor---see below
	Plist_subdivcomp() { ; }
	~Plist_subdivcomp() { ; }

	// Initialize with base Pointlist, and tilesampling in u, v, and t
	// int tilesampling divides up the tiles---0 ==> tilecorners only
	//                                         1 ==> midpoints of edges & faces as well
	//                                         2 ==> 1/4 points along edges & faces
	//                                         3 ==> 1/8 points along edges & faces
	// NOTE: because of bisection of (u, v) coordinates, tilesampling <= subdivdepth will give
	//    limit points on surface, exactly and not approximately.
	// NOTE:  Pointlist_server2::ComputeSubdivBoundaryTiles() must be called before the
	//    the pList is passed to this class, and likewise below in comparemesh
	void init(Pointlist_server2 *pList1, int tilesampling = 3);

	// Compute nearpoints for tiles adjacent to indicated atom in both
	//    base pList and modified, and compute the average difference
	//    of corresponding elements in the two lists of nearpoints
	double comparemesh(Pointlist_server2 *pList2, int modatom_u, int modatom_v);
};

