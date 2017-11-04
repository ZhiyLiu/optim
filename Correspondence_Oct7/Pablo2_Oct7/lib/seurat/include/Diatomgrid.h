/**
 * We need to get *rid* of both Xferlist and Diatomgrid and use M3DFigure directly
 * in seurat, however, it's too much of work to try and rectify other's mistakes,
 * so I'll simply include a pointer to M3DFigure in both these classes and
 * henceforth everyone should use functions/data from that pointer for all new things.
 *
 * - Rohit Saboo
 */

/********************************************************************************/
/*																				*/
/*  	File	:  Diatomgrid.H													*/
/*																				*/
/*	Description:  class declaration for grid of Diatoms.  Just an open struct	*/
/*		with constructor and destructor and a few utility functions.  Had been	*/
/*		included in Diatom.H													*/
/*																				*/
/*	Project :  Seurat															*/
/*	Author  :  A. Thall															*/
/*	Date	:  4. June 2000														*/
/*																				*/
/*	Modifications:  															*/
/*		4. June -- added true_badmesh() and true_badatom() functions			*/
/*		4. June -- made readXferlist() and readXferlist as constructor().		*/
/*		29. Jan 01 -- added update_mesh()										*/
/*		24. June 02 -- removed all .plist Parentset code						*/
/********************************************************************************/

/********************************************************************************/
/* Class Diatomgrid -- quick'n'dirty little open struct, but class-encaps'ed	*/
/*	  to define destructor()													*/
/********************************************************************************/

class ::M3DFigure;

class Diatomgrid
{
public:
    Diatom *dlist;		// a rows x cols array of Diatoms
 						// giving a rectangular mesh in column-major order (index = col*rows + row)
    int rows, cols;
	const M3DFigure* figure;	// The original figure from where this diatomgrid comes.

	Diatomgrid() 
	{ 
		dlist = NULL; 
	}

	Diatomgrid(Xferlist* mfig1);

	~Diatomgrid() 
	{ 
		if(dlist != NULL) delete []dlist;
	}

	void readXferlist(Xferlist *mfig1);

	void CopyXferlist(Xferlist *mfig1);

	int idx(int row, int col) 
	{ 
		return col*rows + row; 
	}

	// update (row, col) element of Diatomgrid with the values in the Xferatom
	void update_mesh(XferAtom *thisatom, int modrow, int modcol);
	void update_mesh(Diatom *thisatom, int modrow, int modcol);

	// This rotates EDGE_M Diatoms to set bvec (q frame) to point outward
	//    from the Grid, based on mesh vectors between the bad Diatom and its
	//    two neighbors on the edge.
	// The other function trues a single selected atom
	void true_badmesh();
	void true_badatom(int row, int col);
};

