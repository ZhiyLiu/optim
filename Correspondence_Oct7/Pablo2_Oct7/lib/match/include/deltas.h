#ifndef DELTAS_H
#define DELTAS_H


#include <vector>
#include "M3DPrimitive.h"


/*
 * FIXME: This entire piece of code only works with quads,
 * What should we do for tubes? - rrs
 */
#define DeltaPrimitive	M3DQuadPrimitive
#define DeltaDeltaPrimitive	M3DQuadPrimitive


/*	Functions for computing "differentials" of atoms, as follows:
		dPrimitive(atom0, atom1) returns atom 1 - atom0.
		ddPrimitive(dAtom0, dAtom1) returns dAtom1 - dAtom0.
		incrPrimitive(atom, dAtom) returns atom + dAtom.
		decrPrimitive(atom, dAtom) returns atom - dAtom.

	These functions do not obey the usual algebraic rules (commutativity,
	associativity, etc.) in the case of computations involving r, for
	which the results are only approximate.
*/
DeltaPrimitive dPrimitive(M3DPrimitive & atom0, M3DPrimitive & atom1);
DeltaDeltaPrimitive ddPrimitive(DeltaPrimitive & dAtom0, DeltaPrimitive & dAtom1);
M3DQuadPrimitive incrPrimitive(M3DPrimitive & atom, DeltaPrimitive & dAtom);
M3DQuadPrimitive decrPrimitive(M3DPrimitive & atom, DeltaPrimitive & dAtom);


/*	Functions for computing "differentials" of vectors of atoms.
	These work by applying the functions above on atom lists.  All
	three lists should be the same length.  If they are not, the
	number of atoms processed will equal the length of the shortest
	list, and a warning message will be printed.
*/

void dPrimitive(std::vector<M3DPrimitive *> & alist0, std::vector<M3DPrimitive *> & alist1,
				std::vector<M3DPrimitive *> & dlist);
void ddPrimitive(std::vector<DeltaPrimitive *> & alist0, std::vector<DeltaPrimitive *> & alist1,
				 std::vector<M3DPrimitive *> & ddlist);
void incrPrimitive(std::vector<M3DPrimitive *> & alist0, std::vector<DeltaPrimitive *> & alist1,
				   std::vector<M3DPrimitive *> & list);
void decrPrimitive(std::vector<M3DPrimitive *> & alist0, std::vector<DeltaPrimitive *> & alist1,
				   std::vector<M3DPrimitive *> & list);


#endif	/* DELTAS_H */

