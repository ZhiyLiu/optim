
#include "deltas.h"

using std::vector;
using std::cout;

/* 
	Function to return the difference in two M3DPrimitives, atom1 - atom0.
	The arguments are not commutative.
*/

DeltaPrimitive dPrimitive(M3DPrimitive & atom0, M3DPrimitive & atom1)
{
	DeltaPrimitive diff;

	diff.setX(atom1.getX() - atom0.getX());
	Quat q0 = atom0.getQ();
	Quat q1 = atom1.getQ();
	q0 = q0.conj();
	diff.setQ(q1*q0);
	diff.setR(atom1.getR()/atom0.getR() - 1.0);	// relative difference
	diff.setTheta(atom1.getTheta() - atom0.getTheta());
	return diff;
}




/* 
	Function to add a difference of M3DPrimitives to an atom to get an
	adjusted atom, atom + dAtom.
*/

M3DQuadPrimitive incrPrimitive(M3DPrimitive & atom, DeltaPrimitive & dAtom)
{
	DeltaPrimitive diff;

	diff.setX(atom.getX() + dAtom.getX());
	Quat q0 = atom.getQ();
	Quat q1 = dAtom.getQ();
	diff.setQ(q1*q0);
	diff.setR(atom.getR()*(1.0 + dAtom.getR()));
	diff.setTheta(atom.getTheta() + dAtom.getTheta());
	return diff;
}



/* 
	Function to subtract a difference of M3DPrimitives from an atom to get
	an adjusted atom, atom - dAtom.
*/

M3DQuadPrimitive decrPrimitive(M3DPrimitive & atom, DeltaPrimitive & dAtom)
{
	DeltaPrimitive diff;

	diff.setX(atom.getX() - dAtom.getX());
	Quat q0 = atom.getQ();
	Quat q1 = dAtom.getQ();
	q1 = q1.conj();
	diff.setQ(q1*q0);
	diff.setR(atom.getR()*(1.0 - dAtom.getR()));
	diff.setTheta(atom.getTheta() - dAtom.getTheta());
	return diff;
}



/* 
	Function to return the difference in two differences of M3DPrimitives,
	dAtom1 - dAtom0.  The arguments are not commutative.
*/

DeltaDeltaPrimitive ddPrimitive(DeltaPrimitive & dAtom0, DeltaPrimitive & dAtom1)
{
	DeltaDeltaPrimitive diff;

	diff.setX(dAtom1.getX() - dAtom0.getX());
	Quat q0 = dAtom0.getQ();
	Quat q1 = dAtom1.getQ();
	q0 = q0.conj();
	diff.setQ(q1*q0);
	diff.setR(dAtom1.getR() - dAtom0.getR());	// simple difference
	diff.setTheta(dAtom1.getTheta() - dAtom0.getTheta());
	return diff;
}


void dPrimitive(vector<M3DPrimitive *> & alist0, vector<M3DPrimitive *> & alist1,
				vector<M3DPrimitive *> & dlist) {
	int n = alist0.size();
	if (alist1.size() != n) {
		cout << "dPrimitive(): lists are not the same length\n";
		if (alist1.size() < n)
			n = alist1.size();
	}
	if (dlist.size() != n) {
		cout << "dPrimitive(): lists are not the same length\n";
		if (dlist.size() < n)
			n = dlist.size();
	}
	for (int i = 0; i < n; i++)
		*dlist[i] = dPrimitive(*alist0[i], *alist1[i]);
}


void ddPrimitive(vector<DeltaPrimitive *> & alist0, vector<DeltaPrimitive *> & alist1,
				 vector<M3DPrimitive *> & ddlist) {
	int n = alist0.size();
	if (alist1.size() != n) {
		cout << "ddPrimitive(): lists are not the same length\n";
		if (alist1.size() < n)
			n = alist1.size();
	}
	if (ddlist.size() != n) {
		cout << "ddPrimitive(): lists are not the same length\n";
		if (ddlist.size() < n)
			n = ddlist.size();
	}
	for (int i = 0; i < n; i++)
		*ddlist[i] = ddPrimitive(*alist0[i], *alist1[i]);
}


void incrPrimitive(vector<M3DPrimitive *> & alist0, vector<DeltaPrimitive *> & alist1,
				   vector<M3DPrimitive *> & list) {
	int n = alist0.size();
	if (alist1.size() != n) {
		cout << "incrPrimitive(): lists are not the same length\n";
		if (alist1.size() < n)
			n = alist1.size();
	}
	if (list.size() != n) {
		cout << "incrPrimitive(): lists are not the same length\n";
		if (list.size() < n)
			n = list.size();
	}
	for (int i = 0; i < n; i++)
		*list[i] = incrPrimitive(*alist0[i], *alist1[i]);
}


void decrPrimitive(vector<M3DPrimitive *> & alist0, vector<DeltaPrimitive *> & alist1,
				   vector<M3DPrimitive *> & list) {
	int n = alist0.size();
	if (alist1.size() != n) {
		cout << "decrPrimitive(): lists are not the same length\n";
		if (alist1.size() < n)
			n = alist1.size();
	}
	if (list.size() != n) {
		cout << "decrPrimitive(): lists are not the same length\n";
		if (list.size() < n)
			n = list.size();
	}
	for (int i = 0; i < n; i++)
		*list[i] = decrPrimitive(*alist0[i], *alist1[i]);
}


