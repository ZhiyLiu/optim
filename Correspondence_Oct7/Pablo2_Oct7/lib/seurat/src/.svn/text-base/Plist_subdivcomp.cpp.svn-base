/********************************************************************************/
/*                                                                              */
/*      File    :  Plist_subdivcomp.cpp											*/
/*                                                                              */
/*      Description:  class functions for new NPC_server, to compare local		*/
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
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <queue>
#include <set>
#define D_PLIST_SUBDIVCOMP
#define D_XFERLIST
#include "Shapedepend.h"

//#define DEBUG		/* Uncomment for verbose output */

using namespace ThallCode;

/********************************************************************************/
/* init() -- set the base mesh and tilesampling, and allocate the meshes for	*/
/*		storing the Bpointlists.												*/
/********************************************************************************/
void Plist_subdivcomp::init(Pointlist_server2 *pList1, int tilesampling)
{
    pListbase = pList1;
	tsamp = tilesampling;
	max_u = (int) pList1->max_u();
	max_v = (int) pList1->max_v();

	// For top&bottom 8 tile region, need 2 * (INTPOW(2, tsamp+1) + 1)*(INTPOW(2, tsampl+1) + 1) 

	// For edge atom, 8 tile region, need     (INTPOW(2, tsamp+2) + 1)*(INTPOW(2, tsampl+1) + 1)

	// for corner atom, 6 tile region, need (INTPOW(2, tsamp+1) + 1)*(INTPOW(2, tsamp+1) + 1)
	//                                       + 2 * (INTPOW(2, tsamp))*(INTPOW(2, tsamp))
}

// typedef enum { u0v0, u0vn, u0vmax, unvmax, umaxvmax, umaxvn, umaxv0, unv0, unvn } atomposition;

/********************************************************************************/	
/* comparemesh() -- Compute nearpoints for tiles adjacent to indicated atom in	*/
/*    both base pList and modified, and compute the average difference of		*/
/*    corresponding elements in the two lists of nearpoints						*/
/********************************************************************************/
double Plist_subdivcomp::comparemesh(Pointlist_server2 *pList2, int modatom_u, int modatom_v)
{
	pListmod = pList2;

	// Have 9 cases:
	//   u0v0, u0vn, u0vmax, unvmax, umaxvmax, umaxvn, umaxv0, unv0, unvn
	// where {unvn} is a central atom, {u0v0, u0vmax, umaxvmax, umaxv0} are corner atoms
	//   and { u0vn, unvmax, umaxvn, unv0 } are edge atoms

	// Determine the atom's mesh position
	atomposition apos;
	if (modatom_u == 0) {
		if (modatom_v == 0)
			apos = u0v0;
		else if (modatom_v == max_v)
			apos = u0vmax;
		else
			apos = u0vn;
	}
	else if (modatom_u == max_u) {
		if (modatom_v == 0)
			apos = umaxv0;
		else if (modatom_v == max_v)
			apos = umaxvmax;
		else
			apos = umaxvn;
	}
	else if (modatom_v == 0)
		apos = unv0;
	else if (modatom_v == max_v)
		apos = unvmax;
	else
		apos = unvn;

	double returnval;

	switch (apos) {
	case unvn:
		returnval = compare_faceatom(modatom_u, modatom_v);
		break;
	case u0v0:
	case u0vmax:
	case umaxvmax:
	case umaxv0:
		returnval = compare_corneratom(apos, modatom_u, modatom_v);
		break;
	case unv0:
	case unvmax:
	case u0vn:
	case umaxvn:
		returnval = compare_edgeatom(apos, modatom_u, modatom_v);
		break;
	default:
		returnval = -1.0;
	}

	return returnval;
}

/********************************************************************************/
/* compare_faceatom() -- for an interior atom, just compare 4-neighborhood of	*/
/*    tiles surrounding it, on top and on bottom.								*/
/********************************************************************************/		
double Plist_subdivcomp::compare_faceatom(int uval, int vval)
{
	// For top&bottom 8 tile region, need 2 * (INTPOW(2, tsamp+1) + 1)*(INTPOW(2, tsampl+1) + 1)
	Bpoint topbase, topmod, bottombase, bottommod;

	int numsamples = INTPOW(2, tsamp+1) + 1;
	double stepsize = 2.0/(numsamples - 1);

	int pointcount = 0;
	double runningtotal = 0;
	DbVector3 oldpnt, newpnt;

	for (int ucnt = 0; ucnt < numsamples; ucnt++) {
		double u_val = (uval - 1) + ucnt*stepsize;
		// clamp u_val at max if nec.
		u_val = (u_val > max_u ? max_u : u_val);

		for (int vcnt = 0; vcnt < numsamples; vcnt++) {
			double v_val = (vval - 1) + vcnt*stepsize;
			// just in case, let's clamp values
			v_val = (v_val > max_v ? max_v : v_val);

			pListbase->subdivBposition(&topbase, u_val, v_val, 1.0);
			pListbase->subdivBposition(&bottombase, u_val, v_val, -1.0);
			pListmod->subdivBposition(&topmod, u_val, v_val, 1.0);
			pListmod->subdivBposition(&bottommod, u_val, v_val, -1.0);

			oldpnt = DbVector3(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
			newpnt = DbVector3(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
			runningtotal += (oldpnt - newpnt).length()*(oldpnt - newpnt).length()/(topbase.rad*topbase.rad);			
			pointcount++;

			oldpnt = DbVector3(bottombase.pnt[0], bottombase.pnt[1], bottombase.pnt[2]);
			newpnt = DbVector3(bottommod.pnt[0], bottommod.pnt[1], bottommod.pnt[2]);
			runningtotal += (oldpnt - newpnt).length()*(oldpnt - newpnt).length()/(topbase.rad*topbase.rad);
			pointcount++;

			glBegin(GL_POINTS);
				glColor3d(1.0, 0.0, 0.0);
				glVertex3d(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
				glVertex3d(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
			glEnd();
			glBegin(GL_POINTS);
				glColor3d(0.0, 1.0, 1.0);
				glVertex3d(bottombase.pnt[0], bottombase.pnt[1], bottombase.pnt[2]);
				glVertex3d(bottommod.pnt[0], bottommod.pnt[1], bottommod.pnt[2]);
			glEnd();
		}

	}
	return runningtotal/pointcount;
}

// typedef enum { u0v0, u0vn, u0vmax, unvmax, umaxvmax, umaxvn, umaxv0, unv0, unvn } atomposition;

/********************************************************************************/
/* compare_edgeatom() -- compare top and bottom two-neighborhood, and the four	*/
/*    side tiles bounding the edge atom											*/	
/********************************************************************************/	
double Plist_subdivcomp::compare_edgeatom(atomposition apos, int uval, int vval)

{
    // For edge atom, 8 tile region, need     (INTPOW(2, tsamp+2) + 1)*(INTPOW(2, tsampl+1) + 1)
	Bpoint topbase, topmod, bottombase, bottommod;

	int numusamples, numvsamples, numtsamples;
	double ustepsize, vstepsize, tstepsize;

	int pointcount = 0;
	double runningtotal = 0;
	DbVector3 oldpnt, newpnt;
	double u_val, v_val, t_val;
	int ucnt, vcnt, tcnt;

	// Two cases:  if along u0 or umax edge
	if (apos == u0vn || apos == umaxvn) {
		numusamples = INTPOW(2, tsamp) + 1;
		ustepsize = 1.0/(numusamples - 1);

		numvsamples = INTPOW(2, tsamp+1) + 1;
		vstepsize = 2.0/(numvsamples - 1);

		// numtsamples and stepsize (all stepsizes should be the same, of course)
		//   is same as for vsamples, but we'll skip first and last row, since
		//   it's already counted as top and bottom edge respectively
		numtsamples = INTPOW(2, tsamp+1) + 1;
		tstepsize = 2.0/(numtsamples - 1);

		// Do top and bottom
		for (ucnt = 0; ucnt < numusamples; ucnt++) {

			if (apos == u0vn)
				u_val = uval + ucnt*ustepsize;
			else
				u_val = uval - ucnt*ustepsize;
			// not necessary to clamp, here, since initial value is integer

			for (vcnt = 0; vcnt < numvsamples; vcnt++) {
				v_val = (vval - 1) + vcnt*vstepsize;
				// just in case, let's clamp values
				v_val = (v_val > max_v ? max_v : v_val);

				pListbase->subdivBposition(&topbase, u_val, v_val, 1.0);
				pListbase->subdivBposition(&bottombase, u_val, v_val, -1.0);
				pListmod->subdivBposition(&topmod, u_val, v_val, 1.0);
				pListmod->subdivBposition(&bottommod, u_val, v_val, -1.0);

				oldpnt = DbVector3(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
				newpnt = DbVector3(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
				runningtotal += (oldpnt - newpnt).length()*(oldpnt - newpnt).length()/(topbase.rad*topbase.rad);
				pointcount++;

				oldpnt = DbVector3(bottombase.pnt[0], bottombase.pnt[1], bottombase.pnt[2]);
				newpnt = DbVector3(bottommod.pnt[0], bottommod.pnt[1], bottommod.pnt[2]);
				runningtotal += (oldpnt - newpnt).length()*(oldpnt - newpnt).length()/(topbase.rad*topbase.rad);
				pointcount++;

				glBegin(GL_POINTS);
				glColor3d(1.0, 0.0, 0.0);
				glVertex3d(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
				glVertex3d(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
				glEnd();
				glBegin(GL_POINTS);
				glColor3d(0.0, 1.0, 1.0);
				glVertex3d(bottombase.pnt[0], bottombase.pnt[1], bottombase.pnt[2]);
				glVertex3d(bottommod.pnt[0], bottommod.pnt[1], bottommod.pnt[2]);
				glEnd();
			}
		}

		// Do sides
		u_val = (apos == u0vn ? 0.0 : max_u);
		for (tcnt = 1; tcnt < numtsamples - 1; tcnt++) {

			t_val = -1.0 + tcnt*tstepsize;
			// Don't need to clamp here either, since skipping first and last

			for (int vcnt = 0; vcnt < numvsamples; vcnt++) {
				v_val = (vval - 1) + vcnt*vstepsize;
				// just in case, let's clamp values
				v_val = (v_val > max_v ? max_v : v_val);

				pListbase->subdivBposition(&topbase, u_val, v_val, t_val);
				pListmod->subdivBposition(&topmod, u_val, v_val, t_val);

				oldpnt = DbVector3(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
				newpnt = DbVector3(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
				runningtotal += (oldpnt - newpnt).length()*(oldpnt - newpnt).length()/(topbase.rad*topbase.rad);
				pointcount++;

				glBegin(GL_POINTS);
				glColor3d(1.0, 0.0, 0.0);
				glVertex3d(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
				glVertex3d(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
				glEnd();
			}
		}
	}
	// Else:  along v0 or vmax edge
	else {
		numvsamples = INTPOW(2, tsamp) + 1;
		vstepsize = 1.0/(numvsamples - 1);

		numusamples = INTPOW(2, tsamp+1) + 1;
		ustepsize = 2.0/(numusamples - 1);

		// numtsamples and stepsize (all stepsizes should be the same, of course)
		//   is same as for usamples, but we'll skip first and last row, since
		//   it's already counted as top and bottom edge respectively
		numtsamples = INTPOW(2, tsamp+1) + 1;
		tstepsize = 2.0/(numtsamples - 1);

		// Do top and bottom
		for (vcnt = 0; vcnt < numvsamples; vcnt++) {

			if (apos == unv0)
				v_val = vval + vcnt*vstepsize;
			else
				v_val = vval - vcnt*vstepsize;
			// not necessary to clamp, here, since initial value is integer

			for (ucnt = 0; ucnt < numusamples; ucnt++) {
				u_val = (uval - 1) + ucnt*ustepsize;
				// just in case, let's clamp values
				u_val = (u_val > max_u ? max_u : u_val);

				pListbase->subdivBposition(&topbase, u_val, v_val, 1.0);
				pListbase->subdivBposition(&bottombase, u_val, v_val, -1.0);
				pListmod->subdivBposition(&topmod, u_val, v_val, 1.0);
				pListmod->subdivBposition(&bottommod, u_val, v_val, -1.0);

				oldpnt = DbVector3(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
				newpnt = DbVector3(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
				runningtotal += (oldpnt - newpnt).length()*(oldpnt - newpnt).length()/(topbase.rad*topbase.rad);
				pointcount++;

				oldpnt = DbVector3(bottombase.pnt[0], bottombase.pnt[1], bottombase.pnt[2]);
				newpnt = DbVector3(bottommod.pnt[0], bottommod.pnt[1], bottommod.pnt[2]);
				runningtotal += (oldpnt - newpnt).length()*(oldpnt - newpnt).length()/(topbase.rad*topbase.rad);
				pointcount++;

				glBegin(GL_POINTS);
				glColor3d(1.0, 0.0, 0.0);
				glVertex3d(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
				glVertex3d(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
				glEnd();
				glBegin(GL_POINTS);
				glColor3d(0.0, 1.0, 1.0);
				glVertex3d(bottombase.pnt[0], bottombase.pnt[1], bottombase.pnt[2]);
				glVertex3d(bottommod.pnt[0], bottommod.pnt[1], bottommod.pnt[2]);
				glEnd();
			}
		}

		// Do sides
		v_val = (apos == unv0 ? 0.0 : max_v);
		for (tcnt = 1; tcnt < numtsamples - 1; tcnt++) {

			t_val = -1.0 + tcnt*tstepsize;
			// Don't need to clamp here either, since skipping first and last
			for (int ucnt = 0; ucnt < numusamples; ucnt++) {

				u_val = (uval - 1) + ucnt*ustepsize;
				// just in case, let's clamp values
				u_val = (u_val > max_u ? max_u : u_val);
#ifdef DEBUG
				fprintf(stderr, "u_val = %f\n", u_val);
#endif

				pListbase->subdivBposition(&topbase, u_val, v_val, t_val);
				pListmod->subdivBposition(&topmod, u_val, v_val, t_val);

#ifdef DEBUG
				if (u_val < (double) uval)
					topbase.print();
#endif

				oldpnt = DbVector3(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
				newpnt = DbVector3(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
				runningtotal += (oldpnt - newpnt).length()*(oldpnt - newpnt).length()/(topbase.rad*topbase.rad);
				pointcount++;

				glBegin(GL_POINTS);
				glColor3d(1.0, 0.0, 0.0);
				glVertex3d(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
				glVertex3d(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
				glEnd();
			}
		}
	}

	return runningtotal/pointcount;
}
/********************************************************************************/
/* compare_corneratom() -- compare top and bottom one-neighborhood and the four	*/
/*    side tiles bounding the corner atom										*/	
/********************************************************************************/	
double Plist_subdivcomp::compare_corneratom(atomposition apos, int uval, int vval)
{
	// for corner atom, 6 tile region, need (INTPOW(2, tsamp+1) + 1)*(INTPOW(2, tsamp+1) + 1)
	//                                       + 2 * (INTPOW(2, tsamp))*(INTPOW(2, tsamp))
	Bpoint topbase, topmod, bottombase, bottommod;

	int numusamples, numvsamples, numtsamples;
	double ustepsize, vstepsize, tstepsize;

	int pointcount = 0;
	double runningtotal = 0;
	DbVector3 oldpnt, newpnt;
	double u_val, v_val, t_val;
	int ucnt, vcnt, tcnt;

	numusamples = INTPOW(2, tsamp) + 1;
	ustepsize = 1.0/(numusamples - 1);

	numvsamples = INTPOW(2, tsamp) + 1;
	vstepsize = 1.0/(numvsamples - 1);

	// numtsamples twice the above---all stepsizes should be the same,
	//   of course---but we'll skip first and last row, since
	//   it's already counted as top and bottom edge respectively
	numtsamples = INTPOW(2, tsamp+1) + 1;
	tstepsize = 2.0/(numtsamples - 1);

	// Do top and bottom
	for (ucnt = 0; ucnt < numusamples; ucnt++) {

		if (apos == u0v0 || apos == u0vmax)
			u_val = uval + ucnt*ustepsize;
		else
			u_val = uval - ucnt*ustepsize;
		// not necessary to clamp, here, since initial value is integer
		//   and incr/decrement should always be to legal value

		for (vcnt = 0; vcnt < numvsamples; vcnt++) {

			if (apos == u0v0 || apos == umaxv0)
				v_val = vval + vcnt*vstepsize;
			else
				v_val = vval - vcnt*vstepsize;
			// not necessary to clamp, here, since initial value is integer
			//   and incr/decrement should always be to legal value

			pListbase->subdivBposition(&topbase, u_val, v_val, 1.0);
			pListbase->subdivBposition(&bottombase, u_val, v_val, -1.0);
			pListmod->subdivBposition(&topmod, u_val, v_val, 1.0);
			pListmod->subdivBposition(&bottommod, u_val, v_val, -1.0);

			oldpnt = DbVector3(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
			newpnt = DbVector3(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
			runningtotal += (oldpnt - newpnt).length()*(oldpnt - newpnt).length()/(topbase.rad*topbase.rad);
			pointcount++;

			oldpnt = DbVector3(bottombase.pnt[0], bottombase.pnt[1], bottombase.pnt[2]);
			newpnt = DbVector3(bottommod.pnt[0], bottommod.pnt[1], bottommod.pnt[2]);
			runningtotal += (oldpnt - newpnt).length()*(oldpnt - newpnt).length()/(topbase.rad*topbase.rad);
			pointcount++;

			glBegin(GL_POINTS);
			glColor3d(1.0, 0.0, 0.0);
			glVertex3d(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
			glVertex3d(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
			glEnd();
			glBegin(GL_POINTS);
			glColor3d(0.0, 1.0, 1.0);
			glVertex3d(bottombase.pnt[0], bottombase.pnt[1], bottombase.pnt[2]);
			glVertex3d(bottommod.pnt[0], bottommod.pnt[1], bottommod.pnt[2]);
			glEnd();
		}
	}

	double u_valCONST, v_valCONST;
	if (apos == u0v0 || apos == u0vmax)
        u_valCONST = 0.0;
	else
		u_valCONST = max_u;

	if (apos == u0v0 || apos == umaxv0)
		v_valCONST = 0.0;
	else
		v_valCONST = max_v;

	for (tcnt = 1; tcnt < numtsamples - 1; tcnt++) {		
		t_val = -1.0 + tcnt*tstepsize;
		// Don't need to clamp here either, since skipping first and last

		// Do u-constant side
		for (int vcnt = 0; vcnt < numvsamples; vcnt++) {
			if (apos == u0v0 || apos == umaxv0)
				v_val = vval + vcnt*vstepsize;
			else
				v_val = vval - vcnt*vstepsize;

			pListbase->subdivBposition(&topbase, u_valCONST, v_val, t_val);
			pListmod->subdivBposition(&topmod, u_valCONST, v_val, t_val);

			oldpnt = DbVector3(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
			newpnt = DbVector3(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
			runningtotal += (oldpnt - newpnt).length()*(oldpnt - newpnt).length()/(topbase.rad*topbase.rad);
			pointcount++;

			glBegin(GL_POINTS);
			glColor3d(1.0, 1.0, 0.0);
			glVertex3d(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
			glVertex3d(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
			glEnd();
		}

		// Do v-constant side (start from index 1, since 0 counted already in u-constant side
		for (int ucnt = 1; ucnt < numusamples; ucnt++) {

			if (apos == u0v0 || apos == u0vmax)
				u_val = uval + ucnt*ustepsize;
			else
				u_val = uval - ucnt*ustepsize;

			pListbase->subdivBposition(&topbase, u_val, v_valCONST, t_val);
			pListmod->subdivBposition(&topmod, u_val, v_valCONST, t_val);

			oldpnt = DbVector3(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
			newpnt = DbVector3(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
			runningtotal += (oldpnt - newpnt).length()*(oldpnt - newpnt).length()/(topbase.rad*topbase.rad);
			pointcount++;

			glBegin(GL_POINTS);
			glColor3d(0.0, 1.0, 0.0);
			glVertex3d(topbase.pnt[0], topbase.pnt[1], topbase.pnt[2]);
			glVertex3d(topmod.pnt[0], topmod.pnt[1], topmod.pnt[2]);
			glEnd();
		}
	}

	return runningtotal/pointcount;
}


