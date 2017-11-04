#include <math.h>
#include <stdarg.h>
#include <assert.h>
#include <vector>

#include "M3DPrimitive.h"
#include "Image3D.h"
#include "NormalAlignedShapeSpace.h"

using namespace std;

const int INVALID_PRIMITIVE_ID = -1;

#ifndef PRODUCTION_VERSION
Image3D * M3DPrimitive::worldImage;
#endif

M3DPrimitive::M3DPrimitive() : x(), hinge(false), selected(0) {}

M3DPrimitive::M3DPrimitive(double x0, double x1, double x2) 
: x(x0,x1,x2), hinge(false), selected(0) {}

M3DPrimitive::M3DPrimitive(const Vector3D& X) 
: x(X), hinge(false), selected(0) {}

M3DPrimitive::M3DPrimitive(const M3DPrimitive & prim)
: x(prim.x), hinge(prim.hinge), selected(prim.selected) {}

// Assignment operator
M3DPrimitive & M3DPrimitive::operator=(const M3DPrimitive & prim)
{
    copy(prim);
    return *this;
}

#ifndef PRODUCTION_VERSION
void M3DPrimitive::setWorld(Image3D * image)
{
	worldImage = image;
}
#endif

void M3DPrimitive::print(ostream & out, char * prefix, bool marked) const
{
    std::string p;
    if (prefix) p = prefix;

    out << p << "x = " << x;
    out << (selected & 0x1 ? "Selected" : "        ") << "  ";
    out << (!(selected & 0x2) ? "SelectedForRegularity" : "        ") << "  ";
    out << (marked ? "Marked" : "      ");
    out << (hinge ? "  (Hinge Atom)\n" : "\n");
}

#if 1
void M3DPrimitive::setR(double r)
{
	// cout << __FILE__ << ":" << __LINE__ << " is at wrong place.\n";
	cout << __FILE__ << ":" << __LINE__ << " shouldn't get called.\n";
}
#endif

Vector3D M3DPrimitive::getBPerp() const
{
    Vector3D b = getB();
    Vector3D n = getN();
    Vector3D bPerp = b.cross(n);
    return bPerp;
}

void M3DPrimitive::writePrimitiveOld(Registry& registry, const char * regStr, ...) const
{
	cout << __FILE__ << ":" << __LINE__ << " is deprecated.\n";
	assert(false);
    char newStr[1024];

    va_list val;
    va_start(val, regStr);
    vsprintf(newStr, regStr, val);
    va_end(val);

    strcat(newStr, ".%s");

    if(type() == M3D_STANDARD_PRIMITIVE)
        registry.setStringValue(newStr, (char*)M3D_STANDARD_PRIMITIVE_STR, "type");
    else if(type() == M3D_END_PRIMITIVE)
    {
        registry.setStringValue(newStr, (char*)M3D_END_PRIMITIVE_STR, "type");
        registry.setDoubleValue(newStr, dynamic_cast<const M3DEndPrimitive*>(this)->getElongation(), "elongation");
    }

    registry.setDoubleValue(newStr, getX().getX(), "x");
    registry.setDoubleValue(newStr, getX().getY(), "y");
    registry.setDoubleValue(newStr, getX().getZ(), "z");
    registry.setDoubleValue(newStr, getR(), "r");
    registry.setDoubleValue(newStr, getQ().getX(), "qx");
    registry.setDoubleValue(newStr, getQ().getY(), "qy");
    registry.setDoubleValue(newStr, getQ().getZ(), "qz");
    registry.setDoubleValue(newStr, getQ().getW(), "qw");
    registry.setDoubleValue(newStr, getTheta() * R_RADIANS_TO_DEGREES, "theta");
    registry.setIntValue(newStr, getSelectionFlags(), "selected");
    if(dynamic_cast<const M3DTubePrimitive*>(this) != NULL) {
        const M3DTubePrimitive* tubePtr = dynamic_cast<const M3DTubePrimitive*>(this);

        registry.setBooleanValue(newStr, tubePtr->isBaseAtom(), "baseAtom");

        strcat(newStr, "[%d]");
        for( int i = 0; i != tubePtr->getNumberOfSpokes(); ++i ) {
            //registry.setDoubleValue(newStr, tubePtr->dr[i], "dr", i);
        }
    }
    else if(dynamic_cast<const M3DQuadPrimitive*>(this) != NULL ) {
    }
    else {
        // Who created a new figure type?
        assert(false);
    }
}

/**
 * TODO: Compute average of atoms.
 */

#define MAX_ITERATION_NUMBER 200
const double EPSILON = 1.0e-03;
bool M3DPrimitive::atomsAlignment(int nPrims, M3DPrimitive** mSet1, M3DPrimitive** mSet2,
	Vector3D &tran, Quat &q, double &s, bool alignAll)
{
	cout << __FILE__ << ":" << __LINE__ << " is deprecated.\n";
	int				i;
	M3DPrimitive	**mSetTmp = new M3DPrimitive* [nPrims];

	Vector3D	cogs[2];

	for(i=0; i<2; i++)
		cogs[i].set(0, 0, 0);

	for(i=0; i<nPrims; i++) {
		if(mSet1[i]->type() == M3D_END_PRIMITIVE)
			mSetTmp[i] = new M3DQuadEndPrimitive();
		else
			mSetTmp[i] = new M3DQuadPrimitive();

		*mSetTmp[i] = *mSet1[i];

		cogs[0] += mSet1[i]->getX();
		cogs[1] += mSet2[i]->getX();
	}

	for(i=0; i<2; i++)
		cogs[i] /= nPrims;

	tran = cogs[1] - cogs[0];

	int				iterNum = 0;
	double			prevDis,
					dis,
					disPlus, disMinus;

	// no need to optimize translation
	// opt on rotation (+scale) only
	VectorND		f,
					deriv;

	if(alignAll) {
		f.setSize(7);
		deriv.setSize(7);
	}
	else {
		f.setSize(6);
		deriv.setSize(6);
	}

	f.set(0.0);
	deriv.set(0.0);

	// translation is already determined by the difference 
	// between the two COGs

	// the rest (rotation/scaling) is calculated by minimizing
	// the summed square geodesic distance between 2 sets of 
	// atoms, using gradient descent
	f.set(0, tran);

	dis = atomSetSquareDistance(nPrims, mSet1, mSet2);

	do {
		prevDis = dis;
		for(i=3; i<f.getSize(); i++) {
			f.set(i, f.get(i) + EPSILON);
			applyVector(f, nPrims, mSet1, mSetTmp);
			disPlus = atomSetSquareDistance(nPrims, mSetTmp, mSet2);

			f.set(i, f.get(i) - EPSILON*2.0);
			applyVector(f, nPrims, mSet1, mSetTmp);
			disMinus = atomSetSquareDistance(nPrims, mSetTmp, mSet2);

			f.set(i, f.get(i) + EPSILON);

			deriv.set(i, (disPlus-disMinus) / 2.0 / EPSILON);
		}

		double	stepSize = 20.0;
		int		trialNum = 0;
		f -= deriv * stepSize;
		applyVector(f, nPrims, mSet1, mSetTmp);
		dis = atomSetSquareDistance(nPrims, mSetTmp, mSet2);

		while(dis >= prevDis && trialNum < 40) {
			// recover the old "f" first
			f += deriv * stepSize;

			// cut the step size to half
			stepSize /= 2.0;
			f -= deriv * stepSize;
			applyVector(f, nPrims, mSet1, mSetTmp);
			dis = atomSetSquareDistance(nPrims, mSetTmp, mSet2);
			trialNum ++;
		}

		if(dis >= prevDis) {
			dis = prevDis;
			f  += deriv * stepSize;
		}

		iterNum++;
	}
	while(prevDis-dis > 1.0e-18 && 
		  iterNum < MAX_ITERATION_NUMBER);

	if(f.getSize() <= 0)
		return false;
	else if(f.getSize() == 1) {
		s = exp(f.get(0));
	}
	else if(f.getSize() == 3) {
		q = q.expMap( Vector3D(f.get(0), f.get(1), f.get(2)) );
	}
	else if(f.getSize() == 4) {
		q = q.expMap( Vector3D(f.get(0), f.get(1), f.get(2)) );
		s = exp(f.get(3));
	}
	else if(f.getSize() == 6) {
		tran.set(f.get(0), f.get(1), f.get(2));
		q = q.expMap( Vector3D(f.get(3), f.get(4), f.get(5)) );
		s = 1.0;
	}
	else if(f.getSize() == 7) {
		tran.set(f.get(0), f.get(1), f.get(2));
		q = q.expMap( Vector3D(f.get(3), f.get(4), f.get(5)) );
		s = exp(f.get(6));
	}
	else 
		return false;

	return true;
}


bool M3DPrimitive::	atomsAlignmentByNormal(int nPrims, M3DPrimitive** mSet1, M3DPrimitive** mSet2,
				   Vector3D &tran, Quat &q, double &s, bool alignAll) {
	cout << __FILE__ << ":" << __LINE__ << " is deprecated.\n";
	//align the mean normal of two atom sets
	int				i;
	M3DPrimitive	**mSetTmp = new M3DPrimitive* [nPrims];

	Vector3D		rotC[2];

	rotC[1] =  mSet1[nPrims-1]->getX(); //the aligning atom itself :the last atom in the set
	rotC[2] =  mSet2[nPrims-1]->getX();

	tran = rotC[1] - rotC[0];

    
	//calculate the quaternion q  to align the mean normal1 to mean normal2
	
	Vector3D		meanN[2];

	//mean Normals
    for(i=0; i<2; i++)
		meanN[i].set(0, 0, 0);

	for(i=0; i<nPrims; i++) {
    	meanN[0] += mSet1[i]->getN();
		meanN[1] += mSet2[i]->getN();
	}
	for(i=0; i<2; i++)
		meanN[i] /= nPrims;

    Vector3D crossTemp;//as the rotation axis
	crossTemp = meanN[0].cross(meanN[1]);
	double theta = acos(meanN[0]*meanN[1]);
	q.setAxisAngle(crossTemp,theta);
	

	// no need to optimize translation
	// opt on rotation (+scale) only

	return true;
}

SimilarityTransform3D  M3DPrimitive::atomsAlignmentSimtrans(int nPrims, M3DPrimitive** mSet1, M3DPrimitive** mSet2) {
	//align the mean normal of two atom sets
	SimTransShapeSpace* shapeSpace = 0;
	
	vector<M3DPrimitive*> refPrims(nPrims, (M3DPrimitive*) 0);
	vector<M3DPrimitive*> targetPrims(nPrims, (M3DPrimitive*) 0);
	vector<int> idxV;
	
	for (int idx =0; idx< nPrims;idx++){
		refPrims[idx] = mSet1[idx];
		targetPrims[idx] = mSet2[idx];
		idxV.push_back(idx);
	}
		
	shapeSpace = new NormalAlignedShapeSpace(refPrims, targetPrims, idxV);	
		
    int transWeight = 1;
    int scaleWeight = 1;	
    int quatWeight = 1;
    int thetaWeight = 1;
	
	// constructing GeodesicDistanceFunction
	
	//cout<<"Opt start\n";
	GeodesicDistanceFunction f(refPrims, targetPrims, shapeSpace);
	f.transWeight = transWeight;
	f.scaleWeight = scaleWeight;
	f.quatWeight = quatWeight;
	f.thetaWeight = thetaWeight;
	
	double bestY = 1e10;
	double bestX = 0;
	for (int x = 0; x < 256; x++)  {
		Vector v(1, (64.0 * x * R_PI) / (128.0));
		double y = f.evaluate(v);
		if (y < bestY) {
			//	cout << "***";
			bestY = y;
			bestX = x;
		}//	cout << "Angle: "  << x << " /64PI;  score: " << y << endl;
		
	}
	
	//	Vector best =  f.doOptimization();
	Vector best = Vector(1,(64.0 * bestX * R_PI) / (128.0));
	double v = f.evaluate(best);	



//	cout << "Final result: " << " bestX: " << best(0) << "... bestVal: " << v << endl;
		
	// Copy the transformed models and save the object	
	/*  vector<M3DPrimitive*> candPrims = f.createCandidateAtoms(best);
    M3DObject* bestObj = refModel->clone(); 
	
	  for (int idx = 0; idx<numPrims; idx++) {
	  if (candPrims[idx]) {
	  *(bestObj->getPrimitivePtr(idx)) = *(candPrims[idx]);
	  }
	  }
	  M3DObjectFile writer;
	  writer.write(outputModelName, *bestObj, 0, false, -1, 0, 0, true); //save PGA
	*/
	SimilarityTransform3D xform = f.createSimTrans(best);
		
	return xform;
}

bool M3DPrimitive:: composePrimitives(const M3DPrimitive * m1,  M3DPrimitive *  m2){
    SymPrimitive * m1Sym =  m1->convert2Sym();
    SymPrimitive  * m2Sym = m2->convert2Sym();
    *m1Sym += *m2Sym;  
    
    if (! m1Sym->convert2Lie(m2)){
        cout << "error at convert2Lie" << endl;
        return false;
    }else {
        if (m1->isSelected() &&  !m2->isSelected())
    //If one of the atom is selected, the composed atom should be selected too. 
            m2->select();
        return true;}
}

bool M3DPrimitive::subtractPrimitives(M3DPrimitive * m1, const M3DPrimitive * m2){
	SymPrimitive * m1Sym =  m1->convert2Sym();
	SymPrimitive  * m2Sym =  m2->convert2Sym();
	*m1Sym -= *m2Sym;             
			 
	if (! m1Sym->convert2Lie(m1)){
		cout << "error at convert2Lie" << endl;
        return false;
	}else 
		return true;
	
}

double M3DPrimitive::atomSetSquareDistance(int nPrims, M3DPrimitive **mSet1, 
											M3DPrimitive **mSet2, double *radius) 
{
	double dis = 0;
	for (int i = 0; i < nPrims; i++) {
		if (radius != NULL)
			dis += mSet1[i]->atomSquareDistance(mSet2[i], &radius[i]);
		else
			dis += mSet1[i]->atomSquareDistance(mSet2[i], NULL);
	}
	return dis;
}

bool M3DPrimitive::applyVector(VectorND & f, int nPrims, M3DPrimitive ** mSetBase,
					M3DPrimitive ** mSet) 
{
	int			i;

	Vector3D	tran;
	Quat		q;
	double		s;
	cout << __FILE__ << ":" << __LINE__ << " is deprecated.\n";

	tran.set(0, 0, 0);
	q.set(1, 0, 0, 0);
	s = 1.0;

	if(f.getSize() <= 0)
		return false;
	else if(f.getSize() == 1) {
		s = exp(f.get(0));
	}
	else if(f.getSize() == 3) {
		q.expMap(Vector3D(f.get(0), f.get(1), f.get(2)));
	}
	else if(f.getSize() == 4) {
		q = q.expMap(Vector3D(f.get(0), f.get(1), f.get(2)));
		s = exp(f.get(3));
	}
	else if(f.getSize() == 6) {
		tran.set(f.get(0), f.get(1), f.get(2));
		q = q.expMap(Vector3D(f.get(3), f.get(4), f.get(5)));
	}
	else if(f.getSize() == 7) {
		tran.set(f.get(0), f.get(1), f.get(2));
		q = q.expMap(Vector3D(f.get(3), f.get(4), f.get(5)));
		s = exp(f.get(6));
	}
	else 
		return false;

	Vector3D	cog,
				deltaX;

	cog.set(0, 0, 0);
	for(i=0; i<nPrims; i++) {
		cog += mSetBase[i]->getX();
	}
	cog /= nPrims;

	for(i=0; i<nPrims; i++) {
		*mSet[i] = *mSetBase[i];
		deltaX = mSet[i]->getX() - cog;

		if(f.getSize() > 1) {
			q.rotateVector(deltaX);
			mSet[i]->rotateBy(q);
		}
		if(f.getSize() == 1 ||
		   f.getSize() == 4 ||
		   f.getSize() == 7   ) {
			deltaX *= s;
			mSet[i]->scaleBy(s);
		}

		mSet[i]->setX( cog + deltaX );

		if(f.getSize() == 6 || f.getSize() == 7)
			mSet[i]->translateBy(tran);
	}

	return true;
}

void M3DPrimitive::copy(const M3DPrimitive & prim)
{
    x = prim.x;
    selected = prim.selected;
    hinge = prim.hinge;
}
