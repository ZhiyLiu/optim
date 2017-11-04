#include <math.h>
#include <iostream>
#include "M3DTubeSimilarityPGAProblem.h"
#include "Tuning.h"


extern int globalVerbosity;


using namespace std;


M3DTubeSimilarityPGAProblem::M3DTubeSimilarityPGAProblem() : M3DSimilarityPGAProblem()
{
	tubePhiMode	= false;
}

M3DTubeSimilarityPGAProblem::~M3DTubeSimilarityPGAProblem()
{
}

// Return # of modes (up to FMAHAMAX) + count of fixed components
// (3*trans, 3*rot, scale).
//		AGG: This function probably should take the figureId as an argument.
//		FIXME: This is important when dealing with both tubes and quads at the 
//		same time - rrs.
int M3DTubeSimilarityPGAProblem::parameterCount() const
{
	int maxModes, numModes;

    if (refObject == NULL) 
        return 0;

    M3DFigureTreeNode* treeNode = refObject->getFigureTreeRoot(tree_index);
    if (treeNode == NULL)
        return 0;

	// Need to know figure id to get the right PG set
	int figureId = treeNode->getFigureId();
	int order = pga->findPGSetId(figureId);

	PGData * pgPtr = pga->getPGDataPtr(order);

	maxModes = (int) tuningWt(FigureMaxPGAModes);
	numModes = pgPtr->numPGs;

	if (maxModes < numModes)
		numModes = maxModes;

	if (0 == (int) tuningWt(FigureMahalanobisOnly))
		numModes += 7;

	// For the phi alignment of a tube
	if( dynamic_cast<M3DTubeFigure*>(refObject->getFigurePtr(figureId)) != NULL ) {
		tubePhiMode	= true;
	}
	else {
		tubePhiMode	= false;
	}

	// To accomodate for the tube phi mode
	if( tubePhiMode ) {
		numModes++;
	}

	return numModes;
}

// This function returns a copy of the model after the PGA, initial, and
// similarity transforms are applied.  In Binary Pablo, as a side effect,
// it also sets similarityObject to contain a copy of the model after the
// PGA and initial transforms are applied, but without the similarity
// transform applied.
M3DObject * M3DTubeSimilarityPGAProblem::createTargetObject(const Vector & x,
	int & figureId, bool predict)
{
    int primitiveCount;
    M3DObject * targetObject;
    M3DFigureTreeNode * treeNode;
    M3DFigure * figure;

    if (refObject == NULL || x.size() != parameterCount())
        return NULL;

    treeNode = refObject->getFigureTreeRoot(tree_index);
    if (treeNode == NULL)
        return NULL;

#ifdef BINARY
	figureId = (int) tuningWt(BpFigureId);	// MultiObject
#else
    figureId = treeNode->getFigureId();
#endif

	//JJ: It does not matter what object is assigned
	//    to 'targetObject'.
    targetObject = refObject->assign();

    figure = targetObject->getFigurePtr(figureId);
    if (figure == NULL) {
        delete targetObject;
        return NULL;
    }
    figure->select();
	M3DTubeFigure* tubeFigure	= dynamic_cast<M3DTubeFigure*>(figure);

    primitiveCount = figure->getPrimitiveCount();

	static int iterationCounter = 0;
	if (globalVerbosity >= 1)
		cout << iterationCounter << ":";
	iterationCounter++;

	std::vector<double> vals	= pgaCoefficients(x);

	// Apply single figure statistics
	int order = pga->findPGSetId(figureId);
	pga->doPGADeform(targetObject, vals, order, false);

    // Apply initial transformation
#ifdef UNFLIPPED
	Vector3D center(0.0, 1.0, 0.0);
#else
	Vector3D center(0.0, 0.0, 0.0);
#endif
#ifdef BINARY
    int useMOM = (tuningWt(BpDoMethodOfMoments) && 
		! tuningWt(BpFigureResetModelBefore));	
    const SimilarityTransform3D * xform = useMOM ? &GLOBALMOMTransform 
		: targetObject->getTransformation();
	figure->applySimilarity( *xform, center );
	similarityObject = targetObject->assign();
#else
	figure->applySimilarity( mainFigProblemSimTransform, center );
	//cout << "mainFigProblemSimTransform" << endl;
	//mainFigProblemSimTransform.print();
	//cout << endl;

	// Apply hand placement transform.
	const SimilarityTransform3D handTransform = pga->getHandPlacementTransform();

	figure->applySimilarity( handTransform, center );
	//cout << "handTransform:" << endl;
	//handTransform.print();
	//cout << endl;
#endif

	SimilarityTransform3D trf;

	if ((int) tuningWt(FigureMahalanobisOnly) == 0) {
		Vector3D trans, axis;
		double angle, scale;
		Quat finalRot;

		// Apply similarity transform
		trans.set(x(0), x(1), x(2));
		trans *= tuningWt(FigureTranslationFactor);
		if (globalVerbosity >= 1)
			cout << " Target obj: T(x,y,z): " << trans.getX() << ' ' << trans.getY() << ' ' << trans.getZ();

		axis.set(x(3), x(4), x(5));
		angle = tuningWt(FigureRotationFactor) * axis.normalize();
		finalRot.setAxisAngle(axis, angle);
		if (globalVerbosity >= 1)
			cout << " R(x,y,z,angle): " << axis.getX() << ' ' << axis.getY() << ' ' << axis.getZ() << ' ' << angle;

		scale = exp(tuningWt(FigureScaleFactor) * x(6));
		if (globalVerbosity >= 1)
			cout << " S: " << scale;

		trf.setRotation(finalRot);
		trf.setScale(scale);
		trf.setTranslation(trans);
	}

	if( tubePhiMode ) {
		assert( tubeFigure != NULL );
		int iphi;
		if (0 != (int) tuningWt(FigureMahalanobisOnly))
			iphi = 0;
		else
			iphi = 7;
		const double phi	= x(iphi)*tuningWt(FigureTubePhiFactor);
		trf.setTubePhi(phi);
		if (globalVerbosity >= 1)
			cout << " phi: " << phi;
	}

	if( addToResidueObject != NULL ) {
		const M3DFigure* baseFigure	= addToResidueObject->getFigurePtr(figureId);
		for( int i = 0; i != primitiveCount; ++i ) {
			M3DPrimitive* residuePrim	= figure->getPrimitivePtr(i);
			M3DPrimitive* prim			= baseFigure->getPrimitivePtr(i)->copyPtr();

			//cout << "Residue Prim:\n";
			//residuePrim->print();
			//cout << "\n";
			//cout << "base prim:\n";
			//prim->print();
			//cout << "\n";

			//M3DPrimitive::composePrimitives(residuePrim,prim);

			SymPrimitive * symPrim = prim->convert2Sym();
			SymPrimitive * symResiduePrim = residuePrim->convert2Sym();
			*symResiduePrim += *symPrim;
			// convert2Lie for tubes picks up the old "phi"
			// rotation for the atom and then saves the result
			// into the passed in atom.
			if (!symResiduePrim->convert2Lie(prim)) {
				cerr << "Error while converting back to lie group represntation." << endl;
				assert(false);
			}
			// composePrimitives does not do the correct job
			// and is also not capable of knowing and doing
			// the correct job as described above.
			//M3DPrimitive::composePrimitives( prim, residuePrim );

			//cout << "composed prim:\n";
			//prim->print();
			//cout << "\n";

			figure->setPrimitivePtr(i, prim);
		}
	}
	//cout << "optimization trf:" << endl;
	//trf.print();
	//cout << endl;
	figure->applySimilarityAboutCOG( &trf );
	if (globalVerbosity >= 1)
		cout << endl;

    return targetObject;
}


std::vector<double> M3DTubeSimilarityPGAProblem::pgaCoefficients(const Vector& x)
{
	int start;
	std::vector<double> coeffs;
	if (0 != (int) tuningWt(FigureMahalanobisOnly))
		start = 0;
	else
		start = 7;
	// To accomodate for the tube phi mode.
	if( tubePhiMode ) {
		start++;
	}
	if (globalVerbosity > 0) cout << "   [PGA: ";
	for (int i = start; i < x.size(); i++) {
		coeffs.push_back(x(i)*tuningWt(FigurePGAFactor));
		if (globalVerbosity > 0) cout << x(i)*tuningWt(FigurePGAFactor) << " ";
	}
	if (globalVerbosity > 0) cout << "]\n";
	return coeffs;
}


double M3DTubeSimilarityPGAProblemBoundsFunction::bound( int i )
{
	assert(i >= 0);
	if( problem->tubePhiMode ) {
		Vector bounds;
		if ((int) tuningWt(FigureMahalanobisOnly) == 0) {
			bounds	= Vector(9,
				tuningWt(FigureTranslationBound), tuningWt(FigureTranslationBound), tuningWt(FigureTranslationBound),
				tuningWt(FigureRotationBound)*R_PI/180.0, tuningWt(FigureRotationBound)*R_PI/180.0, tuningWt(FigureRotationBound)*R_PI/180.0,
				log(tuningWt(FigureScaleBound)),
				tuningWt(FigureTubePhiBound)*R_PI/180.0,
				tuningWt(FigurePGABound) );
		}
		else {
			bounds	= Vector(2, 
				tuningWt(FigureTubePhiBound)*R_PI/180.0,
				tuningWt(FigurePGABound) );
		}
		return bounds( (i < bounds.size()) ? i : (bounds.size() - 1) );
	}
	else {
		// else Use the defaults setup by M3DSimilarityPGAProblemBoundsFunction constructor.
		return bounds( (i < bounds.size()) ? i : (bounds.size() - 1) );
	}
}

double M3DTubeSimilarityPGAProblemBoundsFunction::factor( int i )
{
	assert(i >= 0);
	if( problem->tubePhiMode ) {
		Vector factors;
		if ((int) tuningWt(FigureMahalanobisOnly) == 0) {
			factors	= Vector(9,
				tuningWt(FigureTranslationFactor),
				tuningWt(FigureTranslationFactor),
				tuningWt(FigureTranslationFactor),
				tuningWt(FigureRotationFactor),
				tuningWt(FigureRotationFactor),
				tuningWt(FigureRotationFactor),
				tuningWt(FigureScaleFactor),
				tuningWt(FigureTubePhiFactor),
				tuningWt(FigurePGAFactor));
		}
		else {
			factors	= Vector(2,
				tuningWt(FigureTubePhiFactor),
				tuningWt(FigurePGAFactor));
		}
		return factors( (i < factors.size()) ? i : (factors.size() - 1) );
	}
	else {
		// else Use the defaults setup by M3DSimilarityPGAProblemBoundsFunction constructor.
		return factors( (i < factors.size()) ? i : (factors.size() - 1) );
	}
}

