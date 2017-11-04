#ifndef SIMILARITY_TRANSFORM_3D_H
#define SIMILARITY_TRANSFORM_3D_H

#include "Vector3D.h"
#include "Quat.h"



class SimilarityTransform3D
{
public:
    SimilarityTransform3D();
    SimilarityTransform3D(const Vector3D & trans, const Quat & rot, double s, double p = 0.0 );
    SimilarityTransform3D(const SimilarityTransform3D & trans);
	explicit SimilarityTransform3D(const char * filename);

    SimilarityTransform3D & operator = (const SimilarityTransform3D & trans);

    Vector3D getTranslation() const { return translation; }
    Quat getRotation()        const { return rotation; }
    double getScale()         const { return scale; }
	double getTubePhi()       const { return tubePhi; }


	// Initialize the similarity transform
    void setTranslation(const Vector3D & trans) { translation = trans; }
    void setRotation(const Quat & rot) { rotation = rot; }
    void setScale(double s) { scale = s; }
	void setTubePhi(double p) { tubePhi = p; }

    void setToIdentity();

    // Change the similarity transform
    void scaleBy(double s);
    void translateBy(const Vector3D & v);
    void rotateBy(const Quat & q);
    void invert();
	SimilarityTransform3D inverse() const {
		SimilarityTransform3D i	= *this;
		i.invert();
		return i;
	}
	void rotateAlongAxisBy(double p);

    // Left-multiplication
    void multiplyBy(const SimilarityTransform3D & trans);

	// Apply the similarity transformation to a vector
    void transformVector(Vector3D & v) const;

	bool readSimilarity(const char * filename);
	bool writeSimilarity(const char * filename, bool asMatrix = false);
	// Specify the registry to use for the next I/O operation
	void tie(void * registryPtr) { registry = registryPtr; }
	// This returns true if the transform is the same as in the input file
	bool wasRead() const { return input; }

	void print() const;

private:

	void * registry;
    Vector3D translation;
    Quat rotation;
    double scale;
	//
	// As a similarity ttransform is used for handling correspondence, and so is the tube phi parameter,
	// we need to keep tube phi in here. Maybe this class should be called Correspodencexxx - rohit
	//
	double tubePhi;

	bool input;
};



#endif

