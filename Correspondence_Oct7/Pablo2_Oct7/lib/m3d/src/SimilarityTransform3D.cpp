#include <iostream>
#include "Registry.h"
#include "Matrix4D.h"
#include "SimilarityTransform3D.h"

#include <string.h>

using namespace std;


SimilarityTransform3D::SimilarityTransform3D() :
    translation(0.0, 0.0, 0.0), rotation(0.0, 0.0, 0.0, 1.0), scale(1.0), tubePhi(0.0)
{
	registry = NULL;
	input = false;
}

SimilarityTransform3D::SimilarityTransform3D(const Vector3D & trans,
    const Quat & rot, double s, double p)
{
    translation = trans;
    rotation = rot;
    scale = s;
	tubePhi = p;
	registry = NULL;
	input = false;
}

SimilarityTransform3D::SimilarityTransform3D(
    const SimilarityTransform3D & trans)
{
    translation = trans.getTranslation();
    rotation = trans.getRotation();
    scale = trans.getScale();
	tubePhi = trans.getTubePhi();
	registry = trans.registry;
	input = false;
}

SimilarityTransform3D::SimilarityTransform3D(const char * filename)
{
	registry = NULL;
	readSimilarity(filename);
}

SimilarityTransform3D & SimilarityTransform3D::operator = (
    const SimilarityTransform3D & trans)
{
    translation = trans.getTranslation();
    rotation = trans.getRotation();
    scale = trans.getScale();
	tubePhi = trans.getTubePhi();
	registry = trans.registry;
	input = trans.input;

    return (*this);
}

void SimilarityTransform3D::setToIdentity()
{
    translation.set(0.0, 0.0, 0.0);
    rotation.set(0.0, 0.0, 0.0, 1.0);
    scale = 1.0;
	tubePhi	= 0.0;
	input = false;
}

void SimilarityTransform3D::scaleBy(double s)
{
    translation *= s;
    scale *= s;
	input = false;
}

void SimilarityTransform3D::translateBy(const Vector3D & v)
{
    translation += v;
	input = false;
}

void SimilarityTransform3D::rotateBy(const Quat & q)
{
    rotation = q * rotation;
    rotation.normalize();
    q.rotateVector(translation);
	input = false;
}

void SimilarityTransform3D::rotateAlongAxisBy(double p)
{
    tubePhi	+= p;
	if( tubePhi >= 2*R_PI ) {
		tubePhi	-= 2*R_PI;
	}
	input = false;
}


void SimilarityTransform3D::invert()
{
  SimilarityTransform3D ret;	// start with identity

  // apply inverse operations in reverse order
  ret.translateBy(-translation);
  ret.scaleBy    (1. / scale);
  Quat q = rotation.conj();
  ret.rotateBy   (q);
  ret.rotateAlongAxisBy(-tubePhi);

  translation = ret.translation;	// overwrite *this
  scale       = ret.scale;
  rotation    = ret.rotation;
  tubePhi	  = ret.tubePhi;

}

// Left-multiplication
void SimilarityTransform3D::multiplyBy(const SimilarityTransform3D & trans)
{
    Quat q = trans.getRotation();
    Vector3D newTranslation = translation;
    q.rotateVector(newTranslation);

    translation = newTranslation * trans.getScale() + trans.getTranslation();
    scale *= trans.getScale();
    rotation = trans.getRotation() * rotation;
	tubePhi += trans.tubePhi;
	input = false;
}

void SimilarityTransform3D::transformVector(Vector3D & v) const
{
    rotation.rotateVector(v);
    v *= scale;
    v += translation;
}

bool SimilarityTransform3D::readSimilarity(const char * filename)
{
    char newStr[1024];
	Registry * r;
	bool created, ret;
	double * a;
	int len;

	if (registry == NULL) {
		if (filename == NULL)
			return false;
		r = new Registry;
		registry = r;
		created = true;
	}
	else {
		r = (Registry *) registry;
		created = false;
	}

	if (created) {
		try {
			r->readFromFile(filename);
		}
		catch (RException excp) {
			// Check for "Not a valid file."
			const char * str = excp.message().c_str();
			if (str[0] != 'N' || str[1] != 'o' || str[2] != 't')
				excp.print(std::cout);
			setToIdentity();
			if (created)
				delete (Registry *) registry;
			registry = NULL;
			return false;
		}
	}

	char * base = "model.transformation";
	strcpy(newStr, base);
    strcat(newStr, ".%s");

	// The similarity transform is stored either as a matrix of in 3 parts
	ret = true;
	a = r->getDoubleArray(newStr, &len, "matrix");
	if (a ==  NULL) {
		if (r->hasKey(base)) {
			scale = r->getDoubleValue(newStr, 1.0, "scale");
			translation.setX(r->getDoubleValue(newStr, 0.0, "translation.x"));
			translation.setY(r->getDoubleValue(newStr, 0.0, "translation.y"));
			translation.setZ(r->getDoubleValue(newStr, 0.0, "translation.z"));
			rotation.setX(r->getDoubleValue(newStr, 0.0, "rotation.x"));
			rotation.setY(r->getDoubleValue(newStr, 0.0, "rotation.y"));
			rotation.setZ(r->getDoubleValue(newStr, 0.0, "rotation.z"));
			rotation.setW(r->getDoubleValue(newStr, 1.0, "rotation.w"));
			tubePhi	= r->getDoubleValue(newStr, 0.0, "tubePhi" );
			input = true;
		}
		else
			input = false;
	}
	else {
		input = true;
		if (len == 16) {
			Matrix4D mat(a);
			delete [] a;
			mat.transpose();
			mat.decompose(rotation, scale, translation);
		}
		else {
			cout << "Invalid similarity transformation matrix: discarding" << endl;
			setToIdentity();
			ret = false;
		}
	}

	if (created)
		delete (Registry *) registry;
	registry = NULL;
	return true;
}

bool SimilarityTransform3D::writeSimilarity(const char * filename, bool asMatrix)
{
    char newStr[1024];
	Registry * r;
	bool created;
	bool ret;

	created = false;
	if (registry == NULL) {
		registry = new Registry;
		created = true;
	}
	r = (Registry *) registry;

	strcpy(newStr, "model.transformation");
    strcat(newStr, ".%s");

	if (asMatrix) {
		Matrix4D mat;

		mat.rotate(rotation);
		mat.scale(getScale());
		mat.translate(translation);
		mat.transpose();

		if (created)	// Stand-alone similarity transform file; not a model file
			r->setArrayWidth(4);
		double * array = new double[16];
		for (int i = 0; i < 16; i++)
			array[i] = mat.getData()[i];
		r->setDoubleArray(newStr, 16, array, "matrix"); 
	}
	else {
		r->setDoubleValue(newStr, getScale(), "scale");
		r->setDoubleValue(newStr, translation.getX(), "translation.x");
		r->setDoubleValue(newStr, translation.getY(), "translation.y");
		r->setDoubleValue(newStr, translation.getZ(), "translation.z");
		r->setDoubleValue(newStr, rotation.getX(), "rotation.x");
		r->setDoubleValue(newStr, rotation.getY(), "rotation.y");
		r->setDoubleValue(newStr, rotation.getZ(), "rotation.z");
		r->setDoubleValue(newStr, rotation.getW(), "rotation.w");
		if( tubePhi != 0.0 ) {
			r->setDoubleValue(newStr, tubePhi, "tubePhi");
		}
	}

	ret = true;
	if (filename != NULL) {
		try {
			if (! r->writeToFile(filename))
				cout << "Warning: similarity transformation file may be corrupted"
					<< endl;
		}
		catch (RException excp) {
			cout << excp.message() << endl;
			ret = false;
		}
	}

	if (created)
		delete (Registry *) registry;
	registry = NULL;
	return ret;
}

void SimilarityTransform3D::print() const
{
	cout << "Similarity transformation:\n";
	cout << "    Scale       = " << getScale() << '\n';
	cout << "    Translation = ";
	translation.print();
	cout << '\n';
	cout << "    Rotation    = ";
	rotation.print();
	cout << '\n';
	cout << "    Tube Phi    = " << getTubePhi() << '\n';
}





