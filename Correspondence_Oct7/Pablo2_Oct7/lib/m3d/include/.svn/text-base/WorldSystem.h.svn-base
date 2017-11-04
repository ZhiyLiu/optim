#ifndef WORLD_SYSTEM_H
#define WORLD_SYSTEM_H


#include "Vector3D.h"


/*  WorldSystem is a container class used for reading the
	world coordinate block of .m3d files and for doing some
	minimal processing using it.

    See also class Image3D.h and file M3DObjectFile.cpp.
*/


class WorldSystem
{
public:
	WorldSystem();
	virtual ~WorldSystem();
	void clear();
	void genConversions();
	bool status() const { return ok; }

	// Identity comparison of origin, bound and spacing
	bool operator==(const WorldSystem & w);
	bool equals(const Vector3D & orgn, const Vector3D & bnd,
		const Vector3D & spng);

	const Vector3D * getOrigin() const {
		return ok ? &origin : NULL;
	}
	const Vector3D * getBound() const {
		return ok ? &bound : NULL;
	}
	const Vector3D * getSpacing() const {
		return ok ? &spacing : NULL;
	}
	const std::string * getImagePath() const {
		return ok ? &imagePath : NULL;
	}
	const std::string * getImageTime() const {
		return ok ? &imageModTime : NULL;
	}

    // Coordinate conversion routines
    void modelToWorldCoordinates(Vector3D & coord);
    void worldToModelCoordinates(Vector3D & coord);

	void print();

	friend class M3DObjectFile;

private:
	WorldSystem & operator=(const WorldSystem & w);   // Not implemented
	WorldSystem(WorldSystem & w);       // Not implemented

	Vector3D origin;
	Vector3D bound;
	Vector3D spacing;
	Vector3D sign;
	std::string imagePath;
	std::string imageModTime;
	bool ok;
    double modelToWorldScale;
    double worldToModelScale;
    Vector3D originPixelPos;
};


#endif

