#include "WorldSystem.h"


using namespace std;


// The code code in this class must correspond to that in Image3D.cpp


WorldSystem::WorldSystem()
{
	clear();
}

void WorldSystem::genConversions()
{
	double longest;
	Vector3D zero;

	if (spacing == zero)
		cout << "Warning: using old-style model file, which lacks the voxel spacing"
			<< endl; 

	// Model (.m3d) files produced before March 16, 2007 never had
	// negative spacing values.
	if (spacing.getX() < 0.0) {
		spacing.setX(-spacing.getX());
		sign.setX(-1.0);
	}
	if (spacing.getY() < 0.0) {
		spacing.setY(-spacing.getY());
		sign.setY(-1.0);
	}
	if (spacing.getZ() < 0.0) {
		spacing.setZ(-spacing.getZ());
		sign.setZ(-1.0);
	}
	Vector3D shift(spacing.getX()/2.0, spacing.getY()/2.0, spacing.getZ()/2.0);
	originPixelPos = origin - shift;

    // The scaling from model to world coordinates is the maximum extent
    Vector3D extents = bound - origin;
	extents.abs();
    longest = extents.getX();
    if(extents.getY() > longest)
        longest = extents.getY();
    if(extents.getZ() > longest)
        longest = extents.getZ();
	modelToWorldScale = longest;

    // Scaling from world to model is just the inverse
    worldToModelScale = 1.0 / modelToWorldScale;
}

void WorldSystem::clear()
{
	origin.set(0.0, 0.0, 0.0);
	bound.set(0.0, 0.0, 0.0);
	spacing.set(0.0, 0.0, 0.0);
	sign.set(1.0, 1.0, 1.0);
	imagePath.erase();
	imageModTime.erase();
	ok = false;
}

WorldSystem::~WorldSystem()
{
	imagePath.erase();
	imageModTime.erase();
}

bool WorldSystem::operator==(const WorldSystem & w)
{
	if (origin != w.origin)
		return false;
	if (bound != w.bound)
		return false;
	if (spacing != w.spacing)
		return false;
	// Ignore imagePath, imageModTime, and ok
	return true;
}

bool WorldSystem::equals(const Vector3D & orgn, const Vector3D & bnd,
	const Vector3D & spng)
{
	if (origin != orgn)
		return false;
	if (bound != bnd)
		return false;
	if (spacing != spng)
		return false;
	// Ignore imagePath, imageModTime, and ok
	return true;
}

void WorldSystem::modelToWorldCoordinates(Vector3D & coord)
{
	Vector3D ext = sign;
    ext *= modelToWorldScale;
    coord = coord.vprod(ext);
    coord += originPixelPos;
}

void WorldSystem::worldToModelCoordinates(Vector3D & coord)
{
	Vector3D ext = sign;
    ext *= worldToModelScale;
    coord -= originPixelPos;
    coord = coord.vprod(ext);
}

void WorldSystem::print()
{
	cout << "World coordinate system\nOrigin: ";
	origin.print();
	cout << "Bound: ";
	bound.print();
	cout << "Voxel spacings: ";
	spacing.print();
}




