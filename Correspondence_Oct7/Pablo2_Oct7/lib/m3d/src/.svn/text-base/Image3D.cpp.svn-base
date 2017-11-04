#include <iostream>
#include "Image3D.h"
#include "utility.h"


#include <string.h>

using namespace std;


#ifdef WIN32
#if (_MSC_VER < 1400)

// VC++ 6.0 headers don't have max and min properly defined...
template <class T> const T& __max__ ( const T& a, const T& b ) {
    return (b<a)?a:b;
}

template <class T> const T& __min__ ( const T& a, const T& b ) {
    return (b<a)?b:a;
}

#ifndef MAX
#define MAX __max__
#endif
#ifndef MIN
#define MIN __min__
#endif

#else  /* _MSC_VER */

#ifndef MAX
#define MAX std::max
#endif
#ifndef MIN
#define MIN std::min
#endif

#endif  /* _MSC_VER */
#else  /* WIN32 */

#ifndef MAX
#define MAX std::max
#endif
#ifndef MIN
#define MIN std::min
#endif

#endif  /* WIN32 */


// Caution: making changes to the coordinate system mappings in this
// class can require comparable changes to be made in classes
// ImagePlanes, WorldSystem, and BYU.


const GreyValue OUT_OF_BOUNDS_VALUE = 0;
const GreyValue MIN_GREY_VALUE = 0;
const GreyValue MAX_GREY_VALUE = 65535;
const GreyValue MAX_STACKED_GREY_VALUE = 15;
const GreyValue MAX_GREY_VALUE_WRITTEN = MAX_GREY_VALUE/2;

const char * modalityStr[] = {
    "unknown", "CT", "SHIFTED_CT", "MRI", "DQF"
};

const char * modalityString(modality_t m)
{
    return modalityStr[m];
}


Image3D::Image3D()
{
    voxelArray = NULL;
    yIndexLUT = NULL;
    zIndexLUT = NULL;

    xDim = 0;
    yDim = 0;
    zDim = 0;
	nVoxels = 0;

    xSpacing = 0.0;
    ySpacing = 0.0;
    zSpacing = 0.0;

	originPixelPos = Vector3D(0.0, 0.0, 0.0);

    modelToWorldScale = 1.0;
    worldToModelScale = 1.0;

    intensityWinMin = 0.0;
    intensityWinMax = 1.0;
	intensityWinRange = 1.0;

	minIntens = 1;
	maxIntens = 0;
    intensRange = 0;

	scale_factor = 1.0;
	intens_shift = 0;

    boundaryType = EXTENDED_IMAGE_BOUNDARY;
    //= UNKNOWN_MODALITY;

    is_stacked = false;
    stackedMask = MAX_GREY_VALUE;
#ifdef BINARY
	imageIsStackedCount = 0;
#endif
	getVoxel = &Image3D::getRawVoxel;
}

Image3D::Image3D(int X, int Y, int Z)
{
    int index,
        increment,
        i;
	Vector3D origin(0.0, 0.0, 0.0);

    xDim = X;
    yDim = Y;
    zDim = Z;
	nVoxels = X*Y*Z;

    setSpacingAndOrigin(1.0, 1.0, 1.0, &origin);
    modelToWorldScale = 1.0;
    worldToModelScale = 1.0;

    intensityWinMin = 0.0;
    intensityWinMax = 1.0;
	intensityWinRange = 1.0;

	minIntens = 1;
	maxIntens = 0;
    intensRange = 0;

	scale_factor = 1.0;
	intens_shift = 0;

    yIndexLUT = new int[yDim];
    zIndexLUT = new int[zDim];

    index = 0;
    for(i = 0; i < yDim; i++) {
        yIndexLUT[i] = index;
        index += xDim;
    }

    index = 0;
    increment = xDim * yDim;
    for(i = 0; i < zDim; i++) {
        zIndexLUT[i] = index;
        index += increment;
    }

    voxelArray = new GreyValue[nVoxels];

    boundaryType = EXTENDED_IMAGE_BOUNDARY;
    //mode = UNKNOWN_MODALITY;

    is_stacked = false;
    stackedMask = MAX_GREY_VALUE;
#ifdef BINARY
	imageIsStackedCount = 0;
#endif
	getVoxel = &Image3D::getRawVoxel;
}

void Image3D::clear(GreyValue val, bool reset)
{
	for (int i = 0; i < nVoxels; i++)
		voxelArray[i] = val;

	if (! reset)
		return;
	scale_factor = 1.0;
	intens_shift = 0;
    //mode = UNKNOWN_MODALITY;

	minIntens = 1;
	maxIntens = 0;
    intensRange = 0;
}

Image3D::~Image3D()
{
    // "if (voxelArray != NULL)" is a step 
    // already performed by the delete operator

    delete [] voxelArray;
	delete [] zIndexLUT;
	delete [] yIndexLUT;
    
}

// Image modifiers:
void Image3D::extrudeAndPadInMinusZ(int slices, int padding)
{
	int xi,yi,zi, index;

	if( padding > 0 ) {
		int new_zDim	= zDim + padding;
		int new_nVoxels	= xDim*yDim*new_zDim;

		int *new_zIndexLUT	= new int[new_zDim];
		int increment	= xDim * yDim;
		memcpy( new_zIndexLUT, zIndexLUT, sizeof(zIndexLUT[0])*zDim );
		for( int i = zDim, index = zIndexLUT[zDim-1]; i < new_zDim; i++, index += increment ) {
			new_zIndexLUT[i]	= index;
		}

		GreyValue *new_voxelArray	= new GreyValue[new_nVoxels];
		memcpy( new_voxelArray + increment * padding, voxelArray,
			sizeof(GreyValue) * nVoxels );
		memset( new_voxelArray, 0, increment * padding * sizeof(GreyValue) );

		delete[] voxelArray;
		voxelArray	= new_voxelArray;
		delete[] zIndexLUT;
		zIndexLUT	= new_zIndexLUT;
		zDim		= new_zDim;
		nVoxels		= new_nVoxels;
		Vector3D new_origin	= getWorldOrigin();
		new_origin.setZ( new_origin.getZ() - spacing.getZ() * padding );
		setSpacingAndOrigin( spacing.getX(), spacing.getY(), spacing.getZ(), &new_origin );
	}

	if( slices > 0 ) {
		// Now find lowest z slice with any content
		int z_start	= -1;
		for( zi = 0; zi < zDim; ++zi ) {
			for( yi = 0; yi < yDim; ++yi ) {
				for( xi = 0; xi < xDim; ++xi ) {
					index	= xi + yIndexLUT[yi] + zIndexLUT[zi]; 
					if( voxelArray[index] != 0 ) {
						z_start	= zi;
						xi	= xDim;
						yi	= yDim;
						zi	= zDim;
						break;
					}
				}
			}
		}
		// If image is empty then we have nothing to do.
		if( z_start == -1 ) {
			return;
		}
	#define depth	5
		double cx[depth], cy[depth];
		memset(cx, 0, sizeof(double)*depth);
		memset(cy, 0, sizeof(double)*depth);
		int count;
		for( zi = z_start; zi < z_start + depth; ++zi ) {
			count	= 0;
			for( yi = 0; yi < yDim; ++yi ) {
				for( xi = 0; xi < xDim; ++xi ) {
					index	= xi + yIndexLUT[yi] + zIndexLUT[zi];
					if( voxelArray[index] ) {
						count++;
						cx[zi-z_start]	+= xi;
						cy[zi-z_start]	+= yi;
					}
				}
			}
			cx[zi-z_start]	/= count;
			cy[zi-z_start]	/= count;
		}
		//
		// Now find the regression lines:
		// cx = mx * z + ax and cy = my & z + ay
		//
		double sigma_z = 0.0, sigma_z2 = 0.0;
		double sigma_cx = 0.0, sigma_cy = 0.0;
		double sigma_z_cx = 0.0, sigma_z_cy = 0.0;
		for( zi = z_start; zi < z_start + depth; ++zi ) {
			sigma_z		+= zi;
			sigma_z2	+= zi*zi;
			sigma_cx	+= cx[zi-z_start];
			sigma_cy	+= cy[zi-z_start];
			sigma_z_cx	+= zi * cx[zi-z_start];
			sigma_z_cy	+= zi * cy[zi-z_start];
		}
		const double mx	= (depth * sigma_z_cx - sigma_z * sigma_cx ) / (depth * sigma_z2 - sigma_z * sigma_z );
		const double my	= (depth * sigma_z_cy - sigma_z * sigma_cy ) / (depth * sigma_z2 - sigma_z * sigma_z );
		const double ax	= (sigma_cx - mx * sigma_z) / depth;
		const double ay	= (sigma_cy - my * sigma_z) / depth;
		const int start_cx	= int(cx[0]);
		const int start_cy	= int(cy[0]);
		// Now extude last slice by keeping centroid
		// along the regression line.
		for( zi = z_start - 1; zi >= z_start - slices; --zi ) {
			const int cx	= int(mx * zi + ax);
			const int cy	= int(my * zi + ay);
			const int dx	= start_cx - cx;
			const int dy	= start_cy - cy;
			for( yi = MAX(0,dy); yi < yDim + MIN(0,dy); ++yi ) {
				for( xi = MAX(0,dx); xi < xDim + MIN(0,dx); ++xi ) {
					index	= xi + yIndexLUT[yi] + zIndexLUT[zi];
					const int src_index	= (xi+dx) + yIndexLUT[yi+dy] + zIndexLUT[z_start];
					voxelArray[index]	= voxelArray[src_index];
				}
			}
		}
	}

	return;
}

void Image3D::extrudeAndPadInPlusZ(int slices, int padding)
{
	int xi,yi,zi, index;

	if( padding > 0 ) {
		int new_zDim	= zDim + padding;
		int new_nVoxels	= xDim*yDim*new_zDim;

		int *new_zIndexLUT	= new int[new_zDim];
		int increment	= xDim * yDim;
		memcpy( new_zIndexLUT, zIndexLUT, sizeof(zIndexLUT[0])*zDim );
		for( int i = zDim, index = zIndexLUT[zDim-1]; i < new_zDim; i++, index += increment ) {
			new_zIndexLUT[i]	= index;
		}

		GreyValue *new_voxelArray	= new GreyValue[new_nVoxels];
		memcpy( new_voxelArray, voxelArray,
			sizeof(GreyValue) * nVoxels );
		memset( new_voxelArray + nVoxels, 0, increment * padding * sizeof(GreyValue) );

		delete[] voxelArray;
		voxelArray	= new_voxelArray;
		delete[] zIndexLUT;
		zIndexLUT	= new_zIndexLUT;
		zDim		= new_zDim;
		nVoxels		= new_nVoxels;
	}


	if( slices > 0 ) {
		// Now find highest z slice with any content
		int z_start	= -1;
		for( zi = zDim-1; zi >= 0; --zi ) {
			for( yi = 0; yi < yDim; ++yi ) {
				for( xi = 0; xi < xDim; ++xi ) {
					index	= xi + yIndexLUT[yi] + zIndexLUT[zi]; 
					if( voxelArray[index] != 0 ) {
						z_start	= zi;
						xi	= xDim;
						yi	= yDim;
						zi	= -1;
						break;
					}
				}
			}
		}
		// If image is empty then we have nothing to do.
		if( z_start == -1 ) {
			return;
		}
	#define depth	5
		double cx[depth], cy[depth];
		memset(cx, 0, sizeof(double)*depth);
		memset(cy, 0, sizeof(double)*depth);
		int count;
		for( zi = z_start; zi > z_start - depth; --zi ) {
			count	= 0;
			for( yi = 0; yi < yDim; ++yi ) {
				for( xi = 0; xi < xDim; ++xi ) {
					index	= xi + yIndexLUT[yi] + zIndexLUT[zi];
					if( voxelArray[index] ) {
						count++;
						cx[-zi+z_start]	+= xi;
						cy[-zi+z_start]	+= yi;
					}
				}
			}
			cx[-zi+z_start]	/= count;
			cy[-zi+z_start]	/= count;
		}
		//
		// Now find the regression lines:
		// cx = mx * z + ax and cy = my & z + ay
		//
		double sigma_z = 0.0, sigma_z2 = 0.0;
		double sigma_cx = 0.0, sigma_cy = 0.0;
		double sigma_z_cx = 0.0, sigma_z_cy = 0.0;
		for( zi = z_start; zi > z_start - depth; --zi ) {
			sigma_z		+= zi;
			sigma_z2	+= zi*zi;
			sigma_cx	+= cx[-zi+z_start];
			sigma_cy	+= cy[-zi+z_start];
			sigma_z_cx	+= zi * cx[-zi+z_start];
			sigma_z_cy	+= zi * cy[-zi+z_start];
		}
		const double mx	= (depth * sigma_z_cx - sigma_z * sigma_cx ) / (depth * sigma_z2 - sigma_z * sigma_z );
		const double my	= (depth * sigma_z_cy - sigma_z * sigma_cy ) / (depth * sigma_z2 - sigma_z * sigma_z );
		const double ax	= (sigma_cx - mx * sigma_z) / depth;
		const double ay	= (sigma_cy - my * sigma_z) / depth;
		const int start_cx	= int(cx[0]);
		const int start_cy	= int(cy[0]);
		// Now extude last slice by keeping centroid
		// along the regression line.
		for( zi = z_start+1; zi <= z_start + slices; ++zi ) {
			const int cx	= int(mx * zi + ax);
			const int cy	= int(my * zi + ay);
			const int dx	= start_cx - cx;
			const int dy	= start_cy - cy;
			for( yi = MAX(0,dy); yi < yDim + MIN(0,dy); ++yi ) {
				for( xi = MAX(0,dx); xi < xDim + MIN(0,dx); ++xi ) {
					index	= xi + yIndexLUT[yi] + zIndexLUT[zi];
					const int src_index	= (xi+dx) + yIndexLUT[yi+dy] + zIndexLUT[z_start];
					voxelArray[index]	= voxelArray[src_index];
				}
			}
		}
	}

	return;
}

// Model bound, the world coordinates of the far edge of the last voxel
Vector3D Image3D::getModelBound() const
{
	// See modelToWorldCoordinates().
	// Note that modelToWorldScale is not used here, so this will
	// not return the same as modelToWorldCoordinates(1.0, 1.0, 1.0).
	// (Model space usually extends beyond the bound in at least
	// one dimension).
    Vector3D bound(xExtent, yExtent, zExtent);
	bound = bound.vprod(sign);
    bound += originPixelPos;
	return bound;
}

// Real-world origin, at the center of the first voxel
Vector3D Image3D::getWorldOrigin() const
{
	Vector3D shift = 0.5*spacing;
	return originPixelPos + shift;
}

// Real-world bound, the center of the last voxel
Vector3D Image3D::getImageBound() const
{
	Vector3D diff((xDim - 0.5), (yDim - 0.5), (zDim - 0.5));
    diff = diff.vprod(spacing);
	return originPixelPos + diff;
}

double Image3D::getXWorldOrigin() const
{
	return originPixelPos.getX() + 0.5*spacing.getX();
}

double Image3D::getYWorldOrigin() const
{
	return originPixelPos.getY() + 0.5*spacing.getY();
}

double Image3D::getZWorldOrigin() const
{
	return originPixelPos.getZ() + 0.5*spacing.getZ();
}

double Image3D::getXImageBound() const
{
	return originPixelPos.getX() + (xDim - 0.5)*spacing.getX();
}

double Image3D::getYImageBound() const
{
	return originPixelPos.getY() + (yDim - 0.5)*spacing.getY();
}

double Image3D::getZImageBound() const
{
	return originPixelPos.getZ() + (zDim - 0.5)*spacing.getZ();
}

double Image3D::maxExtent() const
{
	double longest;

    longest = xExtent;
    if(yExtent > longest)
        longest = yExtent;
    if(zExtent > longest)
        longest = zExtent;
	return longest;
}

void Image3D::getSign(int * s) const
{
	s[0] = (int) sign.getX();
	s[1] = (int) sign.getY();
	s[2] = (int) sign.getZ();
}

void Image3D::modelToWorldCoordinates(Vector3D & coord) const
{
	// world_pt = model_origin + model_pt*(world_extent/model_extent)
	//          = originPixelPos + model_pt*signed_extents
	//          = originPixelPos + model_pt*modelToWorldScale*sign
	Vector3D ext = sign;
    ext *= modelToWorldScale;
    coord = coord.vprod(ext);
    coord += originPixelPos;
}

void Image3D::worldToModelCoordinates(Vector3D & coord) const
{
	// model_pt = (world_pt - model_origin)/(world_extent/model_extent)
	//          = (world_pt - originPixelPos)/signed_extents
	//          = (world_pt - originPixelPos)*worldToModelScale*sign
	Vector3D ext = sign;
    ext *= worldToModelScale;
    coord -= originPixelPos;
    coord = coord.vprod(ext);
}

// Compute fractional slice indexes from world coordinates
void Image3D::worldToImageCoordinates(Vector3D & coord) const
{
    coord.setX(coord.getX()*worldToImageScale.getX() + worldToImageOffset.getX());
    coord.setY(coord.getY()*worldToImageScale.getY() + worldToImageOffset.getY());
    coord.setZ(coord.getZ()*worldToImageScale.getZ() + worldToImageOffset.getZ());
}

// Compute nearest integer slice indexes from world coordinates
void Image3D::worldToImageCoordinates(const Vector3D & coord, int slices[3]) const
{
    double x, y, z;
    x = coord.getX()*worldToImageScale.getX() + worldToImageOffset.getX();
    y = coord.getY()*worldToImageScale.getY() + worldToImageOffset.getY();
    z = coord.getZ()*worldToImageScale.getZ() + worldToImageOffset.getZ();
    slices[0] = (int) (x + (x < 0.0 ? -0.5 : 0.5));
    slices[1] = (int) (y + (y < 0.0 ? -0.5 : 0.5));
    slices[2] = (int) (z + (z < 0.0 ? -0.5 : 0.5));
}

// Compute fractional slice indexes from model coordinates
void Image3D::modelToImageCoordinates(Vector3D & coord) const
{
    coord.setX(coord.getX()*modelToImageScale.getX() + modelToImageOffset);
    coord.setY(coord.getY()*modelToImageScale.getY() + modelToImageOffset);
    coord.setZ(coord.getZ()*modelToImageScale.getZ() + modelToImageOffset);
}

// Compute nearest integer slice indexes from model coordinates
void Image3D::modelToImageCoordinates(const Vector3D & coord, int slices[3]) const
{
    double x, y, z;
    x = coord.getX()*modelToImageScale.getX() + modelToImageOffset;
    y = coord.getY()*modelToImageScale.getY() + modelToImageOffset;
    z = coord.getZ()*modelToImageScale.getZ() + modelToImageOffset;
    slices[0] = (int) (x + (x < 0.0 ? -0.5 : 0.5));
    slices[1] = (int) (y + (y < 0.0 ? -0.5 : 0.5));
    slices[2] = (int) (z + (z < 0.0 ? -0.5 : 0.5));
}

double Image3D::imageToModelDistance(Vector3D & pt0, Vector3D & pt1) const
{
	Vector3D a = pt0;
	imageToModelCoordinates(a);
	Vector3D b = pt1;
	imageToModelCoordinates(b);
	a -= b;
	return a.norm();
}

double Image3D::imageToWorldDistance(Vector3D & pt0, Vector3D & pt1) const
{
	Vector3D a = pt0;
	imageToWorldCoordinates(a);
	Vector3D b = pt1;
	imageToWorldCoordinates(b);
	a -= b;
	return a.norm();
}

double Image3D::modelToWorldDistance(Vector3D & pt0, Vector3D & pt1) const
{
	Vector3D a(modelXToWorld(pt0.getX()), modelYToWorld(pt0.getY()),
		modelZToWorld(pt0.getZ()));
	Vector3D b(modelXToWorld(pt1.getX()), modelYToWorld(pt1.getY()),
		modelZToWorld(pt1.getZ()));
	a -= b;
	return a.norm();
}

double Image3D::modelToImageDistance(Vector3D & pt0, Vector3D & pt1) const
{
	Vector3D a = pt0.vprod(modelToImageScale);
	Vector3D b = pt1.vprod(modelToImageScale);
	a -= b;
	return a.norm();
}

double Image3D::worldToModelDistance(Vector3D & pt0, Vector3D & pt1) const
{
	Vector3D a(worldXToModel(pt0.getX()), worldYToModel(pt0.getY()),
		worldZToModel(pt0.getZ()));
	Vector3D b(worldXToModel(pt1.getX()), worldYToModel(pt1.getY()),
		worldZToModel(pt1.getZ()));
	a -= b;
	return a.norm();
}

double Image3D::worldToImageDistance(Vector3D & pt0, Vector3D & pt1) const
{
	Vector3D a = pt0.vprod(worldToImageScale);
	a += worldToImageOffset;
	Vector3D b = pt1.vprod(worldToImageScale);
	b += worldToImageOffset;
	a -= b;
	return a.norm();
}

bool Image3D::clipToImage(int slices[3]) const
{
	bool inside = true;

    if (slices[0] < 0) { slices[0] = 0; inside = false; }
    if (slices[0] >= xDim) { slices[0] = xDim - 1; inside = false; }
    if (slices[1] < 0) { slices[1] = 0; inside = false; }
    if (slices[1] >= yDim) { slices[1] = yDim - 1; inside = false; }
    if (slices[2] < 0) { slices[2] = 0; inside = false; }
    if (slices[2] >= zDim) { slices[2] = zDim - 1; inside = false; }

	return inside;
}

bool Image3D::clipToImage(Vector3D & coord) const
{
	bool inside = true;

    if (coord.getX() < 0.0) { coord.setX(0.0); inside = false; }
    if (coord.getX() > xDim - 1) { coord.setX(xDim - 1); inside = false; }
    if (coord.getY() < 0.0) { coord.setY(0.0); inside = false; }
    if (coord.getY() > yDim - 1) { coord.setY(yDim - 1); inside = false; }
    if (coord.getZ() < 0.0) { coord.setZ(0.0); inside = false; }
    if (coord.getZ() > zDim - 1) { coord.setZ(zDim - 1); inside = false; }

	return inside;
}

void Image3D::imageToModelCoordinates(Vector3D & coord) const
{
    coord.setX((coord.getX() - modelToImageOffset)/modelToImageScale.getX());
    coord.setY((coord.getY() - modelToImageOffset)/modelToImageScale.getY());
    coord.setZ((coord.getZ() - modelToImageOffset)/modelToImageScale.getZ());
}

void Image3D::imageToModelCoordinates(const int slices[3], Vector3D & coord) const
{
    coord.setX((slices[0] - modelToImageOffset)/modelToImageScale.getX());
    coord.setY((slices[1] - modelToImageOffset)/modelToImageScale.getY());
    coord.setZ((slices[2] - modelToImageOffset)/modelToImageScale.getZ());
}

void Image3D::imageToWorldCoordinates(Vector3D & coord) const
{
	coord.setX(spacing.getX()*(coord.getX() - modelToImageOffset) + originPixelPos.getX());
	coord.setY(spacing.getY()*(coord.getY() - modelToImageOffset) + originPixelPos.getY());
	coord.setZ(spacing.getZ()*(coord.getZ() - modelToImageOffset) + originPixelPos.getZ());
}

void Image3D::imageToWorldCoordinates(const int slices[3], Vector3D & coord) const
{
	coord.setX(spacing.getX()*(slices[0] - modelToImageOffset) + originPixelPos.getX());
	coord.setY(spacing.getY()*(slices[1] - modelToImageOffset) + originPixelPos.getY());
	coord.setZ(spacing.getZ()*(slices[2] - modelToImageOffset) + originPixelPos.getZ());
}

double Image3D::getXModelSpacing() const
{
    return 1.0 / modelToImageScale.getX();
}

double Image3D::getYModelSpacing() const
{
    return 1.0 / modelToImageScale.getY();
}

double Image3D::getZModelSpacing() const
{
    return 1.0 / modelToImageScale.getZ();
}

Vector3D Image3D::getModelSpacing() const
{
    return Vector3D(1.0, 1.0, 1.0).vdiv(modelToImageScale);
}

void Image3D::print(const char * filename) const
{
	cout << "\nImage Properties\n----------------\n";
	if (filename)
		cout << "File: " << filename << '\n';
	cout << "Dimensions: " << xDim << " x " << yDim << " x " << zDim << '\n';
	printWorld();
	cout << "Intensity range: [" << minIntens << ", " << maxIntens << "]\n";
	cout << "Actual (file) intensity range: [" << mapRelativeToActual(0) << ", "
		<< mapRelativeToActual(1) << "]\n";
	if (is_stacked)
		cout << "Image is stacked\n";
    if (modality() != UNKNOWN_MODALITY) {
        const char * s = modalityString();
        cout << "Modality is " << s << '\n';
    }
	cout << endl;
}

void Image3D::printWorld() const
{
	cout << "World coordinate system of image:\n  Origin: ";
	getWorldOrigin().print();
/*
    cout << "  Center: ";
    Vector3D center(0.5, 0.5, 0.5);
    modelToWorldCoordinates(center);
    center.print();
*/
	cout << "  Bound: ";
	getModelBound().print();
	cout << "  Voxel spacings: ";
	getSpacing().print();
    cout << "  Extents: " << getXExtent() << " x " << getYExtent()
        << " x " << getZExtent() << '\n';
}

// Note: this function does no flipping in Y or Z.  The dimensions
// must all be greater than zero, unless voxels is NULL.  However,
// this function can appear to cause a flip, because it deletes the
// Y and Z index lookup tables, into which a previous flip may have
// been encoded.  The final argument should normally not be used.
void Image3D::setVoxels(GreyValue * voxels, int numX, int numY, int numZ,
	bool ignore)
{
    int index, increment, i;

    if (voxelArray != NULL)
        delete [] voxelArray;

    if (yIndexLUT != NULL)
        delete [] yIndexLUT;

    if (zIndexLUT != NULL)
        delete [] zIndexLUT;

    voxelArray = voxels;
    xDim = numX;
    yDim = numY;
    zDim = numZ;

    if (voxelArray != NULL)
    {
		nVoxels = xDim*yDim*zDim;
		if (! ignore)
			calc_intensity_range(voxels, nVoxels, minIntens, maxIntens);
        intensRange = maxIntens - minIntens;	// AGG: Should this be in the samve block as the above stmt?

        yIndexLUT = new int[yDim];
        zIndexLUT = new int[zDim];

		index = 0;
		for (i = 0; i < yDim; i++) {
			yIndexLUT[i] = index;
			index += xDim;
		}

        index = 0;
        increment = xDim * yDim;
        for (i = 0; i < zDim; i++) {
            zIndexLUT[i] = index;
            index += increment;
        }
    }
    else {
		nVoxels = 0;
		minIntens = 1;
		maxIntens = 0;
        intensRange = 0;

        yIndexLUT = NULL;
        zIndexLUT = NULL;
    }
}

// This function replaces the voxels without changing the lookup tables.
// See setVoxels() above.
bool Image3D::replaceVoxels(GreyValue * voxels)
{
    if (voxelArray == NULL) {
		if (voxels != NULL)
			cout << "Error: attempting to replace nonexistent voxels" << endl;
		return false;
	}
	if (voxels == NULL) {
		cout << "Error: attempting to replace voxels from an empty array" << endl;
		return false;
	}

    delete [] voxelArray;
    voxelArray = voxels;

	calc_intensity_range(voxels, nVoxels, minIntens, maxIntens);
    intensRange = maxIntens - minIntens;
	return true;
}

void Image3D::setVoxel(int x, int y, int z, GreyValue val)
{
    if (voxelArray == NULL ||
		x < 0 || x >= xDim || y < 0 || y >= yDim || z < 0 || z >= zDim)
    {
        //cout << "Error: voxel values cannot be set off of the image\n";
        return;
    }
    voxelArray[x + yIndexLUT[y] + zIndexLUT[z]] = val;
}

void Image3D::adjustRange() {
    if (nVoxels == 0)
        return;
	calc_intensity_range(voxelArray, nVoxels, minIntens, maxIntens);
    intensRange = maxIntens - minIntens;
}

/*
	Function setSpacingAndOrigin() must be called whenever a new voxel
	array is provided, so that the relationship between image (voxel)
	coordinates and patient (world) coordinates, as well as to model
	space, can be established.

	The first three arguments are the signed widths of each voxel.  The
	signs indicate the directions of the corresponding axes in world
	coordinate space.  Thus, a negative value of xSpace indicates that to
	step from the first voxel of a row to the second voxel on that row,
	one must subtract |xSpace| from the current X world coordinate to get
	the new X world coordinate.  A negative sign for an axis indicates
	that the image and world coordinate axes along that dimension are
	in opposite directions.

	The origin argument should contain the world coordinates (centimeters)
	of the center of the first voxel of the image data.


	Explanation of Pablo's coordinate systems.

	We assume that the parameters provided in input files are correct.
	This means that the number of slices, distance between voxel centers,
	and world coordinate origin values are known.  It is essential to
	understand that the origin specified to this function must be the
	center of the first voxel in each of the three Cartesian directions.
	This definition of the origin corresponds to that in .raw3 files, but
	may differ from that used in other file formats.

	We define our model space using the longest axis in world space.
	Coordinate 0.0 is assigned to the left edge of the first voxel and 1.0
	to the right edge of the last voxel.  This simplifies computations
	and works best with OpenGL.  The model space definitions for the
	shorter axes are reduced accordingly on the assumption that the world
	coordinate system uses the same units (e.g. cm) for each axis.

	Image data is stored in the same order as it is in input files.
	Image space is defined so the indicies correspond to voxel centers.
	This means there is a half-voxel shift that must be applied to all
	conversions to image coordinates.  Most voxels can be thought of as a
	closed-to-open interval in each cardinal direction.  The first voxel,
	normally thought of as having image coordinate 0, is an exception,
	occupying a open-open interval in image space.  It can be considered
	to have all-floating point image coordinates in the range (-0.5, 0.5),
	which get rounded to 0 by the functions that convert to image
	coordinates.  The adjacent voxel has coordinates in the range
	[0.5, 1.5).  Integer image coordinates are always rounded up (down,
	if negative) when they are converted from floating point to integer.
	This means the outer edge of the last voxel on any axis is actually
	off of the image; converting from model or world coordinates to it
	will give an invalid image index.

	For example, voxel (1, 0) in the figure below occupies the world
	space interval [4.0625, 4.125) in the X-direction.  Converting any
	point in this range to image coordinates will give the range [0.5,
	1.5).  Converting to the integer image coordinate will give 1.  The
	final voxel in the X-direction below occupies world coordinate
	interval [4.4375, 4.5625).  Converting this range to image coordinates
	will give [3.5, 4.5), or rounded integer value 4.  Note that X world
	coordinate 4.5625 will convert to integer coordinate 5.  Likewise,
	converting model coordinate of 1.0 in the longest direction (X) to
	image coordinates, will produce a voxel index of 5, which is equal to
	the dimension (e.g. slice count) in that direction, an invalid image
	index.

	In order to allow negative indexes of in image space, standard
	rounding is used when converting to voxel indexes.  That is, negative
	indexes are rounded down.  Consequently, voxel 0 covers an open
	interval.  This leads to the unavoidable property that the model space
	origin is infinitesimally off the image and belongs to voxel
	(-1, -1, -1).  In practice, this should not matter, since off-image
	voxels don't really exist.

	The world coordinate origin stored inside this class (in variable
	originPixelPos) corresponds to the model space origin.  This is
	slightly different from the origin passed to this function, which is
	the world coordinate value of the center of the first voxel.  The
	difference is a shift of -0.5*(xSpacing, ySpacing, zSpacing).

	EXAMPLE:

	Here is a 2D illustration of these ideas using a 4 x 5 image with
	voxel spacings of 0.125 in each direction and (.m3d file) origin of
	(4.0, 2.0).

                                 X
                               ---->
                                                     Model   World   Image
      (3.9375,      0     1     2     3     4        =====   =====   =====
         1.9375) .-----------------------------. <--  0.0   1.9375   -0.5
                 |     |     |     |     |     |
             0   |  x  |     |     |     |  x  | <----------  2.0     0.0
                 |     |     |     |     |     |
                 +-----------------------------+
                 |     |     |     |     |     |
        |    1   |     |     |     |     |     |
        |        |     |     |     |     |     |
      Y |        +-----------------------------+ <--  0.4    2.1875   1.5
        |        |     |     |     |     |     |
        V    2   |     |     |     |     |     |
                 |     |     |     |     |     |
                 +-----------------------------+
                 |     |     |     |     |     |
             3   |  x  |     |     |     |  x  | <----------  2.375   3.0
                 |     |     |     |     |     |
                 '-----------------------------' <--  0.8 (NOTE: Not 1.0)

                 ^  ^  ^        ^           ^  ^
                 |  |  |        |           |  |
                 |  |  |        |           |  |
                    |  |                    |
     Model:    0.0             0.5             1.0
     World: 3.9375 4.0 4.0625  4.25        4.5 4.5625
     Image:   -0.5 0.0 0.5     2.0         4.0 4.5 (NOTE: Rounds to 5)


	For this image, the full-image X-extent is 0.625 and the Y-extent
	is 0.5 in world coordinates.  The origin stored in the class will be
	(3.9375, 1.9375).  Then the bottom right corner of the diagram above
	corresponds to world coordinate point (4.5625, 2.4375).

	The scale factors for conversions between the systems are:

		modelToWorldScale = 0.625
		worldToModelScale = 1.0/0.625 = 1.6

		modelToImageScale = (0.625/0.125, 0.625/0.125) = (5.0, 5.0)
		modelToImageOffset = (-0.5, -0.5, -0.5)

		worldToImageScale = (1/0.125, 1.0/0.125) = (8.0, 8.0)
		worldToImageOffset = (-4.0*8.0, -2.0*8.0) = (-32.0, -16.0)
	
	NOTICE: In Pablo before March 2007, images having a negative spacing
	for any axis were inverted along that axis when read.  Since this
	caused the slice coordinates for a given point in world space to differ
	from values used in clinical practice, the flipping was eliminated.
	This meant that old models could not be used without some compensating
	flip, and this flip is now done by reordering the contents of the
	lookup table for Y-slices, yIndexLUT, which is used in function
	getVoxels().  All code to effect this flip (see also ImagePlanes.cpp)
	is contained within #ifndef UNFLIPPED statements.  At UNC-CH, this
	variable is never defined, so that the flip occurs.  If a version of
	Pablo is used that does not do the flip, then it must be used to create
	new models (.m3d files), or else old ones must be flipped before they
	can be used in optimizations.  The flipping is permitted in Y and Z, but
	not X, because no lookup table is used for X in the getVoxel() function.
*/
void Image3D::setSpacingAndOrigin(double xSpace, double ySpace, double zSpace,
	const Vector3D * origin)
{
	double longest;
#ifndef UNFLIPPED
	int i;
	int index;
	int increment;
#endif

    spacing = Vector3D(xSpace, ySpace, zSpace);  // A vector; values may be < 0
    if (xSpace < 0.0) {

        xSpacing = -xSpace;
        sign.setX(-1.0);
		// AGG: Special case: must replace getVoxelValue() so it does "x = xDim -1 - x;"
		// before accessing voxelArray.  The ImagePlanes class will also have to be changed.
		cout << "Error: This program cannot be used on this image because it contains a\n";
		cout << "    negative voxel size in the X direction.  Doing this is presently\n";
		cout << "    unsupported.  While the program may appear to run, it will not yield\n";
		cout << "    correct results.  One possible alternative, is to negate the input\n";
		cout << "    file's current 3 voxel spacing values and then rerun.  The program\n";
		cout << "    will then operate without error, although the results will be\n";
		cout << "    inverted in all three cardinal directions.\n";
    }
    else {
        xSpacing = xSpace;
        sign.setX(1.0);
    }
    xExtent = xSpacing * (double) xDim;

    if (ySpace < 0.0) {
        ySpacing = -ySpace;
        sign.setY(-1.0);

#ifndef UNFLIPPED
		// Recompute the y-axis lookup table
		if (nVoxels > 0) {
			index = xDim*(yDim - 1);
			for (i = 0; i < yDim; i++) {
				yIndexLUT[i] = index;
				index -= xDim;
			}
		}
#endif
    }
    else {
        ySpacing = ySpace;
        sign.setY(1.0);
    }
    yExtent = ySpacing * (double) yDim;

    if (zSpace < 0.0) {
        zSpacing = -zSpace;
        sign.setZ(-1.0);

#ifndef UNFLIPPED
		// Recompute the z-axis lookup table
		if (nVoxels > 0) {
			index = xDim*yDim*(zDim - 1);
			increment = xDim * yDim;
			for (i = 0; i < zDim; i++) {
				zIndexLUT[i] = index;
				index -= increment;
			}
		}
#endif
    }
    else {
        zSpacing = zSpace;
        sign.setZ(1.0);
    }
    zExtent = zSpacing * (double) zDim;

	Vector3D shift(0.5*spacing);
	if (origin != NULL)
		originPixelPos = *origin - shift;
	else
		originPixelPos = -shift;

    // The scaling from model to world coordinates is the maximum extent
    longest = xExtent;
    if(yExtent > longest)
        longest = yExtent;
    if(zExtent > longest)
        longest = zExtent;
	modelToWorldScale = longest;	// Positive-valued; the sign is
    // taken into account in the conversion functions

    // Scaling from world to model is just the inverse
    worldToModelScale = 1.0 / modelToWorldScale;	// Positive-valued;
    // the sign is taken into account in the conversion functions

    // Scaling from model to image coordinates
    modelToImageScale.setX(longest / xSpacing);     // All positive-valued
    modelToImageScale.setY(longest / ySpacing);
    modelToImageScale.setZ(longest / zSpacing);

    // Scaling from world to image coordinates.
	// The values in worldToImageScale may be negative.
    worldToImageScale.set(Vector3D(1.0, 1.0, 1.0).vdiv(spacing));

    worldToImageOffset.setX(-getXWorldOrigin()*worldToImageScale.getX());
    worldToImageOffset.setY(-getYWorldOrigin()*worldToImageScale.getY());
    worldToImageOffset.setZ(-getZWorldOrigin()*worldToImageScale.getZ());
}

// Set all voxels >= threshold to value; otherwise set them to 0.
// It is assume that at least one voxel will be above threshold.
void Image3D::makeBinary(GreyValue threshold, GreyValue value)
{
	for (int i = 0; i < nVoxels; i++)
		voxelArray[i] = (voxelArray[i] >= threshold) ? value : (GreyValue) 0;
    minIntens = 0;
    maxIntens = value;
	intensRange = value;
	is_stacked = false;
	scale_factor = 1.0;
	intens_shift = 0;
    //mode = UNKNOWN_MODALITY;
}

GreyValue Image3D::getRawVoxel(int index) const
{
    return voxelArray[index];
}

#ifdef AE2_BUILD

GreyValue Image3D::getMaskedVoxel(int index) const
{
    return IMAGE_BITS & voxelArray[index];
}

// ConStruct must use masking when segmentation is done in it.
// When done in a backgrounded Pablo2 job, masking should not be
// used.  This function is called from ConStruct as needed.
void Image3D::setMasking(bool yesNo)
{
	if (yesNo)
		getVoxel = &Image3D::getMaskedVoxel;
	else
		getVoxel = &Image3D::getRawVoxel;
}

#endif	/* AE2_BUILD */

GreyValue Image3D::getVoxelValue(int x, int y, int z)
{
    if (voxelArray == NULL || yIndexLUT == NULL || zIndexLUT == NULL)
        return OUT_OF_BOUNDS_VALUE;

    if ((boundaryType == CONSTANT_IMAGE_BOUNDARY) && 
       (x < 0 || x >= xDim || y < 0 || y >= yDim || z < 0 || z >= zDim))
    {
		cout <<
			"Error: boundaryType is CONSTANT_IMAGE_BOUNDARY in call to Image3D::getVoxelValue\n"
			<< endl;
        return OUT_OF_BOUNDS_VALUE;
    }
    else if(boundaryType == EXTENDED_IMAGE_BOUNDARY)
    {
        if(x < 0)
            x = 0;
        else if(x >= xDim)
            x = xDim - 1;
        if(y < 0)
            y = 0;
        else if(y >= yDim)
            y = yDim - 1;
        if(z < 0)
            z = 0;
        else if(z >= zDim)
            z = zDim - 1;
    }

	// If the image is stacked binary, then mask the value returned
	// with stackedMask.
	if (is_stacked) 
	{
		int ind = x + yIndexLUT[y] + zIndexLUT[z];
		return (GreyValue) ((voxelArray[ind] & stackedMask) != 0);
	}

#ifdef AE2_BUILD
    // This is for use with AE2.  It cannot be used for builds of Pablo.
    // If the image is scaled on input (Pablo's default) or has a range
    // larger than 12 bits, this will cause loss of data.
    return IMAGE_BITS & voxelArray[x + yIndexLUT[y] + zIndexLUT[z]];
#else
    return voxelArray[x + yIndexLUT[y] + zIndexLUT[z]];
#endif
}

double Image3D::getInterpolatedVoxelValue(double x, double y, double z)
{
    int xIndex,
        yIndex,
        zIndex;

    double xt, yt, zt;

    double v[8];
    double val;

    if(voxelArray == NULL)
       return OUT_OF_BOUNDS_VALUE;

    xIndex = (int) x;
    yIndex = (int) y;
    zIndex = (int) z;

    xt = (double) x - (double) xIndex;
    yt = (double) y - (double) yIndex;
    zt = (double) z - (double) zIndex;

    v[0] = getVoxelValue(xIndex, yIndex, zIndex);
    v[1] = getVoxelValue(xIndex + 1, yIndex, zIndex);
    v[2] = getVoxelValue(xIndex, yIndex + 1, zIndex);
    v[3] = getVoxelValue(xIndex + 1, yIndex + 1, zIndex);
    v[4] = getVoxelValue(xIndex, yIndex, zIndex + 1);
    v[5] = getVoxelValue(xIndex + 1, yIndex, zIndex + 1);
    v[6] = getVoxelValue(xIndex, yIndex + 1, zIndex + 1);
    v[7] = getVoxelValue(xIndex + 1, yIndex + 1, zIndex + 1);

	// This works on voxel centers:  if xt, yt and zt are 0, no
	// interpolation occurs.
    val = (1 - zt) * ((1 - yt) * (v[0] * (1 - xt) + v[1] * xt)
        + yt * (v[2] * (1 - xt) + v[3] * xt))
        + zt * ((1 - yt) * (v[4] * (1 - xt) + v[5] * xt)
        + yt * (v[6] * (1 - xt) + v[7] * xt));

    return val;
}

GreyValue Image3D::getWindowedVoxelValue(int x, int y, int z)
{
    if(voxelArray == NULL || yIndexLUT == NULL || zIndexLUT == NULL)
        return OUT_OF_BOUNDS_VALUE;

	// Variant of mapDisplayToRelative()
	double initVal = (getVoxelValue(x, y, z) - minIntens)/intensRange;
    double finalVal;

    if (initVal <= intensityWinMin)
        finalVal = 0.0;
    else if (initVal >= intensityWinMax)
        finalVal = 1.0;
    else
        finalVal = (initVal - intensityWinMin)/intensityWinRange;

    return mapRelativeToDisplay(finalVal);
}

double Image3D::getWindowedInterpolatedVoxelValue(double x, double y, double z)
{
    int xIndex,
        yIndex,
        zIndex;

    double xt, yt, zt;

    double v[8];
    double val;

    if(voxelArray == NULL)
       return OUT_OF_BOUNDS_VALUE;

    xIndex = (int) x;
    yIndex = (int) y;
    zIndex = (int) z;

    xt = (double) x - (double) xIndex;
    yt = (double) y - (double) yIndex;
    zt = (double) z - (double) zIndex;

    v[0] = getVoxelValue(xIndex, yIndex, zIndex);
    v[1] = getVoxelValue(xIndex + 1, yIndex, zIndex);
    v[2] = getVoxelValue(xIndex, yIndex + 1, zIndex);
    v[3] = getVoxelValue(xIndex + 1, yIndex + 1, zIndex);
    v[4] = getVoxelValue(xIndex, yIndex, zIndex + 1);
    v[5] = getVoxelValue(xIndex + 1, yIndex, zIndex + 1);
    v[6] = getVoxelValue(xIndex, yIndex + 1, zIndex + 1);
    v[7] = getVoxelValue(xIndex + 1, yIndex + 1, zIndex + 1);

	// This works on voxel centers:  if xt, yt and zt are 0, no
	// interpolation occurs.
    val = (1 - zt) * ((1 - yt) * (v[0] * (1 - xt) + v[1] * xt)
        + yt * (v[2] * (1 - xt) + v[3] * xt))
        + zt * ((1 - yt) * (v[4] * (1 - xt) + v[5] * xt)
        + yt * (v[6] * (1 - xt) + v[7] * xt));

	// Variant of mapDisplayToRelative()
	double initVal = (val - minIntens)/intensRange;
    double finalVal;

    if (initVal <= intensityWinMin)
        finalVal = 0.0;
    else if (initVal >= intensityWinMax)
        finalVal = 1.0;
    else
        finalVal = (initVal - intensityWinMin)/intensityWinRange;

	return minIntens + finalVal*intensRange;	// Variant of mapRelativeToDisplay()
}

#ifdef BINARY

void Image3D::pushImageIsStacked(bool isStacked, GreyValue mask)
{
	if (imageIsStackedCount == sizeof(imageIsStackedMaskHistory)/sizeof(imageIsStackedMaskHistory[0])) {
		std::cout << "Error: Image stack is full" << std::endl;
		return;
	}
	// Save current state onto stack
	//imageIsStackedIsStackedHistory[imageIsStackedCount] = getIsImageStacked();
	//imageIsStackedMaskHistory[imageIsStackedCount] = getStackedMask();

	// Set current values to passed values
	setIsImageStacked(isStacked);
	setStackedMask(mask);

	imageIsStackedCount++;
}

void Image3D::popImageIsStacked()
{
	if (imageIsStackedCount == 0) {
		std::cout << "Error: Image stack is empty" << std::endl;
		return;
	}
	imageIsStackedCount--;
	setIsImageStacked(imageIsStackedIsStackedHistory[imageIsStackedCount]);
	setStackedMask(imageIsStackedMaskHistory[imageIsStackedCount]);
}

#endif


/*	Test of coordinate systems.
	To use this, call it from pablo.cpp.

	Coordinate system conversions used in the tests are shown below for two
	cases.

  Case 1.  image spacing = (0.12, -0.125, 0.13).

	X.  Model:  0.0         0.2         0.4   0.5   0.6         0.8         1.0
		World:  3.94  4.0         4.12        4.24        4.36       4.48   4.54
		         |-----+-----|-----+-----|-----+-----|-----+-----|-----+-----|
		Image: -0.5    0    0.5    1    1.5    2    2.5    3    3.5    4    4.5

	Y.  Model:  0.0    0.10416_ 0.20833_       0.4166_ 0.52083_ 0.625           0.8333_
		World:  2.0625  2.0            1.875            1.75            1.625   1.5625
		         |-------+-------|-------+-------|-------+-------|-------+-------|
		Image: -0.5      0      0.5      1      1.5      2      2.5      3      3.5

	Z.  Model:  0.0     0.1083_ 0.2166          0.4333  0.5416  0.65
		World: -0.065   0.0             0.13            0.26    0.325
		         |-------+-------|-------+-------|-------+-------|
		Image: -0.5      0      0.5      1      1.5      2      2.5

  Case 2.  image spacing = (0.125, -0.125, 0.125).

	X.  Model:  0.0         0.2         0.4   0.5   0.6         0.8         1.0
		World: 3.9375 4.0         4.125       4.25        4.375       4.5   4.5625
		         |-----+-----|-----+-----|-----+-----|-----+-----|-----+-----|
		Image: -0.5    0    0.5    1    1.5    2    2.5    3    3.5    4    4.5

	Y.  Model:  0.0           0.2           0.4    0.5    0.6           0.8
		World:  2.0625 2.0           1.875         1.75         1.625   1.5625
		         |------+------|------+------|------+------|------+------|
		Image: -0.5     0     0.5     1     1.5     2     2.5     3     3.5

	Z.  Model:  0.0           0.2           0.4    0.5    0.6
		World: -0.0625 0.0           0.125         0.25   0.3125
		         |------+------|------+------|------+------|
		Image: -0.5     0     0.5     1     1.5     2     2.5
*/
/*
void testImage3D()
{
	int xDim, yDim, zDim;
	double xSpacing, ySpacing, zSpacing;
	Vector3D zero(0.0, 0.0, 0.0);
	Vector3D underHalf(0.499995, 0.499995, 0.499995);
	Vector3D half(0.5, 0.5, 0.5);
	Vector3D overHalf(0.500005, 0.500005, 0.500005);
	Vector3D one(1.0, 1.0, 1.0);
	Vector3D v, w;
	int slices[3];
	double len;

	// Parameters
	xDim = 5;
	yDim = 4;
	zDim = 3;
	xSpacing = 0.12;
	ySpacing = -0.125;
	zSpacing = 0.13;
//	xSpacing = 0.125;
//	ySpacing = -0.125;
//	zSpacing = 0.125;
	Vector3D origin(4.0, 2.0, 0.0);

	// Initialization
	Image3D im(xDim, yDim, zDim);
	im.clear();
	im.setSpacingAndOrigin(xSpacing, ySpacing, zSpacing, &origin);

	// Report
	cout << "World origin: (" << im.getXWorldOrigin() << ", " << im.getYWorldOrigin()
		<< ", " << im.getZWorldOrigin() << ")\n";
	cout << "Dimensions: (" << im.getXDim() << ", " << im.getYDim()
		<< ", " << im.getZDim() << ")\n";
	cout << "Spacings (cm): (" << im.getXSpacing() << ", "
		<< im.getYSpacing() << ", " << im.getZSpacing() << ")\n";
	cout << "Spacings (model-space): ";
    im.getModelSpacing().print();
	cout << "Extents (cm): (" << im.getXExtent() << ", " << im.getYExtent()
		<< ", " << im.getZExtent() << ")\n";
	cout << "Model extents: (" << im.getXModelExtent() << ", " << im.getYModelExtent()
		<< ", " << im.getZModelExtent() << ")\n";
	cout << "Model origin: (" << im.getXModelOrigin() << ", "
		<< im.getYModelOrigin() << ", " << im.getZModelOrigin() << ")\n";
	cout << "Center of the last voxel: ";
	im.getImageBound().print();
	cout << "Model bound in world coordinates: ";
	im.getModelBound().print();

	// Test 0
	cout << "\nTest 0\n";
	w = zero;
	im.modelToImageCoordinates(w);
	cout << "Model (0.0, 0.0, 0.0) is at image coordinates ";
	w.print();

	w = zero;
	w.set(im.imageXToWorld((int) w.getX()), im.imageYToWorld((int) w.getY()),
		im.imageZToWorld((int) w.getZ()));	// Note: not really truncating
	cout << "These image coordinates correspond to world coordinates ";
	w.print();

	w = zero;
	im.modelToImageCoordinates(w, slices);
	cout << "These model coordinates correspond to slices ";
	cout << '(' << slices[0] << ", " << slices[1] << ", " << slices[2]
		<< ")\n";

	v = zero;
	im.modelToWorldCoordinates(v);
	cout << "Model (0.0, 0.0, 0.0) is at world coordinates ";
	v.print();

	w = v;
	im.worldToModelCoordinates(v);
	cout << "These world coordinates invert to model coordinates ";
	v.print();

	im.worldToImageCoordinates(w);
	cout << "These world coordinates correspond to image coordinates ";
	w.print();

	v = w;
    im.imageToWorldCoordinates(w);
	cout << "These image coordinates correspond to world coordinates\n\t";
	w.print();

	im.imageToModelCoordinates(v);
	cout << "These image coordinates correspond to model coordinates ";
	v.print();
	cout << '\n';

	// Before-image test
	cout << "Before-image test\n";
	w = Vector3D(-0.0001, -0.0001, -0.0001);
	v = w;
	cout << "Model (" << w.getX() << ", " <<  w.getY() << ", "
		<<  w.getZ() << ") is at image coordinates\n\t";
	im.modelToImageCoordinates(w);
	w.print();

	cout << "These model coordinates correspond to image coordinates ";
	im.modelToImageCoordinates(v, slices);
	cout << '(' << slices[0] << ", " << slices[1] << ", " << slices[2]
		<< ")\n";

	w = Vector3D(-1.0, -1.0, -1.0);
	cout << "Image (" << w.getX() << ", " <<  w.getY() << ", "
		<<  w.getZ() << ") is at model coordinates\n\t";
	im.imageToModelCoordinates(w);
	w.print();
	v = w;

	cout << "These model coordinates correspond to image coordinates\n\t";
	im.modelToImageCoordinates(w);
	w.print();

	cout << "These model coordinates correspond to image coordinates\n\t";
	im.modelToImageCoordinates(v, slices);
	cout << '(' << slices[0] << ", " << slices[1] << ", " << slices[2]
		<< ")\n";

	w = Vector3D(-0.5005, -0.50048, -0.50046153846);
	v = w;
	cout << "Image (" << w.getX() << ", " <<  w.getY() << ", "
		<<  w.getZ() << ") is at model coordinates\n\t";
	im.imageToModelCoordinates(w);
	w.print();

	cout << "These image coordinates correspond to world coordinates\n\t";
	im.imageToWorldCoordinates(v);
	v.print();

	cout << "These world coordinates correspond to model coordinates\n\t";
	im.worldToModelCoordinates(v);
	v.print();

	// Origin test
	cout << "\nOrigin test\n";
	w = Vector3D(0.2, 0.8, 1.2);
	cout << "Image coordinates (" << w.getX() << ", " << w.getY() << ", "
		<< w.getZ() << ")\n\tare at world coordinates (";
	im.imageToWorldCoordinates(w);
	cout << w.getX() << ", " << w.getY() << ", " << w.getZ() << ")\n";

	im.worldToImageCoordinates(w, slices);
	cout << "These world coordinates invert to model coordinates ";
	cout << '(' << slices[0] << ", " << slices[1] << ", " << slices[2]
		<< ")\n";

	w = Vector3D(-0.2, -0.8, -1.2);
	cout << "Image coordinates (" << w.getX() << ", " << w.getY() << ", "
		<< w.getZ() << ")\n\tare at world coordinates (";
	im.imageToWorldCoordinates(w);
	cout << w.getX() << ", " << w.getY() << ", " << w.getZ() << ")\n";

	im.worldToImageCoordinates(w, slices);
	cout << "These world coordinates invert to model coordinates ";
	cout << '(' << slices[0] << ", " << slices[1] << ", " << slices[2]
		<< ")\n";

	w = Vector3D(-0.6, -0.5, -0.4);
	cout << "Image coordinates (" << w.getX() << ", " << w.getY() << ", "
		<< w.getZ() << ")\n\tare at world coordinates (";
	im.imageToWorldCoordinates(w);
	cout << w.getX() << ", " << w.getY() << ", " << w.getZ() << ")\n";

	im.worldToImageCoordinates(w, slices);
	cout << "These world coordinates invert to model coordinates ";
	cout << '(' << slices[0] << ", " << slices[1] << ", " << slices[2]
		<< ")\n";

	// Test 1st voxel center
	cout << "\nTest 1st voxel center\n";
	cout << "Center of first voxel is at world coordinates ";
	origin.print();
	w = origin;
	im.worldToModelCoordinates(w);
	cout << "This corresponds to model coordinates ";
	w.print();

	v = w;
	im.modelToImageCoordinates(w);
	cout << "These model coordinates are at image coordinates ";
	w.print();

	w = v;
	im.modelToImageCoordinates(w, slices);
	cout << "These model coordinates correspond to slices ";
	cout << '(' << slices[0] << ", " << slices[1] << ", " << slices[2]
		<< ")\n";

	im.modelToWorldCoordinates(v);
	cout << "These model coordinates are at world coordinates ";
	v.print();

	w = v;
	im.worldToModelCoordinates(v);
	cout << "These world coordinates invert to model coordinates ";
	v.print();

	im.worldToImageCoordinates(w);
	cout << "These world coordinates correspond to image coordinates ";
	w.print();

	v = w;
    im.imageToWorldCoordinates(w);
	cout << "These image coordinates correspond to world coordinates ";
	w.print();

	im.imageToModelCoordinates(v);
	cout << "These image coordinates correspond to model coordinates ";
	v.print();
	cout << '\n';

	// Test 0.5
	cout << "Test 0.5\n";
	w = half;
	im.modelToImageCoordinates(w);
	cout << "Model (0.5, 0.5, 0.5) is at image coordinates ";
	w.print();

	w = half;
	im.modelToImageCoordinates(w, slices);
	cout << "These model coordinates correspond to slices ";
	cout << '(' << slices[0] << ", " << slices[1] << ", " << slices[2]
		<< ")\n";

	v = half;
	im.modelToWorldCoordinates(v);
	cout << "Model (0.5, 0.5, 0.5) is at world coordinates ";
	v.print();

	w = v;
	im.worldToModelCoordinates(v);
	cout << "These world coordinates invert to model coordinates ";
	v.print();

	im.worldToImageCoordinates(w);
	cout << "These world coordinates correspond to image coordinates ";
	w.print();

	v = w;
    im.imageToWorldCoordinates(w);
	cout << "These image coordinates correspond to world coordinates ";
	w.print();

	im.imageToModelCoordinates(v);
	cout << "These image coordinates correspond to model coordinates ";
	v.print();
	cout << '\n';

	// Test < 0.5
	cout << "Test < 0.5\n";
	w = underHalf;
	im.modelToImageCoordinates(w);
	cout << "Model (0.499995, 0.499995, 0.499995) is at image coordinates\n\t";
	w.print();

	w = underHalf;
	im.modelToImageCoordinates(w, slices);
	cout << "These model coordinates correspond to slices ";
	cout << '(' << slices[0] << ", " << slices[1] << ", " << slices[2] << ")\n";

	v = underHalf;
	im.modelToWorldCoordinates(v);
	cout << "Model (0.499995, 0.499995, 0.499995) is at world coordinates\n\t`";
	v.print();

	w = v;
	im.worldToModelCoordinates(v);
	cout << "These world coordinates invert to model coordinates\n\t";
	v.print();

	im.worldToImageCoordinates(w);
	cout << "These world coordinates correspond to image coordinates\n\t";
	w.print();

	v = w;
    im.imageToWorldCoordinates(w);
	cout << "These image coordinates correspond to world coordinates\n\t";
	w.print();

	im.imageToModelCoordinates(v);
	cout << "These image coordinates correspond to model coordinates\n\t";
	v.print();
	cout << '\n';

	// Test > 0.5
	cout << "Test > 0.5\n";
	w = overHalf;
	im.modelToImageCoordinates(w);
	cout << "Model (0.500005, 0.500005, 0.500005) is at image coordinates\n\t";
	w.print();

	w = overHalf;
	im.modelToImageCoordinates(w, slices);
	cout << "These model coordinates correspond to slices ";
	cout << '(' << slices[0] << ", " << slices[1] << ", " << slices[2]
		<< ")\n";

	v = overHalf;
	im.modelToWorldCoordinates(v);
	cout << "Model (0.500005, 0.500005, 0.500005) is at world coordinates\n\t";
	v.print();

	w = v;
	im.worldToModelCoordinates(v);
	cout << "These world coordinates invert to model coordinates\n\t";
	v.print();

	im.worldToImageCoordinates(w);
	cout << "These world coordinates correspond to image coordinates\n\t";
	w.print();

	v = w;
    im.imageToWorldCoordinates(w);
	cout << "These image coordinates correspond to world coordinates\n\t";
	w.print();

	im.imageToModelCoordinates(v);
	cout << "These image coordinates correspond to model coordinates\n\t";
	v.print();
	cout << '\n';

	// Test 1
	cout << "Test 1\n";
	v = Vector3D(im.getXModelOrigin() + im.getXExtentVec(),
		im.getYModelOrigin() + im.getYExtentVec(),
		im.getZModelOrigin() + im.getZExtentVec());
	cout << "The world coordinate limit (ends of the last voxel) is\n\tat ";
	v.print();

	w = v;
	im.worldToImageCoordinates(w);
	cout << "This corresponds to image coordinates ";
	w.print();

	im.worldToModelCoordinates(v);
	cout << "The world coordinate limit corresponds to model coordinates ";
	v.print();
	w = v;

	im.modelToWorldCoordinates(v);
	cout << "These model coordinates invert to world coordinates ";
	v.print();

	im.modelToImageCoordinates(w, slices);
	cout << "These model coordinates correspond to slices ";
	cout << '(' << slices[0] << ", " << slices[1] << ", " << slices[2]
		<< ")\n";

	im.modelToImageCoordinates(w);
	cout << "These model coordinates correspond to image coordinates ";
	w.print();

	w.set(im.getXDim() - 1, im.getYDim() - 1, im.getZDim() - 1);
	cout << "The last image voxel has image coordinates ";
	w.print();

	w.set(im.imageXToWorld((int) w.getX()), im.imageYToWorld((int) w.getY()),
		im.imageZToWorld((int) w.getZ()));	// Note: not really truncating
	cout << "It has world coordinates ";
	w.print();

	// Cardinal distance tests
	cout << "\nCardinal distance tests\n";
	cout << "The world extents of (" << im.getXExtent() << ", "
		<< im.getYExtent() << ", " << im.getZExtent()
		<< ") have cardinal distances in image\n\tspace of: ";
	v = Vector3D(im.worldXToImageDistance(im.getXExtent()),
		im.worldYToImageDistance(im.getYExtent()), 
		im.worldZToImageDistance(im.getZExtent()));
	v.print();

	cout << "These image distances convert to world space distances as: ";
	w = Vector3D(im.imageXToWorldDistance(v.getX()),
		im.imageYToWorldDistance(v.getY()), 
		im.imageZToWorldDistance(v.getZ()));
	w.print();

	cout << "These image distances convert to model space distances as: ";
	w = Vector3D(im.imageXToModelDistance(v.getX()),
		im.imageYToModelDistance(v.getY()), 
		im.imageZToModelDistance(v.getZ()));
	w.print();

	cout << "The world extents in model space are: ";
	v = Vector3D(im.worldToModelDistance(im.getXExtent()),
		im.worldToModelDistance(im.getYExtent()), 
		im.worldToModelDistance(im.getZExtent()));
	v.print();

	cout << "These model distances convert to image space distances as: ";
	w = Vector3D(im.modelXToImageDistance(v.getX()),
		im.modelYToImageDistance(v.getY()), 
		im.modelZToImageDistance(v.getZ()));
	w.print();

	cout << "These model distances convert to world space distances as: ";
	w = Vector3D(im.modelToWorldDistance(v.getX()),
		im.modelToWorldDistance(v.getY()), 
		im.modelToWorldDistance(v.getZ()));
	w.print();

	// General distance tests
	cout << "\nGeneral distance tests\n";
	len = im.modelToWorldDistance(zero, half);
	cout << "The (0, 0, 0) to (0.5, 0.5, 0.5) model-space diagonal has a length of\n\t"
		<< len << " cm\n";
	len = im.modelToImageDistance(zero, half);
	cout << "The same diagonal has a length of " << len << " voxels\n";
	w = half;
	im.modelToWorldCoordinates(w);
	len = im.worldToImageDistance(im.getModelOrigin(), w);
	cout << "The distance between corresponding world coordinates is "
		<< len << " voxels\n";

	v = im.getModelOrigin();
	len = im.worldToModelDistance(w, v);
	cout << "The diagonal world distance corresponds to a length of " << len
		<< " in model space\n";

	v -= w;
	cout << "The same diagonal has a length in world space of " << v.norm() << " cm\n";

	v = zero;
	im.modelToImageCoordinates(v);
	w = half;
	im.modelToImageCoordinates(w);
	len = im.imageToModelDistance(v, w);
	cout << "The distance between the corresponding image points in model space is\n\t" << len << '\n';
	len = im.imageToWorldDistance(v, w);
	cout << "The same distance in world space is\n\t" << len << '\n';

	cout << "\nImage point (1.0, 1.0, 1.0) is at model (" << im.imageXToModel(1.0)
		<< ", " << im.imageYToModel(1.0) << ", " << im.imageZToModel(1.0) << ")\n";
	v = Vector3D(2.50005, 0.5, 1.499995);
	w = v;
	im.imageToModelCoordinates(w);
	cout << "Image point (" << v.getX() << ", " << v.getY() << ", " << v.getZ()
		<< ") is at model (" << w.getX() << ", " << w.getY() << ", " << w.getZ()
		<< ")\n";
	len = im.imageToModelDistance(one, v);
	cout << "The distance between these points is " << len << " (model)\n";
	v = one;
	im.imageToModelCoordinates(v);
	cout << "This distance should be "
		<< sqrt((v.getX() - w.getX())*(v.getX() - w.getX())
		+ (v.getY() - w.getY())*(v.getY() - w.getY())
		+ (v.getZ() - w.getZ())*(v.getZ() - w.getZ())) << endl;

	cout << "Image point (1.0, 1.0, 1.0) is at world (" << im.imageXToWorld(1.0)
		<< ", " << im.imageYToWorld(1.0) << ", " << im.imageZToWorld(1.0) << ")\n";
	v = Vector3D(2.50005, 0.5, 1.499995);
	w = v;
	im.imageToWorldCoordinates(w);
	cout << "Image point (" << v.getX() << ", " << v.getY() << ", " << v.getZ()
		<< ") is at world (" << w.getX() << ", " << w.getY() << ", " << w.getZ()
		<< ")\n";
	len = im.imageToWorldDistance(one, v);
	cout << "The distance between these points is " << len << " (world)\n";
	v = one;
	im.imageToWorldCoordinates(v);
	cout << "This distance should be "
		<< sqrt((v.getX() - w.getX())*(v.getX() - w.getX())
		+ (v.getY() - w.getY())*(v.getY() - w.getY())
		+ (v.getZ() - w.getZ())*(v.getZ() - w.getZ())) << endl;

}
*/

