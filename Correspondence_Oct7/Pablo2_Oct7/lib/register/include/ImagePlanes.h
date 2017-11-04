#ifndef IMAGE_PLANES_H
#define IMAGE_PLANES_H

#include "Image3D.h"


/* Class ImagePlanes

    This class is used to store a pointer to an image and copies of
	displayable (intensity mapped) 2D image planes derived from that
    image.  The slices stored here may be raw or smoothed by trilinear
    interpolation.

    The image provided to this class at construction is considered
	to be owned by this class.  In particular, it will be deleted
	when the destructor is used.

    Access to the planes is gained by the get*CutPlane() functions.
	These are designed to be called repeatedly, as when a slider is
	moved with continuous updating.  If between successive calls, the
	coordinates or intensity range of a planes changes, then the next
	call to the get*CutPlane() function will return a coarse image
	quickly made without interpolation.  However, once the parameters
	defining the plane are unchanged for two successive calls to the
	get*CutPlane() function, then an interpolated version of it will
	be generated and returned.
*/


class ImagePlanes
{

public:

    ImagePlanes();
    ImagePlanes(Image3D * im, int size = 512);
    void setImagePtr(Image3D * im, int size = 512, bool external = false);

    ~ImagePlanes();

	void setSmoothing(bool smoothImages);

    Image3D * getImagePtr() const { return image; }
    int getImagePlaneSize() const { return imagePlaneSize; }

	// Get the current planes
    unsigned char * getXCutPlane(bool final);
    unsigned char * getYCutPlane(bool final);
    unsigned char * getZCutPlane(bool final);
    unsigned char * getArbCutPlane(bool final);

	// World-space positions of the image slices
    double getXCutPlaneWorldPos() const { return xCutPlaneWorldPos; }
    double getYCutPlaneWorldPos() const { return yCutPlaneWorldPos; }
    double getZCutPlaneWorldPos() const { return zCutPlaneWorldPos; }

	// Model-space positions of the current image planes
    double getXCutPlaneModelPos() const;
    double getYCutPlaneModelPos() const;
    double getZCutPlaneModelPos() const;

	bool getCutPlaneModelPositions(double & x, double & y,
        double & z) const;

	// Image-space ([0.0, 1.0]) coordinate (slice number) of current plane
    int getXCutPlaneIndex() const;
    int getYCutPlaneIndex() const;
    int getZCutPlaneIndex() const;

    Vector3D & getArbCutPlanePos() { return arbCutPlaneWorldPos; }
    Vector3D & getArbCutPlaneXDir() { return arbCutPlaneXDir; }
    Vector3D & getArbCutPlaneYDir() { return arbCutPlaneYDir; }

    double getArbCutPlaneExtent() const { return arbCutPlaneExtent; }

	// Change the current plane
	void setXCutPlaneSlice(int sliceIndex);
	void setYCutPlaneSlice(int sliceIndex);
	void setZCutPlaneSlice(int sliceIndex);

	// Change the current plane by specifying world coordinates
    void setXCutPlaneWorldPos(double pos) { xCutPlaneWorldPos = pos; }
    void setYCutPlaneWorldPos(double pos) { yCutPlaneWorldPos = pos; }
    void setZCutPlaneWorldPos(double pos) { zCutPlaneWorldPos = pos; }

    void setArbCutPlaneWorldPos(const Vector3D & pos);
    void setArbCutPlaneXDir(const Vector3D & xdir);
    void setArbCutPlaneYDir(const Vector3D & ydir);

    void setArbCutPlaneExtent(double extent);

	// Set the display intensities in the "color" map (lookup table)
	// used for filling the image planes.  Arguments min and max
	// should be in relative intensities, [0, 1].
    void setIntensityWindow(double min, double max);


private:

    Image3D * image;
	bool smooth;
	GreyValue lut_max;

    // Pixel size (1D) of an image plane (image planes are square)
    int imagePlaneSize;

    // Image planes data
    unsigned char * xCutPlane,
                  * yCutPlane,
                  * zCutPlane,
                  * arbCutPlane;

	// Direction of axes
	bool yFlip, zFlip;

    // World coordinate positions of cut planes
    double xCutPlaneWorldPos, yCutPlaneWorldPos, zCutPlaneWorldPos;

    // Position and directional vectors to define arbitrary plane
    Vector3D arbCutPlaneWorldPos;
    Vector3D arbCutPlaneXDir;
    Vector3D arbCutPlaneYDir;

    // Size of arbitrary plane in model coordinates (image is square)
    double arbCutPlaneExtent;

    unsigned char * colorLUT;

	void init(Image3D * im, int size, bool smoothness);

    // Functions to directly fill the pixels of the cut planes with image data
    void fastFillXCutPlane();
    void fastFillYCutPlane();
    void fastFillZCutPlane();
	void fastFillArbCutPlane();

    // Functions to fill the pixels of the cut planes with trilerped image data
    void fillXCutPlane();
    void fillYCutPlane();
    void fillZCutPlane();
    void fillArbCutPlane();

    // Memory allocate and clean up functions
    void allocatePlaneMemory();
    void freePlaneMemory();
};


#endif


