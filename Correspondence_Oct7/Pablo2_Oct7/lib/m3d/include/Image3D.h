#ifndef IMAGE3D_H
#define IMAGE3D_H

#include <stdio.h>
#include <math.h>
#include "Vector3D.h"


/*
    Class Image3D stores voxels and handles coordinate transformations
    between three 3D spaces, image, model and world.

    The relationship between image (voxel) and world (physical) coordinates
    (in cm in all three directions) is specified in the image files.  Model
    space maps the world coordinates uniformly so that the longest axis
    (largest extent) has length 1.0.  The other two axes must then be
    contained within the model-space range [0.0, 1.0].  Image space counts
    voxel centers beginning with zero at the first voxel.  Unlike the other
    two spaces, the origin of world coordinate space may begin anywhere.
    The model and image space axes point in corresponding directions; the
    world coordinate axes will parallel the other axes, but may arbitrarily
    be reversed.  In particular, the image and model spaces are right-handed,
    while world space may be left-handed or right-handed.

    This class takes into account the fact that medical scanners report
    the world coordinates of voxels at their centers.  This is reflected
    in the two sets of functions provided below for returning the origin.
    This also means that getInterpolatedVoxelValue() called with integer
    slice indexes will return the intensity of a single voxel, instead of
    a value interpolated from several voxels.  Additional details about
    the coordinate systems are provided in file Image3D.cpp.

    The intensity range of the pixels stored by this class does not change
    after construction, and so it is not recorded in the class.  For
    distance maps, the range is [0, MAX_GREY_VALUE_WRITTEN], although it is
    recorded in the output file to give an automatic windowing on input as
    [MIN_GREY_VALUE, MAX_GREY_VALUE].  Images exported from a model are in
    the range [0, MAX_GREY_VALUE_WRITTEN].  For all other images, the range
    is [MIN_GREY_VALUE, MAX_GREY_VALUE].  See file RAWImageFile.cpp and
    functions exportToImage() and exportDistanceMap() in P3DControl.cpp.
    Also see class AllImageIO.

    Stacked images are created by combining multiple images into different
    bits of the image.  Thus each bit may represent a different figure.  Such
    images are used for computing a multi-object image match with a
    different image slice for each figure.

    See also class WorldSystem.
*/


typedef unsigned short GreyValue;

extern const GreyValue MIN_GREY_VALUE;
extern const GreyValue MAX_GREY_VALUE;
extern const GreyValue MAX_GREY_VALUE_WRITTEN;
extern const GreyValue MAX_STACKED_GREY_VALUE;

#ifdef AE2_BUILD
#ifndef IMAGE_BITS
#define IMAGE_BITS  0xfff
#endif
#endif

// This tells how to handle out-of-bounds image values
enum ImageBoundaryType
{
    CONSTANT_IMAGE_BOUNDARY,
    EXTENDED_IMAGE_BOUNDARY
};

/*  Enumeration of medical image modalities.  Images having a particular
    modality may be treated in a specific manner when read or written.
    The modality is not use in class Image3D, although it is recorded there.

    A modality of CT implies that the intensity values fall within the
    calibrated range [-1024, 3071].  With PlanIm-format files, it is common to
    add 1024 to the intensities.  To handle this, modality SHIFTED_CT has been
    defined for the intensity range [0, 4095].  MRI images do not have any
    specific intensity range.  DQF intensities are not pixels, and are handled
    exclusively by class DQFImage.

    If this enumeration is changed, corresponding changes may be needed in the
    ImageIO files.  See file ImageStruct.h.
*/
enum modality_t { UNKNOWN_MODALITY = 0, CT, SHIFTED_CT, MRI, DQF };
const char * modalityString(modality_t m);  // Printable name of the modality

class Image3D;
typedef GreyValue (Image3D::* access_voxel_t)(int index) const;

class Image3D
{
public:

    Image3D();
    Image3D(int X, int Y, int Z);
    virtual ~Image3D();

#ifdef AE2_BUILD
    void setMasking(bool yesNo);
#endif

    // Function to specify the image's coordinate systems.  The first three
    // arguments are the inter-voxel spacings in world coordinates, any of
    // which may be negative.  Thus the world system may be left-handed.  If
    // no origin is supplied, the origin (center for the first voxel) will be
    // assumed to be at (0.0, 0.0, 0.0).  If an origin is provided, it should
    // be the world coordinates of the center of the first voxel.  This
    // function should only be used after the dimensions are set.
    void setSpacingAndOrigin(double xSpace, double ySpace, double zSpace,
        const Vector3D * origin = NULL);

    // Specify the actual intensity range of the image voxels.  This range
    // is set when an image is read, but is not used by this class.  It is
    // only here for use by the ImagePlanes class and for other purposes.
    // This range is not used when writing images to disk.
    void setRange(GreyValue min, GreyValue max) {
        minIntens = min;
        maxIntens = max;
        intensRange = max - min;
    }

    // Image modifiers:
    void extrudeAndPadInMinusZ(int slices, int padding);
    void extrudeAndPadInPlusZ(int slices, int padding);
    
    // Range of intensities of voxels stored in this object.  This is
    // provided mainly for the mappings to display intensities (below).
    // It is not used when writing images to disk.
    void range(GreyValue & min, GreyValue & max) const {
        min = minIntens;
        max = maxIntens;
    }

    // Function to specify how the image was mapped when input.
    void setIntensityMapping(double scale, int shift);
    void getIntensityMapping(double & scale, int & shift);

    // Function to replace specific image voxels.  Intensities outside
    // the prior intensity range should not be used.  However, if they
    // are, function adjustRange() may be used to recompute the stored
    // intensity range.  If voxels outside the original range are used,
    // the mappings to actual intensities below will give erroneous
    // values for those voxels.
    void setVoxel(int x, int y, int z, GreyValue val);

    void adjustRange();

    // Function to replace all image voxels.  Calling this function does
    // not change the world coordinate system, voxel access, or intensity
    // mappings.  It is normally necessary to call setSpacingAndOrigin()
    // after using this function.  Failing to do so may result in the
    // image being upside down in Y or Z.  The last argument should
    // normally be specified; if used, setRange() should be called.
    void setVoxels(GreyValue * voxels, int numX, int numY, int numZ,
        bool ignore = false);
    // Function to replace all image voxels without changing the image
    // dimensions or spacings.  The intensity range will be recomputed.
    bool replaceVoxels(GreyValue * voxels);

    // Function to convert a grey image into a binary image.  Voxels
    // having intensities equal to or above the threshold will be assigned
    // the supplied value.  All others will be set to 0.
    void makeBinary(GreyValue threshold, GreyValue value);

    int imageIndex(int x, int y, int z) const {
        return x + yIndexLUT[y] + zIndexLUT[z];
    }
    int imageIndex(int slices[3]) const {
        return slices[0] + yIndexLUT[slices[1]] + zIndexLUT[slices[2]];
    }
    // Set the image intensities to val.  If reset is true, the modality
    // will be set to unknown and the intensity mapping variables will
    // be reset.
    void clear(GreyValue val = 0, bool reset = true);

    bool clipToImage(int slices[3]) const;
    bool clipToImage(Vector3D & coord) const;
    bool imageCoordinatesInBounds(int x, int y, int z) const {
        return ! (x < 0 || y < 0 || z < 0 ||
            x >= xDim || y >= yDim || z >= zDim);
    }
    bool imageCoordinatesInBounds(int slices[3]) const {
        return ! (slices[0] < 0 || slices[1] < 0 || slices[2] < 0 ||
            slices[0] >= xDim || slices[1] >= yDim || slices[2] >= zDim);
    }

    // Actual intensities (for example, CT values) are those stored as shorts
    // in image files (.raw3).  Display voxel values are the internal voxel
    // values of this class stored as unsigned shorts (GreyValue) with a
    // maximum range of [0, 65535].  These are mapped through the lookup
    // table in class ImagePlanes before being displayed.  They are often
    // used directly for model fitting.  Relative voxel values are in the
    // range [0.0, 1.0] and are used for slice selection in the ImagePlanes
    // class.

    // The functions below may be used to convert intensities among these three
    // systems.  In conversion to actual voxel values, the mapping may not be
    // exact, since the voxels are rounded in the original input conversion
    // (see RAWImageFile::read()).  However, images with known modalities (e.g.
    // CT) are stored with predefined intensity ranges, often without remapping,
    // and the extreme values of the intensity sliders may fall outside the
    // actual range of image intensities.
    int mapRelativeToActual(double intens) const;
    double mapActualToRelative(int intens) const;
    GreyValue mapRelativeToDisplay(double intens) const;
    double mapDisplayToRelative(GreyValue intens) const;
    GreyValue mapActualToDisplay(int intens) const;
    int mapDisplayToActual(GreyValue intens) const;
    // Extra versions for higher precision
    double mapDisplayToRelative(double intens) const;
    double mapDisplayToActual(double intens) const;

    // Number of slices per Cartesian direction
    int getXDim() const { return xDim; }
    int getYDim() const { return yDim; }
    int getZDim() const { return zDim; }
    int getVoxelCount() const { return nVoxels; }

    // Signed (vector) physical (world) spacing between slices
    double getXSpacing() const { return spacing.getX(); }
    double getYSpacing() const { return spacing.getY(); }
    double getZSpacing() const { return spacing.getZ(); }
    Vector3D getSpacing() const { return spacing; }
    // Physical (world) spacing (distance) between slices
    double getAbsXSpacing() const { return xSpacing; }
    double getAbsYSpacing() const { return ySpacing; }
    double getAbsZSpacing() const { return zSpacing; }
    Vector3D getAbsSpacing() const {
        return Vector3D(xSpacing, ySpacing, zSpacing);
    }
    // Directions of the Cartesian axes
    Vector3D getSign() const { return sign; }
    void getSign(int * s) const;    // Array of size 3

    // Model-space (non-negative) distance between slices
    double getXModelSpacing() const;
    double getYModelSpacing() const;
    double getZModelSpacing() const;
    Vector3D getModelSpacing() const;

    // Physical (world) lengths of the axes
    double getXExtent() const { return xExtent; }
    double getYExtent() const { return yExtent; }
    double getZExtent() const { return zExtent; }
    double maxExtent() const;
    // Physical (world) axes (edges of the space) as vectors
    double getXExtentVec() const { return xExtent*sign.getX(); }
    double getYExtentVec() const { return yExtent*sign.getY(); }
    double getZExtentVec() const { return zExtent*sign.getZ(); }

    // Model space lengths (non-negative) of the axes
    double getXModelExtent() const { return xExtent*worldToModelScale; }
    double getYModelExtent() const { return yExtent*worldToModelScale; }
    double getZModelExtent() const { return zExtent*worldToModelScale; }

    // World coordinate position of the model space origin.
    // This origin is a half voxel before the origin in the image file.
    Vector3D getModelOrigin() const { return originPixelPos; }
    double getXModelOrigin() const { return originPixelPos.getX(); }
    double getYModelOrigin() const { return originPixelPos.getY(); }
    double getZModelOrigin() const { return originPixelPos.getZ(); }
    // The world coordinate position of the model-space bound, which
    // is the far edge of the last voxel.
    Vector3D getModelBound() const;

    // World coordinate origin of the image.  This origin is at the
    // center of the first voxel and equals the world-coordinate
    // origin recorded in the header of the image file.
    Vector3D getWorldOrigin() const;
    double getXWorldOrigin() const;
    double getYWorldOrigin() const;
    double getZWorldOrigin() const;
    // The world's bound corresponds to the model bound (above).

    // The image-space bound is the center of the last voxel.  These
    // functions return it in world coordinates.  This is a half
    // voxel less than the model-space bound (above).
    Vector3D getImageBound() const;
    double getXImageBound() const;
    double getYImageBound() const;
    double getZImageBound() const;

    // Axial coordinate conversion routines
    double imageXToModel(int xSlice) const;
    double imageYToModel(int ySlice) const;
    double imageZToModel(int zSlice) const;

    double imageXToWorld(int xSlice) const;
    double imageYToWorld(int ySlice) const;
    double imageZToWorld(int zSlice) const;

    // These three functions assume their arguments are >= 0.  Also, since
    // model space may exceed image space in up to two directions, these
    // may return invalid slice numbers.
    int modelXToImage(double xPos) const;
    int modelYToImage(double yPos) const;
    int modelZToImage(double zPos) const;

    double modelXToWorld(double xPos) const;
    double modelYToWorld(double yPos) const;
    double modelZToWorld(double zPos) const;

    // These three functions assume their arguments are >= the world origin
    int worldXToImage(double xPos) const;
    int worldYToImage(double yPos) const;
    int worldZToImage(double zPos) const;

    double worldXToModel(double xPos) const;
    double worldYToModel(double yPos) const;
    double worldZToModel(double zPos) const;

    // General coordinate conversion routines
    void imageToModelCoordinates(Vector3D & coord) const;
    void imageToModelCoordinates(const int slices[3], Vector3D & coord) const;

    void imageToWorldCoordinates(Vector3D & coord) const;
    void imageToWorldCoordinates(const int slices[3], Vector3D & coord) const;

    void modelToImageCoordinates(Vector3D & coord) const;
    void modelToImageCoordinates(const Vector3D & coord, int slices[3]) const;

    void modelToWorldCoordinates(Vector3D & coord) const;

    void worldToImageCoordinates(Vector3D & coord) const;
    void worldToImageCoordinates(const Vector3D & coord, int slices[3]) const; 

    void worldToModelCoordinates(Vector3D & coord) const;

    // The next function is deprecated - do not use it
    double getModelToWorldScale() const { return modelToWorldScale; }

    // Distance (cardinal direction) conversion routines.  These functions
    // convert a distance in the first space to a distance in the second
    // along the indicated cardinal direction.  Distances are considered
    // to be positive, scalar quantities, so a distance in one coordinate
    // system will alw1ays convert to a positive distance in another,
    // regardless of the handedness of the coordinate systems used.
    double imageXToModelDistance(double xLen) const;
    double imageYToModelDistance(double yLen) const;
    double imageZToModelDistance(double zLen) const;

    double imageXToWorldDistance(double xLen) const;
    double imageYToWorldDistance(double yLen) const;
    double imageZToWorldDistance(double zLen) const;

    // Caution: model space may exceed image space in up to two directions
    double modelXToImageDistance(double xLen) const;
    double modelYToImageDistance(double yLen) const;
    double modelZToImageDistance(double zLen) const;

    double worldXToImageDistance(double xLen) const;
    double worldYToImageDistance(double yLen) const;
    double worldZToImageDistance(double zLen) const;

    // For these the cardinal directions are all equivalent
    double modelToWorldDistance(double len) const;
    double worldToModelDistance(double len) const;
    // It is impossible to convert a voxel distance to world or
    // model coordinates, without knowing the angle at which the
    // distance was measured.

    // General distance (arbitrary direction) mappings.  These
    // functions convert the specified points between the coordinate
    // systems and then compute and return the distance between them.
    double imageToModelDistance(Vector3D & pt0, Vector3D & pt1) const;
    double imageToWorldDistance(Vector3D & pt0, Vector3D & pt1) const;
    double modelToWorldDistance(Vector3D & pt0, Vector3D & pt1) const;
    double modelToImageDistance(Vector3D & pt0, Vector3D & pt1) const;
    double worldToModelDistance(Vector3D & pt0, Vector3D & pt1) const;
    double worldToImageDistance(Vector3D & pt0, Vector3D & pt1) const;

    // Obtain voxels by slice index
    GreyValue getVoxelValue(int x, int y, int z);
    double getInterpolatedVoxelValue(double x, double y, double z);

    // Obtain voxels by array index.  No error checking is performed.
    GreyValue getVoxelValue(int index) const;   // Do not use for stacked images

    // Obtain intensity-windowed voxels by slice index
    GreyValue getWindowedVoxelValue(int x, int y, int z);
    double getWindowedInterpolatedVoxelValue(double x, double y, double z);

    // This function should not be used to bypass the voxel-access
    // functions above.  Directly indexing into the array may yield
    // incorrect intensity values, depending on the image file loaded.
    // Only use this function when all voxels of the image are treated
    // identically.
    GreyValue * getVoxels() const { return voxelArray; }

    // Functions used for interactive adjustment of image display levels,
    // as returned from the getWindowed*() functions above.
    void intensityWindow(double min, double max) {
        intensityWinMin = min;
        intensityWinMax = max;
        intensityWinRange = max - min;
    }
    double getMinIntensityWindow() const { return intensityWinMin; }
    double getMaxIntensityWindow() const { return intensityWinMax; }

    void print(const char * filename = NULL) const;
    void printWorld() const;

    // Medical image modality
    void setModality(modality_t modality) { /*mode = modality;*/ }
    modality_t modality() const { return mode; }
    const char * modalityString() const { return ::modalityString(mode); }

    // Stacked image information.
    void setIsImageStacked(bool yesno) { is_stacked = yesno; }
    bool getIsImageStacked() { return is_stacked; }
    void setStackedMask(GreyValue mask) { stackedMask = mask; }
    GreyValue getStackedMask() { return stackedMask; }

#ifdef BINARY
    // MultiObject: change image access values (isStacked, mask)
    // especially suited to bitStacked images -- for example, MOM
    // needs to know where a particular set of object(s) is/are,
    // but penetration needs to know where other objects are.
    // "Push" saves the specified isStacked and mask on the stack
    // (as in FILO queue) and sets the current values to the
    // passed values. "Pop" replaces the current values with the
    // top values from the stack (which are popped). These make it
    // easier to use values temporarily without having the
    // application record the current values for later
    // replacement.
    void pushImageIsStacked(bool isStacked = false, GreyValue mask = 0xFFFF);
    void popImageIsStacked();
#endif

protected:

    GreyValue * voxelArray;

    // Number of voxels in x, y, and z (image coordinate bounds)
    int xDim, yDim, zDim;
    int nVoxels;

    // Spacing between voxels in world coordinates.  Positive valued.
    double xSpacing, ySpacing, zSpacing;
    Vector3D spacing;  // Signed voxel spacings in world coordinates
    Vector3D sign;  // Sign (+/-1.0) to multiply by axial spacings

    // Extent of each full dimension (model-space 0.0 to 1.0) in world
    // coordinates (world coordinate bounds).  Positive valued.
    double xExtent, yExtent, zExtent;

    // Position of model-space origin in world coordinates
    Vector3D originPixelPos;    // Signed value

    double modelToWorldScale;   // Positive valued; sign must be considered in use
    double worldToModelScale;   // Positive valued; sign must be considered in use
    Vector3D worldToImageScale; // Signed value
    Vector3D modelToImageScale; // Positive valued
    Vector3D worldToImageOffset;    // Signed value

    // Look-up tables for y and z indexes
    int * yIndexLUT;
    int * zIndexLUT;

    double intensityWinMin, intensityWinMax, intensityWinRange;
    GreyValue minIntens, maxIntens, intensRange;    // Range of stored voxels

    double scale_factor;
    int intens_shift;

    ImageBoundaryType boundaryType;

    bool is_stacked;
    GreyValue stackedMask;

private:

    Image3D(Image3D & image);   // Not implemented
    Image3D operator=(Image3D & image); // Not implemented

    access_voxel_t getVoxel;
    GreyValue getRawVoxel(int index) const;
#ifdef AE2_BUILD
    GreyValue getMaskedVoxel(int index) const;
#endif

#ifdef BINARY
    // See pushImageIsStacked(), popImageIsStacked()
    int imageIsStackedCount;                    // # of items on stack
    GreyValue imageIsStackedMaskHistory[10];    // stack
    bool imageIsStackedIsStackedHistory[10];
#endif

    modality_t mode;
};


// Constant corresponding to the half-voxel shifts between origins
const double modelToImageOffset = -0.5;


// If the arguments of the next three functions are outside the world and
// sufficiently negative, then the returned values will be rounded incorrectly.
inline int Image3D::worldXToImage(double xPos) const
{
    return (int) (xPos*worldToImageScale.getX() + worldToImageOffset.getX() + 0.5);
}

inline int Image3D::worldYToImage(double yPos) const
{
    return (int) (yPos*worldToImageScale.getY() + worldToImageOffset.getY() + 0.5);
}

inline int Image3D::worldZToImage(double zPos) const
{
    return (int) (zPos*worldToImageScale.getZ() + worldToImageOffset.getZ() + 0.5);
}

// If the arguments of the next three functions are < 0, then
// the returned values may be improperly rounded and incorrect.
inline int Image3D::modelXToImage(double xPos) const
{
    return (int) (xPos*modelToImageScale.getX() + modelToImageOffset + 0.5);
}

inline int Image3D::modelYToImage(double yPos) const
{
    return (int) (yPos*modelToImageScale.getY() + modelToImageOffset + 0.5);
}

inline int Image3D::modelZToImage(double zPos) const
{
    return (int) (zPos*modelToImageScale.getZ() + modelToImageOffset + 0.5);
}

inline double Image3D::worldXToModel(double xPos) const
{
    return sign.getX()*(xPos - originPixelPos.getX())*worldToModelScale;
}

inline double Image3D::worldYToModel(double yPos) const
{
    return sign.getY()*(yPos - originPixelPos.getY())*worldToModelScale;
}

inline double Image3D::worldZToModel(double zPos) const
{
    return sign.getZ()*(zPos - originPixelPos.getZ())*worldToModelScale;
}

inline double Image3D::modelXToWorld(double xPos) const
{
    return sign.getX()*xPos*modelToWorldScale + originPixelPos.getX();
}

inline double Image3D::modelYToWorld(double yPos) const
{
    return sign.getY()*yPos*modelToWorldScale + originPixelPos.getY();
}

inline double Image3D::modelZToWorld(double zPos) const
{
    return sign.getZ()*zPos*modelToWorldScale + originPixelPos.getZ();
}

inline double Image3D::imageXToWorld(int xSlice) const
{
    return spacing.getX()*(xSlice - modelToImageOffset) + originPixelPos.getX();
}

inline double Image3D::imageYToWorld(int ySlice) const
{
    return spacing.getY()*(ySlice - modelToImageOffset) + originPixelPos.getY();
}

inline double Image3D::imageZToWorld(int zSlice) const
{
    return spacing.getZ()*(zSlice - modelToImageOffset) + originPixelPos.getZ();
}

inline double Image3D::imageXToModel(int xSlice) const
{
    return (xSlice - modelToImageOffset) / modelToImageScale.getX();
}

inline double Image3D::imageYToModel(int ySlice) const
{
    return (ySlice - modelToImageOffset) / modelToImageScale.getY();
}

inline double Image3D::imageZToModel(int zSlice) const
{
    return (zSlice - modelToImageOffset) / modelToImageScale.getZ();
}

inline double Image3D::worldXToImageDistance(double xLen) const
{
    return xLen*fabs(worldToImageScale.getX());
}

inline double Image3D::worldYToImageDistance(double yLen) const
{
    return yLen*fabs(worldToImageScale.getY());
}

inline double Image3D::worldZToImageDistance(double zLen) const
{
    return zLen*fabs(worldToImageScale.getZ());
}

inline double Image3D::modelXToImageDistance(double xLen) const
{
    return xLen*modelToImageScale.getX();
}

inline double Image3D::modelYToImageDistance(double yLen) const
{
    return yLen*modelToImageScale.getY();
}

inline double Image3D::modelZToImageDistance(double zLen) const
{
    return zLen*modelToImageScale.getZ();
}

inline double Image3D::imageXToModelDistance(double xLen) const
{
    return xLen/modelToImageScale.getX();
}

inline double Image3D::imageYToModelDistance(double yLen) const
{
    return yLen/modelToImageScale.getY();
}

inline double Image3D::imageZToModelDistance(double zLen) const
{
    return zLen/modelToImageScale.getZ();
}

inline double Image3D::imageXToWorldDistance(double xLen) const
{
    return xLen*fabs(spacing.getX());
}

inline double Image3D::imageYToWorldDistance(double yLen) const
{
    return yLen*fabs(spacing.getY());
}

inline double Image3D::imageZToWorldDistance(double zLen) const
{
    return zLen*fabs(spacing.getZ());
}

inline double Image3D::worldToModelDistance(double len) const
{
    return len*worldToModelScale;
}

inline double Image3D::modelToWorldDistance(double len) const
{
    return modelToWorldScale*len;
}

inline GreyValue Image3D::getVoxelValue(int index) const
{
    return (this->*getVoxel)(index);
}

inline void Image3D::setIntensityMapping(double scale, int shift)
{
    scale_factor = scale;
    intens_shift = shift;
}

inline void Image3D::getIntensityMapping(double & scale, int & shift)
{
    scale = scale_factor;
    shift = intens_shift;
}

/*  The rationale for the intensity mapping calculations is as follows.
    Images are expected to be read from files, and the input classes
    make the mapping found below in mapActualToDisplay().  Function
    mapDisplayToActual() is its inverse.  Next, the relative intensity
    range is defined to encompass the intensity range of the stored
    voxels, so that a relative intensity of 0.0 corresponds to display
    intensity minIntensity, and 1.0 corresponds to maxIntensity.  This
    determines the computations in functions mapRelativeToDisplay() and
    mapDisplayToRelative().  Finally, functions  mapRelativeToActual()
    and mapActualToRelative() are determined by composing the other
    functions.  Note that with some modalities (e.g. CT) the actual
    intensity range may fall within the relative range.  This will only
    have a minor effect on the operation of intensity sliders.

    Note: a couple of these functions are reimplemented in Image3D.cpp.
    It they are changed, changes also may be needed there.
*/

inline int Image3D::mapRelativeToActual(double intens) const
{
    double i = (intens*intensRange + minIntens)/scale_factor + intens_shift;
    return (int) ((i < 0.0) ? i - 0.5 : i + 0.5);
}

inline double Image3D::mapActualToRelative(int intens) const
{
    return ((intens - intens_shift)*scale_factor - minIntens)/intensRange;
}

inline GreyValue Image3D::mapRelativeToDisplay(double intens) const
{
    return minIntens + (GreyValue) (intens*intensRange + 0.5);
}

inline double Image3D::mapDisplayToRelative(GreyValue intens) const
{
    return ((double) (intens - minIntens))/intensRange;
}

inline double Image3D::mapDisplayToRelative(double intens) const
{
    return (intens - minIntens)/intensRange;
}

inline GreyValue Image3D::mapActualToDisplay(int intens) const
{
    return (GreyValue) (0.5 + (intens - intens_shift)*scale_factor);
}

inline int Image3D::mapDisplayToActual(GreyValue intens) const
{
    double i = intens/scale_factor + intens_shift;
    return (int) ((i < 0.0) ? i - 0.5 : i + 0.5);
}

inline double Image3D::mapDisplayToActual(double intens) const
{
    return intens/scale_factor + intens_shift;
}


#endif

