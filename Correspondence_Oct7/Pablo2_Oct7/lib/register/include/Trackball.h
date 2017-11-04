#ifndef TRACKBALL_H
#define TRACKBALL_H

#include <math.h>

#include "Quat.h"
#include "M3DObject.h"
#include "P3DControl.h"

#define TRACKBALLSIZE 0.8

// Bits to represent which keys are pressed
#define NO_KEY    0
#define CTRL_KEY  1
#define SHIFT_KEY 2
#define ALT_KEY   4
#define META_KEY  8

// Maximum amount we are allowed to scale with the trackball
#define MAX_SCALE    100.0f
#define MIN_SCALE    0.001f
#define SCALE_FACTOR 5.0f

// Type of tracking (transform view or object)
enum TrackMode {
    NONE, 
    VIEW_ROTATE, VIEW_TRANSLATE, VIEW_SCALE,
    OBJECT_ROTATE, OBJECT_PIROUETTE, OBJECT_TRANSLATE, OBJECT_SCALE, OBJECT_SCALE_WIDTH
};

// Type of rotation (3D or 2D)
enum RotationType {
    ROTATE_3D, ROTATE_2D
};

// Type of mouse event
enum MouseEventType {
    MOUSE_PUSH, MOUSE_DRAG, MOUSE_RELEASE
};

// Which orthographic view we are using
enum ViewType {
//    XY_ORTHO,
//    YZ_ORTHO,
//    XZ_ORTHO,
    PERSPECTIVE
};

// Structure to represent a mouse event
struct MouseEventStruct
{
    MouseEventType type;            // Type of mouse event
    int            button;          // Which button was pressed
    int            x, y;            // Coordinates in window
    int            width, height;   // Dimensions of window
    int            modifiers;       // What keys are pressed (shift, ctrl, ...)
    ViewType       viewType;
};


#ifdef AE2_BUILD

// Type of kb transformation option
enum KBTransformationType
{
    KB_NONE, KB_ROTATE, KB_TRANSLATE, KB_SCALE
};

// Structure to represent a keyboard event
struct KBTransformationStruct
{
    KBTransformationType type;      // Type of keyboard event
    int            x, y;            // Coordinates in window
    int            width, height;   // Dimensions of window
    ViewType       viewType;
};

#endif  /* AE2_BUILD */


class Trackball
{
public:
    Trackball(double width = 0.0, double height = 0.0);
    ~Trackball() {}

    void reset(float * transMat);

    void setControl(P3DControl * cntrl) { control = cntrl; }

    void setViewDims(double width, double height)
    {
        viewWidth = width;
        viewHeight = height;
    }

    void setMode(TrackMode mode);
    TrackMode getMode() { return trackMode; }
    bool getModeClass() { return trackMode >= OBJECT_ROTATE; }

    float getCurrScale() { return currScale; }
    void getCurrTrans(Vector3D &trans)
    {
        trans.setX(transformMatrix[3][0]);
        trans.setY(transformMatrix[3][1]);
        trans.setZ(transformMatrix[3][2]);
    }
	float* getCurrMatrix() {return &currMatrix[0][0];}
	float* getMatrix() {return &transformMatrix[0][0];}

    void mouseEvent(const MouseEventStruct & mouse);
#ifdef AE2_BUILD
	// to simulate mouse based model transformations
	// via kb
	void kbEvent(const KBTransformationStruct & kb);
#endif  /* AE2_BUILD */

    void doViewTransform();

    // Functions to manually change transform matrix
    void translate(double x, double y, double z);
    void rotate(double angle, double x, double y, double z);
    void multMatrix(double * mat);

    // Translation of the viewport by a vector specified in object space coordinates
    void translateView(double x, double y, double z);

    // Rotation by angle of the viewport about a vector specified in object space coordinates
    void rotateView(double angle, double x, double y, double z);

    // Don't allow the trackball to rotate (view rotation)
    void disableRotation() { canRotate = false; }

    void setRotationType(RotationType t) { rotationType = t; }

    // Returns true if the view trackball is currently in use
	// to change the view.
    bool isTracking()
    {
        if((trackMode & 0x0F) == VIEW_ROTATE || (trackMode & 0x0F) == VIEW_TRANSLATE ||
            (trackMode & 0x0F) == VIEW_SCALE)
        {
            return true;
        }

        return false;
    }

    // Sets the current matrix to our base matrix
    void saveMatrix()
    {
        int i, j;
        for(i = 0; i < 4; i++)
            for(j = 0; j < 4; j++)
                transformMatrix[i][j] = currMatrix[i][j];
    }

    void printMatrix() const;			// Debugging function
    const char * getMouseMode() const;	// Debugging function

private:
    void determineMode(const MouseEventStruct & mouse);

    void rotate3D(Quat &q, double p1x, double p1y, double p2x, double p2y);
    void rotate2D(Quat &q, double deltax);

    void translateObject(const Vector3D &trans);
    void rotateObject(Quat objQuat);

    void translateTileSet(const Vector3D &trans);
    void rotateTileSet(Quat tileQuat);
    void pirouetteTileSet(Quat tileQuat);

    TrackMode trackMode;

    RotationType rotationType;

    P3DControl * control;

    Quat rotQuat;

    float currScale;
    float dScale;

    Vector3D dTrans;

    // The column-major, view transformation matrix
    float transformMatrix[4][4];

    // The column-major, current transformation matrix
	// May include other transformations that are not done by the trackball
    float currMatrix[4][4];

    // Position of mouse in current and previous operations
    int currMouseX, currMouseY,
        lastMouseX, lastMouseY;

    // Dimensions of view in world coordinates
    double viewWidth,
           viewHeight;

    // Let us disable rotations (mainly for orthogonal views)
    bool canRotate;
};

#endif

