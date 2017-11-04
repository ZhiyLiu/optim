#include "Trackball.h"

//#define DEBUG	1

Trackball::Trackball(double width, double height) : rotQuat()
{
	control = NULL;
	trackMode = NONE;

	dTrans.set(0.0, 0.0, 0.0);
	currScale = 1.0f;
	dScale = 0.0f;

	viewWidth = width;
	viewHeight = height;

	canRotate = true;
	rotationType = ROTATE_3D;
}

void Trackball::setMode(TrackMode mode)
{
#ifdef DEBUG
	cout << "Trackball mode set to " << mode << endl;
#endif
	trackMode = mode;
}

// Reset the view transformation matrix
void Trackball::reset(float * transMat)
{
	int i, j;
	for(i = 0; i < 4; i++)
	{
		for(j = 0; j < 4; j++)
		{
			transformMatrix[i][j] = *transMat;
			currMatrix[i][j] = *transMat;
			transMat++;
		}
	}

	rotQuat = Quat();

    // Find out the scale of this matrix (NOTE: Assuming a similarity transform)
    Vector3D test(transformMatrix[0][0], transformMatrix[1][0], transformMatrix[2][0]);
    currScale = (float) test.norm();

	dScale = 0.0f;
	dTrans.set(0.0, 0.0, 0.0);
	trackMode = NONE;
}

void Trackball::mouseEvent(const MouseEventStruct & mouse)
{
	Vector3D objTrans;
	double objScale;
	Quat objQuat;
	Vector3D objCOG;
	M3DObject * object;
	TileSet *tileSet;

	// Set what mode we are in
	if (mouse.type == MOUSE_PUSH)
	{
		determineMode(mouse);

		if (control != NULL)
		{
			if (trackMode == OBJECT_ROTATE || trackMode == OBJECT_TRANSLATE ||
			   trackMode == OBJECT_SCALE || trackMode == OBJECT_PIROUETTE)
				control->startOperation();

		}

		currMouseX = lastMouseX = mouse.x;
		currMouseY = lastMouseY = mouse.y;

		rotQuat = Quat();
		dTrans.set(0.0, 0.0, 0.0);
		dScale = 0.0f;
#ifdef DEBUG
		cout << "mouseEvent: push - mode = " << trackMode << endl;
#endif
	}

	else if (mouse.type == MOUSE_DRAG)
	{
		currMouseX = mouse.x;
		currMouseY = mouse.y;

		switch (trackMode)
		{
			case VIEW_ROTATE:
				rotate3D(rotQuat, 
					(2.0 * lastMouseX / mouse.width) - 1.0,
					1.0 - (2.0 * lastMouseY / mouse.height),
					(2.0 * currMouseX / mouse.width) - 1.0,
					1.0 - (2.0 * currMouseY / mouse.height));

				lastMouseX = currMouseX;
				lastMouseY = currMouseY;

				break;

			case VIEW_TRANSLATE:
				dTrans.set(viewWidth * (currMouseX - lastMouseX) / mouse.width,
					-viewHeight * (currMouseY - lastMouseY) / mouse.height,
					0.0);

				lastMouseX = currMouseX;
				lastMouseY = currMouseY;

				break;

			case VIEW_SCALE:
				dScale = SCALE_FACTOR * (currMouseY - lastMouseY) / mouse.height;

				lastMouseY = currMouseY;

				break;

			case OBJECT_ROTATE:
				if (control == NULL)
					break;

				object = control->getObjectPtr();
				tileSet = control->getTileSetPtr();
				if (object || tileSet) {
					if (rotationType == ROTATE_3D) {
						objCOG = object ? object->getCOG() : tileSet->getCOG();
						objCOG.multByMatrix4(&(transformMatrix[0][0]));
						
						rotate3D(objQuat,
							(2.0 * lastMouseX / mouse.width) - 1.0 - objCOG.getX(),
							1.0 - (2.0 * lastMouseY / mouse.height) - objCOG.getY(),
							(2.0 * currMouseX / mouse.width) - 1.0 - objCOG.getX(),
							1.0 - (2.0 * currMouseY / mouse.height) - objCOG.getY());
					}
					else {
						rotate2D(objQuat,
							(double)(currMouseX - lastMouseX) / (double) mouse.width);
					}
					
					rotateObject(objQuat.conj());
					rotateTileSet(objQuat.conj());
					
					lastMouseX = currMouseX;
					lastMouseY = currMouseY;
				}

				break;

			case OBJECT_PIROUETTE:	// only tileSet's can pirouette
				if (control == NULL)
					break;

				tileSet = control->getTileSetPtr();
				if (tileSet) {
					if(rotationType == ROTATE_3D) {
						objCOG = tileSet->getCOG();	// all tile-figures share a rotation
						objCOG.multByMatrix4(&(transformMatrix[0][0]));
						
						rotate3D(objQuat,
							(2.0 * lastMouseX / mouse.width) - 1.0 - objCOG.getX(),
							1.0 - (2.0 * lastMouseY / mouse.height) - objCOG.getY(),
							(2.0 * currMouseX / mouse.width) - 1.0 - objCOG.getX(),
							1.0 - (2.0 * currMouseY / mouse.height) - objCOG.getY());
					}
					else {
						rotate2D(objQuat,
							(double)(currMouseX - lastMouseX) / (double) mouse.width);
					}
					
					pirouetteTileSet(objQuat.conj());
					
					lastMouseX = currMouseX;
					lastMouseY = currMouseY;
				}

				break;


			case OBJECT_TRANSLATE:
				if (control == NULL)
					break;

				objTrans.set(viewWidth * (currMouseX - lastMouseX) / mouse.width,
					-viewHeight * (currMouseY - lastMouseY) / mouse.height,
					0.0);

				translateObject(objTrans);

				lastMouseX = currMouseX;
				lastMouseY = currMouseY;

				break;

			case OBJECT_SCALE:
				if (control == NULL)
					break;

				objScale = SCALE_FACTOR * (currMouseY - lastMouseY) / mouse.height;
				control->scale(1.0 + objScale);

				lastMouseY = currMouseY;

				break;

			case OBJECT_SCALE_WIDTH:
				if (control == NULL)
					break;

				objScale = SCALE_FACTOR * (currMouseY - lastMouseY) / mouse.height;
				control->scaleWidth(1.0 + objScale);

				lastMouseY = currMouseY;

				break;

			default:
				break;
		};
#ifdef BINARY
#ifdef DEBUG
		object = control->getObjectPtr();
		// Calc & print model's regularity for each mouse event
		if (object)
			cout << "REGU="
			 << object->dist2FromAveOfNeighbors(
				  (enum PrimNeighborhoodDefn) (int) tuningWt(BpAtomNeighborhood),
				  0,
				  (enum DistanceType) (int) tuningWt(BpAtomDistanceType))
			 <<  endl;
#endif
#endif

#ifdef DEBUG
		cout << "mouseEvent: drag - mode = " << trackMode << endl;
#endif
	}

	else if (mouse.type == MOUSE_RELEASE)
	{
#ifdef DEBUG
		cout << "mouseEvent: release - mode = " << trackMode << endl;
#endif
		if (control != NULL) {
			if (trackMode == OBJECT_ROTATE || trackMode == OBJECT_TRANSLATE
				|| trackMode == OBJECT_SCALE)
			{
				control->endOperation();
			}
		}
		trackMode = NONE;
	}
}

#ifdef AE2_BUILD

void Trackball::kbEvent(const KBTransformationStruct & kb)
{
	Vector3D objTrans;
	double objScale;
	Quat objQuat;
	Vector3D objCOG;
	M3DObject * object;
	TileSet *tileSet;

	if (control == NULL)
		return;

	// kb.x/kb.y are supposed to store deltaX and deltaY
	currMouseX = lastMouseX + kb.x;
	currMouseY = lastMouseY + kb.y;
	
	switch (kb.type)
	{
	    case KB_ROTATE:
		    object = control->getObjectPtr();
		    tileSet = control->getTileSetPtr();
		    if (object || tileSet) {
			    if (rotationType == ROTATE_3D) {
				    objCOG = object ? object->getCOG() : tileSet->getCOG();
				    objCOG.multByMatrix4(&(transformMatrix[0][0]));
				    
				    rotate3D(objQuat,
					    (2.0 * lastMouseX / kb.width) - 1.0 - objCOG.getX(),
					    1.0 - (2.0 * lastMouseY / kb.height) - objCOG.getY(),
					    (2.0 * currMouseX / kb.width) - 1.0 - objCOG.getX(),
					    1.0 - (2.0 * currMouseY / kb.height) - objCOG.getY());
			    }
			    else {
				    rotate2D(objQuat,
					    (double)(currMouseX - lastMouseX) / (double) kb.width);
			    }

			    rotateObject(objQuat.conj());
			    rotateTileSet(objQuat.conj());

			    lastMouseX = currMouseX;
			    lastMouseY = currMouseY;
		    }

		    break;

	    case KB_TRANSLATE:
		    objTrans.set(viewWidth * (currMouseX - lastMouseX) / kb.width,
			    -viewHeight * (currMouseY - lastMouseY) / kb.height,
			    0.0);

		    translateObject(objTrans);

		    lastMouseX = currMouseX;
		    lastMouseY = currMouseY;

		    break;
		    
	    case KB_SCALE:

		    objScale = SCALE_FACTOR * (currMouseY - lastMouseY) / kb.height;
		    control->scale(1.0 + objScale);

		    lastMouseY = currMouseY;

		    break;
	    default:
		    break;
	};
#ifdef DEBUG
	cout << "kbEvent: kb event type = " << kb.type << endl;
#endif
}

#endif  /* AE2_BUILD */

void Trackball::doViewTransform()
{
	double rotMatrix[4][4];

	if(trackMode == VIEW_SCALE)
	{
		float lastScale = currScale;
		currScale *= (1.0f + dScale);

		// Make sure we don't go out of range
		if(currScale > MAX_SCALE)
		{
			// Amount of change in scale needed to get to MAX_SCALE
			dScale = (MAX_SCALE / lastScale) - 1.0f;
			currScale = MAX_SCALE;
		}

		else if(currScale < MIN_SCALE)
		{
			// Amount of change in scale needed to get to MIN_SCALE
			dScale = (MIN_SCALE / lastScale) - 1.0f;
			currScale = MIN_SCALE;
		}

		glScalef(1.0f + dScale, 1.0f + dScale, 1.0f + dScale);
	}

	else if(trackMode == VIEW_TRANSLATE)
		glTranslated(dTrans.getX(), dTrans.getY(), dTrans.getZ());

	else if(trackMode == VIEW_ROTATE)
	{
		rotQuat.buildRotMatrix(&(rotMatrix[0][0])); 
		glMultMatrixd(&(rotMatrix[0][0]));
	}

	glMultMatrixf(&(transformMatrix[0][0]));

	if(trackMode != NONE)
		glGetFloatv(GL_MODELVIEW_MATRIX, &(transformMatrix[0][0]));

	glGetFloatv(GL_MODELVIEW_MATRIX, &(currMatrix[0][0]));
}

// Open GL translation function						AGG: discard?
void Trackball::translate(double x, double y, double z)
{
	glTranslated(x, y, z);
	glGetFloatv(GL_MODELVIEW_MATRIX, &(currMatrix[0][0]));
}

// Open GL rotation function						AGG: discard?
void Trackball::rotate(double angle, double x, double y, double z)
{
	glRotated(angle, x, y, z);
	glGetFloatv(GL_MODELVIEW_MATRIX, &(currMatrix[0][0]));
}

// Open GL matrix composition function						AGG: discard?
void Trackball::multMatrix(double * mat)
{
	glMultMatrixd(mat);
	glGetFloatv(GL_MODELVIEW_MATRIX, &(currMatrix[0][0]));
}

// Translation of the viewport by a vector specified in object space coordinates
void Trackball::translateView(double x, double y, double z)
{
	glTranslated(x, y, z);
	glGetFloatv(GL_MODELVIEW_MATRIX, &(transformMatrix[0][0]));
}

// Rotation by angle of the viewport about a vector specified in object space coordinates
void Trackball::rotateView(double angle, double x, double y, double z)
{
	glRotated(angle, x, y, z);
	glGetFloatv(GL_MODELVIEW_MATRIX, &(transformMatrix[0][0]));
}

void Trackball::determineMode(const MouseEventStruct & mouse)
{
#ifdef DEBUG
	cout << "Mouse modifiers: " << mouse.modifiers << " = "
		<< ((mouse.modifiers & ALT_KEY) ? "ALT " : "")
		<< ((mouse.modifiers & META_KEY) ? "META " : "")
		<< ((mouse.modifiers & CTRL_KEY) ? "CTRL " : "")
		<< ((mouse.modifiers & SHIFT_KEY) ? "SHIFT " : "") << endl;
#endif

	// Left mouse w/ no keys pressed is a view rotation
	if((mouse.button == 1) && (mouse.modifiers == NO_KEY) && canRotate)
		trackMode = VIEW_ROTATE;

	// Middle mouse w/ no keys pressed is a view translate
	else if((mouse.button == 2) && (mouse.modifiers == NO_KEY))
		trackMode = VIEW_TRANSLATE;

	// Left mouse w/ shift key simulates middle mouse (view translate)
	else if((mouse.button == 1) && (mouse.modifiers == SHIFT_KEY))
		trackMode = VIEW_TRANSLATE;

	// Right mouse w/ no keys pressed is a view scale
	else if((mouse.button == 3) && (mouse.modifiers == NO_KEY))
		trackMode = VIEW_SCALE;

	// Left mouse w/ Alt key is an object rotation
	else if((mouse.button == 1) && (mouse.modifiers == ALT_KEY || mouse.modifiers == META_KEY
		|| mouse.modifiers == (SHIFT_KEY | CTRL_KEY)))
			trackMode = OBJECT_ROTATE;

    // Left mouse w/ Alt + shift key is an object pirouette
	else if((mouse.button == 1) && (mouse.modifiers == (SHIFT_KEY | ALT_KEY)
		|| mouse.modifiers == (SHIFT_KEY | META_KEY)))
            trackMode = OBJECT_PIROUETTE;

	// Middle mouse w/ Alt key is an object translation
	else if((mouse.button == 2) && (mouse.modifiers == ALT_KEY || mouse.modifiers == META_KEY
		|| mouse.modifiers == (SHIFT_KEY | CTRL_KEY)))
			trackMode = OBJECT_TRANSLATE;

	// Left mouse w/ shift key simulates middle mouse (+ ALT gives object translate)
	else if((mouse.button == 1) && (mouse.modifiers == (SHIFT_KEY | ALT_KEY)
		|| mouse.modifiers == (SHIFT_KEY | ALT_KEY)))
			trackMode = OBJECT_TRANSLATE;

	// Right mouse w/ Alt key is an object scale
	else if((mouse.button == 3) && (mouse.modifiers == ALT_KEY || mouse.modifiers == META_KEY
		|| mouse.modifiers == (SHIFT_KEY | CTRL_KEY)))
			trackMode = OBJECT_SCALE;

    // Right mouse w/ Alt + shift key is an object scale width
	else if((mouse.button == 3) && (mouse.modifiers == (SHIFT_KEY | ALT_KEY)
		|| mouse.modifiers == (SHIFT_KEY | META_KEY)))
            trackMode = OBJECT_SCALE_WIDTH;

	else
		trackMode = NONE;
}

// Project an x,y pair onto a sphere of radius r OR a hyperbolic sheet
// if we are away from the center of the sphere.
static double projectToSphere(double r, double x, double y)
{
	double distSquared, z;

	distSquared = (x * x) + (y * y);

	// inside a sphere
	if(distSquared < ((r * r) / 2.0))
		z = sqrt((r * r) - distSquared);

	// On a hyperbola
    else
        z = (r * r) / (2.0 * sqrt(distSquared));

	return z;
}

// Sets the Quat to a rotation projected to a sphere
void Trackball::rotate3D(Quat &q, double p1x, double p1y, double p2x, double p2y)
{
	// Zero rotation
	if((p1x == p2x) && (p1y == p2y))
	{
		q = Quat();
		return;
	}

	Vector3D axis;  // Axis of rotation
	double   theta; // how much to rotate

	// Here are our two points projected to the sphere
	Vector3D p1 (p1x, p1y, projectToSphere(TRACKBALLSIZE, p1x, p1y)),
		     p2 (p2x, p2y, projectToSphere(TRACKBALLSIZE, p2x, p2y));

	axis = p2.cross(p1);

	// Figure out how much to rotate around that axis.
	double t = (p1 - p2).norm() / (2.0 * TRACKBALLSIZE);

	// Avoid problems with out-of-control values...
	if(t > 1.0)
		t = 1.0;
	else if(t < -1.0)
		t = -1.0;

	theta = 2.0 * asin(t);

	// Create the rotation quaternion
	q.setAxisAngle(axis, theta);
}

// Sets the Quat to a rotation in the plane of the current view
void Trackball::rotate2D(Quat &q, double deltax)
{
	if (deltax == 0.0)
		return;

	// Inverse of the current view scale
	double invScale = 1.0 / currScale;
	double scaleFact = invScale * invScale;

	// Find inverse of view rotation matrix (upper-left 3x3 of view matrix)
	// Since a rotation matrix is orthogonal, we just divide out
	// by the scale and then take the transpose
	// The top two rows of this matrix are our new X and Y axis
	double invRotMat[2][3];

	invRotMat[0][0] = currMatrix[0][0] * scaleFact;
	invRotMat[0][1] = currMatrix[1][0] * scaleFact;
	invRotMat[0][2] = currMatrix[2][0] * scaleFact;

	invRotMat[1][0] = currMatrix[0][1] * scaleFact;
	invRotMat[1][1] = currMatrix[1][1] * scaleFact;
	invRotMat[1][2] = currMatrix[2][1] * scaleFact;

	Vector3D newX (&(invRotMat[0][0])),
		     newY (&(invRotMat[1][0]));

	// Axis of rotation is the normal to the viewing plane
	Vector3D normal = newX.cross(newY);

	if (deltax > 1.0)
		deltax = 1.0;
	else if (deltax < -1.0)
		deltax = -1.0;

	q.setAxisAngle(normal, deltax * R_TWO_PI);
}

// Translates the object with a translation vector
// The vector is given relative to the view coordinates (just X and Y matter)
// and must be transformed into world coordinates by the inverse view matrix
void Trackball::translateObject(const Vector3D &trans)
{
	if(control == NULL)
		return;

	// Inverse of the current view scale
	double invScale = 1.0 / currScale;
	double scaleFact = invScale * invScale;

	// Find inverse of view rotation matrix (upper-left 3x3 of view matrix)
	// Since a rotation matrix is orthogonal, we just divide out
	// by the scale and then take the transpose
	// The top two rows of this matrix are our new X and Y axis
	double invRotMat[2][3];

	invRotMat[0][0] = currMatrix[0][0] * scaleFact;
	invRotMat[0][1] = currMatrix[1][0] * scaleFact;
	invRotMat[0][2] = currMatrix[2][0] * scaleFact;

	invRotMat[1][0] = currMatrix[0][1] * scaleFact;
	invRotMat[1][1] = currMatrix[1][1] * scaleFact;
	invRotMat[1][2] = currMatrix[2][1] * scaleFact;

	Vector3D newX (&(invRotMat[0][0])),
		     newY (&(invRotMat[1][0]));

	Vector3D finalTrans = newX * trans.getX() + newY * trans.getY();

	// Translate the object with the transformed translation vector
	control->translate(finalTrans);
}

// Rotates the object by the given Quat, which is transformed by the view matrix
void Trackball::rotateObject(Quat objQuat)
{
	if(control == NULL)
		return;

	double n = objQuat.norm();

	double invScale = 1.0 / currScale;
	Vector3D axis = objQuat.getVector();

	// Break out if we have a zero vector
	if(axis * axis == 0)
		return;

	axis.normalize();

	Vector3D oldAxis = axis;

	// Transform the quaternion's axis by the view matrix
	axis.setX(transformMatrix[0][0] * oldAxis.getX() +
			  transformMatrix[0][1] * oldAxis.getY() +
			  transformMatrix[0][2] * oldAxis.getZ());

	axis.setY(transformMatrix[1][0] * oldAxis.getX() +
			  transformMatrix[1][1] * oldAxis.getY() +
			  transformMatrix[1][2] * oldAxis.getZ());

	axis.setZ(transformMatrix[2][0] * oldAxis.getX() +
			  transformMatrix[2][1] * oldAxis.getY() +
			  transformMatrix[2][2] * oldAxis.getZ());

	// Re-normalize the quaternion's length
	axis *= invScale;

	double halfAngle = acos(objQuat.getW());
	double sinHalfAngle = sin(halfAngle);

	objQuat.setX(axis.getX() * sinHalfAngle);
	objQuat.setY(axis.getY() * sinHalfAngle);
	objQuat.setZ(axis.getZ() * sinHalfAngle);

	// Perform the rotation on the object with the transformed quat
	control->rotate(objQuat);
}



// Rotates the object by the given Quat, which is transformed by the view matrix
void Trackball::rotateTileSet(Quat objQuat)
{
	if(control == NULL)
		return;

	double n = objQuat.norm();

	double invScale = 1.0 / currScale;
	Vector3D axis = objQuat.getVector();

	// Break out if we have a zero vector
	if(axis * axis == 0)
		return;

	axis.normalize();

	Vector3D oldAxis = axis;

	// Transform the quaternion's axis by the view matrix
	axis.setX(transformMatrix[0][0] * oldAxis.getX() +
			  transformMatrix[0][1] * oldAxis.getY() +
			  transformMatrix[0][2] * oldAxis.getZ());

	axis.setY(transformMatrix[1][0] * oldAxis.getX() +
			  transformMatrix[1][1] * oldAxis.getY() +
			  transformMatrix[1][2] * oldAxis.getZ());

	axis.setZ(transformMatrix[2][0] * oldAxis.getX() +
			  transformMatrix[2][1] * oldAxis.getY() +
			  transformMatrix[2][2] * oldAxis.getZ());

	// Re-normalize the quaternion's length
	axis *= invScale;

	double halfAngle = acos(objQuat.getW());
	double sinHalfAngle = sin(halfAngle);

	objQuat.setX(axis.getX() * sinHalfAngle);
	objQuat.setY(axis.getY() * sinHalfAngle);
	objQuat.setZ(axis.getZ() * sinHalfAngle);

	// Perform the rotation on the object with the transformed quat
	control->rotateTileSet(objQuat);
}



// Rotates the object by the given Quat, which is transformed by the view matrix
void Trackball::pirouetteTileSet(Quat objQuat)
{
	if(control == NULL)
		return;

	double n = objQuat.norm();

	double invScale = 1.0 / currScale;
	Vector3D axis = objQuat.getVector();

	// Break out if we have a zero vector
	if(axis * axis == 0)
		return;

	axis.normalize();

	Vector3D oldAxis = axis;

	// Transform the quaternion's axis by the view matrix
	axis.setX(transformMatrix[0][0] * oldAxis.getX() +
			  transformMatrix[0][1] * oldAxis.getY() +
			  transformMatrix[0][2] * oldAxis.getZ());

	axis.setY(transformMatrix[1][0] * oldAxis.getX() +
			  transformMatrix[1][1] * oldAxis.getY() +
			  transformMatrix[1][2] * oldAxis.getZ());

	axis.setZ(transformMatrix[2][0] * oldAxis.getX() +
			  transformMatrix[2][1] * oldAxis.getY() +
			  transformMatrix[2][2] * oldAxis.getZ());

	// Re-normalize the quaternion's length
	axis *= invScale;

	double halfAngle = acos(objQuat.getW());
	double sinHalfAngle = sin(halfAngle);

	objQuat.setX(axis.getX() * sinHalfAngle);
	objQuat.setY(axis.getY() * sinHalfAngle);
	objQuat.setZ(axis.getZ() * sinHalfAngle);

	// Perform the rotation on the object with the transformed quat
	control->pirouetteTileSet(objQuat);
}

void Trackball::printMatrix() const
{
    int i, j;
    for(i = 0; i < 4; i++) {
        for(j = 0; j < 3; j++)
            cout << transformMatrix[j][i] << ", ";
        cout << transformMatrix[3][i] << '\n';
    }
	cout << flush;
}

static const char *names[] = {
	"NONE", "VIEW_ROTATE", "VIEW_TRANSLATE", "VIEW_SCALE",
    "OBJECT_ROTATE", "OBJECT_TRANSLATE", "OBJECT_SCALE",
	"OBJECT_SCALE_WIDTH", "TILE_PIROUETTE"
};

const char * Trackball::getMouseMode() const
{
	return names[(int) trackMode];
}


