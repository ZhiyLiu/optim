#include <math.h>
#include <sys/stat.h>
#include <stdio.h>
#include <assert.h>

#include "M3DPGA.h"


const int NUM_PGA_PARAMS = 9;

#ifndef S_ISREG
#define S_ISREG(m)        (((S_IFREG & m) != 0) ? true : false)
#endif

M3DPGA::M3DPGA(M3DObject * object) : pgaData()
{
    if(object != NULL)
        initializeMean(object);
}

void M3DPGA::initializeMean(M3DObject * object)
{
    if(pgaData.mean != NULL)
    {
        delete pgaData.mean;
        pgaData.mean = NULL;
    }

    if(object != NULL)
        pgaData.mean = object->assign();
}

#ifndef BINARY
bool M3DPGA::readPGAFile(const char * pgaFilename)
{
    struct stat buf;
    FILE * fp;
    int ret;
    int i;


    if(pgaFilename == NULL)
        return false;

    if(stat(pgaFilename, &buf) != 0 || ! S_ISREG(buf.st_mode))
        return false;

    fp = fopen(pgaFilename, "r");
    ret = fscanf(fp, "%d\n%d\n", &(pgaData.numVectors), &(pgaData.vectorSize));
    if(ret != 2)
        return false;

    if(pgaData.v != NULL)
        delete pgaData.v;

    pgaData.v = new double[pgaData.vectorSize * pgaData.numVectors];
    for(i = 0; i < pgaData.vectorSize * pgaData.numVectors; i++)
    {
        ret = fscanf(fp, "%lf\n", &(pgaData.v[i]));
        if(ret != 1)
        {
            delete pgaData.v;
            pgaData.v = NULL;
            return false;
        }
    }
    fclose(fp);

    return true;
}
#endif

bool M3DPGA::pgaDeform(M3DObject * object, std::vector<double> & vals)
{
    int primitiveCount;
    int i, j;

    M3DPrimitive * prim,
                 * meanPrim;
    Quat q, q0, q1;
    Vector3D axis;
    double theta;
    Vector3D pos, n0, n1;
    int numParams;


    if(object == NULL || pgaData.mean == NULL || pgaData.v == NULL)
        return false;

    primitiveCount = object->getPrimitiveCount();

    if(NUM_PGA_PARAMS*primitiveCount != pgaData.vectorSize ||
		primitiveCount != pgaData.mean->getPrimitiveCount()) {
		printf("[FATAL] Incorrect number of PGA parameters\n");
        return false;
	}

    if(pgaData.numVectors > vals.size())
        numParams = vals.size();
    else
        numParams = pgaData.numVectors;

    double * v1 = new double[pgaData.vectorSize];
    for(i = 0; i < pgaData.vectorSize; i++)
        v1[i] = 0.0;

    for(i = 0; i < numParams; i++)
        for(j = 0; j < pgaData.vectorSize; j++)
            v1[j] += vals[i] * pgaData.v[j + i * pgaData.vectorSize];


    for(i = 0; i < primitiveCount; i++)
    {
        prim = object->getPrimitivePtr(i);
        meanPrim = pgaData.mean->getPrimitivePtr(i);

        if(prim == NULL || meanPrim == NULL)
            continue;

        pos.set(v1[i*NUM_PGA_PARAMS], v1[i*NUM_PGA_PARAMS + 1], v1[i*NUM_PGA_PARAMS + 2]);

        prim->setX(meanPrim->getX() + pos);
        prim->setR(meanPrim->getR() * exp(v1[i*NUM_PGA_PARAMS + 3]));

        n0 = meanPrim->getNormalizedY0();
        n1 = meanPrim->getNormalizedY1();
        q0 = ShapeSpace::S2::rotationFromOrigin(n0);
        q1 = ShapeSpace::S2::rotationFromOrigin(n1);

        n0 = ShapeSpace::S2::Exp(Vector2D(v1[i*NUM_PGA_PARAMS + 6], v1[i*NUM_PGA_PARAMS + 7]));
        n1 = ShapeSpace::S2::Exp(Vector2D(v1[i*NUM_PGA_PARAMS + 4], v1[i*NUM_PGA_PARAMS + 5]));
        q0.rotateVector(n0);
        q1.rotateVector(n1);

        symToLieAtom(q, theta, n1, n0);
        prim->setQ(q);
        prim->setTheta(theta);

        if(prim->type() == M3D_END_PRIMITIVE)
            (dynamic_cast<M3DQuadEndPrimitive *>(prim))->setElongation((dynamic_cast<M3DQuadEndPrimitive *>(meanPrim))->getElongation() * exp(v1[i*NUM_PGA_PARAMS + 8]));
    }

    delete [] v1;

    return true;
}


/**
 * This function is deprecated. Avoid it's use.
 */
void M3DPGA::symToLieAtom(Quat & q, double & theta, Vector3D n0, Vector3D n1)
{
//	printf("M3DPGA::symToLieAtom called. This function is deprecated.\n");
//	assert(false);

	Vector3D b, n, bPerp;

    // Choose n to be the difference of the two spokes
	// (n is in the same direction as Y1, not Y0?)

    if (n1 == n0) {
		printf("Two spokes identical in M3DPGA::symToLieAtom!\n");

        // Handle the case where n0 = n1
        if ( fabs(n0.getX()) < fabs(n0.getY()) )
            if ( fabs(n0.getX()) < fabs(n0.getZ()) )
                n.set(1.0, 0.0, 0.0);
            else
                n.set(0.0, 0.0, 1.0);
        else
            if ( fabs(n0.getY()) > fabs(n0.getZ()) )
                n.set(0.0, 0.0, 1.0);
            else
                n.set(0.0, 1.0, 0.0);
		n = n - (n * n0) * n0;
    }
    else
        n = n0 - n1;
    n.normalize();

    b = n0 + n1;
    // Handle the case where n0 = -n1
    if (b.norm() == 0) {
		printf("Two spokes back to back in M3DPGA::symToLieAtom!\n");

		// If n0 and n1 are back to back,
		// pick a set of local frames
        if ( fabs(n0.getX()) < fabs(n0.getY()) )
            if ( fabs(n0.getX()) < fabs(n0.getZ()) )
                b.set(1.0, 0.0, 0.0);
            else
                b.set(0.0, 0.0, 1.0);
        else
            if ( fabs(n0.getY()) > fabs(n0.getZ()) )
                b.set(0.0, 0.0, 1.0);
            else
                b.set(0.0, 1.0, 0.0);
    }
    b = b - (b * n) * n;
    b.normalize();

    bPerp = b.cross(n);

    theta = acos(b * n0);
    q = Quat(b, n, bPerp);
}


