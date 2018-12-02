/*
 *  Zhiyuan Liu
 * Feb 22, 2018
 * invoke optimizer of NEWUOA to optimize angles as well as length
 * It also contains experiment of random initialization vs. sophisticated initialization
 *
 */
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "M3DNEWUOAOptimizer.h"
//#include "itkImage.h"
//#include "itkImageFileReader.h"
#include "ImageDistanceMap.h"
#include "RAWImageFile.h"
#include "M3DObjectFile.h"
#include "WorldSystem.h"
#include "AllImageIO.h"
#include <time.h>

M3DObject* readModel(const char * fileName)
{
        M3DObjectFile objectFile;
        SimilarityTransform3D * xform;
        M3DObject * obj;
        WorldSystem * world;

        if (fileName == NULL || strlen(fileName) == 0)
                return NULL;

        xform = new SimilarityTransform3D;
        world = new WorldSystem;
    int markedPrimitiveId = 0;
        obj = objectFile.read(fileName, markedPrimitiveId,
                xform, world, NULL);

        // cout << "Model file read successfully" << endl ;

        if (obj != NULL) {
        obj->setWorld(world, true);	// Sets world in object and object->that
        }
        else {
        delete xform;
                delete world;
                return NULL;
        }

        if (xform->wasRead()) {
#ifdef DEBUG
                cout << "Similarity transformation loaded" << endl;
#endif
        }
        else {			// The transformation in object is set to identity
                cout << "Similarity transformation not loaded" << endl;
        }
    obj->setTransformation(xform);
        obj->loadedObject()->setTransformation(xform);

        return obj;
}

Image3D* readImage(const char* fileName)
{
    AllImageIO file;

        if (fileName == NULL || strlen(fileName) == 0)
                return NULL;

        Image3D* image = file.read(fileName, false, false);
    image->setIsImageStacked(false);
    return image;
}

void randomizeModel(M3DObject* model)
{
    /* initialize random seed: */
    srand(time(NULL));
    M3DFigure* figure = model->getFigurePtr(0);
    int primitiveCount = figure->getPrimitiveCount();
    int spokeIndex = 0;
    for(int i = 0; i < primitiveCount; ++i)
    {
        // The i-th primitive
        M3DPrimitive* currentPrimitive = figure->getPrimitivePtr(i);
        M3DQuadEndPrimitive* endPrimitive = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);

        if(endPrimitive == NULL)
        {
            // it is an standard primitive
            M3DQuadPrimitive* quadAtom = dynamic_cast<M3DQuadPrimitive*>(currentPrimitive);

            // update U0 and U1 in quadAtom
//            Vector3D vRandom1(rand() % 100, rand() % 100, rand() % 100);
//            Vector3D vRandom2(-vRandom1);
//            vRandom1.normalize();
//            vRandom2.normalize();
//            quadAtom->setU0(vRandom1);
//            quadAtom->setU1(vRandom2);

            // update R0 and R1
            double r0 = (rand() % 100) / (double)100;
            double r1 = (rand() % 100) / (double)100;
            quadAtom->setR0(r0);
            quadAtom->setR1(r1);

            // update spokeIndex
            spokeIndex += 2;
        }
        else
        {
            M3DQuadPrimitive* standardVersion = dynamic_cast<M3DQuadPrimitive*> (currentPrimitive);

            // compute new U0 and U1
//            Vector3D vRandom1(rand() % 100, rand() % 100, rand() % 100);
//            Vector3D vRandom2(rand() % 100, rand() % 100, rand() % 100);
//            vRandom1.normalize();
//            vRandom2.normalize();

//            // update U0 and U1 in quadAtom
//            standardVersion->setU0(vRandom1);
//            standardVersion->setU1(vRandom2);
            // update R0 and R1
            double r0 = (rand() % 100) / (double)100;
            double r1 = (rand() % 100) / (double)100;
            double r2 = (rand() % 100) / (double)100;
            standardVersion->setR0(r0);
            standardVersion->setR1(r1);

            M3DQuadEndPrimitive* endVersion = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);

            // update UEnd
//            Vector3D vRandom3(rand() % 100, rand() % 100, rand() % 100);
//            vRandom3.normalize();
//            endVersion->setUEnd(vRandom3);
            // update REnd
            endVersion->setREnd(r2);

            spokeIndex += 3;
        }
    }

}

int main(int argc, char** argv)
{
    if(argc < 10)
    {
        std::cout << "[Usage:] optimizeSrep [modelName] [imagePath] [wtImageMatch] [wtNormalPenalty] [wtSradPenalty] [stepSize] [endCriterion] [maxIterations] [outputFileName]" << std::endl;
        return -1;
    }

    char* modelName = argv[1];
    char* imagePath = argv[2];
    char* sWtImageMatch = argv[3];
    char* sWtNormalMatch = argv[4];
    char* sWtSradPenalty = argv[5];
    char* sStepSize = argv[6];
    char* sEndCriterion = argv[7];
    char* sMaxIterations = argv[8];
    char* sOutputFileName = argv[9];
    char* sRandomize = NULL;

    if(argc == 11)
    {
        sRandomize = argv[10];
    }

    // read image
    Image3D* image3D = readImage(imagePath);
    ImageDistanceMap* binaryDistanceMap = new ImageDistanceMap(image3D);
    binaryDistanceMap->fromImage3D(image3D);
    M3DNEWUOAOptimizer optimizer;

    // read model
    M3DObject* object = readModel(modelName);

    // if randomize the initial model
    if(sRandomize != NULL && atoi(sRandomize) == 1)
    {
        randomizeModel(object);
        M3DObjectFile outputRandomFile;
        outputRandomFile.write("random_model.m3d", *object);
    }

    optimizer.setImage(binaryDistanceMap);
    optimizer.setObject(object);

    M3DObject* modelAfterOptimizer = new M3DObject;
    double wtImageMatch = atof(sWtImageMatch);
    double wtNormalMatch = atof(sWtNormalMatch);
    double wtSrepMatch = atof(sWtSradPenalty);
    double stepSize = atof(sStepSize);
    double endCriterion = atof(sEndCriterion);
    int maxIterations = atoi(sMaxIterations);
    std::cout << "Parameter wtImageMatch:" << wtImageMatch << ", wtNormalMatch:" << wtNormalMatch << ", wtSrepMatch:"
              << wtSrepMatch << ", stepSize:" <<stepSize << ", endCriterion:" << endCriterion << ", maxIterations:"
              << maxIterations << std::endl;
    time_t lt = time(NULL);
    FILE *logFile = fopen("newuoa.log", "w+");
    optimizer.setLogFilePtr(logFile);
    fprintf(logFile, "===========start to optimize at time:%s===========\n", asctime(localtime(&lt)));
    optimizer.perform(modelAfterOptimizer, wtImageMatch, wtNormalMatch,
                      wtSrepMatch, stepSize, endCriterion, maxIterations,sOutputFileName);
    fprintf(logFile, "\noutput to file: %s\n", sOutputFileName);
    fprintf(logFile, "===========end of optimize at time:%s===========\n", asctime(localtime(&lt)));
    return 0;
}
