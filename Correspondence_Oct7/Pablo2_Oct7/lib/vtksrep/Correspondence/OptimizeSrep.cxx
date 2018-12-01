/*
 *  Zhiyuan Liu
 * Feb 22, 2018
 * invoke optimizer of NEWUOA to optimize angles
 */
#include <iostream>
#include "M3DNEWUOAOptimizer.h"
//#include "itkImage.h"
//#include "itkImageFileReader.h"
#include "ImageDistanceMap.h"
#include "RAWImageFile.h"
#include "M3DObjectFile.h"
#include "WorldSystem.h"
#include "AllImageIO.h"

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

    // read image
    Image3D* image3D = readImage(imagePath);
    ImageDistanceMap* binaryDistanceMap = new ImageDistanceMap(image3D);
    binaryDistanceMap->fromImage3D(image3D);
    M3DNEWUOAOptimizer optimizer;

    // read model
    M3DObject* object = readModel(modelName);

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

    optimizer.perform(modelAfterOptimizer, wtImageMatch, wtNormalMatch,
                      wtSrepMatch, stepSize, endCriterion, maxIterations,sOutputFileName);

    return 0;
}
