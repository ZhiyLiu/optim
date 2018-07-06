#include <iostream>
#include "M3DSpokeLengthOptimizer.h"
#include "M3DSpokeLengthOptimizer2D.h"
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
		return false;

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
    if(argc < 3)
    {
        std::cout << "[Usage:] optimizeLength <model file> <image file> <dimension(optional)>" << std::endl;
        return -1;
    }

    char* modelName = argv[1];
    char* imagePath = argv[2];

    // read image
    Image3D* image3D = readImage(imagePath);
    ImageDistanceMap* binaryDistanceMap = new ImageDistanceMap(image3D);
    binaryDistanceMap->fromImage3D(image3D);
    if(argc == 4)
    {
        char* mode = argv[3];
        if(strcmp(mode, "2") == 0)
        {
            // 2D optimization
            M3DSpokeLengthOptimizer2D optimizer;

            // read model
            M3DObject* object = readModel(modelName);

            optimizer.setImage(binaryDistanceMap);
            optimizer.setObject(object);

            M3DObject* modelAfterOptimizer = new M3DObject;
            optimizer.perform(modelAfterOptimizer);

        }
    }
    else
    {
        M3DSpokeLengthOptimizer optimizer;

        // read model
        M3DObject* object = readModel(modelName);

        optimizer.setImage(binaryDistanceMap);
        optimizer.setObject(object);

        M3DObject* modelAfterOptimizer = new M3DObject;
        optimizer.perform(modelAfterOptimizer);

    }

    return 0;
}