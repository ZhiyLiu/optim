#ifndef METHOD_OF_MOMENTS
#define METHOD_OF_MOMENTS

// Initialization code for finding scales, centroids, and rotations 
//    in voxel images and m-rep single quad figure objects

class MethodOfMoments
{
public:
	#ifdef BINARY
	static void objectMoments(M3DObject * object, double & volume, Vector3D & center, Quat & rotation);
	#endif
	static void figureMoments(M3DFigure * figure, double & volume, Vector3D & center, Quat & rotation);
	static void imageMoments(Image3D * image, double & volume, Vector3D & center, Quat & rotation);

	// Helpers
	static void covarianceToQuat(Matrix covariance, Quat & quat);
	static void leastAngleRotation(Quat rotation1, Quat rotation2, Quat & bestRotation);

	#ifdef BINARY
	static void positionObjectOnImage(M3DObject * object, Image3D * image, bool verbose = false);
	#endif
};

#endif

