/* The purpose of this class is to apply tps on new srep representation
 *
 * Zhiyuan Liu
 * 2018.05
 */
#include <iostream>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkSphereSource.h>
#include <vtkMetaImageWriter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataReader.h>
#include "thinplatesplinepdmtosrep.h"
//#include "itkTransformFactoryBase.h"
//#include "itkTransformFactory.h"
//#include "itkTransformFileReader.h"
#include "itkThinPlateSplineExtended.h"

int main()
{
    // Get the filename from the command line
//    std::string inputFilename = argv[1];

    // Get all data from the file up spokes
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName("../up.vtp");//inputFilename.c_str());
    reader->Update();
    vtkPolyData* upSpokes = reader->GetOutput();

    // Get down spokes
    vtkSmartPointer<vtkXMLPolyDataReader> reader2 = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader2->SetFileName("../down.vtp");//inputFilename.c_str());
    reader2->Update();
    vtkPolyData* downSpokes = reader2->GetOutput();

    // Get crest
    vtkSmartPointer<vtkXMLPolyDataReader> reader3 = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader3->SetFileName("../crest.vtp");//inputFilename.c_str());
    reader3->Update();
    vtkPolyData* crestSpokes = reader3->GetOutput();

    typedef double CoordinateRepType;
	typedef itkThinPlateSplineExtended TransformType;
	typedef itk::Point< CoordinateRepType, 3 > PointType;
	typedef std::vector< PointType > PointArrayType;
	typedef TransformType::PointSetType PointSetType;
	typedef PointSetType::Pointer PointSetPointer;
	typedef PointSetType::PointIdentifier PointIdType;


	thinplatesplinepdmtosrep obj;
	std::string templateModelName=argv[1];
	std::string transformFileName=argv[2];
	std::string sourceLandMarkFileName=argv[3];
	std::string outputPrefix=argv[4];

	ifstream inFile;
	inFile.open(transformFileName.c_str());
	if(!inFile) {
		std::cerr << "Unable to open the file: " << transformFileName << std::endl;
		return EXIT_FAILURE;
	}
	itkThinPlateSplineExtended::DMatrixType D;
	itkThinPlateSplineExtended::AMatrixType A;
	itkThinPlateSplineExtended::BMatrixType B;
	// first read in the size
	string buffer;
	std::getline(inFile,buffer);
	std::istringstream ss1(buffer);
	int nRows = 0;
	int nCols = 0;
	char tmp;
	ss1 >> nRows >> tmp >> nCols;
	D.set_size(nRows,nCols);
	for(int i = 0; i < nRows; i++) {
		for(int j = 0; j < nCols; j++) {
			string buffer2;
			std::getline(inFile, buffer2, ',');
			double entry = atof(buffer2.c_str());
			D(i,j) = entry;
		}
	}
	buffer.clear();
	std::getline(inFile, buffer);
	std::getline(inFile, buffer);
	for(int i = 0; i < A.rows(); i++) {
		for(int j = 0; j < A.cols(); j++) {
			string buffer2;
			std::getline(inFile, buffer2,',');
			double entry = atof(buffer2.c_str());
			A(i,j) = entry;
		}
	}
	buffer.clear();
	std::getline(inFile, buffer);
	std::getline(inFile, buffer);
	for(int i = 0; i < B.size(); i++) {
		string buffer2;
		std::getline(inFile, buffer2, ',');
		double entry = atof(buffer2.c_str());
		B(i) = entry;
	}
	TransformType::Pointer tps = TransformType::New();
	tps->setDMatrix(D);
	tps->setAMatrix(A);
	tps->setBVector(B);

	PointSetType::Pointer sourceLandMarks = PointSetType::New();
	PointSetType::PointsContainer::Pointer sourceLandMarkContainer = sourceLandMarks->GetPoints();
	vtkSmartPointer<vtkPolyDataReader> reader_source = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkSmartPointer<vtkPolyData> polyData_source = vtkSmartPointer<vtkPolyData>::New();
	reader_source->SetFileName(sourceLandMarkFileName.c_str());
	reader_source->Update();
	polyData_source = reader_source->GetOutput();
	PointIdType id_s = itk::NumericTraits< PointIdType >::Zero;
	PointType p1;
	// Read in the source points set
	for(unsigned int i = 0; i < polyData_source->GetNumberOfPoints(); i += 10){
		double p[3];
		polyData_source->GetPoint(i,p);
		// This is identical to:
		// polydata->GetPoints()->GetPoint(i,p);
		p1[0] = p[0];
		p1[1] = p[1];
		p1[2] = p[2];
		sourceLandMarkContainer->InsertElement(id_s, p1);
		id_s++;
	}
	tps->SetSourceLandmarks(sourceLandMarks);

    // iterate up spokes
    vtkPointData* ptData = upSpokes.GetPointData();
    for(int i = 0; i < upSpokes.GetNumberOfPoints(); ++i)
    {

    }
    return 1;
}