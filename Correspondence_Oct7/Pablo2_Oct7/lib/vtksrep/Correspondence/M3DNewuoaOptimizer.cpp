/* This is an wrapper of optimizer using NEWUOA
 *
 * Zhiyuan Liu
 * 2017.11
 */
#include "M3DNewuoaOptimizer.h"

#include <iostream>
#include <exception>
#include "newuoa.h"
#include "toolsfunc.h"
#include "M3DQuadFigure.h"
#include "movespokes.h"
#include "M3DAtomPredictorQuad.h"
#include "SimilarityComputer.h"
#include "M3DInterpolater.h"

M3DNewuoaOptimizer::M3DNewuoaOptimizer()
    :mSreps(NULL){}

M3DNewuoaOptimizer::M3DNewuoaOptimizer(M3DObject* sreps)
{
    if(sreps != NULL)
    {
        mSreps = sreps;
    }
}

double M3DNewuoaOptimizer::computeSradPenalty()
{
    double sradPenalty = 0.0;
    // assume figure is quad figure
    M3DAtomPredictorQuad *atomPredictor = new M3DAtomPredictorQuad();
    M3DFigure *figure = mSreps->getFigurePtr(mFigureIndex);
    sradPenalty = atomPredictor->getFigureRSradPenalty(figure, -1, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold));
    return sradPenalty;
}

void M3DNewuoaOptimizer::interpolateSRep(std::vector<M3DSpoke> *outputSpokes)
{
    if(outputSpokes == NULL)
    {
        return;
    }
    outputSpokes->clear();
    // interpolation level from config file, currently use level 2
    int interpolationLevel = 2;
    std::vector<double> steps;
    steps.push_back(.25);
    steps.push_back(.75);

    M3DFigure *figure = mSreps->getFigurePtr(mFigureIndex);
    M3DInterpolater interpolater(figure);

    M3DQuadFigure* quadFig = dynamic_cast<M3DQuadFigure*>( figure );

    int umax = quadFig->getRowCount() - 1;
	int vmax = quadFig->getColumnCount() - 1;

	int u,v,i, j;
	for (u = 0; u < umax; u++)
    {

		for (v = 0; v < vmax; v++)
		{
			M3DQuadPrimitive* quadAtom = dynamic_cast<M3DQuadPrimitive*>( quadFig->getPrimitivePtr(u, v) );
			Vector3D X = quadAtom->getX();
			Vector3D U0 = quadAtom->getU0();
			Vector3D U1 = quadAtom->getU1();
			double R0 = quadAtom->getR0();
			double R1 = quadAtom->getR1();

			M3DSpoke top = M3DSpoke(X, U0, R0);
			M3DSpoke bottom = M3DSpoke(X, U1, R1);

			outputSpokes->push_back(top);
			outputSpokes->push_back(bottom);

			for (i = 0; i < steps.size(); i++)
			{
				for(j = 0; j < steps.size(); j++)
				{
					outputSpokes->push_back(*(interpolater.interpolateSpoke(figure, u + steps[i], v + steps[j], 0)));
					outputSpokes->push_back(*(interpolater.interpolateSpoke(figure, u + steps[i], v + steps[j], 1)));

					outputSpokes->push_back(*(interpolater.interpolateSpoke(figure, u, v + steps[j], 0)));
					outputSpokes->push_back(*(interpolater.interpolateSpoke(figure, u, v + steps[j], 0)));

					outputSpokes->push_back(*(interpolater.interpolateSpoke(figure, u + steps[i], v, 1)));
				}
            }
		}
    }

    interpolateCrestSpokes(interpolationLevel, outputSpokes);
 }

void M3DNewuoaOptimizer::generateVtkSrep(M3DQuadFigure* quadfig, vtkSmartPointer<vtkSRep>& srepfig)
{
    if(quadfig == NULL)
    {
        std::cout << "quadfig = nullptr in M3DNewuoaOptimizer::generateVtkSrep" << std::endl;
        return;
    }

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();

    vtkSRep::VectorSRepIdsType pointsIds;

    vtkSRep::RadiusVectorType allradius;
    vtkSRep::SpokesVectorType allspokes;

    vtkIdType vtkAtomId = -1;

    for(int u = 0; u < quadfig->getRowCount(); u++){
        pointsIds.push_back(vtkSRep::VectorIdsType());
        for(int v = 0; v < quadfig->getColumnCount(); v++){

            M3DQuadPrimitive* prim0 = dynamic_cast<M3DQuadPrimitive*>(quadfig->getPrimitivePtr(u, v));

            Vector3D x = prim0->getX();
            Vector3D u0 = prim0->getU0();
            Vector3D u1 = prim0->getU1();

            vtkSRep::VectorVNLType vnlspokes;
            vtkSRep::VNLType s(3);
            s[0] = u0.getX();
            s[1] = u0.getY();
            s[2] = u0.getZ();
            vnlspokes.push_back(s);

            s[0] = u1.getX();
            s[1] = u1.getY();
            s[2] = u1.getZ();
            vnlspokes.push_back(s);

            vtkSRep::VectorDoubleType radius;
            radius.push_back(prim0->getR0());
            radius.push_back(prim0->getR1());

            if(u == 0 || u == quadfig->getRowCount() - 1 || v == 0 || v == quadfig->getColumnCount() - 1){

                M3DQuadEndPrimitive* prim0 = dynamic_cast<M3DQuadEndPrimitive*>(quadfig->getPrimitivePtr(u, v));
                Vector3D uend = prim0->getUEnd();

                s[0] = uend.getX();
                s[1] = uend.getY();
                s[2] = uend.getZ();

                vnlspokes.push_back(s);

                radius.push_back(prim0->getREnd());

            }

            int numCols = quadfig->getColumnCount();
            pointsIds[u].push_back(hubpos->InsertNextPoint(x.getX(), x.getY(), x.getZ()));

            allspokes.push_back(vnlspokes);
            allradius.push_back(radius);
        }
    }


    for(unsigned i = 0; i < pointsIds.size() - 1; i++){
         for(unsigned j = 0; j < pointsIds[i].size() - 1; j++){

             vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
             quad->GetPointIds()->SetId(0, pointsIds[i][j]);
             quad->GetPointIds()->SetId(1, pointsIds[i+1][j]);
             quad->GetPointIds()->SetId(2, pointsIds[i+1][j+1]);
             quad->GetPointIds()->SetId(3, pointsIds[i][j+1]);

             //quad->Print(cout);

             cellarray->InsertNextCell(quad);

         }
     }

    srepfig->SetPoints(hubpos);
    srepfig->SetPolys(cellarray);
    srepfig->SetAllSpokes(allspokes);
    srepfig->SetAllRadius(allradius);
//    srepfig->SetGridTopolgyIds(pointsIds);

    for(unsigned i = 0; i < pointsIds.size(); i++){
        pointsIds[i].clear();
    }
    pointsIds.clear();
    for(unsigned i = 0; i < allradius.size(); i++){
        allradius[i].clear();
        allspokes[i].clear();
    }
    allradius.clear();
    allspokes.clear();
}
void M3DNewuoaOptimizer::interpolateCrestSpokes(int interpolationLevel, std::vector<M3DSpoke>* crestSpokes)
{
    // step 1. generate vtk srep
    if(mSreps == NULL)
    {
        std::cout << "mSreps = nullptr in M3DNewuoaOptimizer::interpolateCrestSpokes" << std::endl;
        return;
    }
    vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
    M3DFigure * figure = mSreps->getFigurePtr(mFigureIndex) ;
    M3DQuadFigure* quadfig = dynamic_cast<M3DQuadFigure*>(figure);
    generateVtkSrep(quadfig, srepfig);
    
    // step 2. interpolate using vtk interface
    vtkIdType vtkAtomId = -1;
    int spokeId = -1;
    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
    interpolatecrestspokes->SetInput(srepfig);
    interpolatecrestspokes->SetInterpolationLevel(interpolationLevel);
//    interpolatecrestspokes->SetAtomId(vtkAtomId);

//    interpolatecrestspokes->SetGamma_t(1);
//    interpolatecrestspokes->SetGamma_theta(1);

    interpolatecrestspokes->SetSpokeType(vtkSRep::CREST_SPOKE);

    interpolatecrestspokes->Update();
    vtkSRep* srepcrest = interpolatecrestspokes->GetSRepOutput();


    for(unsigned i = 0; i < srepcrest->GetNumberOfPoints(); i++){

        double point[3];

        srepcrest->GetPoint(i, point);

        vtkSRep::VectorVNLType currentspokes = srepcrest->GetSpokes(i);
        vtkSRep::VectorDoubleType radius = srepcrest->GetSpokesRadius(i);

        for(unsigned j = 0; j < currentspokes.size(); j++){
            M3DSpoke spoke = M3DSpoke(point[0], point[1], point[2], currentspokes[j][0], currentspokes[j][1], currentspokes[j][2], radius[j]);

            crestSpokes->push_back(spoke);
        }
    }

}

int iterNum = 0;
/* Update spoke properties after each optimization */
void M3DNewuoaOptimizer::updateFigure(const double *coeff, int figureId)
{
    try
    {
        std::cout << "Updating figure id: " << figureId << ", after the " << iterNum++ << "th iteration" << std::endl;
        
        M3DFigure* figure = mSreps->getFigurePtr(figureId);
        int spokeCount = figure->getSpokeCount();
        for(int iSpoke = 0; iSpoke < spokeCount; ++iSpoke)
        {
            std::cout << "The length of the " << iSpoke + 1 << "th spoke is:" << exp(coeff[iSpoke]) << std::endl;
        }
        int primitiveCount = figure->getPrimitiveCount();
        for(int i = 0; i < primitiveCount; ++i)
        {
            // The i-th primitive
            M3DPrimitive* currentPrimitive = figure->getPrimitivePtr(i);
            M3DQuadEndPrimitive* endPrimitive = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);
            if(endPrimitive == NULL)
            {
                // it is an standard primitive
                M3DQuadPrimitive* primitive = dynamic_cast<M3DQuadPrimitive*>(currentPrimitive);
                double r0 = primitive->getR0();
                primitive->setR0(r0 * exp(coeff[i]));

                double r1 = primitive->getR1();
                primitive->setR1(r1 * exp(coeff[i]));
             }
            else
            {
                M3DQuadPrimitive* standardVersion = dynamic_cast<M3DQuadPrimitive*> (currentPrimitive);
                double r0 = standardVersion->getR0();
                standardVersion->setR0(r0 * exp(coeff[i]));

                double r1 = standardVersion->getR1();
                standardVersion->setR1(r1 * exp(coeff[i]));

                M3DQuadEndPrimitive* endVersion = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);
                double r = endVersion->getREnd();
                endVersion->setREnd(r * exp(coeff[i]));
            }
        }
    }
    catch(std::exception& e)
    {
        std::cout << "[Error]updateFigure has an exception:"<< e.what() << std::endl;
    }
}
/* Compute the total entropy. */
double M3DNewuoaOptimizer::getObjectiveFunctionValue(const double *coeff, double w1, double w2)
{
    // input coeff from newuoa is the coeff of length. it could be positive or negative

    // 0. Update new length to each spoke
    updateFigure(coeff, mFigureIndex);
    
    double objFunctionValue = 0.0;
    // coeff are now lengths of spokes

    // 1. Image match
    // TODO: should read from config file
    double w_ImageMatch = 9.0;
    double w_sradPenalty = 1.0;
    double w_srepModelPenalty = 9.0;  // may be too large, was 10

    // 1.1 Interpolate spokes
    interpolateSRep(&mSpokesAfterInterp);

    // 1.2 measure sum square of distance from implied boundary to expected image boundary
    SimilarityComputer similarityCompter;
    similarityCompter.setSrepModel(mSpokesAfterInterp);
    similarityCompter.setTargetImage(mSignedDistanceImage);

    double imageMatch = 0.0;
    if(similarityCompter.compute(&imageMatch) == false)
    {
        std::cout << "[Error]Error encountered when compute similarity measure" << std::endl;
        return -999.0;
    }

    // 2. measure how far from regular srep model
    double sradPenalty = computeSradPenalty();

    // 3. spoke model penalty
    DistanceType distType = (enum DistanceType) (int) tuningWt(BpSpokeDistanceType);
    double spokeModelpenalty = mSreps->dist2FromObject(mSreps->loadedObject(),
				mFigureIndex, distType);

    // 4. image normal match
//    double imageNormalMatch = 
    objFunctionValue = w_ImageMatch * imageMatch + w_sradPenalty * sradPenalty
        + w_srepModelPenalty * spokeModelpenalty;
    
    return objFunctionValue;
}
double M3DNewuoaOptimizer::operator() (double *coeff)
{
    double cost = 0.0;
    cost = this->getObjectiveFunctionValue(coeff, 15, 1);
    return cost;
}

void M3DNewuoaOptimizer::setImage(ImageDistanceMap* imageDistanceMap)
{
    mSignedDistanceImage = imageDistanceMap;
}

void M3DNewuoaOptimizer::setObject(M3DObject* sreps)
{
    mSreps = sreps;
}
/* Main entry of optimizer */
int M3DNewuoaOptimizer::perform(M3DObject* outputModel)
{
    if(mSreps == NULL || mSignedDistanceImage == NULL)
    {
        std::cout << "[Error]object or signed distance image is empty in perform function" << std::endl;
        return -1;
    }
    
    // optimization
    M3DFigure* figure = mSreps->getFigurePtr(mFigureIndex);
    int spokeCount = figure->getSpokeCount();
    double *coeffOfLength = new double[spokeCount];  // optimize length's coefficients instead of length itself, to avoid negative length output by newuoa

    for(unsigned int i =0;i<spokeCount;i++){
        coeffOfLength[i] = 0.0;
    }
    min_newuoa(spokeCount,coeffOfLength,*this,0.1, 0.000001, 30000);

    // coeffOfLength here is variable at optimum
    outputModel = mSreps;
    delete[] coeffOfLength;
    return 0;

}
