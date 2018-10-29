/* The purpose of this class is to compute similarity between two images
 *
 * Zhiyuan Liu
 * 2017.12
 */
#include "SimilarityComputer.h"
#include "M3DObject.h"
#include "ImageDistanceMap.h"
//#include "toolsfunc.h"
#include "M3DQuadFigure.h"
#include "movespokes.h"
#include "M3DAtomPredictorQuad.h"
#include "M3DObjectFile.h"

SimilarityComputer::SimilarityComputer():
    mDistanceImage(NULL), mFigure(NULL){}

SimilarityComputer::~SimilarityComputer()
{
    if(mFigure != NULL)
    {
        delete mFigure;
        mFigure = NULL;
    }
}
void SimilarityComputer::setTargetImage(ImageDistanceMap* target)
{
    mDistanceImage = target;
}

void SimilarityComputer::setSrepModel(M3DObject* sreps)
{
    mSreps = sreps;
}

void SimilarityComputer::setTargetFigureIndex(int figureId)
{
    mFigureIndex = figureId;
}

void SimilarityComputer::setFigure(M3DQuadFigure* f)
{
    if(mFigure != NULL)
    {
        delete mFigure;
        mFigure = NULL;
    }
    mFigure = new M3DQuadFigure(*f);
}

double SimilarityComputer::getSpokeNormalMatch(std::vector<M3DSpoke*>& spokes)
{
    double result = 0.0;
    for(int i = 0; i < spokes.size(); ++i)
    {
        Vector3D gradDir = mDistanceImage->getGradDistance(spokes[i]->getB());
        Vector3D spokeDir = spokes[i]->getU();
        double spokeLength = spokes[i]->getR();
        spokeDir.normalize();
        gradDir.normalize();
        double dotProduct = spokeDir * gradDir;

//        if(dotProduct < -1)
//        {
//            dotProduct = -1;
//        }
//        else if(dotProduct > 1)
//        {
//            dotProduct = 1;
//        }
        result = result + spokeLength*spokeLength*(1-dotProduct);// sum of r^2(1-UD)

    }
    result /= spokes.size();
    return result;
}

bool SimilarityComputer::compute(double *similarityMeasure, double *normalMatch)
{
    // step 1: validate data
    if(mDistanceImage == NULL)
    {
        std::cout << "[Error]SimilarityComputer requires both target image and reference model" << std::endl;
        return false;
    }
    double match = 0.0;
    double dilationFactor = 0.0;
    //toolsfunc tools;

    // step 2: for each spoke, compute image match nearby interpolated neighbors.
    M3DFigure* figure = mSreps->getFigurePtr(mFigureIndex);
    int spokeCount = figure->getSpokeCount();
    // get all the spokes
    int primitiveCount = figure->getPrimitiveCount();
    double totalDistance = 0.0;
    double totalNormalMatch = 0.0;
    double maxNormalMatch = -100.0;
    for(int i = 0; i < primitiveCount; ++i)
    {
        // The i-th primitive
        M3DPrimitive* currentPrimitive = figure->getPrimitivePtr(i);
        M3DQuadEndPrimitive* endPrimitive = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);
        M3DQuadPrimitive* quadAtom = dynamic_cast<M3DQuadPrimitive*>(currentPrimitive);
        Vector3D X = quadAtom->getX();

        Vector3D U0 = quadAtom->getU0();
        Vector3D U1 = quadAtom->getU1();
        double R0 = quadAtom->getR0();
        double R1 = quadAtom->getR1();

        M3DSpoke top = M3DSpoke(X, U0, R0);
        M3DSpoke bottom = M3DSpoke(X, U1, R1);
        std::vector<M3DSpoke*> topNeighbors, botNeighbors;
        getRelevantSpokes(figure, quadAtom, i, 0, &topNeighbors);
        getRelevantSpokes(figure, quadAtom, i, 1, &botNeighbors);

        double topDistance = getSSD(topNeighbors);
        double botDistance = getSSD(botNeighbors);
        totalDistance += topDistance + botDistance;
//        std::vector<M3DSpoke*> topArray, botArray; // solution that assume angles are correct
//        topArray.push_back(new M3DSpoke(X,U0, R0));
//        botArray.push_back(new M3DSpoke(X,U1, R1));
//        totalDistance = getSSD(topArray) + getSSD(botArray);

        double topNormalMatch = getSpokeNormalMatch(topNeighbors);
        double botNormalMatch = getSpokeNormalMatch(botNeighbors);
        totalNormalMatch += topNormalMatch + botNormalMatch;
        //maxNormalMatch = max(topNormalMatch + botNormalMatch, maxNormalMatch);

        topNeighbors.clear();
        botNeighbors.clear();
        if(endPrimitive != NULL)
        {
            M3DQuadEndPrimitive* endVersion = dynamic_cast<M3DQuadEndPrimitive*>(currentPrimitive);
            double r = endVersion->getREnd();
            Vector3D UEnd = endVersion->getUEnd();
            M3DSpoke endSpoke = M3DSpoke(X, UEnd, r);

            std::vector<M3DSpoke*> endNeighbors;
            getRelevantSpokes(figure, quadAtom, i, 2, &endNeighbors);
            double endDistance = getSSD(endNeighbors);
            // totalDistance = max(totalDistance, endDistance);
            totalDistance+= endDistance;

            totalNormalMatch += getSpokeNormalMatch(endNeighbors);
            maxNormalMatch = max(maxNormalMatch, getSpokeNormalMatch(endNeighbors));
            endNeighbors.clear();
        }

    }

    *similarityMeasure = totalDistance;
    *normalMatch = totalNormalMatch;

    return true;
}
double SimilarityComputer::getSSD(std::vector<M3DSpoke*>& spokes)
{
    double result = 0.0;
    for(int i = 0; i < spokes.size(); ++i)
    {
        double d = fabs(mDistanceImage->getDistance(spokes[i]->getB()));
        result += d * d;
    }

    result /= spokes.size();
    return result;
}
void SimilarityComputer::getRelevantSpokes(M3DFigure* f, M3DQuadPrimitive* quadAtom, int atomId, int spokeId, std::vector<M3DSpoke*>* outSpokes)
{
        M3DQuadFigure* figure = dynamic_cast<M3DQuadFigure*>(f);
        int umax = figure->getRowCount() - 1;
        int vmax = figure->getColumnCount() - 1;
        int numCols = figure->getColumnCount();

        int u = atomId / numCols;
        int v = atomId % numCols;

        Vector3D X = quadAtom->getX();
        Vector3D U;
        double R;
        if (spokeId == 0)
        {
            U = quadAtom->getU0();
            R = quadAtom->getR0();
        }
        else if (spokeId == 1)
        {
            U = quadAtom->getU1();
            R = quadAtom->getR1();
        }
        else
        {
            U = dynamic_cast<M3DQuadEndPrimitive*>(quadAtom)->getUEnd();
            R = dynamic_cast<M3DQuadEndPrimitive*>(quadAtom)->getREnd();
        }

        M3DSpoke* thisSpoke = new M3DSpoke(X, U, R);
        outSpokes->push_back(thisSpoke);

        int i, j;

        vector<double> steps;
        steps.push_back(0.25);
        steps.push_back(.75);
        if (spokeId == 0 || spokeId == 1)
        {

            for (i = 0; i < steps.size(); i++)
            {
                for(j = 0; j < steps.size(); j++)
                {
                    if (u == 0)
                    {
                        if (v == 0)
                        {
                            outSpokes->push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u + steps[i], v, spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u, v + steps[i], spokeId));
                        }
                        else if (v == vmax)
                        {
                            outSpokes->push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], spokeId));
                            //spokes.push_back(interpolateSpoke(figure, u + steps[i], v, spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u, v - steps[i], spokeId));
                        }
                        else
                        {
                            outSpokes->push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u + steps[i], v, spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u, v + steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u, v - steps[j], spokeId));
                        }
                    }
                    else if (u == umax)
                    {
                        if (v == 0)
                        {
                            outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v, spokeId));
                            //outSpokes->push_back(interpolateSpoke(figure, u, v + steps[i], spokeId));

                        }
                        else if (v == vmax)
                        {
                            outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], spokeId));
                            //outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v, spokeId));
                            //outSpokes->push_back(interpolateSpoke(figure, u, v - steps[i], spokeId));
                        }
                        else
                        {
                            outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v, spokeId));
                            //outSpokes->push_back(interpolateSpoke(figure, u, v + steps[j], spokeId));
                            //outSpokes->push_back(interpolateSpoke(figure, u, v - steps[j], spokeId));
                        }
                    }
                    else
                    {
                        if (v==0)
                        {
                            outSpokes->push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u + steps[i], v, spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v, spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u, v + steps[i], spokeId));

                        }
                        else if (v == vmax)
                        {
                            outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], spokeId));
                            //outSpokes->push_back(interpolateSpoke(figure, u + steps[i], v, spokeId));
                            //outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v, spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u, v - steps[i], spokeId));

                        }
                        else
                        {
                            outSpokes->push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u + steps[i], v, spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u - steps[i], v, spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u, v + steps[i], spokeId));
                            outSpokes->push_back(interpolateSpoke(figure, u, v - steps[i], spokeId));
                        }
                    }
                }
            }
        }

        if(spokeId == 2){
            vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
            vtkIdType vtkAtomId = GetVtkSrepFig(figure, srepfig);
            GetCrestSpokes(srepfig, 2, *outSpokes,false, vtkAtomId, spokeId);
        }
}

void SimilarityComputer::GetCrestSpokes(vtkSmartPointer<vtkSRep> srepfig, int level, vector<M3DSpoke*>& spokes,bool istube, vtkIdType vtkAtomId, int spokeId){

    if(false){//Set this to true to disable the crest interpolation, without crest interpolation, it should return only the base crest spokes
        level = 0;
    }

    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
    interpolatecrestspokes->SetInput(srepfig);
    interpolatecrestspokes->SetInterpolationLevel(level);
    interpolatecrestspokes->SetAtomId(vtkAtomId);
    if(vtkAtomId != -1){
        interpolatecrestspokes->SetGamma_t(0.5);
        interpolatecrestspokes->SetGamma_theta(0.5);
    }else{
        interpolatecrestspokes->SetGamma_t(1);
        interpolatecrestspokes->SetGamma_theta(1);
    }
    if(istube){        
        interpolatecrestspokes->SetCyclicCurve(false);
        interpolatecrestspokes->SetCyclicSpokes(true);
        interpolatecrestspokes->SetSpokeType(spokeId);
    }else{
        interpolatecrestspokes->SetCyclicCurve(true);
        interpolatecrestspokes->SetCyclicSpokes(false);
        if(spokeId == 0){
            interpolatecrestspokes->SetSpokeType(vtkSRep::TOP_SPOKE);
        }else if(spokeId == 1){
            interpolatecrestspokes->SetSpokeType(vtkSRep::BOTTOM_SPOKE);
        }else if(spokeId == 2){
            interpolatecrestspokes->SetSpokeType(vtkSRep::CREST_SPOKE);
        }
    }

    interpolatecrestspokes->Update();


    vtkSRep* srepcrest = interpolatecrestspokes->GetSRepOutput();


    for(unsigned i = 0; i < srepcrest->GetNumberOfPoints(); i++){

        double point[3];

        srepcrest->GetPoint(i, point);

        vtkSRep::VectorVNLType currentspokes = srepcrest->GetSpokes(i);        

        for(unsigned j = 0; j < currentspokes.size(); j++){
            vtkSRep::VNLType s = currentspokes[j];
            double radius = s.magnitude();
            s.normalize();
            M3DSpoke* spoke = new M3DSpoke(point[0], point[1], point[2], s[0], s[1], s[2], radius);

            spokes.push_back(spoke);
        }
    }

}

vtkIdType SimilarityComputer::GetVtkSrepFig(M3DFigure* figure, vtkSmartPointer<vtkSRep>& srepfig, int atomId){
    M3DQuadFigure* quadfig = dynamic_cast<M3DQuadFigure*>( figure );

    vtkIdType vtkAtomId = -1;

    if(quadfig){

        vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
        vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();

        vtkSRep::VectorSRepIdsType pointsIds;

        vtkSRep::RadiusVectorType allradius;
        vtkSRep::SpokesVectorType allspokes;

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
                if(atomId != -1 && atomId / numCols == u && atomId % numCols == v){
                    vtkAtomId = hubpos->InsertNextPoint(x.getX(), x.getY(), x.getZ());
                    pointsIds[u].push_back(vtkAtomId);
                }else{
                    pointsIds[u].push_back(hubpos->InsertNextPoint(x.getX(), x.getY(), x.getZ()));
                }


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
        //srepfig->SetGridTopolgyIds(pointsIds);

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

    }else{

        M3DTubeFigure* tubefig = dynamic_cast<M3DTubeFigure*>(figure);

            vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
            vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();

            vtkSRep::VectorSRepIdsType pointsIds;

            vtkSRep::RadiusVectorType allradius;
            vtkSRep::SpokesVectorType allspokes;



            pointsIds.push_back(vtkSRep::VectorIdsType());

            for(int u = 0; u < (int)tubefig->getPrimitiveCount(); u++){

                M3DTubePrimitive* prim0 = dynamic_cast<M3DTubePrimitive*>(tubefig->getPrimitivePtr(u));
                Vector3D x = prim0->getX();

                vtkIdType id = hubpos->InsertNextPoint(x.getX(), x.getY(), x.getZ());

                if(atomId == u){
                    vtkAtomId = id;
                }

                pointsIds[0].push_back(id);

                vtkSRep::VectorVNLType vnlspokes;
                vtkSRep::VectorDoubleType radius;

                for(int v = 0; v < (int)tubefig->getNumberOfSpokes(); v++){


                    Vector3D u0 = prim0->getYN(v);


                    vtkSRep::VNLType s(3);
                    s[0] = u0.getX();
                    s[1] = u0.getY();
                    s[2] = u0.getZ();

                    radius.push_back(s.magnitude());
                    vnlspokes.push_back(s.normalize());


                }

                allspokes.push_back(vnlspokes);
                allradius.push_back(radius);
            }


            for(unsigned i = 0; i < pointsIds.size(); i++){
                 for(unsigned j = 0; j < pointsIds[i].size() - 1; j++){

                     vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                     line->GetPointIds()->SetId(0, pointsIds[i][j]);
                     line->GetPointIds()->SetId(1, pointsIds[i][j+1]);

                     //quad->Print(cout);

                     cellarray->InsertNextCell(line);

                 }
             }

            srepfig->SetPoints(hubpos);
            srepfig->SetLines(cellarray);
            srepfig->SetAllSpokes(allspokes);
            srepfig->SetAllRadius(allradius);

    }

    return vtkAtomId;

}

M3DSpoke* SimilarityComputer::interpolateSpoke(M3DFigure *figure, double u, double v, int side)
{
    int ubase = (int)floor(u);
	int vbase = (int)floor(v);

	u = u - ubase;
	v = v - vbase;

	M3DQuadFigure *tempFigure = dynamic_cast<M3DQuadFigure*>(figure);

	// Get four corner atoms

	M3DPrimitive *atom11 = tempFigure->getPrimitivePtr(ubase,vbase);
	M3DQuadPrimitive* quadAtom11 = dynamic_cast<M3DQuadPrimitive*>( atom11 );

	M3DPrimitive *atom21 = tempFigure->getPrimitivePtr(ubase+1,vbase);
	M3DQuadPrimitive* quadAtom21 = dynamic_cast<M3DQuadPrimitive*>( atom21 );

	M3DPrimitive *atom12 = tempFigure->getPrimitivePtr(ubase,vbase+1);
	M3DQuadPrimitive* quadAtom12 = dynamic_cast<M3DQuadPrimitive*>( atom12 );

	M3DPrimitive *atom22 = tempFigure->getPrimitivePtr(ubase+1,vbase+1);
	M3DQuadPrimitive* quadAtom22 = dynamic_cast<M3DQuadPrimitive*>( atom22 );

//	cout << "Positions: " << endl;
//	cout << quadAtom11->getX().getX() << ", " << quadAtom11->getX().getY() << ", " << quadAtom11->getX().getZ() << endl;
//	cout << quadAtom21->getX().getX() << ", " << quadAtom21->getX().getY() << ", " << quadAtom21->getX().getZ() << endl;
//	cout << quadAtom12->getX().getX() << ", " << quadAtom12->getX().getY() << ", " << quadAtom12->getX().getZ() << endl;
//	cout << quadAtom22->getX().getX() << ", " << quadAtom22->getX().getY() << ", " << quadAtom22->getX().getZ() << endl;

//	cout << "Spokes: " << endl;
//	cout << quadAtom11->getU0().getX() << ", " << quadAtom11->getU0().getY() << ", " << quadAtom11->getU0().getZ() << endl;
//	cout << quadAtom21->getU0().getX() << ", " << quadAtom21->getU0().getY() << ", " << quadAtom21->getU0().getZ() << endl;
//	cout << quadAtom12->getU0().getX() << ", " << quadAtom12->getU0().getY() << ", " << quadAtom12->getU0().getZ() << endl;
//	cout << quadAtom22->getU0().getX() << ", " << quadAtom22->getU0().getY() << ", " << quadAtom22->getU0().getZ() << endl;

	// Get neighboring atoms for use in derivatives
	// If ubase/vbase is 0 or ubase+1/vbase+1 is max, cant use neighboring, use current instead
	// Currently using non-central differences for ease, need to fix later
	/*M3DPrimitive *atom10 = tempFigure->getPrimitivePtr(1,0);
	M3DQuadPrimitive* quadAtom10 = dynamic_cast<M3DQuadPrimitive*>( atom10 );

	M3DPrimitive *atom20 = tempFigure->getPrimitivePtr(2,0);
	M3DQuadPrimitive* quadAtom20 = dynamic_cast<M3DQuadPrimitive*>( atom20 );

	M3DPrimitive *atom31 = tempFigure->getPrimitivePtr(3,1);
	M3DQuadPrimitive* quadAtom31 = dynamic_cast<M3DQuadPrimitive*>( atom31 );

	M3DPrimitive *atom32 = tempFigure->getPrimitivePtr(3,2);
	M3DQuadPrimitive* quadAtom32 = dynamic_cast<M3DQuadPrimitive*>( atom32 );

	M3DPrimitive *atom23 = tempFigure->getPrimitivePtr(2,3);
	M3DQuadPrimitive* quadAtom23 = dynamic_cast<M3DQuadPrimitive*>( atom23 );

	M3DPrimitive *atom13 = tempFigure->getPrimitivePtr(1,3);
	M3DQuadPrimitive* quadAtom13 = dynamic_cast<M3DQuadPrimitive*>( atom13 );

	M3DPrimitive *atom02 = tempFigure->getPrimitivePtr(0,2);
	M3DQuadPrimitive* quadAtom02 = dynamic_cast<M3DQuadPrimitive*>( atom02 );

	M3DPrimitive *atom01 = tempFigure->getPrimitivePtr(0,1);
	M3DQuadPrimitive* quadAtom01 = dynamic_cast<M3DQuadPrimitive*>( atom01 );*/

	// Take the positions of the atoms and do simple finite difference derivative calculations
	/*Vector3D x10 = quadAtom10->getX();
	Vector3D x20 = quadAtom20->getX();
	Vector3D x31 = quadAtom31->getX();
	Vector3D x32 = quadAtom32->getX();
	Vector3D x23 = quadAtom23->getX();
	Vector3D x13 = quadAtom13->getX();
	Vector3D x02 = quadAtom02->getX();
	Vector3D x01 = quadAtom01->getX();*/

	Vector3D x11 = quadAtom11->getX();
	Vector3D x21 = quadAtom21->getX();
	Vector3D x12 = quadAtom12->getX();
	Vector3D x22 = quadAtom22->getX();

	Vector3D U0_11, U0_12, U0_21, U0_22;

	if (side == 0)
	{
		U0_11 = quadAtom11->getU0();
		U0_21 = quadAtom21->getU0();
		U0_12 = quadAtom12->getU0();
		U0_22 = quadAtom22->getU0();
	}
	else
	{
		U0_11 = quadAtom11->getU1();
		U0_21 = quadAtom21->getU1();
		U0_12 = quadAtom12->getU1();
		U0_22 = quadAtom22->getU1();
	}

	Vector3D du11, du21, du12, du22, dv11, dv21, dv12, dv22;

	Vector3D dU0du11, dU0dv11, dU0du21, dU0dv21, dU0du12, dU0dv12, dU0du22, dU0dv22;

	// Get derivatives for each corner, starting at 11, then 21, 12, and 22

	// Corner 11
	// du

	if (ubase == 0)
	{
		du11 = x21 - x11;
		dU0du11 = U0_21 - U0_11;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase-1, vbase));
		du11 = (x21 - tempAtom->getX()) / 2;

		if (side == 0)
			dU0du11 = (U0_21 - tempAtom->getU0()) / 2;
		else
			dU0du11 = (U0_21 - tempAtom->getU1()) / 2;
	}

	//dv

	if (vbase == 0)
	{
		dv11 = x12 - x11;
		dU0dv11 = U0_12 - U0_11;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase, vbase-1));
		dv11 = (x12 - tempAtom->getX()) / 2;

		if (side == 0)
			dU0dv11 = (U0_12 - tempAtom->getU0()) / 2;
		else
			dU0dv11 = (U0_12 - tempAtom->getU1()) / 2;
	}

	// Corner 21

	// du

	if (ubase + 1 == tempFigure->getRowCount() - 1)
	{
		du21 = x21 - x11;
		dU0du21 = U0_21 - U0_11;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase+2, vbase));
		du21 = (tempAtom->getX() - x11) / 2;

		if (side == 0)
			dU0du21 = (tempAtom->getU0() - U0_11) / 2;
		else
			dU0du21 = (tempAtom->getU1() - U0_11) / 2;
	}

	//dv

	if (vbase == 0)
	{
		dv21 = x22 - x21;
		dU0dv21 = U0_22 - U0_21;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase+1, vbase-1));
		dv21 = (x22 - tempAtom->getX()) / 2;

		if (side == 0)
			dU0dv21 = (U0_22 - tempAtom->getU0()) / 2;
		else
			dU0dv21 = (U0_22 - tempAtom->getU1()) / 2;
	}

	// Corner 12

	// du

	if (ubase == 0)
	{
		du12 = x22 - x12;
		dU0du12 = U0_22 - U0_12;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase-1,vbase+1));
		du12 = (x22 - tempAtom->getX()) / 2;

		if (side == 0)
			dU0du12 = (U0_22 - tempAtom->getU0()) / 2;
		else
			dU0du12 = (U0_22 - tempAtom->getU1()) / 2;
	}

	if (vbase + 1 == tempFigure->getColumnCount() - 1)
	{
		dv12 = x12 - x11;
		dU0dv12 = U0_12 - U0_11;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase, vbase+2));
		dv12 = (tempAtom->getX() - x11) / 2;

		if (side == 0)
			dU0dv12 = (tempAtom->getU0() - U0_11) / 2;
		else
			dU0dv12 = (tempAtom->getU1() - U0_11) / 2;
	}

	// Corner 22
	
	if (ubase + 1 == tempFigure->getRowCount() - 1)
	{
		du22 = x22 - x12;
		dU0du22 = U0_22 - U0_12;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase+2,vbase+1));
		du22 = (tempAtom->getX() - x12) / 2;

		if (side == 0)
			dU0du22 = (tempAtom->getU0() - U0_12) / 2;
		else
			dU0du22 = (tempAtom->getU1() - U0_12) / 2;
	}

	if (vbase + 1 == tempFigure->getColumnCount() - 1)
	{
		dv22 = x22 - x21;
		dU0dv22 = U0_22 - U0_21;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase+1,vbase+2));
		dv22 = (tempAtom->getX() - x21) / 2;

		if (side == 0)
			dU0dv22 = (tempAtom->getU0() - U0_21) / 2;
		else
			dU0dv22 = (tempAtom->getU1() - U0_21) / 2;
	}

//	cout << endl << "U derivatives: " << endl;
//	cout << du11.getX() << ", " << du11.getY() << ", " << du11.getZ() << endl;
//	cout << du21.getX() << ", " << du21.getY() << ", " << du21.getZ() << endl;
//	cout << du12.getX() << ", " << du12.getY() << ", " << du12.getZ() << endl;
//	cout << du22.getX() << ", " << du22.getY() << ", " << du22.getZ() << endl;

//	cout << endl << "V derivatives: " << endl;
//	cout << dv11.getX() << ", " << dv11.getY() << ", " << dv11.getZ() << endl;
//	cout << dv21.getX() << ", " << dv21.getY() << ", " << dv21.getZ() << endl;
//	cout << dv12.getX() << ", " << dv12.getY() << ", " << dv12.getZ() << endl;
//	cout << dv22.getX() << ", " << dv22.getY() << ", " << dv22.getZ() << endl;


	/*Vector3D du11 = x21 - x11;
	Vector3D du21 = x21 - x11;
	Vector3D du12 = x22 - x12;
	Vector3D du22 = x22 - x12;

	Vector3D dv11 = x12 - x11;
	Vector3D dv21 = x22 - x21;
	Vector3D dv12 = x12 - x11;
	Vector3D dv22 = x22 - x21;*/

	

	/*Vector3D du11 = (x21 - x01) / 2;
	Vector3D du21 = (x31 - x11) / 2;
	Vector3D du12 = (x22 - x02) / 2;
	Vector3D du22 = (x32 - x12) / 2;

	Vector3D dv11 = (x12 - x10) / 2;
	Vector3D dv21 = (x22 - x20) / 2;
	Vector3D dv12 = (x13 - x11) / 2;
	Vector3D dv22 = (x23 - x21) / 2;*/

	// Get unit normals for the four corners and project derivatives on to the tangent plane
	Vector3D n11 = quadAtom11->getU0() - quadAtom11->getU1();
	n11.normalize();
	Vector3D n21 = quadAtom21->getU0() - quadAtom21->getU1();
	n21.normalize();
	Vector3D n12 = quadAtom12->getU0() - quadAtom22->getU1();
	n12.normalize();
	Vector3D n22 = quadAtom22->getU0() - quadAtom22->getU1();
	n22.normalize();

	Vector3D du11t = du11 - (du11 * n11) * n11;
	Vector3D du21t = du21 - (du21 * n21) * n21;
	Vector3D du12t = du12 - (du12 * n12) * n12;
	Vector3D du22t = du22 - (du22 * n22) * n22;

	Vector3D dv11t = dv11 - (dv11 * n11) * n11;
	Vector3D dv21t = dv21 - (dv21 * n21) * n21;
	Vector3D dv12t = dv12 - (dv12 * n12) * n12;
	Vector3D dv22t = dv22 - (dv22 * n22) * n22;

	// Build matrices for hermite interpolation of medial sheet
	double hx[16] = { x11.getX(), x21.getX(), du11t.getX(), du21t.getX(), 
		x12.getX(), x22.getX(), du12t.getX(), du22t.getX(), dv11t.getX(), 
		dv21t.getX(), 0, 0, dv21t.getX(), dv22t.getX(), 0, 0};

	double hy[16] = { x11.getY(), x21.getY(), du11t.getY(), du21t.getY(), 
		x12.getY(), x22.getY(), du12t.getY(), du22t.getY(), dv11t.getY(), 
		dv21t.getY(), 0, 0, dv21t.getY(), dv22t.getY(), 0, 0};

	double hz[16] = { x11.getZ(), x21.getZ(), du11t.getZ(), du21t.getZ(), 
		x12.getZ(), x22.getZ(),du12t.getZ(), du22t.getZ(), dv11t.getZ(), 
		dv21t.getZ(), 0, 0, dv21t.getZ(), dv22t.getZ(), 0, 0};

	Matrix hxmat = Matrix(4, 4, hx, true);
	Matrix hymat = Matrix(4, 4, hy, true);
	Matrix hzmat = Matrix(4, 4, hz, true);

	double hu[4] = { h1(u), h2(u), h3(u), h4(u) };
	double hv[4] = { h1(v), h2(v), h3(v), h4(v) };
	Matrix humat = Matrix(1, 4, hu, true);
	Matrix hvmat = Matrix(4, 1, hv, true);

	Matrix xn = humat * hxmat * hvmat;
	Matrix yn = humat * hymat * hvmat;
	Matrix zn = humat * hzmat * hvmat;

	// Calculate sRad matrices using finite differences
	/*Vector3D dU0du11 = (quadAtom21->getU0() - quadAtom01->getU0()) / 2;
	Vector3D dU0dv11 = (quadAtom12->getU0() - quadAtom10->getU0()) / 2;

	Vector3D dU0du21 = (quadAtom31->getU0() - quadAtom11->getU0()) / 2;
	Vector3D dU0dv21 = (quadAtom22->getU0() - quadAtom20->getU0()) / 2;

	Vector3D dU0du12 = (quadAtom22->getU0() - quadAtom02->getU0()) / 2;
	Vector3D dU0dv12 = (quadAtom12->getU0() - quadAtom10->getU0()) / 2;

	Vector3D dU0du22 = (quadAtom32->getU0() - quadAtom12->getU0()) / 2;
	Vector3D dU0dv22 = (quadAtom23->getU0() - quadAtom21->getU0()) / 2;*/

	/*Vector3D dU0du11, dU0dv11, dU0du21, dU0dv21, dU0du12, dU0dv12, dU0du22, dU0dv22;

	if (side == 0)
	{
		dU0du11 = quadAtom21->getU0() - quadAtom11->getU0();
		dU0dv11 = quadAtom12->getU0() - quadAtom11->getU0();

		dU0du21 = quadAtom21->getU0() - quadAtom11->getU0();
		dU0dv21 = quadAtom22->getU0() - quadAtom21->getU0();

		dU0du12 = quadAtom22->getU0() - quadAtom12->getU0();
		dU0dv12 = quadAtom12->getU0() - quadAtom11->getU0();

		dU0du22 = quadAtom22->getU0() - quadAtom12->getU0();
		dU0dv22 = quadAtom22->getU0() - quadAtom21->getU0();
	}
	else
	{
		dU0du11 = quadAtom21->getU1() - quadAtom11->getU1();
		dU0dv11 = quadAtom12->getU1() - quadAtom11->getU1();

		dU0du21 = quadAtom21->getU1() - quadAtom11->getU1();
		dU0dv21 = quadAtom22->getU1() - quadAtom21->getU1();

		dU0du12 = quadAtom22->getU1() - quadAtom12->getU1();
		dU0dv12 = quadAtom12->getU1() - quadAtom11->getU1();

		dU0du22 = quadAtom22->getU1() - quadAtom12->getU1();
		dU0dv22 = quadAtom22->getU1() - quadAtom21->getU1();
	}*/

	// These form the non-orthogonal medial coordinate system for the spoke derivatives

	//Vector3D U0_11, U0_12, U0_21, U0_22;

	//if (side == 0)
	//{
	//	U0_11 = quadAtom11->getU0();
	//	U0_21 = quadAtom21->getU0();
	//	U0_12 = quadAtom12->getU0();
	//	U0_22 = quadAtom22->getU0();
	//}
	//else
	//{
	//	U0_11 = quadAtom11->getU1();
	//	U0_21 = quadAtom21->getU1();
	//	U0_12 = quadAtom12->getU1();
	//	U0_22 = quadAtom22->getU1();
	//}

	double basis11[9] = { du11.getX(), du11.getY(), du11.getZ(), dv11.getX(), dv11.getY(), dv11.getZ(), 
		U0_11.getX(), U0_11.getY(), U0_11.getZ() };

	double basis21[9] = { du21.getX(), du21.getY(), du21.getZ(), dv21.getX(), dv21.getY(), dv21.getZ(), 
		U0_21.getX(), U0_21.getY(), U0_21.getZ() };

	double basis12[9] = { du12.getX(), du12.getY(), du12.getZ(), dv12.getX(), dv12.getY(), dv12.getZ(), 
		U0_12.getX(), U0_12.getY(), U0_12.getZ() };

	double basis22[9] = { du22.getX(), du22.getY(), du22.getZ(), dv22.getX(), dv22.getY(), dv22.getZ(), 
		U0_22.getX(), U0_22.getY(), U0_22.getZ() };
	


	Matrix C11 = Matrix(3, 3, basis11, true);
	Matrix C21 = Matrix(3, 3, basis21, true);
	Matrix C12 = Matrix(3, 3, basis12, true);
	Matrix C22 = Matrix(3, 3, basis22, true);

	Matrix Ci11, Ci21, Ci12, Ci22;
	
	C11.inverse(Ci11);
	C21.inverse(Ci21);
	C12.inverse(Ci12);
	C22.inverse(Ci22);

	double dU0du11_coeffs[3] = { dU0du11.getX(), dU0du11.getY(), dU0du11.getZ() };
	double dU0dv11_coeffs[3] = { dU0dv11.getX(), dU0dv11.getY(), dU0dv11.getZ() };
	double dU0du21_coeffs[3] = { dU0du21.getX(), dU0du21.getY(), dU0du21.getZ() };
	double dU0dv21_coeffs[3] = { dU0dv21.getX(), dU0dv21.getY(), dU0dv21.getZ() };
	double dU0du12_coeffs[3] = { dU0du12.getX(), dU0du12.getY(), dU0du12.getZ() };
	double dU0dv12_coeffs[3] = { dU0dv12.getX(), dU0dv12.getY(), dU0dv12.getZ() };
	double dU0du22_coeffs[3] = { dU0du22.getX(), dU0du22.getY(), dU0du22.getZ() };
	double dU0dv22_coeffs[3] = { dU0dv22.getX(), dU0dv22.getY(), dU0dv22.getZ() };

	Matrix Au11 = Matrix(3,1,dU0du11_coeffs,true);
	Matrix Av11 = Matrix(3,1,dU0dv11_coeffs,true);
	Matrix Au21 = Matrix(3,1,dU0dv21_coeffs,true);
	Matrix Av21 = Matrix(3,1,dU0dv21_coeffs,true);
	Matrix Au12 = Matrix(3,1,dU0du12_coeffs,true);
	Matrix Av12 = Matrix(3,1,dU0dv12_coeffs,true);
	Matrix Au22 = Matrix(3,1,dU0du22_coeffs,true);
	Matrix Av22 = Matrix(3,1,dU0dv22_coeffs,true);

	Matrix Bu11 = Ci11 * Au11;
	Matrix Bv11 = Ci11 * Av11;
	Matrix Bu21 = Ci21 * Au21;
	Matrix Bv21 = Ci21 * Av21;
	Matrix Bu12 = Ci12 * Au12;
	Matrix Bv12 = Ci12 * Av12;
	Matrix Bu22 = Ci22 * Au22;
	Matrix Bv22 = Ci22 * Av22;

	double b11[4] = { -1*Bu11(0,0), -1*Bu11(1,0), -1*Bv11(0,0), -1*Bv11(1,0) };
	double b21[4] = { -1*Bu21(0,0), -1*Bu21(1,0), -1*Bv21(0,0), -1*Bv21(1,0) };
	double b12[4] = { -1*Bu12(0,0), -1*Bu12(1,0), -1*Bv12(0,0), -1*Bv12(1,0) };
	double b22[4] = { -1*Bu22(0,0), -1*Bu22(1,0), -1*Bv22(0,0), -1*Bv22(1,0) };

	Matrix Srad11 = Matrix(2,2,b11,true);
	Matrix Srad21 = Matrix(2,2,b21,true);
	Matrix Srad12 = Matrix(2,2,b12,true);
	Matrix Srad22 = Matrix(2,2,b22,true);

	Matrix rSrad11 = Srad11 * quadAtom11->getR();
	Matrix rSrad21 = Srad21 * quadAtom21->getR();
	Matrix rSrad12 = Srad12 * quadAtom12->getR();
	Matrix rSrad22 = Srad22 * quadAtom22->getR();

//	cout << endl << "rSrad: " << endl;
//	cout << rSrad11(0,0) << ", " << rSrad11(0,1) << ", " << rSrad11(1,0) << ", " << rSrad11(1,1) << endl;
//	cout << rSrad21(0,0) << ", " << rSrad21(0,1) << ", " << rSrad21(1,0) << ", " << rSrad21(1,1) << endl;
//	cout << rSrad12(0,0) << ", " << rSrad12(0,1) << ", " << rSrad12(1,0) << ", " << rSrad12(1,1) << endl;
//	cout << rSrad22(0,0) << ", " << rSrad22(0,1) << ", " << rSrad22(1,0) << ", " << rSrad22(1,1) << endl;

	Vector L11, L21, L12, L22;
	Matrix V11, V21, V12, V22;

	rSrad11.factorEV(L11, V11, NON_SYM);
	rSrad21.factorEV(L21, V21, NON_SYM);
	rSrad12.factorEV(L12, V12, NON_SYM);
	rSrad22.factorEV(L22, V22, NON_SYM);

	//cout << "HERE" << endl;

	double Lam1_11, Lam2_11, Lam1_21, Lam2_21, Lam1_12, Lam2_12, Lam1_22, Lam2_22;
	Vector e1_11, e2_11, e1_21, e2_21, e1_12, e2_12, e1_22, e2_22;

	// Corner 1
	if (L11(0) > L11(1))
	{
		Lam1_11 = L11(0);
		Lam2_11 = L11(1);

		e1_11 = V11.getColumn(0);
		e2_11 = V11.getColumn(1);
	}
	else
	{
		Lam1_11 = L11(1);
		Lam2_11 = L11(0);

		e1_11 = V11.getColumn(1);
		e2_11 = V11.getColumn(0);
	}

	// Corner 2
	if (L21(0) > L21(1))
	{
		Lam1_21 = L21(0);
		Lam2_21 = L21(1);

		e1_21 = V21.getColumn(0);
		e2_21 = V21.getColumn(1);
	}
	else
	{
		Lam1_21 = L21(1);
		Lam2_21 = L21(0);

		e1_21 = V21.getColumn(1);
		e2_21 = V21.getColumn(0);
	}

	//Corner 3
	if (L12(0) > L12(1))
	{
		Lam1_12 = L12(0);
		Lam2_12 = L12(1);

		e1_12 = V12.getColumn(0);
		e2_12 = V12.getColumn(1);
	}
	else
	{
		Lam1_12 = L12(1);
		Lam2_12 = L12(0);

		e1_12 = V12.getColumn(1);
		e2_12 = V12.getColumn(0);
	}

	//Corner 4
	if (L22(0) > L22(1))
	{
		Lam1_22 = L22(0);
		Lam2_22 = L22(1);

		e1_22 = V22.getColumn(0);
		e2_22 = V22.getColumn(1);
	}
	else
	{
		Lam1_22 = L22(1);
		Lam2_22 = L22(0);

		e1_22 = V22.getColumn(1);
		e2_22 = V22.getColumn(0);
	}

//	cout << endl << "Eigenvectors: " << endl;
//	cout << "(" << e1_11(0) << ", " << e1_11(1) << "), (" << e2_11(0) << ", " << e2_11(1) << ")" << endl;
//	cout << "(" << e1_21(0) << ", " << e1_21(1) << "), (" << e2_21(0) << ", " << e2_21(1) << ")" << endl;
//	cout << "(" << e1_12(0) << ", " << e1_12(1) << "), (" << e2_12(0) << ", " << e2_12(1) << ")" << endl;
//	cout << "(" << e1_22(0) << ", " << e1_22(1) << "), (" << e2_22(0) << ", " << e2_22(1) << ")" << endl;

	double testlambda[4] = { Lam1_11, 0, 0, Lam2_11 };

	
	
	Matrix testL = Matrix(2,2,testlambda,true);
	Matrix testinv;

	Matrix testsrad = V11 * testL * V11.inverse(testinv);

//	cout << endl << endl;
//	cout << testsrad(0,0) << ", " << testsrad(0,1) << ", " << testsrad(1,0) << ", " << testsrad(1,1) << endl;
//	cout << endl << endl;

//	cout << endl << "Eigenvalues: " << endl;
//	cout << Lam1_11 << ", " << Lam2_11 << endl;
//	cout << Lam1_21 << ", " << Lam2_21 << endl;
//	cout << Lam1_12 << ", " << Lam2_12 << endl;
//	cout << Lam1_22 << ", " << Lam2_22 << endl;


	double logLam1_11 = log(1 - Lam1_11);
	double logLam2_11 = log(1 - Lam2_11);
	double logLam1_21 = log(1 - Lam1_21);
	double logLam2_21 = log(1 - Lam2_21);
	double logLam1_12 = log(1 - Lam1_12);
	double logLam2_12 = log(1 - Lam2_12);
	double logLam1_22 = log(1 - Lam1_22);
	double logLam2_22 = log(1 - Lam2_22);

	double logAvg1 = (1-u)*(1-v)*logLam1_11 + (u)*(1-v)*logLam1_21 + (1-u)*(v)*logLam1_12 + (u)*(v)*logLam1_22;
	double logAvg2 = (1-u)*(1-v)*logLam2_11 + (u)*(1-v)*logLam2_21 + (1-u)*(v)*logLam2_12 + (u)*(v)*logLam2_22;

	double Lam1 = 1 - exp(logAvg1);
	double Lam2 = 1 - exp(logAvg2);

	// Uses a form of bilinear interpolation to get the eigenvectors: This should probably be a Frechet mean of the thetas 
	// on the unit circle in the future

	/*double theta1_11 = atan2(e1_11(1), e1_11(0));
	double theta2_11 = atan2(e2_11(1), e2_11(0));
	double theta1_21 = atan2(e1_21(1), e1_21(0));
	double theta2_21 = atan2(e2_21(1), e2_21(0));
	double theta1_12 = atan2(e1_12(1), e1_12(0));
	double theta2_12 = atan2(e2_12(1), e2_12(0));
	double theta1_22 = atan2(e1_22(1), e1_22(0));
	double theta2_22 = atan2(e2_22(1), e2_22(0));*/

	double PI = 3.14159265;

	double avgx1, avgy1, avgx2, avgy2, avgtheta1, avgtheta2;//, at1_11_21, at2_11_21, at1_12_22, at2_12_22;

	avgx1 = (1-u)*(1-v)*e1_11(0) + (u)*(1-v)*e1_21(0) + (1-u)*(v)*e1_12(0) + (u)*(v)*e1_22(0);
	avgy1 = (1-u)*(1-v)*e1_11(1) + (u)*(1-v)*e1_21(1) + (1-u)*(v)*e1_12(1) + (u)*(v)*e1_22(1);
	avgtheta1 = atan2(avgy1, avgx1);

	avgx2 = (1-u)*(1-v)*e2_11(0) + (u)*(1-v)*e2_21(0) + (1-u)*(v)*e2_12(0) + (u)*(v)*e2_22(0);
	avgy2 = (1-u)*(1-v)*e2_11(1) + (u)*(1-v)*e2_21(1) + (1-u)*(v)*e2_12(1) + (u)*(v)*e2_22(1);
	avgtheta2 = atan2(avgy2, avgx2);

	// Go from 11 to 21 (u direction)
	//double thetadist1_11_21 = theta1_21 - theta1_11;
	//if (thetadist1_11_21 > PI)
	//{
	//	theta1_21 -= 2*PI;
	//	thetadist1_11_21 = theta1_21 - theta1_11;
	//}
	//else if (thetadist1_11_21 < -PI)
	//{
	//	theta1_21 += 2*PI;
	//	thetadist1_11_21 = theta1_21 - theta1_11;
	//}
	//at1_11_21 = theta1_11 + u*thetadist1_11_21;

	//double thetadist2_11_21 = theta2_21 - theta2_11;
	//if (thetadist2_11_21 > PI)
	//{
	//	theta2_21 -= 2*PI;
	//	thetadist2_11_21 = theta2_21 - theta2_11;
	//}
	//else if (thetadist2_11_21 < -PI)
	//{
	//	theta2_21 += 2*PI;
	//	thetadist2_11_21 = theta2_21 - theta2_11;
	//}
	//at2_11_21 = theta2_11 + u*thetadist2_11_21;

	//// Go from 12 to 22 (u direction)
	//double thetadist1_12_22 = theta1_22 - theta1_12;
	//if (thetadist1_12_22 > PI)
	//{
	//	theta1_22 -= 2*PI;
	//	thetadist1_12_22 = theta1_22 - theta1_12;
	//}
	//else if (thetadist1_12_22 < -PI)
	//{
	//	theta1_22 += 2*PI;
	//	thetadist1_12_22 = theta1_22 - theta1_12;
	//}
	//at1_12_22 = theta1_12 + u*thetadist1_12_22;

	//double thetadist2_12_22 = theta2_22 - theta2_12;
	//if (thetadist2_12_22 > PI)
	//{
	//	theta2_22 -= 2*PI;
	//	thetadist2_12_22 = theta2_22 - theta2_12;
	//}
	//else if (thetadist2_12_22 < -PI)
	//{
	//	theta2_22 += 2*PI;
	//	thetadist2_12_22 = theta2_22 - theta2_12;
	//}
	//at2_12_22 = theta2_12 + u*thetadist2_12_22;

	//// Now go v direction
	//double thetadist1_v = at1_12_22 - at1_11_21;
	//if (thetadist1_v > PI)
	//{
	//	at1_12_22 -= 2*PI;
	//	thetadist1_v = at1_12_22 - at1_11_21;
	//}
	//else if (thetadist1_v < -PI)
	//{
	//	at1_12_22 += 2*PI;
	//	thetadist1_v = at1_12_22 - at1_11_21;
	//}
	//avgtheta1 = at1_11_21 + v*thetadist1_v;

	//double thetadist2_v = at2_12_22 - at2_11_21;
	//if (thetadist2_v > PI)
	//{
	//	at2_12_22 -= 2*PI;
	//	thetadist2_v = at2_12_22 - at2_11_21;
	//}
	//else if (thetadist2_v < -PI)
	//{
	//	at2_12_22 += 2*PI;
	//	thetadist2_v = at2_12_22 - at2_11_21;
	//}
	//avgtheta2 = at2_11_21 + v*thetadist2_v;

	double neweigenv[4] = { cos(avgtheta1), sin(avgtheta1), cos(avgtheta2), sin(avgtheta2) };

	Matrix NewV = Matrix(2,2,neweigenv,true);

	double newlambda[4] = { Lam1, 0, 0, Lam2 };

	Matrix NewL = Matrix(2,2,newlambda,true);
	
	/*cout << V11(0,0) << ", " << V11(1,0)  << ", " << V11(0,1) << ", " << V11(1,1) << endl;
	cout << V21(0,0) << ", " << V21(1,0)  << ", " << V21(0,1) << ", " << V21(1,1) << endl;
	cout << V12(0,0) << ", " << V12(1,0)  << ", " << V12(0,1) << ", " << V12(1,1) << endl;
	cout << V22(0,0) << ", " << V22(1,0)  << ", " << V22(0,1) << ", " << V22(1,1) << endl;
	cout << NewV(0,0) << ", " << NewV(1,0) << ", " << NewV(0,1) << ", " << NewV(1,1) << endl;*/
	
	Matrix NewVi;
	NewV.inverse(NewVi);

	Matrix NewrSrad = NewV * NewL * NewVi;

	double curru = 0;
	double currv = 0;
	double du = u / 10;
	double dv = v / 10;

	if (side == 0)
		U0_11 = quadAtom11->getU0();
	else
		U0_11 = quadAtom11->getU1();
	
	double Pdata[6] = { du11.getX(), dv11.getX(), du11.getY(), dv11.getY(), du11.getZ(), dv11.getZ() };
	Matrix P = Matrix(2,3,Pdata,true);

	double Udata[3] = { U0_11.getX(), U0_11.getY(), U0_11.getZ() };
	Matrix U = Matrix(1,3,Udata,true);

	double ru11 = du11 * U0_11;
	double rv11 = dv11 * U0_11;

	double ru21 = du21 * U0_21;
	double rv21 = dv21 * U0_21;

	double ru12 = du12 * U0_12;
	double rv12 = dv12 * U0_12;

	double ru22 = du22 * U0_22;
	double rv22 = dv22 * U0_22;

//	cout << endl << "Ru & Rv: " << endl;
//	cout << ru11 << ", " << rv11 << endl;
//	cout << ru21 << ", " << rv21 << endl;
//	cout << ru12 << ", " << rv12 << endl;
//	cout << ru22 << ", " << rv22 << endl;

	double Idata[9] = {1,0,0,0,1,0,0,0,1};
	Matrix I = Matrix(3,3,Idata,true);

	Matrix Q = P * ( (U.t() * U) - I );

	Matrix R = -1 * P * U.t();

	Matrix dS = ( rSrad11.t() * Q ) + ( R * U );

	Vector dSdu = dS.getRow(0);
	Vector dSdv = dS.getRow(1);

	Vector3D dSu = Vector3D(dSdu(0), dSdu(1), dSdu(2));
	Vector3D dSv = Vector3D(dSdv(0), dSdv(1), dSdv(2));

	Vector3D S0;

	if (side == 0)
		S0 = U0_11 * quadAtom11->getR0();
	else 
		S0 = U0_11 * quadAtom11->getR1();

	Vector3D NewSpoke = S0 + (du * dSu) + (dv * dSv);

	curru = curru + du;
	currv = currv + dv;

	hu[0] = h1(curru); hu[1] = h2(curru); hu[2] = h3(curru); hu[3] = h4(curru); 
	hv[0] = h1(currv); hv[1] = h2(currv); hv[2] = h3(currv); hv[3] = h4(currv);	

	humat = Matrix(1,4,hu,true);
	hvmat = Matrix(4,1,hv,true);

	xn = humat * hxmat * hvmat;
	yn = humat * hymat * hvmat;
	zn = humat * hzmat * hvmat;

	Vector3D newpos = Vector3D(xn(0,0), yn(0,0), zn(0,0));

	for (int i = 0; i < 10; i++)
	{
		// To find rSrad, need to calculate lambdas and eigenvectors
		logAvg1 = (1-curru)*(1-currv)*logLam1_11 + (curru)*(1-currv)*logLam1_21 + (1-curru)*(currv)*logLam1_12 + (curru)*(currv)*logLam1_22;
		logAvg2 = (1-curru)*(1-currv)*logLam2_11 + (curru)*(1-currv)*logLam2_21 + (1-curru)*(currv)*logLam2_12 + (curru)*(currv)*logLam2_22;
	
		Lam1 = 1 - exp(logAvg1);
		Lam2 = 1 - exp(logAvg2);


		avgx1 = (1-curru)*(1-currv)*e1_11(0) + (curru)*(1-currv)*e1_21(0) + (1-curru)*(currv)*e1_12(0) + (curru)*(currv)*e1_22(0);
		avgy1 = (1-curru)*(1-currv)*e1_11(1) + (curru)*(1-currv)*e1_21(1) + (1-curru)*(currv)*e1_12(1) + (curru)*(currv)*e1_22(1);
		avgtheta1 = atan2(avgy1, avgx1);

		avgx2 = (1-u)*(1-v)*e2_11(0) + (u)*(1-v)*e2_21(0) + (1-u)*(v)*e2_12(0) + (u)*(v)*e2_22(0);
		avgy2 = (1-u)*(1-v)*e2_11(1) + (u)*(1-v)*e2_21(1) + (1-u)*(v)*e2_12(1) + (u)*(v)*e2_22(1);
		avgtheta2 = atan2(avgy2, avgx2);

		neweigenv[0] = cos(avgtheta1); neweigenv[1] = sin(avgtheta1); neweigenv[2] = cos(avgtheta2); neweigenv[3] = sin(avgtheta2);
		NewV = Matrix(2,2,neweigenv,true);
		NewV.inverse(NewVi);

		newlambda[0] = Lam1; newlambda[1] = 0; newlambda[2] = 0; newlambda[3] = Lam2;
		NewL = Matrix(2,2,newlambda,true);

		NewrSrad = NewV * NewL * NewVi;

		double hu[4] = { h1(curru), h2(curru), h3(curru), h4(curru) };
		double hv[4] = { h1(currv), h2(currv), h3(currv), h4(currv) };
		Matrix humat = Matrix(1,4,hu,true);
		Matrix hvmat = Matrix(4,1,hv,true);

		double hup[4] = { h1p(curru), h2p(curru), h3p(curru), h4p(curru) };
		double hvp[4] = { h1p(currv), h2p(currv), h3p(currv), h4p(currv) };
		Matrix hupmat = Matrix(1,4,hup,true);
		Matrix hvpmat = Matrix(4,1,hvp,true);

		Matrix pux = hupmat * hxmat * hvmat;
		Matrix puy = hupmat * hymat * hvmat;
		Matrix puz = hupmat * hzmat * hvmat;

		Matrix pvx = humat * hxmat * hvpmat;
		Matrix pvy = humat * hymat * hvpmat;
		Matrix pvz = humat * hzmat * hvpmat;
	
		double pdata[6] = { pux(0,0), pvx(0,0), puy(0,0), pvy(0,0), puz(0,0), pvz(0,0) };
		P = Matrix(2,3,pdata,true);

		double r = NewSpoke.normalize();

		double udata[3] = { NewSpoke.getX(), NewSpoke.getY(), NewSpoke.getZ() };
		U = Matrix(1,3,udata,true);

		Q = P * ( (U.t() * U) - I );
		
		R = -1 * P * U.t();
		
		dS = ( rSrad11.t() * Q ) + ( R * U );

		dSdu = dS.getRow(0);
		dSdv = dS.getRow(1);

		dSu = Vector3D(dSdu(0), dSdu(1), dSdu(2));
		dSv = Vector3D(dSdv(0), dSdv(1), dSdv(2));

		Vector3D TempSpoke = r * NewSpoke;

		NewSpoke = TempSpoke + (du * dSu) + (dv * dSv);

		curru = curru + du;
		currv = currv + dv;

	}

	hu[0] = h1(u); hu[1] = h2(u); hu[2] = h3(u); hu[3] = h4(u); 
	hv[0] = h1(v); hv[1] = h2(v); hv[2] = h3(v); hv[3] = h4(v);	

	humat = Matrix(1,4,hu,true);
	hvmat = Matrix(4,1,hv,true);

	xn = humat * hxmat * hvmat;
	yn = humat * hymat * hvmat;
	zn = humat * hzmat * hvmat;

	newpos = Vector3D(xn(0,0), yn(0,0), zn(0,0));

	//cout << "{" << quadAtom11->getX().getX() + (quadAtom11->getR0() * quadAtom11->getU0().getX()) << ", " << quadAtom11->getX().getY() + (quadAtom11->getR0() * quadAtom11->getU0().getY()) << ", " << quadAtom11->getX().getZ() + (quadAtom11->getR0() * quadAtom11->getU0().getZ()) << "}" << endl;
	//cout << "{" << quadAtom21->getX().getX() + (quadAtom21->getR0() * quadAtom21->getU0().getX()) << ", " << quadAtom21->getX().getY() + (quadAtom21->getR0() * quadAtom21->getU0().getY()) << ", " << quadAtom21->getX().getZ() + (quadAtom21->getR0() * quadAtom21->getU0().getZ()) << "}" << endl;
	//cout << "{" << quadAtom12->getX().getX() + (quadAtom12->getR0() * quadAtom12->getU0().getX()) << ", " << quadAtom12->getX().getY() + (quadAtom12->getR0() * quadAtom12->getU0().getY()) << ", " << quadAtom12->getX().getZ() + (quadAtom12->getR0() * quadAtom12->getU0().getZ()) << "}" << endl;
	//cout << "{" << quadAtom22->getX().getX() + (quadAtom22->getR0() * quadAtom22->getU0().getX()) << ", " << quadAtom22->getX().getY() + (quadAtom22->getR0() * quadAtom22->getU0().getY()) << ", " << quadAtom22->getX().getZ() + (quadAtom22->getR0() * quadAtom22->getU0().getZ()) << "}" << endl;
	//cout << "{" << newpos.getX() + NewSpoke.getX() << ", " << newpos.getY() + NewSpoke.getY() << ", " << newpos.getZ() + NewSpoke.getZ() << endl;

	//cout << quadAtom11->getR0() << endl;
	//cout << quadAtom21->getR0() << endl;
	//cout << quadAtom12->getR0() << endl;
	//cout << quadAtom22->getR0() << endl;
	//cout << NewSpoke.norm() << endl;

	//cout << quadAtom11->getX().getX() << ", " << quadAtom11->getX().getY() << ", " << quadAtom11->getX().getZ() << endl;
	//cout << quadAtom21->getX().getX() << ", " << quadAtom21->getX().getY() << ", " << quadAtom21->getX().getZ() << endl;
	//cout << quadAtom12->getX().getX() << ", " << quadAtom12->getX().getY() << ", " << quadAtom12->getX().getZ() << endl;
	//cout << quadAtom22->getX().getX() << ", " << quadAtom22->getX().getY() << ", " << quadAtom22->getX().getZ() << endl;

	double NewR = NewSpoke.normalize();

	//cout << newpos.getX() << ", " << newpos.getY() << ", " << newpos.getZ() << endl;

	//cout << newpos.getX() << ", " << NewR << ", " << NewSpoke.getX() << endl;

	M3DSpoke* spoke = new M3DSpoke(newpos, NewSpoke, NewR);
	

	return spoke;

}
