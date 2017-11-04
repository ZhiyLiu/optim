/* Get tps transformation from source and target PDM (.vtk points), apply this transform onto target srep,
 * we get the source s-rep's spoke info, re-generate this s-rep (target s-rep).
 *
 * Liyun Tu
 * May 16, 2014
*/



#include "thinplatesplinepdmtosrep.h"

thinplatesplinepdmtosrep::thinplatesplinepdmtosrep()
{
}


/* Get the tps transformation from the two input .vtk file (PDM).
 * Apply tps transformation to source s-rep, generate a s-rep for target image.
*/
thinplatesplinepdmtosrep::TransformType::Pointer thinplatesplinepdmtosrep::get_tps(const char * sourcefilename, const char * targetfilename){

    // Landmarks correspondances may be associated with the SplineKernelTransforms
    // via Point Set containers. Let us define containers for the landmarks.
    // Define container for landmarks
    PointSetType::Pointer sourceLandMarks = PointSetType::New();
    PointSetType::Pointer targetLandMarks = PointSetType::New();
    PointType p1; PointType p2; // same as double p1[3];
    PointSetType::PointsContainer::Pointer sourceLandMarkContainer = sourceLandMarks->GetPoints();
    PointSetType::PointsContainer::Pointer targetLandMarkContainer = targetLandMarks->GetPoints();

    PointIdType id_s = itk::NumericTraits< PointIdType >::Zero;
    PointIdType id_t = itk::NumericTraits< PointIdType >::Zero;

    // Read in the source points set
    vtkSmartPointer<vtkPolyData> polyData_source; // same as vtkPolyData* polyData_source;
    vtkSmartPointer<vtkPolyData> polyData_target;

    vtkSmartPointer<vtkPolyDataReader> reader_source = vtkSmartPointer<vtkPolyDataReader>::New();
    vtkSmartPointer<vtkPolyDataReader> reader_target = vtkSmartPointer<vtkPolyDataReader>::New();

    reader_source->SetFileName(sourcefilename);
    reader_source->Update();
    polyData_source = reader_source->GetOutput();
    cout<<"----sourcefilename: "<<sourcefilename<<endl;

    reader_target->SetFileName(targetfilename);
    reader_target->Update();
    polyData_target = reader_target->GetOutput();
    cout<<"----targetfilename: "<<targetfilename<<endl;

    // Read in the source points set
    for(unsigned int i = 0; i < polyData_source->GetNumberOfPoints(); i++){
        double p[3];
        polyData_source->GetPoint(i,p);
        // This is identical to:
        // polydata->GetPoints()->GetPoint(i,p);
        p1[0] = p[0];
        p1[1] = p[1];
        p1[2] = p[2];
        //std::cout << "Point " << id_s << " : (" << p1[0] << " " << p1[1] << " " << p1[2] << ")" << std::endl;
        sourceLandMarkContainer->InsertElement(id_s, p1);
        id_s++;
    }

    // Read in the target points set
    for(unsigned int i = 0; i < polyData_target->GetNumberOfPoints(); i++){
        double p[3];
        polyData_target->GetPoint(i,p);
        p2[0] = p[0];
        p2[1] = p[1];
        p2[2] = p[2];
        //std::cout << "Point " << id_t << " : (" << p2[0] << " " << p2[1] << " " << p2[2] << ")" << std::endl;
        targetLandMarkContainer->InsertElement(id_t, p2);
        id_t++;
    }

    TransformType::Pointer tps = TransformType::New();
    tps->SetSourceLandmarks(sourceLandMarks);
    tps->SetTargetLandmarks(targetLandMarks);

    cout<<"Computing W Matrix... "<<endl;
    tps->ComputeWMatrix();
    cout<<"Compute W Matrix finished!"<<endl;

    // *******************verifying transform*******************
/*    // Read in the source points set
    vtkPoints * sourcePoints = polyData_source->GetPoints();
    for(unsigned int i = 0; i < polyData_source->GetNumberOfPoints(); i++){
        double p[3];
        polyData_source->GetPoint(i,p);
        PointType pos;
        for (int i = 0; i < 3; i++) {
            pos[i] = p[i];
        }
        pos = tps->TransformPoint( pos );
        sourcePoints->SetPoint(i,pos[0],pos[1],pos[2]);
    }
    polyData_source->SetPoints(sourcePoints);
    vtkSmartPointer<vtkPolyDataWriter> polyDataWriter_Source = vtkSmartPointer<vtkPolyDataWriter>::New();
    polyDataWriter_Source->SetInput(polyData_source);
    std::string outputFile = std::string(sourcefilename) + "_transform.vtk"; //This should similar to target image (test shows exactly the same!?)

    polyDataWriter_Source->SetFileName(outputFile.c_str());
    polyDataWriter_Source->Update();
*/
    return tps;
}



thinplatesplinepdmtosrep::PointType thinplatesplinepdmtosrep::Vector3DtoPointType(Vector3D point){
    PointType pos;
    pos[0] = point.getX();
    pos[1] = point.getY();
    pos[2] = point.getZ();

    return pos;
}


double thinplatesplinepdmtosrep::calculateSpokeLength(PointType tail, PointType tip){
    PointType spokeVector;
    for(unsigned int dim = 0; dim < 3; dim++){
        spokeVector[dim] = tip[dim] - tail[dim];
    }

    // Compute the spoke length
    double spokeRadiu=0;
    for(unsigned int dim=0; dim<3; dim++){
        spokeRadiu += spokeVector[dim] * spokeVector[dim];
    }

    return sqrt(spokeRadiu);
}

Vector3D thinplatesplinepdmtosrep::calculateSpokeDirection(PointType tail, PointType tip){
    Vector3D spokeVector;

    // Compute the spoke length
    double spokeRadiu = calculateSpokeLength(tail, tip);

    // Normalize the spoke's direction to a unit vector
    spokeVector.setX((tip[0] - tail[0]) / spokeRadiu);
    spokeVector.setY((tip[1] - tail[1]) / spokeRadiu);
    spokeVector.setZ((tip[2] - tail[2]) / spokeRadiu);

    return spokeVector;
}


/* sourcefilename: template surface vtk points
 * targetfilename: the new surface we are generating the s-rep for.
 * sourcesrepfilename: the template s-rep (m3d)./tps_srep <tamplate PDM, target PDM, tamplate s-rep and target s-rep folders>
*/
int thinplatesplinepdmtosrep::tps_to_srep(const char * sourcefilename, const char * targetfilename, const char* sourcesrepfilename,
                                          const char* targetsrepFolder){

    // Get tps transform info
    TransformType::Pointer tps = get_tps(sourcefilename, targetfilename);

    // Read in source s-rep.
    M3DQuadFigure* quadFig = tls.GetQuadFigure(sourcesrepfilename);
    int rowNums = quadFig->getRowCount();
    int colNums = quadFig->getColumnCount();

    // Loop each atom
    for(unsigned int u = 0; u < rowNums; u++){
        for(unsigned int v = 0; v < colNums; v++){
            // For up and down spokes
            M3DQuadPrimitive* prim = dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v));
            // Skeletal point
            Vector3D tail = prim->getX();
            // Up boundary point
            Vector3D up_tip = prim->getX() + prim->getR0()*prim->getU0();
            // Down boundary point
            Vector3D down_tip = prim->getX() + prim->getR1()*prim->getU1();

            // Performance tps on these three points
            PointType transformed_tail = tps->TransformPoint(Vector3DtoPointType(tail));
            PointType transformed_up_tip = tps->TransformPoint(Vector3DtoPointType(up_tip));
            PointType transformed_down_tip = tps->TransformPoint(Vector3DtoPointType(down_tip));

            // Get transformed spoke radius and direction
            double upSpokeRadiu = calculateSpokeLength(transformed_tail, transformed_up_tip);
            double downSpokeRadiu = calculateSpokeLength(transformed_tail, transformed_down_tip);

            Vector3D upSpokeDir = calculateSpokeDirection(transformed_tail, transformed_up_tip);
            Vector3D downSpokeDir = calculateSpokeDirection(transformed_tail, transformed_down_tip);

            // Set transformed spoke radius and direction to this atom
            prim->setR0(upSpokeRadiu);
            prim->setR1(downSpokeRadiu);
            prim->setU0(upSpokeDir);
            prim->setU1(downSpokeDir);

            // For crest spokes
            if(u==0 || u==rowNums-1 || v==0 || v==colNums-1) {
                M3DQuadEndPrimitive* endPrim = dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(u, v));
                Vector3D crest_tip = endPrim->getX() + endPrim->getREnd()*endPrim->getUEnd();

                PointType transform_crest_tip = tps->TransformPoint(Vector3DtoPointType(crest_tip));

                // Get transformed spoke radius and direction
                double crestSpokeRadiu = calculateSpokeLength(transformed_tail, transform_crest_tip);
                Vector3D crestSpokeDir = calculateSpokeDirection(transformed_tail, transform_crest_tip);

                // Set transformed spoke radius and direction to this atom
                endPrim->setREnd(crestSpokeRadiu);
                endPrim->setUEnd(crestSpokeDir);
            }

            prim->setX(transformed_tail[0], transformed_tail[1], transformed_tail[2]);
        }
    }

    // Save this quadFig to m3d file
    saveModel(quadFig, targetfilename, sourcesrepfilename, targetsrepFolder);

    return EXIT_SUCCESS;
}



void thinplatesplinepdmtosrep::saveModel(M3DQuadFigure* quadFig, const char * targetfilename, const char* sourcesrepfilename,
                                         const char* targetsrepFolder){
    Registry registry;
    registry.readFromFile(sourcesrepfilename, false);

    // Write the changed quadFig to m3d file.
    char figureStr[1024];
    //sprintf(figureStr, "model.figure[%d]", figNum);
    sprintf(figureStr, "model.figure[%d]", 0);//In m3d file, every figure is identified as figure[0], cannot change the index, why?
    //cout<<"Message from movespokes.cpp; model.figure[%d] is: "<<figureStr<<endl;
    quadFig->writeFigure(figureStr,registry); //write the figure into registry.

    string filename = string(targetsrepFolder) + tls.getBaseFileName(targetfilename) + ".m3d";

    cout<<"----Write srep to: "<< filename << endl;
    //overwrite the input file to save the changes for next loop.
    registry.writeToFile(filename.c_str());

    cout<<"----Finished!" << endl;
}



/* sourcefilename: template surface vtk points
 * targetPDM: the file pathname of the new surface we are generating the s-rep for.
 * sourcesrepfilename: the template s-rep (m3d)./tps_srep <tamplate PDM, target PDM, tamplate s-rep and target s-rep folders>
*/
int thinplatesplinepdmtosrep::multiple_template_tps(const char * tamplatePDMFileList, const char * targetPDM, const char* tamplateSrepFileList,
                                          const char* targetSrepFolder){
    // get template PDMs
    vector< std::string > templatePDMs = readTempalte(tamplatePDMFileList);

    // get template s-reps
    vector< std::string > templateSreps = readTempalte(tamplateSrepFileList);

    int templateNum = templatePDMs.size();
    if(templateNum != templateSreps.size()){
        std::cerr << "Msg from thinplatesplinepdmtosrep::multiple_template_tps: something wrong with the input templates!!!" << std::endl;
        return EXIT_FAILURE;
    }

    // loop each target PDM, get the tps transform from each template
    vector<TransformType::Pointer > tpsList;
    for(unsigned int i = 0; i < templateNum; i++){
        // Get tps transform from each template to this targetPDM
        TransformType::Pointer tps = get_tps(templatePDMs[i].c_str(), targetPDM);

        tpsList.push_back(tps);
    }

    // read template sreps into a quad list, avoid repeat use "GetQuadFigure" which has memory leak
    vector<M3DQuadFigure* > quadFigList;
    for(unsigned int i = 0; i < templateNum; i++){
        // Read in source s-rep.
        M3DQuadFigure* quadFig = tls.GetQuadFigure(templateSreps[i].c_str());

        quadFigList.push_back(quadFig);
    }

    // apply each tps to its tamplate s-rep, and compute the mean of the transformed s-reps
    M3DQuadFigure* quadFig = computeMeanSrep(tpsList, true, quadFigList);

    // save the mean s-rep to m3d file
    string targetfilename = string(targetSrepFolder) + "\/" + tls.getBaseFileName(targetPDM) + ".m3d";
    string anysrepfilename = templateSreps[0]; //any srep is ok, just used to initialize the register object.
    saveModel(quadFig, targetfilename.c_str(), anysrepfilename.c_str(), targetSrepFolder);

    return EXIT_SUCCESS;
}


M3DQuadFigure* thinplatesplinepdmtosrep::computeMeanSrep(vector<TransformType::Pointer > tpsList, bool useTPS, vector<M3DQuadFigure* > quadFigList){
    // initial a quadfig, can be read from any srep
    M3DQuadFigure* meanQuadFig = quadFigList[0];

    int rowNums = meanQuadFig->getRowCount();
    int colNums = meanQuadFig->getColumnCount();

    // Loop each atom, get the mean point
    for(unsigned int u = 0; u < rowNums; u++){
        for(unsigned int v = 0; v < colNums; v++){
            // skeletal point (1)
            PointType tail = meanPoint(tpsList, useTPS, quadFigList, u, v, 1);

            // up boundary point (3)
            PointType up_tip = meanPoint(tpsList, useTPS, quadFigList, u, v, 3);

            // down boundary point (5)
            PointType down_tip = meanPoint(tpsList, useTPS, quadFigList, u, v, 5);

            // Get transformed spoke radius and direction
            double upSpokeRadiu = calculateSpokeLength(tail, up_tip);
            double downSpokeRadiu = calculateSpokeLength(tail, down_tip);

            Vector3D upSpokeDir = calculateSpokeDirection(tail, up_tip);
            Vector3D downSpokeDir = calculateSpokeDirection(tail, down_tip);

            M3DQuadPrimitive* prim = dynamic_cast<M3DQuadPrimitive*>(meanQuadFig->getPrimitivePtr(u,v));

            // set transformed spoke radius and direction to this atom
            prim->setR0(upSpokeRadiu);
            prim->setR1(downSpokeRadiu);
            prim->setU0(upSpokeDir);
            prim->setU1(downSpokeDir);

            // for crest atoms
            if(u==0 || u==rowNums-1 || v==0 || v==colNums-1) {
                // crest boundary point (7)
                PointType crest_tip = meanPoint(tpsList, useTPS, quadFigList, u, v, 7);

                M3DQuadEndPrimitive* endPrim = dynamic_cast<M3DQuadEndPrimitive*>(meanQuadFig->getPrimitivePtr(u, v));

                // Get transformed spoke radius and direction
                double crestSpokeRadiu = calculateSpokeLength(tail, crest_tip);
                Vector3D crestSpokeDir = calculateSpokeDirection(tail, crest_tip);

                // Set transformed spoke radius and direction to this atom
                endPrim->setREnd(crestSpokeRadiu);
                endPrim->setUEnd(crestSpokeDir);
            }

            // The hub position must be set at the end, if not, the quadFigList[0]'s crest spokes will use a transformed tails. Or you can
            // instead "M3DQuadFigure* meanQuadFig = quadFigList[0];" with a clone qiadFig, which will not change quadFigList[0].
            // set transformed skeletal point position
            prim->setX(tail[0], tail[1], tail[2]);
        }
    }

    return meanQuadFig;
}


/* Avarage the spoke direction, the spoke length, and the skeletal position.*/
M3DQuadFigure* thinplatesplinepdmtosrep::compute_mean_srep_new(vector<M3DQuadFigure* > quadFigList){
    // initial a quadfig, can be read from any srep
    M3DQuadFigure* meanQuadFig = quadFigList[0];

    int rowNums = meanQuadFig->getRowCount();
    int colNums = meanQuadFig->getColumnCount();

    // Loop each atom, get the mean point
    for(unsigned int u = 0; u < rowNums; u++){
        for(unsigned int v = 0; v < colNums; v++){
            // compute mean spoke
            Vector3D up_spoke_dir = mean_direction(quadFigList, u, v, 0);
            Vector3D down_spoke_dir = mean_direction(quadFigList, u, v, 1);

            double up_spoke_radius = mean_radius(quadFigList, u, v, 0);
            double down_spoke_radius = mean_radius(quadFigList, u, v, 1);

            M3DQuadPrimitive* prim = dynamic_cast<M3DQuadPrimitive*>(meanQuadFig->getPrimitivePtr(u,v));

            // set transformed spoke radius and direction to this atom
            prim->setR0(up_spoke_radius);
            prim->setR1(down_spoke_radius);
            prim->setU0(up_spoke_dir);
            prim->setU1(down_spoke_dir);

            // for crest atoms
            if(u==0 || u==rowNums-1 || v==0 || v==colNums-1) {
                // compute mean spoke
                Vector3D crest_spoke_dir = mean_direction(quadFigList, u, v, 2);
                double crest_spoke_radius = mean_radius(quadFigList, u, v, 2);

                M3DQuadEndPrimitive* endPrim = dynamic_cast<M3DQuadEndPrimitive*>(meanQuadFig->getPrimitivePtr(u, v));

                // Set transformed spoke radius and direction to this atom
                endPrim->setREnd(crest_spoke_radius);
                endPrim->setUEnd(crest_spoke_dir);
            }

            // skeletal point
            PointType skeletal_position = mean_point(quadFigList, u, v);
            prim->setX(skeletal_position[0], skeletal_position[1], skeletal_position[2]);
        }
    }

    return meanQuadFig;
}


/* compute the mean spoke direction of the atom (u,v)
 * useTPS: if true, use tps to each point; if false, do not use tps, just compute the avarage point.
 * if useTPS = false, tpsList can be pasted as NULL.
*/
Vector3D thinplatesplinepdmtosrep::mean_direction(vector<M3DQuadFigure* > quadFigList, int u, int v, int side){
    PointType mean_dir;
    mean_dir[0] = 0.0;
    mean_dir[1] = 0.0;
    mean_dir[2] = 0.0;

    int templateNum = quadFigList.size();

    // loop each template srep
    for(unsigned int i = 0; i < templateNum; i++){
        M3DQuadFigure* quadFig_i = quadFigList[i];
        M3DQuadPrimitive* prim = dynamic_cast<M3DQuadPrimitive*>(quadFig_i->getPrimitivePtr(u,v));

        switch(side){
        case 0:
        {
            // Up spoke
            Vector3D spoke_dir = prim->getU0();

            mean_dir = sumPoint(mean_dir, Vector3DtoPointType(spoke_dir));
        }
            break;
        case 1:
        {
            // Down spoke
            Vector3D spoke_dir = prim->getU1();

            mean_dir = sumPoint(mean_dir, Vector3DtoPointType(spoke_dir));
        }
            break;
        case 2:
        {
            // Crest spoke
            M3DQuadEndPrimitive* endPrim = dynamic_cast<M3DQuadEndPrimitive*>(quadFig_i->getPrimitivePtr(u, v));
            Vector3D spoke_dir = endPrim->getUEnd();

            mean_dir = sumPoint(mean_dir, Vector3DtoPointType(spoke_dir));
        }
            break;
        }
    }

    // compute the mean spoke direction
    mean_dir[0] /= templateNum;
    mean_dir[1] /= templateNum;
    mean_dir[2] /= templateNum;

    // normalize the direction to be a unit vector
    double distance = 0.0;
    for(unsigned int dim = 0; dim < 3; dim++){
        distance += pow(mean_dir[dim], 2);
    }

    double spokelength = sqrt(distance);

    Vector3D norMeanDir;
    norMeanDir.setX(mean_dir[0] / spokelength);
    norMeanDir.setY(mean_dir[1] / spokelength);
    norMeanDir.setZ(mean_dir[2] / spokelength);

    return norMeanDir;
}


/* compute the mean spoke radius of the atom (u,v)
*/
double thinplatesplinepdmtosrep::mean_radius(vector<M3DQuadFigure* > quadFigList, int u, int v, int side){
    double mean_radius = 0.0;

    int templateNum = quadFigList.size();

    // loop each template srep
    for(unsigned int i = 0; i < templateNum; i++){
        M3DQuadFigure* quadFig_i = quadFigList[i];
        M3DQuadPrimitive* prim = dynamic_cast<M3DQuadPrimitive*>(quadFig_i->getPrimitivePtr(u,v));

        switch(side){
        case 0:
            mean_radius += prim->getR0();
            break;
        case 1: // Down spoke
            mean_radius += prim->getR1();
            break;
        case 2: // Crest spoke
            M3DQuadEndPrimitive* endPrim = dynamic_cast<M3DQuadEndPrimitive*>(quadFig_i->getPrimitivePtr(u, v));

            mean_radius += endPrim->getREnd();
            break;
        }
    }

    // compute the mean spoke radius
    mean_radius /= templateNum;

    return mean_radius;
}


/* compute the mean spoke of atom (u,v, side)
 * useTPS: if true, use tps to each point; if false, do not use tps, just compute the avarage point.
 * if useTPS = false, tpsList can be pasted as NULL.
*/
thinplatesplinepdmtosrep::PointType thinplatesplinepdmtosrep::mean_point(vector<M3DQuadFigure* > quadFigList, int u, int v){
    PointType meanPoint;
    meanPoint[0] = 0.0;
    meanPoint[1] = 0.0;
    meanPoint[2] = 0.0;

    int templateNum = quadFigList.size();

    // loop each template srep
    for(unsigned int i = 0; i < templateNum; i++){
        M3DQuadFigure* quadFig_i = quadFigList[i];
        M3DQuadPrimitive* prim = dynamic_cast<M3DQuadPrimitive*>(quadFig_i->getPrimitivePtr(u,v));

        // Skeletal point
        Vector3D tail = prim->getX();

        meanPoint = sumPoint(meanPoint, Vector3DtoPointType(tail));
    }

    // compute the mean point
    for(unsigned dim = 0;  dim < 3; dim++){
        meanPoint[dim] /= templateNum;
    }

    return meanPoint;
}


/* compute the mean atom
 * useTPS: if true, use tps to each point; if false, do not use tps, just compute the avarage point.
 * if useTPS = false, tpsList can be pasted as NULL.
*/
thinplatesplinepdmtosrep::PointType thinplatesplinepdmtosrep::meanPoint(vector<TransformType::Pointer > tpsList, bool useTPS,
                                                                        vector<M3DQuadFigure* > quadFigList, int u, int v, int pointSide){
    PointType point;
    point[0] = 0.0;
    point[1] = 0.0;
    point[2] = 0.0;

    int templateNum = quadFigList.size();

    // loop each template srep
    for(unsigned int i = 0; i < templateNum; i++){
        M3DQuadFigure* quadFig_i = quadFigList[i];
        M3DQuadPrimitive* prim = dynamic_cast<M3DQuadPrimitive*>(quadFig_i->getPrimitivePtr(u,v));

        // Skeletal point
        Vector3D tail = prim->getX();

        switch(pointSide){
        case 1:
        {
            // skeletal point
            if(useTPS) {                
                PointType transformed_tail = tpsList[i]->TransformPoint(Vector3DtoPointType(tail));
                point = sumPoint(point, transformed_tail);
            }
            else {
                point = sumPoint(point, Vector3DtoPointType(tail));                
            }
        }
            break;
        case 3:
        {
            // Up boundary point
            Vector3D up_tip = prim->getX() + prim->getR0()*prim->getU0();

            if(useTPS) {                
                PointType transformed_up_tip = tpsList[i]->TransformPoint(Vector3DtoPointType(up_tip));
                point = sumPoint(point, transformed_up_tip);
            }
            else {
                point = sumPoint(point, Vector3DtoPointType(up_tip));                
            }
        }
            break;
        case 5:
        {
            // Down boundary point
            Vector3D down_tip = prim->getX() + prim->getR1()*prim->getU1();

            if(useTPS){
                PointType transformed_down_tip = tpsList[i]->TransformPoint(Vector3DtoPointType(down_tip));
                point = sumPoint(point, transformed_down_tip);
            }
            else {
                point = sumPoint(point, Vector3DtoPointType(down_tip));
            }
        }
            break;
        case 7:
        {
            M3DQuadEndPrimitive* endPrim = dynamic_cast<M3DQuadEndPrimitive*>(quadFig_i->getPrimitivePtr(u, v));
            Vector3D crest_tip = endPrim->getX() + endPrim->getREnd()*endPrim->getUEnd();
            //std::cout<<"====endPrim->getX(): "<<endPrim->getX()<<endl; // equals to tail

            if(useTPS){
                PointType transformed_crest_tip = tpsList[i]->TransformPoint(Vector3DtoPointType(crest_tip));
                point = sumPoint(point, transformed_crest_tip);
            }
            else {
                point = sumPoint(point, Vector3DtoPointType(crest_tip));
            }
        }
            break;
        }
    }

    // compute the mean point
    PointType meanPoint;
    meanPoint[0] = point[0] / templateNum;
    meanPoint[1] = point[1] / templateNum;
    meanPoint[2] = point[2] / templateNum;

    return meanPoint;
}



/* listName: each line stores a pathname
*/
vector< std::string > thinplatesplinepdmtosrep::readTempalte( const char* listName) {
    vector< std::string > filePathNames;

    std::ifstream listFile;
    listFile.open( listName, std::ios_base::in );
    if (!listFile) {
      std::cerr << "Unable to open shape list \"" << listName << "\"!" << std::endl;
      EXIT_FAILURE;
    }

    // parse the list file
    std::string currentFileName;
    while (!listFile.eof()) {
      std::getline( listFile, currentFileName );
      if (currentFileName.empty() || currentFileName[0]=='#') { continue; }

      filePathNames.push_back(currentFileName);
    }

    listFile.close();

    return filePathNames;
}


thinplatesplinepdmtosrep::PointType thinplatesplinepdmtosrep::sumPoint(PointType point1, PointType point2){
    PointType point;
    for(unsigned dim = 0; dim < 3; dim++){
        point[dim] = point1[dim] + point2[dim];
    }

    return point;
}



/* This mean srep is simply the avarage of the input sreps (Arithmetic mean).
 * used to validate the multiply template tps by compute each srep from the single template, and then compute the mean of the target sreps,
 * to check if the average of single template equals to the multiply template tps.
 * it should be the same.
*/
void thinplatesplinepdmtosrep::avarage_srep(vector< std::string > srepPathNames, const char* outFileName) {
    int srepNum = srepPathNames.size();

    // read template sreps into a quad list
    vector<M3DQuadFigure* > quadFigList;
    for(unsigned int i = 0; i < srepNum; i++){
        // Read in source s-rep.
        M3DQuadFigure* quadFig = tls.GetQuadFigure(srepPathNames[i].c_str());

        quadFigList.push_back(quadFig);
    }

    vector<TransformType::Pointer > nullVal;
    M3DQuadFigure* avarageQuadFig = computeMeanSrep(nullVal, false, quadFigList);

    // save the mean s-rep to m3d file
    Registry registry;
    registry.readFromFile(srepPathNames[0].c_str(), false);

    // Write the changed quadFig to m3d file.
    char figureStr[1024];
    sprintf(figureStr, "model.figure[%d]", 0);
    avarageQuadFig->writeFigure(figureStr,registry); //write the figure into registry.

    cout<<"----Write srep to: "<< outFileName << endl;
    registry.writeToFile(outFileName);
    cout<<"----Finished!" << endl;
}



/* This mean srep is computed by: avarage the skeletal position, avarage the spoke direction, avarage the spoke radii.
 * SINCE AVARAGE THE BOUNDARY POINTS as done in avarage_srep DOES NOT MAKE SENCE.
*/
void thinplatesplinepdmtosrep::avarage_srep_new(vector< std::string > srepPathNames, const char* outFileName) {
    int srepNum = srepPathNames.size();

    // read template sreps into a quad list
    vector<M3DQuadFigure* > quadFigList;
    for(unsigned int i = 0; i < srepNum; i++){
        // Read in source s-rep.
        M3DQuadFigure* quadFig = tls.GetQuadFigure(srepPathNames[i].c_str());

        quadFigList.push_back(quadFig);
    }

    M3DQuadFigure* avarageQuadFig = compute_mean_srep_new(quadFigList);

    // save the mean s-rep to m3d file
    Registry registry;
    registry.readFromFile(srepPathNames[0].c_str(), false);

    // Write the changed quadFig to m3d file.
    char figureStr[1024];
    sprintf(figureStr, "model.figure[%d]", 0);
    avarageQuadFig->writeFigure(figureStr,registry); //write the figure into registry.

    cout<<"----Write srep to: "<< outFileName << endl;
    registry.writeToFile(outFileName);
    cout<<"----Finished!" << endl;
}


/* This mean srep is computed by calling the matlab CPNS code.
 * The spokes are Euclideanized to the composite spaece and compute the mean, and then Euclideanized back to the s-rep space.
*/
void thinplatesplinepdmtosrep::cpns_avarage_srep(vector< std::string > srepPathNames, const char* outFileName) {
    int srepNum = srepPathNames.size();

    // read template sreps into a quad list
    vector<M3DQuadFigure* > quadFigList;
    for(unsigned int i = 0; i < srepNum; i++){
        // Read in source s-rep.
        M3DQuadFigure* quadFig = tls.GetQuadFigure(srepPathNames[i].c_str());

        quadFigList.push_back(quadFig);
    }

    M3DQuadFigure* avarageQuadFig = compute_mean_srep_new(quadFigList);

    // save the mean s-rep to m3d file
    Registry registry;
    registry.readFromFile(srepPathNames[0].c_str(), false);

    // Write the changed quadFig to m3d file.
    char figureStr[1024];
    sprintf(figureStr, "model.figure[%d]", 0);
    avarageQuadFig->writeFigure(figureStr,registry); //write the figure into registry.

    cout<<"----Write srep to: "<< outFileName << endl;
    registry.writeToFile(outFileName);
    cout<<"----Finished!" << endl;
}









