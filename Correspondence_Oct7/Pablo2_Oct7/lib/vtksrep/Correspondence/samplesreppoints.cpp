#include "samplesreppoints.h"

using namespace std;

samplesreppoints::samplesreppoints()
{
}


/* This function sample the s-rep points with the shared skeletal points counted for each spoke.
 *
*/
MatrixType samplesreppoints::samplePointsWithRepeatedSkeletal(M3DQuadFigure* quadFig){

    int rowNum = quadFig->getRowCount();
    int colNum = quadFig->getColumnCount();
    int crestAtomNum = rowNum*2 + colNum*2 - 4;
    int interiorAtomNum = rowNum*colNum - crestAtomNum;

    int spokeNum = interiorAtomNum*2 + crestAtomNum*3;//106
    int pointNum = spokeNum;// * 2; //212

    MatrixType points(3, pointNum);

    M3DQuadPrimitive* prim;
    M3DQuadEndPrimitive* endPrim;

    Vector3D spokeInfo;

    // Loop each atom of the quadfig
    int pointIndex = 0;
    for(unsigned u = 0; u < rowNum; u++){
        for(unsigned v = 0; v < colNum; v++){

            prim = dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v));

            // Get up spoke tail and tip
            /*spokeInfo = prim->getX();
            points[0][pointIndex] = spokeInfo.getX();
            points[1][pointIndex] = spokeInfo.getY();
            points[2][pointIndex] = spokeInfo.getZ();
            pointIndex++;*/

            spokeInfo = prim->getX() + prim->getR0()*prim->getU0();
            points[0][pointIndex] = spokeInfo.getX();
            points[1][pointIndex] = spokeInfo.getY();
            points[2][pointIndex] = spokeInfo.getZ();
            pointIndex++;

            // Get down spoke tail and tip
            /*spokeInfo = prim->getX();
            points[0][pointIndex] = spokeInfo.getX();
            points[1][pointIndex] = spokeInfo.getY();
            points[2][pointIndex] = spokeInfo.getZ();
            pointIndex++;*/

            spokeInfo = prim->getX() + prim->getR1()*prim->getU1();
            points[0][pointIndex] = spokeInfo.getX();
            points[1][pointIndex] = spokeInfo.getY();
            points[2][pointIndex] = spokeInfo.getZ();
            pointIndex++;

            // Get crest spoke tail and tip
            if(u==0 || u==rowNum-1 || v==0 || v==colNum-1){
                endPrim = dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(u, v));

                /*spokeInfo = endPrim->getX();
                points[0][pointIndex] = spokeInfo.getX();
                points[1][pointIndex] = spokeInfo.getY();
                points[2][pointIndex] = spokeInfo.getZ();
                pointIndex++;*/

                spokeInfo = endPrim->getX() + endPrim->getREnd()*endPrim->getUEnd();
                points[0][pointIndex] = spokeInfo.getX();
                points[1][pointIndex] = spokeInfo.getY();
                points[2][pointIndex] = spokeInfo.getZ();
                pointIndex++;
            }
        }
    }

    return points;
}



/* This function sample the s-rep points with each skeletal point count once.
 *
*/
/*void samplesreppoints::noRepeatSkeletalPoints(){

}*/



///* repeatSkeletalPoints: if true, sample points with repeated spoke tails;
// * if false: sample points with no repeated spoke tails.
//*/
//vector<MatrixType> samplesreppoints::sampleSrepPoints(const char* srepFolder, const char* outputFolder, bool repeatSkeletalPoints){
//    // Get the srep filenames from given folder
//    vector< std::string > inputSreps;
//    inputSreps = tls.getModelname(string(srepFolder));

//    // Read in each sreps
//    vector<MatrixType> srepPointsOfTrainingSet; // srep points of the training set.
//    MatrixType srepPoints;
//    for(unsigned int i = 0; i < inputSreps.size(); i++ ){
//        const char* filename = inputSreps[i].c_str();
//        cout<<"---filename: "<<filename<<endl;
//        M3DQuadFigure* quadfig = tls.GetQuadFigure(filename);

//        if(repeatSkeletalPoints) {
//            // sample repeat skeletal points
//            srepPoints = samplePointsWithRepeatedSkeletal(quadfig);
//        }
//        else {
//            // sample no repeat skeletal points
//            //srepPoints = noRepeatSkeletalPoints(quadfig);
//        }

//        // Save to srepPoints
//        srepPointsOfTrainingSet.push_back(srepPoints);

//        // Save shapes to vtk file
//        string bfname = tls.getBaseFileName(inputSreps[i]);
//        string outputFile = std::string(outputFolder) +  bfname + ".vtk";
//        saveToVTK(srepPoints, outputFile);
//    }

//    // Save the whole training set to one file, row is a s-rep, each column is a feature.
//    string outFile = std::string(outputFolder) +  "srep_PDMs.txt";
//    saveToText(srepPointsOfTrainingSet, outFile.c_str());

//    return srepPointsOfTrainingSet;
//}





/* Sample the optimized s-rep point. The skeletal sheet split to three part: up side, down side and fold curve.
*/
void samplesreppoints::sampleSrepPoints(const char* srepFolder, vector< std::string > filename, const char* outputFolder, const char* coeffFolder,
                                                             int rowNum, int colNum, int type, int interpolationLevel){
    int crestAtomNum = rowNum*2 + colNum*2 - 4; // 28. crest spoke number, each spoke has one variable.
    int interiorAtomNum = rowNum*colNum - crestAtomNum; // 11
    int varStand = crestAtomNum*1 - 4 + interiorAtomNum*2; // 46. variables for standard spokes (up or down spokes).

    int srepNum = filename.size();

    // reserve space for storing the sampled points of each s-rep in the training set
    vector<vtkSmartPointer<vtkPoints> > sampledPoints;
    for(unsigned int i = 0; i < srepNum; i++){
        vtkSmartPointer<vtkPoints> temp = vtkSmartPointer< vtkPoints >::New();
        sampledPoints.push_back(temp);
    }

    // get the files storing the coeff string
    string up_str_file, down_str_file, crest_str_file;
    if(!strcmp ( "0", coeffFolder)){ // If set to sample the srep before correspondence
        up_str_file = "0";
        down_str_file = "0";
        crest_str_file = "0";
    }
    else {
        // /NIRAL/work/ltu/WorkSpace/Test/30AllAtomsFromAtomStage/after_optimized_coeffs
        up_str_file = string(coeffFolder) + "/up_str.txt";
        down_str_file = string(coeffFolder) + "/down_str.txt";
        crest_str_file = string(coeffFolder) + "/crest_str.txt";
    }

    // Sample points
    sampleShiftedSrep(srepFolder, filename, 0, up_str_file, varStand, sampledPoints, type, interpolationLevel); // up side
    sampleShiftedSrep(srepFolder, filename, 1, down_str_file, varStand, sampledPoints, type, interpolationLevel); //down side
    sampleShiftedSrep(srepFolder, filename, 2, crest_str_file, crestAtomNum, sampledPoints, type, interpolationLevel); //crest side

    // Save sampled points of each s-rep into a .vtk file
    for(unsigned int i = 0; i < srepNum; i++){
        string bfname = tls.getBaseFileName(filename[i]);
        string outputFile = std::string(outputFolder) +  bfname + ".vtk";
        //cout<<"---------------------sampledPoints size: "<<sampledPoints[i]->GetNumberOfPoints()<<endl;
        // save the ith srep's points.
        saveVtkPointsToVTK(sampledPoints[i], outputFile);
        cout<<"----Finish saving the "<<i<<"th srep's  to: "<<outputFile<<endl;
    }

    // save the whole traing set points to file, for later use (input to PNS main to Euclidenize).
    //string outFile = std::string(outputFolder) +  "srep_PDMs_afterCorrespondence.txt";
    //saveVtkPointsVectorToText(sampledPoints, outFile.c_str());
    //cout<<"------Successfully saved the srep points of the training set into file: "<<outFile<<endl;
}



/* The sampled point will be stored in X.
*/
void samplesreppoints::sampleShiftedSrep(const char* srepFolder, vector< std::string > filename, int side, string vars_str_file, int varNums,
                                         vector<vtkSmartPointer<vtkPoints> > &sampledPoints,int type, int interpolationLevel){
    // sample the srep before correspondence
    if(!strcmp ( "0", vars_str_file.c_str())){
        vector<double> vars;
        // coeffs all set to 0.
        for(unsigned int i =0 ; i < varNums; i++){
            vars.push_back(0);
        }

        for(unsigned int q = 0; q < filename.size(); q++){ //Loop each srep
            // Get this srep
            string fn = string(srepFolder) + filename[q];
            M3DQuadFigure* quadfig = tls.GetQuadFigure(fn.c_str());

            // Shift the srep and get the shifted surface points
            vtkSmartPointer< vtkPoints > surfacePos = vtkSmartPointer< vtkPoints >::New();
            getSurfacePoints(quadfig, side, vars, varNums, surfacePos, type, interpolationLevel);

            // save the points to corresponding srep entry.
            for(unsigned int m = 0; m < surfacePos->GetNumberOfPoints(); m++){ // for each point
                double p[3];
                surfacePos->GetPoint(m, p);
                sampledPoints[q]->InsertNextPoint(p[0], p[1], p[2]);
            }
        }
    }

    // sample the srep after correspondence
    else {
        //stringstream ss(vars_str);
        ifstream ss(vars_str_file.c_str());

        if (! ss.is_open()){
            std::cerr << "Msg from samplesreppoints::sampleShiftedSrep: Unable to open shape list \"" << vars_str_file << "\"!" << std::endl;
            EXIT_FAILURE;
        }

        vector<double> vars;

        for(unsigned int q = 0; q < filename.size(); q++){ //Loop each srep
            // Get this srep
            string fn = string(srepFolder) + filename[q];
            M3DQuadFigure* quadfig = tls.GetQuadFigure(fn.c_str());

            // Get up side points
            string token;
            int i = 0;
            while (i < varNums) {
                ss >> token;
                vars.push_back(atof(token.c_str()));
                i++;
            }
            //cout<<"---------------------------var size: "<< vars.size()<<endl;

            // Shift the srep and get the shifted surface points
            vtkSmartPointer< vtkPoints > surfacePos = vtkSmartPointer< vtkPoints >::New();
            getSurfacePoints(quadfig, side, vars, varNums, surfacePos, type, interpolationLevel);

            // Clear for next
            vars.clear();

            // save the points to corresponding srep entry.
            for(unsigned int m = 0; m < surfacePos->GetNumberOfPoints(); m++){ // for each point
                double p[3];
                surfacePos->GetPoint(m, p);
                sampledPoints[q]->InsertNextPoint(p[0], p[1], p[2]);
            }
        }
    }
}


/* Get surface point and save to the reference passing para surfacePoints.
 * This function is diff from getPoints() by getting the surface point under specific interpolation level.
*/
void samplesreppoints::getSurfacePoints(M3DQuadFigure* quadFig, int side, vector<double> vars, int varNums,
                                        vtkSmartPointer< vtkPoints > surfacePos, int type, int interpolationLevel){
    int rowNum = quadFig->getRowCount();
    int colNum = quadFig->getColumnCount();
    // Convert coeffs from vector to pointer
    double *coeff = new double [varNums]; // storing up or down spokes variables.
    for(unsigned int i = 0; i < varNums; i++){
        coeff[i] = vars[i];
    }

    M3DQuadFigure* tempQuadFig = dynamic_cast<M3DQuadFigure*>(quadFig->clone()); //Move all base on original srep.

    switch (side){
    case 0 : // up or down
    case 1 :
    {        
        slidestandardspokes mStandardSpoke;
        // Shift spokes by given vars.
        mStandardSpoke.updateSpokeInfo(tempQuadFig, coeff, side);
    }
        break;
    case 2 : // crest
    {        
        slidecrestspokes mCrestSpoke;
        // Shift crest spokes
        visualizecrest vcrest;
        vtkSmartPointer<vtkSRep> srepFig = vcrest.getSrepFig(tempQuadFig);
        mCrestSpoke.updateCrestSpokeInfo(srepFig, tempQuadFig, coeff); //Move base on original srep.
    }
        break;
    }

    getSrepPoints(tempQuadFig, side, surfacePos, type, interpolationLevel);

    delete coeff;
}


/* type: 1 (only collect skeletal point), 3 (only boundary point), 5 (skeletal point + boundary point)
 * interpolationLevel: 0 (sample the original s-rep spoke's), 1,...n (sample the points under interpolation level i)
*/
void samplesreppoints::getSrepPoints(M3DQuadFigure* quadFig, int side, vtkSmartPointer< vtkPoints > surfacePos, int type, int interpolationLevel){
    if(side == 2){        
        // Get crest side surface points
        getCrestSidePoints(quadFig, surfacePos, interpolationLevel,type);
    }
    else if(side == 0 || side==1){
        // Get up or down side surface points
        getUpOrDownSidePoints(quadFig, surfacePos, interpolationLevel, side, type);
    }
    else {
        cout<<"Message from samplesreppoints::getSurfacePoints:: Invaild side type!!"<<endl;        
        exit(1);
    }
}


/* The sampled points will be store in pos and passed to previous function.
*/
void samplesreppoints::getUpOrDownSidePoints(M3DQuadFigure* quadFig, vtkSmartPointer< vtkPoints > pos, int interpolationLevel, int side,
                                                     int type){
    int rowNums = quadFig->getRowCount();
    int colNums = quadFig->getColumnCount();
    int step = pow((double)2, (double)interpolationLevel);

    // Get interpolated u & v coordinate, storing in interpolatedU & interpolatedV
    vector<double> interpolatedU;
    vector<double> interpolatedV;
    tls.getUVCoordinate(interpolatedU, interpolatedV, rowNums, colNums, step);

    // Get all the interpolated boundary points (on the surface)
    M3DQuadInterpolater *tpm = new M3DQuadInterpolater(quadFig);
    Vector3D point;
    for(unsigned int i = 0; i < interpolatedU.size(); i++){
        if(type==1 || type==5) { // Get skeletal points
            point = tpm->interpolateQuadSpoke(quadFig,interpolatedU[i],interpolatedV[i], side)->getX();
            pos->InsertNextPoint(point.getX(),point.getY(),point.getZ());
        }
        if(type==3 || type==5) { // Get boundary points
            point = tpm->interpolateQuadSpoke(quadFig,interpolatedU[i],interpolatedV[i], side)->getB();
            pos->InsertNextPoint(point.getX(),point.getY(),point.getZ());
        }
    }

    delete tpm;
}


/* For those crest spokes' share a same tail, that (interpolated)skeletal point only count once. */
void samplesreppoints::getCrestSidePoints(M3DQuadFigure* quadFig, vtkSmartPointer< vtkPoints > pos, int interpolationLevel, int type){

    // New a crest interplator....
    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
    visualizecrest visCrest;
    vtkSmartPointer<vtkSRep> srepfig = visCrest.getSrepFig(quadFig);
    interpolatecrestspokes->SetInput(srepfig);
    interpolatecrestspokes->SetInterpolationLevel(interpolationLevel);
    interpolatecrestspokes->Update();

    int rowNums = quadFig->getRowCount();
    int colNums = quadFig->getColumnCount();
    int crestQuadNum = rowNums*2 + (colNums-2)*2;
    int step = pow((double)2, (double)interpolationLevel);

    vector<double> subpoints; //store the v coordinate of all the sub-point along spoke tip connection curve.

    // Along the curve between two correspondence up and down spokes.(v direction in crest interpolate method)
    tls.splitLine(0,1,step, subpoints); // Up spoke tip is consider as 0, down spoke 1. The medial crest is 0.5.

    // Loop each crest quad.
    //int counter= 0;
    for(unsigned int i = 0; i< crestQuadNum; i++){
        for(unsigned int m = 0; m < subpoints.size()-1; m++){ //loop each sub points in v direction
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s1 = interpolatecrestspokes->GetInterpolatedPoint(i,subpoints[m]);

            // Get the skeletal points.
            if(type==1 || type==5) {
                pos->InsertNextPoint(s1[0], s1[1], s1[2]);
            }

            // Get the boundary points
            if(type==3 || type==5) {
                if(interpolationLevel==0){ // If true, add the fold one spoke tip.
                    // Get each interpolated boundary spoke direction at this point.
                    vtkSRep::VNLType p0 = interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints[m],0.5); //here subpoints[m] = 0;
                    pos->InsertNextPoint(s1[0] + p0[0], s1[1] + p0[1], s1[2] + p0[2]);
                }
                else {
                    for(unsigned int n = 0; n < subpoints.size(); n++){ //loop each sub points in t direction, without discard the up and down tip pos.
                        // Get each interpolated boundary spoke direction at this point.
                        vtkSRep::VNLType p0 = interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints[m],subpoints[n]);
                        pos->InsertNextPoint(s1[0] + p0[0], s1[1] + p0[1], s1[2] + p0[2]);
                        //counter++;
                        //cout<<"----"<<counter<<endl;
                    }
                }
            }
        }
    }
}



void samplesreppoints::saveToText(vector<MatrixType> data, const char* filename){
    std::ofstream fout;
    fout.open(filename);

    if(fout)  {
        for(int i =0; i< data.size();i++){
            MatrixType X = data[i];
            for(unsigned int c = 0; c < X.columns(); c++){
                for(unsigned int r = 0; r < 3; r++){
                    fout<< X[r][c] <<" ";
                }
            }
            fout<<endl;
        }

        cout<<"Successfully saved data to: "<<filename<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();
}




/* Save the srep points to a vtk file.
 * Matrix X is: 3-by-n
*/
void samplesreppoints::saveToVTK(MatrixType X, string outputFile){

    // Save each new shape points to vtk file
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New(); // same as vtkPolyData* polyData;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    for(unsigned int c = 0; c < X.columns(); c++){
        points->InsertNextPoint(X[0][c], X[1][c], X[2][c]);
    }

    polyData->SetPoints(points);
    vtkSmartPointer<vtkPolyDataWriter> polyDataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    polyDataWriter->SetInput(polyData);

    polyDataWriter->SetFileName(outputFile.c_str());
    polyDataWriter->Update();
}


/* save vtkPoints vector to a txt file. Each row is a s-rep, each column is a feature.*/
void samplesreppoints::saveVtkPointsVectorToText(vector<vtkSmartPointer<vtkPoints> > data, const char* filename){
    std::ofstream fout;
    fout.open(filename);

    if(fout)  {
        // for each s-rep
        for(unsigned int i = 0; i < data.size();i++){
            vtkSmartPointer<vtkPoints> points = data[i];
            // for each point
            for(unsigned int j = 0; j < points->GetNumberOfPoints(); j++){
                double p[3];
                points->GetPoint(j, p);
                fout << p[0] << " " << p[1]<<" " << p[2] << " ";
            }
            fout<<endl;
        }

        cout<<"Successfully saved data to: "<<filename<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();
}


/* Save the srep points to a vtk file.
 * points: store all the points of a srep.
*/
void samplesreppoints::saveVtkPointsToVTK(vtkSmartPointer<vtkPoints> points, string outputFile){

    // Save each new shape points to vtk file
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New(); // same as vtkPolyData* polyData;

    polyData->SetPoints(points);
    vtkSmartPointer<vtkPolyDataWriter> polyDataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    polyDataWriter->SetInput(polyData);

    polyDataWriter->SetFileName(outputFile.c_str());
    polyDataWriter->Update();
}


/* Save the srep points to a vtk file.
 * Matrix X is: 3-by-n, holding the s-rep boundary points in sequence:
 * 39 up points, 39 down points, 28 fold points.
*/
/*void samplesreppoints::saveToVTK_2(MatrixType X, string outputFile){

    // Save each new shape points to vtk file
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New(); // same as vtkPolyData* polyData;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    for(unsigned int c = 0; c < X.columns(); c++){
        points->InsertNextPoint(X[0][c], X[1][c], X[2][c]);
    }

    polyData->SetPoints(points);

    int upSpokeNum = row*col;
    int foldSpokeNum = 2*(row+col-2);

    //create triangles on the points in the polydata
    // Up side triangles
    for(unsigned int i = 0; i < col-1; i++){
        for(unsigned int j = 0; j < row-1; j++){
            vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
            triangle1->GetPointIds()->SetId ( 0, i + col*j);
            triangle1->GetPointIds()->SetId ( 1, i + col*(j+1));
            triangle1->GetPointIds()->SetId ( 2, i + 1 + col*(j+1));

            vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
            triangle2->GetPointIds()->SetId ( 0, i + col*j);
            triangle2->GetPointIds()->SetId ( 1, i + 1 + col*j);
            triangle2->GetPointIds()->SetId ( 2, i + 1 + col*(j+1));
        }
    }

    // Down side triangle
    for(unsigned int i = 0; i < col-1; i++){
        for(unsigned int j = 0; j < row-1; j++){
            vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
            triangle1->GetPointIds()->SetId ( 0, upSpokeNum + i + col*j);
            triangle1->GetPointIds()->SetId ( 1, upSpokeNum + i + col*(j+1));
            triangle1->GetPointIds()->SetId ( 2, upSpokeNum + i + 1 + col*(j+1));

            vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
            triangle2->GetPointIds()->SetId ( 0, upSpokeNum + i + col*j);
            triangle2->GetPointIds()->SetId ( 1, upSpokeNum + i + 1 + col*j);
            triangle2->GetPointIds()->SetId ( 2, upSpokeNum + i + 1 + col*(j+1));
        }
    }

    // Fold triangle: up crest tip + fold spoke tip
    // The first col-2 fold spokes
    for(unsigned int i = 0; i < col-1; i++){
        vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
        triangle1->GetPointIds()->SetId ( 0, i);
        triangle1->GetPointIds()->SetId ( 1, i + 1);
        triangle1->GetPointIds()->SetId ( 2, 2*upSpokeNum + i + 1);

        vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
        triangle2->GetPointIds()->SetId ( 0, i);
        triangle2->GetPointIds()->SetId ( 1, 2*upSpokeNum + i);
        triangle2->GetPointIds()->SetId ( 2, 2*upSpokeNum + i + 1);
    }
    // The v=col-1
    for(unsigned int r = 0; r < row-1; r++){
        int atomIndex = col-1+r*col;
        vtkSmartPointer<vtkTriangle> triangle_col0 = vtkSmartPointer<vtkTriangle>::New();
        triangle_col0->GetPointIds()->SetId ( 0, atomIndex);
        triangle_col0->GetPointIds()->SetId ( 1, 2*upSpokeNum + col-1 + 2*(r+1));
        triangle_col0->GetPointIds()->SetId ( 2, atomIndex + col);

        vtkSmartPointer<vtkTriangle> triangle_col1 = vtkSmartPointer<vtkTriangle>::New();
        triangle_col1->GetPointIds()->SetId ( 0, atomIndex);
        triangle_col1->GetPointIds()->SetId ( 1, 2*upSpokeNum + col-1 + 2*(r+1));
        triangle_col1->GetPointIds()->SetId ( 2, 2*upSpokeNum + col-1 + 2*r);
    }
    // The u=row-1, v=col-1,...,1
    int crestIndex = col + (row-1-1)*2 + col;
    for(unsigned int v = col-1; v>0; v--){
        vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
        triangle1->GetPointIds()->SetId ( 0, v);
        triangle1->GetPointIds()->SetId ( 1, v - 1);
        triangle1->GetPointIds()->SetId ( 2, 2*upSpokeNum + crestIndex - 1);

        vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
        triangle2->GetPointIds()->SetId ( 0, v);
        triangle2->GetPointIds()->SetId ( 1, 2*upSpokeNum + crestIndex);
        triangle2->GetPointIds()->SetId ( 2, 2*upSpokeNum + crestIndex - 1);

        crestIndex--;
    }
    // The v=0, u=row-1,...0
    for(unsigned int u = row-1; u > 0; u--){
        int atomIndex = u*col;
        crestInce = col + u*1;
        vtkSmartPointer<vtkTriangle> triangle_col0 = vtkSmartPointer<vtkTriangle>::New();
        triangle_col0->GetPointIds()->SetId ( 0, u*col);
        triangle_col0->GetPointIds()->SetId ( 1, u*col - col);
        triangle_col0->GetPointIds()->SetId ( 2, 2*upSpokeNum + col);

        vtkSmartPointer<vtkTriangle> triangle_col1 = vtkSmartPointer<vtkTriangle>::New();
        triangle_col1->GetPointIds()->SetId ( 0, u*col);
        triangle_col1->GetPointIds()->SetId ( 1, 2*upSpokeNum + (u-1)*2 + col);
        triangle_col1->GetPointIds()->SetId ( 2, 2*upSpokeNum + col);
    }




    // Get the up side crest tips



      vtkSmartPointer<vtkTriangle> triangle2 =
        vtkSmartPointer<vtkTriangle>::New();
      triangle2->GetPointIds()->SetId ( 0, 2 );
      triangle2->GetPointIds()->SetId ( 1, 3 );
      triangle2->GetPointIds()->SetId ( 2, 0 );



    //add the triangles to the list of triangles
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    triangles->InsertNextCell(triangle1);
    triangles->InsertNextCell(triangle2);



    vtkSmartPointer<vtkPolyDataWriter> polyDataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    polyDataWriter->SetInput(polyData);

    polyDataWriter->SetFileName(outputFile.c_str());
    polyDataWriter->Update();
}*/

