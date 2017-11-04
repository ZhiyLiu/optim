#include "shiftedsrep.h"


typedef vnl_matrix<double> MatrixType;

shiftedsrep::shiftedsrep()
{
}

shiftedsrep::shiftedsrep(const char* rootDir, vector< std::string > srepNames, double rX, double rY, double rZ, double opacity)
{
    this->rootDir = rootDir;
    this->srepNames = srepNames;

    this->crestAtomNum = rowNum * 2 + colNum * 2 - 4; // 28. crest spoke number, each spoke has one variable.
    this->interiorAtomNum = rowNum * colNum - this->crestAtomNum; // 11
    this->varStand = this->crestAtomNum * 1 - 4 + this->interiorAtomNum * 2; // 46. variables for standard spokes (up or down spokes).

    this->rX = rX;
    this->rY = rY;
    this->rZ = rZ;
    this->opacity = opacity;

    this->srepNum = srepNames.size();

    // Read in all the training s-reps, store in member variable "quadFigList".
    readTrainingSreps();
}


shiftedsrep::shiftedsrep(double rX, double rY, double rZ, double opacity){
    this->rX = rX;
    this->rY = rY;
    this->rZ = rZ;
    this->opacity = opacity;
}


/* Compute the mean srep. The mean of each side is compute seperately.
*/
void shiftedsrep::computeMeanFig(string up_str, string down_str, string crest_str){

    // Move and update each side of each s-rep
    computeMeanSide(up_str, 0, this->varStand);
    computeMeanSide(down_str, 1, this->varStand);
    computeMeanSide(crest_str, 2, this->crestAtomNum);
}

/* Return the mean srep. up or down side has 39 points, crest has 28. total 106 points for each srep's boundary.
 * This function only compute the mean srep, if want higher level surface point,
 * can call getMeanSurfacePoints(interpolationLevel), use the mean srep and reinterpolate.
*/
void shiftedsrep::computeMeanSide(string vars_str, int side, int varNums){
    stringstream ss(vars_str);

    // Move and update each side of each s-rep
    int srepNum = this->quadFigList.size();
    vector<MatrixType> boundaryPoints, skeletalPoints;
    vector<double> vars;

    for(unsigned int q = 0; q < srepNum; q++){ //Loop each srep
        // Get up side points
        string token;
        int i = 0;
        while (i < varNums) {
            ss >> token;
            vars.push_back(atof(token.c_str()));
            //cout<<"---------------------------var "<< vars[i]<<endl;
            i++;
        }
        //cout<<"---------------------------var size: "<< vars.size()<<endl;
        getPoints(this->quadFigList[q], side, vars, varNums, boundaryPoints, skeletalPoints);

        // Clear for next
        vars.clear();
    }

    // Calculate the mean point across all the sreps
    int pointsNum = boundaryPoints[0].columns();
    //cout<<"--------------boundaryPoints: "<<pointsNum<<endl; //28 for crest, 39 for up or down.
    MatrixType mean_bpos(3,pointsNum);
    MatrixType mean_spos(3,pointsNum);
    mean_bpos = meanSidePointPair(boundaryPoints);
    mean_spos = meanSidePointPair(skeletalPoints);

    // Calculate the mean side srep.
    meanSide(mean_bpos, mean_spos, side);
}


/* points: store n srep's * side points.
 * Loop each of the n srep. Sum the correspondence xyz of each sreps, get a matrix storing 3*pointsNum double.
*/
MatrixType shiftedsrep::meanSidePointPair(vector<MatrixType> points){
    int srepNum = points.size();
    int pointsNum = points[0].columns();
    MatrixType mean(3,pointsNum);
    mean.fill(0);

    // Loop each srep
    for(unsigned int k = 0; k< srepNum; k++){
        MatrixType X = points[k];

        //cout<<"------------------------the points used to compute mean is: "<<endl;
        for(unsigned int m = 0; m< pointsNum; m++){
            for(unsigned int dim=0; dim<3; dim++){
                mean[dim][m] += X[dim][m];
            }
        }
    }

    // Get the mean xyz.
    for(unsigned int m = 0; m < pointsNum; m++){
        for(unsigned int dim = 0; dim < 3; dim++){
            mean[dim][m] /= srepNum;
        }
    }

    return mean;
}



/* mean_bpos: storing all the mean boundary points
 * mean_spos: storing all the mean skeletal points. Correspondence to each other.
 * This function first copy a quadFig as a template, then reset the spoke info that will be the new quadFig.
*/
void shiftedsrep::meanSide(MatrixType mean_bpos, MatrixType mean_spos, int side){
    // Copy a arbitray quadFig as a template
    M3DQuadFigure* tempFig = this->quadFigList[0];

    MatrixType spokeVector = mean_bpos - mean_spos;

    // Set spoke tail (hub position)
    switch (side){
    case 0 : // up or down
    case 1 :
    {
        // Get spoke info
        M3DQuadPrimitive* prim;
        // Loop each atoms
        int atomIndex = 0; // points index
        for(unsigned u = 0; u<rowNum; u++){
            for(unsigned v = 0; v<colNum; v++){
                prim = dynamic_cast<M3DQuadPrimitive*>(tempFig->getPrimitivePtr(u,v));

                // Set hub position
                prim->setX(mean_spos[0][atomIndex], mean_spos[1][atomIndex],mean_spos[2][atomIndex]);

                // Set spoke radius
                double spokeRadiu=0;
                for(unsigned int dim=0; dim<3; dim++){
                    spokeRadiu += spokeVector[dim][atomIndex] * spokeVector[dim][atomIndex];
                }

                // Set spoke direction
                double sqrtRadiu = sqrt(spokeRadiu);
                Vector3D spokeDirection;
                spokeDirection.setX(spokeVector[0][atomIndex] / sqrtRadiu);
                spokeDirection.setY(spokeVector[1][atomIndex] / sqrtRadiu);
                spokeDirection.setZ(spokeVector[2][atomIndex] / sqrtRadiu);

                if(side == 1){
                    prim->setR1(sqrtRadiu);
                    prim->setU1(spokeDirection);                    
                }
                else{
                    prim->setR0(sqrtRadiu);
                    prim->setU0(spokeDirection);
                }

                atomIndex++;
            }
        }

        if(side == 1)
            this->downMeanQuadFig = tempFig;
        else
            this->upMeanQuadFig = tempFig;
    }
        break;
    case 2 : // crest
    {
        // Get spoke info of the moved crest
        int atomIndex = 0;
        M3DQuadEndPrimitive* primEnd;
        for(unsigned int u=0; u< rowNum; u++){
            for(unsigned int v=0; v< colNum; v++){

                if(u==0 || u==(rowNum-1) || v==0 || v==(colNum-1)){

                    primEnd = dynamic_cast<M3DQuadEndPrimitive*>(tempFig->getPrimitivePtr(u,v));

                    if(!primEnd){ //If point to nothing
                        cout<<"Message from: shiftedsrep::getPoints:: the primEnd pointer point is NULL!! EXIT_FAILURE!!"<<endl;
                        exit(1);
                    }

                    // Set hub position
                    primEnd->setX(mean_spos[0][atomIndex], mean_spos[1][atomIndex],mean_spos[2][atomIndex]);

                    // Set spoke radius
                    double spokeRadiu=0;
                    for(unsigned int dim=0; dim<3; dim++){
                        spokeRadiu += spokeVector[dim][atomIndex] * spokeVector[dim][atomIndex];
                    }

                    // Set spoke direction
                    double sqrtRadiu = sqrt(spokeRadiu);
                    Vector3D spokeDirection;
                    spokeDirection.setX(spokeVector[0][atomIndex] / sqrtRadiu);
                    spokeDirection.setY(spokeVector[1][atomIndex] / sqrtRadiu);
                    spokeDirection.setZ(spokeVector[2][atomIndex] / sqrtRadiu);

                    primEnd->setREnd(sqrtRadiu);
                    primEnd->setUEnd(spokeDirection);

                    atomIndex++;
                }
            }
        }

        this->crestMeanQuadFig = tempFig;

    }
        break;
    }
}



void shiftedsrep::readTrainingSreps(){
    cout<<"-------------Input : "<<this->srepNames.size()<<"sreps."<<endl;
    for(int q = 0; q < this->srepNames.size(); q++){

        string filename = this->rootDir + string("/models/compute_mean/") + this->srepNames[q];

        M3DQuadFigure* qFig = this->tls.GetQuadFigure(filename.c_str());
        this->quadFigList.push_back(qFig);
    }
}



/* Get surface point and save to the reference passing para surfacePoints
*/
void shiftedsrep::getPoints(M3DQuadFigure* quadFig, int side, vector<double> vars, int varNums, vector<MatrixType> &boundaryPos,
                            vector<MatrixType> &skeletalPos){
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
        int atomNum = this->crestAtomNum + this->interiorAtomNum;
        MatrixType bp(3, atomNum);
        MatrixType sp(3, atomNum);

        slidestandardspokes mStandardSpoke;
        // Shift spokes by given vars.
        mStandardSpoke.updateSpokeInfo(tempQuadFig, coeff, side);

        // Get spoke info
        M3DQuadPrimitive* prim;
        // Loop the interior atoms
        int atomIndex = 0;
        for(unsigned u = 0; u<rowNum; u++){
            for(unsigned v = 0; v<colNum; v++){
                prim = dynamic_cast<M3DQuadPrimitive*>(tempQuadFig->getPrimitivePtr(u,v));

                // Skeletal points.
                Vector3D point = prim->getX();
                sp[0][atomIndex] = point.getX();
                sp[1][atomIndex] = point.getY();
                sp[2][atomIndex] = point.getZ();

                // Boundary points.
                if(side==1){
                    point = prim->getX() + prim->getR1()*prim->getU1();
                }
                else{
                    point = prim->getX() + prim->getR0()*prim->getU0();
                }
                bp[0][atomIndex] = point.getX();
                bp[1][atomIndex] = point.getY();
                bp[2][atomIndex] = point.getZ();

                atomIndex++;
            }
        }

        boundaryPos.push_back(bp);
        skeletalPos.push_back(sp);
    }
        break;
    case 2 : // crest
    {
        MatrixType bp(3, this->crestAtomNum);
        MatrixType sp(3, this->crestAtomNum);

        slidecrestspokes mCrestSpoke;
        // Shift crest spokes
        visualizecrest vcrest;
        vtkSmartPointer<vtkSRep> srepFig = vcrest.getSrepFig(tempQuadFig);
        mCrestSpoke.updateCrestSpokeInfo(srepFig, tempQuadFig, coeff); //Move base on original srep.

        // Get spoke info of the moved crest
        int atomIndex = 0;
        M3DQuadEndPrimitive* primEnd;
        for(unsigned int u=0; u< rowNum; u++){
            for(unsigned int v=0; v< colNum; v++){

                if(u==0 || u==(rowNum-1) || v==0 || v==(colNum-1)){

                    primEnd = dynamic_cast<M3DQuadEndPrimitive*>(tempQuadFig->getPrimitivePtr(u,v));

                    if(!primEnd){ //If point to nothing
                        cout<<"Message from: shiftedsrep::getPoints:: the primEnd pointer point is NULL!! EXIT_FAILURE!!"<<endl;
                        exit(1);
                    }

                    // Boundary points.
                    Vector3D point = primEnd->getX() + primEnd->getREnd()*primEnd->getUEnd();
                    bp[0][atomIndex] = point.getX();
                    bp[1][atomIndex] = point.getY();
                    bp[2][atomIndex] = point.getZ();

                    // Skeletal points.
                    point = primEnd->getX();
                    sp[0][atomIndex] = point.getX();
                    sp[1][atomIndex] = point.getY();
                    sp[2][atomIndex] = point.getZ();

                    atomIndex++;
                }
            }
        }
        boundaryPos.push_back(bp);
        skeletalPos.push_back(sp);
    }
        break;
    }

    delete coeff;
}



void shiftedsrep::drawMeanSrep(){

}



/* Draw correspondence spokes. Use point index as index search in the color map. */
void shiftedsrep::drawSpokes(bool showup, bool showdown, bool showcrest, int interpolationLevel, vtkSmartPointer<vtkRenderer> renderer){

    visualizesrepsurface vss(this->rX, this->rY, this->rZ, this->opacity);
    int quadNum = (rowNum-1)*(colNum-1);//quadNum means the quad numbers on skeletal sheet.
    int step = pow((double)2, (double)interpolationLevel);
    int subQuadPointsNum = (step+1)*(step+1);

    // Get interpolated u & v coordinate, storing in interpolatedU & interpolatedV
    VectorQuadPoint interpolatedU;
    VectorQuadPoint interpolatedV;
    vss.getInterpolateUVCoordinate(interpolatedU, interpolatedV, rowNum , colNum, step);

    if(showup){
        vtkSmartPointer< vtkPoints > hubpos_medialsheet = vtkSmartPointer< vtkPoints >::New();

        // Get all the interpolated boundary points
        M3DQuadInterpolater *tpm = new M3DQuadInterpolater(this->upMeanQuadFig);
        Vector3D point_s, point_b;
        // Get up spokes
        vtkSmartPointer<vtkCellArray> cellarraypointsline_up = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer< vtkPoints > hubpos_up = vtkSmartPointer< vtkPoints >::New();
        for(unsigned int i = 0; i < quadNum; i++){
            for(unsigned int j = 0; j < subQuadPointsNum; j++){
                point_b = tpm->interpolateQuadSpoke(this->upMeanQuadFig,interpolatedU[i][j],interpolatedV[i][j],0)->getB();
                point_s = tpm->interpolateQuadSpoke(this->upMeanQuadFig,interpolatedU[i][j],interpolatedV[i][j],0)->getX();

                hubpos_medialsheet->InsertNextPoint(point_s.getX(),point_s.getY(),point_s.getZ());

                vtkIdType id0 = hubpos_up->InsertNextPoint(point_s.getX(),point_s.getY(),point_s.getZ());
                vtkIdType id1 = hubpos_up->InsertNextPoint(point_b.getX(),point_b.getY(),point_b.getZ());

                vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

                medialsheetline->GetPointIds()->SetId(0, id0);
                medialsheetline->GetPointIds()->SetId(1, id1);

                cellarraypointsline_up->InsertNextCell(medialsheetline);
            }
        }

        double quadColor[3] = {1,0.1,1};//gray //{0,0.8,0};//Green.
        vss.addSpokesToRender_2(hubpos_up, cellarraypointsline_up, renderer, quadColor);

        // Draw the up medial sheet.
        vss.drawQuads(hubpos_medialsheet, quadNum, subQuadPointsNum, step, renderer);

        delete tpm;
    }

    if(showdown){
        vtkSmartPointer< vtkPoints > hubpos_medialsheet = vtkSmartPointer< vtkPoints >::New();

        // Get all the interpolated boundary points
        M3DQuadInterpolater *tpm = new M3DQuadInterpolater(this->upMeanQuadFig);
        Vector3D point_s, point_b;

        // Get down spokes
        vtkSmartPointer<vtkCellArray> cellarraypointsline_down = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer< vtkPoints > hubpos_down = vtkSmartPointer< vtkPoints >::New();
        for(unsigned int i = 0; i < quadNum; i++){
            for(unsigned int j = 0; j < subQuadPointsNum; j++){
                point_b = tpm->interpolateQuadSpoke(this->downMeanQuadFig,interpolatedU[i][j],interpolatedV[i][j],1)->getB();
                point_s = tpm->interpolateQuadSpoke(this->downMeanQuadFig,interpolatedU[i][j],interpolatedV[i][j],1)->getX();

                hubpos_medialsheet->InsertNextPoint(point_s.getX(),point_s.getY(),point_s.getZ());

                vtkIdType id0 = hubpos_down->InsertNextPoint(point_s.getX(),point_s.getY(),point_s.getZ());
                vtkIdType id1 = hubpos_down->InsertNextPoint(point_b.getX(),point_b.getY(),point_b.getZ());

                vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

                medialsheetline->GetPointIds()->SetId(0, id0);
                medialsheetline->GetPointIds()->SetId(1, id1);

                cellarraypointsline_down->InsertNextCell(medialsheetline);
            }
        }

        double quadColor[3] = {0.3,0.3,0.3};//gray //{0,0.8,0};//Green.
        vss.addSpokesToRender_2(hubpos_down, cellarraypointsline_down, renderer, quadColor);

        // Draw the up medial sheet.
        vss.drawQuads(hubpos_medialsheet, quadNum, subQuadPointsNum, step, renderer);

        delete tpm;
    }

    if(showcrest){
        // Draw crest spokes
        double quadColor[3] = {0.3,0.3,0.3};//gray //{0,0.8,0};//Green.
        visualizecrest visualObject_crest(this->crestMeanQuadFig, interpolationLevel, quadColor, renderer,this->rX, this->rY, this->rZ, this->opacity);
        visualObject_crest.drawCrestSpokes();
    }
}



vtkSmartPointer< vtkPoints > shiftedsrep::getMeanSurfacePoints(int interpolationLevel){
    visualizesrepsurface vss(this->rX, this->rY, this->rZ, this->opacity);

    vtkSmartPointer< vtkPoints > surfacePoints = vtkSmartPointer< vtkPoints >::New();

    // Get up and down surface points
    vss.getUpOrDownSurfacePoints(this->upMeanQuadFig, surfacePoints, interpolationLevel, 0);
    vss.getUpOrDownSurfacePoints(this->downMeanQuadFig, surfacePoints, interpolationLevel, 1);

    // Get crest surface points
    vss.getCrestSurfacePoints(this->crestMeanQuadFig, surfacePoints, interpolationLevel);

    return surfacePoints;

}


void shiftedsrep::drawPointsSetSurfaceOfTheMean(int interpolationLevel, vtkSmartPointer<vtkRenderer> renderer){
    visualizesrepsurface vss(this->rX, this->rY, this->rZ, this->opacity);

    vtkSmartPointer< vtkPoints > surfacePoints = getMeanSurfacePoints(interpolationLevel);

    vss.drawPointsSetSurface_zcolor(surfacePoints, renderer);
}


/* Return the point number of the surface.
*/
/*int shiftedsrep::getPointsNum(int interpolationLevel){
    vtkSmartPointer< vtkPoints > surfacePoints = vtkSmartPointer< vtkPoints >::New();

    visualizesrepsurface vss(this->rX, this->rY, this->rZ);

    vss.getUpOrDownSurfacePoints(this->quadFigList[0], surfacePoints, interpolationLevel, 0);
    int pointNum = surfacePoints->GetNumberOfPoints();
    cout<<"------------up-----------------pointNum is: "<<pointNum<<endl;
    vss.getUpOrDownSurfacePoints(this->quadFigList[0], surfacePoints, interpolationLevel, 1);
    pointNum = surfacePoints->GetNumberOfPoints();
    cout<<"--------------+down---------------pointNum is: "<<pointNum<<endl;
    vss.getCrestSurfacePoints(this->quadFigList[0], surfacePoints, interpolationLevel);
    pointNum = surfacePoints->GetNumberOfPoints();
    cout<<"------------------+ crest-----------pointNum is: "<<pointNum<<endl;

    return pointNum;
}*/


/* Return the point number of the surface. These point are not repeat. For the crest point here discard the up and down tip, only
 * keep the interpolated point between them, because the up or down tip has already counted in the standard side.
 * level = 0, there are 78 surface point;
 * level = 1, there are 306 surface point;
 * level = 2, there are 1218 surface point;
*/
int shiftedsrep::calPointsNum(int interpolationLevel){

    int step = pow((double)2, (double)interpolationLevel);

    // Standard side
    int standQuads = (rowNum-1)*(colNum-1);
    int posEachStandQuad = step * step;
    int posStandSide = standQuads * posEachStandQuad + (colNum-1)*step + (rowNum-1)*step + 1;

    int crestQuads = rowNum*2 + (colNum-2)*2;
    int posEachCrestQuad = (step-1)*step;
    int posCrestSide = crestQuads * posEachCrestQuad;

    return posStandSide * 2 + posCrestSide;
}


/* Return the point number of the up or down side's surface.
 * level = 0, there are 39 standard surface point each side;
 * level = 1, there are 125 standard surface point each side;
 * level = 2, there are 441 standard surface point each side;
*/
int shiftedsrep::getSSPN(int interpolationLevel){ // Get the standard side point number.

    int step = pow((double)2, (double)interpolationLevel);

    // Standard side
    int standQuads = (rowNum-1)*(colNum-1);
    int posEachStandQuad = step * step;
    int posStandSide = standQuads * posEachStandQuad + (colNum-1)*step + (rowNum-1)*step + 1;

    return posStandSide;
}



/* Deformation the reconstructed mean surface along the first three eigenmodes.
 * The training sample is stored into a 3*n by N matrix, where n is point number, N is srep number.
 * Here we have 3*n variables, N observations.
*/
void shiftedsrep::deformOriginalSurfaceMean( int interpolationLevel){
    // Step 1: Get surface points of the training set. 3*n by 30 matrix.
    visualizesrepsurface vss(this->rX, this->rY, this->rZ, this->opacity);
    int pointNum = calPointsNum(interpolationLevel);

    // Define a matrix(3*n by 30), holding all the surface points of the 30 sreps.
    MatrixType X(3*pointNum, this->srepNum); // Each row is feature, each column is a s-rep,
    X.fill(0);
    for(unsigned int i = 0; i < this->srepNum; i++){
        // Get surface point of each srep
        vtkSmartPointer< vtkPoints > surfacePos = vtkSmartPointer< vtkPoints >::New();
        vss.getUpOrDownSurfacePoints(this->quadFigList[i], surfacePos, interpolationLevel, 0);
        vss.getUpOrDownSurfacePoints(this->quadFigList[i], surfacePos, interpolationLevel, 1);
        vss.getCrestSurfacePoints(this->quadFigList[i], surfacePos, interpolationLevel);

        // Save points to matrix
        for(unsigned int j = 0; j < pointNum; j++){
            double p[3];
            surfacePos->GetPoint(j, p);
            for(unsigned int dim = 0; dim < 3; dim++){
                X(j*3+dim, i) = p[dim];
            }
        }
    }
    saveMatrix("/NIRAL/work/ltu/WorkSpace/May15/Pablo2/lib/vtksrep/Correspondence/deform_mean/X.txt", X);

    deformAlongEigMode(X, pointNum, this->srepNum);
}



void shiftedsrep::deformAlongEigMode(MatrixType X, int pointNum, int sampleNum){
    int n = 3*pointNum; // variables number (dimesion)
    // Step 2: Compute the mean vector of each variable. (sample mean, mean shape, n by 1 vector)
    MatrixType XBar(n, 1);
    XBar.fill(0);
    for(unsigned int r = 0; r < n; r++){
        for(unsigned int c = 0; c < sampleNum; c++){
            XBar[r][0] += X[r][c];
        }

        XBar[r][0] /= sampleNum; //XBar[r][0] is the rth variable's mean over all the s-reps.
    }
    saveMatrix("/NIRAL/work/ltu/WorkSpace/May15/Pablo2/lib/vtksrep/Correspondence/deform_mean/XBar.txt", XBar);

    // Step 3: Center each variable, still store in X.
    for(unsigned int r = 0; r < n; r++){
        for(unsigned int c = 0; c < sampleNum; c++){
            X[r][c] -= XBar[r][0];
        }
    }
    saveMatrix("/NIRAL/work/ltu/WorkSpace/May15/Pablo2/lib/vtksrep/Correspondence/deform_mean/centered_X.txt", X);

    // Step 4: Compute the covariance matrix of centered X.
    MatrixType covMatrix = (X * X.transpose()) / (sampleNum-1);
    //saveMatrix("/NIRAL/work/ltu/WorkSpace/May15/Pablo2/lib/vtksrep/Correspondence/deform_mean/corMatrix.txt", covMatrix);

    // Step 5: Find the eigenvectors(PC, column vector, each column is a PC) and eigenvalues(diagonal of V) of corMatrix.
    vnl_symmetric_eigensystem<double> eig(covMatrix);
    vnl_diag_matrix<double> eigenValues(n, n);
    eigenValues = eig.D; //eig.D is default by increasing order.

    /*cout<<"--------eigenValues is: "<<endl;
    for(unsigned r=0; r<eigenValues.rows();r++){
        cout<<eigenValues(r,r)<<endl;
    }*/
    MatrixType eigenVals(n,1);
    for(unsigned r=0; r<eigenValues.rows();r++){
        eigenVals[r][0] = eigenValues(r,r);
    }
    saveMatrix("/NIRAL/work/ltu/WorkSpace/May15/Pablo2/lib/vtksrep/Correspondence/deform_mean/eigenValues.txt", eigenVals);

    MatrixType eigenVector(n,n);
    eigenVector = eig.V;
    /*cout<<"--------eigenVector is: "<<endl;
    for(unsigned int r = 0; r < n; r++){
        for(unsigned int c = 0; c < n; c++){
            cout<<eigenVector[r][c]<<"  ";
        }
        cout<<endl;
    }*/
    //saveMatrix("/NIRAL/work/ltu/WorkSpace/May15/Pablo2/lib/vtksrep/Correspondence/deform_mean/eigenVector.txt", eigenVector);

    // Step 6: Compute the standard deviation, which is the square root of the eigenvalue.
    double sigma[3];
    // Choose the first three dominant eigenvalue and its eigenvector.
    for(unsigned int dim=0; dim<3; dim++){
        sigma[dim] = sqrt(eigenValues(n-1-dim,n-1-dim)); // The last three is the bigger three.
        //cout<<"-----eigenValues("<<dim<<","<<dim<<") is: "<<eigenValues(n-1-dim,n-1-dim)<<endl;
        //cout<<"-----sigma["<<dim<<"] is: "<<sigma[dim]<<endl;
    }
    // Save the three correspondence eigenvector. sigma is in descent order, so do the corrEigenVector.
    MatrixType corrEigenVector(n, 3);
    for(unsigned int i = 0; i<n; i++){
        for(unsigned int dim=0; dim<3; dim++){
            corrEigenVector[i][dim] = eigenVector[i][n-1-dim];
        }
    }
    //saveMatrix("/NIRAL/work/ltu/WorkSpace/May15/Pablo2/lib/vtksrep/Correspondence/deform_mean/corrEigenVector.txt", corrEigenVector);

    // Step 7: Calculate the deformation vector and deform the mean along the correspondence eigenvector.
    // Newshape = mean + eigenVector*weighted*sigma
    MatrixType newShape(n, 3); // Each column is a srep.
    int sdWeight = -3;//3;
    for(unsigned int mode=0; mode<3; mode++){
        double w = sdWeight * sigma[mode];
        for(unsigned int r = 0; r < n; r++){
            corrEigenVector[r][mode] *= w;
        }
    }    

    for(unsigned int r = 0; r < n; r++){
        for(unsigned int mode=0; mode<3; mode++){
            newShape[r][mode] = XBar[r][0] + corrEigenVector[r][mode];            
        }
    }

    // Save new shapes to vtk file.
    for(unsigned int mode=0; mode<3; mode++){
        char num[20];
        sprintf(num, "%d", mode);
        std::string outputFile = std::string("/NIRAL/work/ltu/WorkSpace/May15/Pablo2/lib/vtksrep/Correspondence/deform_mean/newShape_") +  num + ".vtk";
        saveToVTK(pointNum, newShape, outputFile, mode);
    }

    // Save the mean shape to vtk file
    string outputFile = string("/NIRAL/work/ltu/WorkSpace/May15/Pablo2/lib/vtksrep/Correspondence/deform_mean/meanShape.vtk");
    saveToVTK(pointNum, XBar, outputFile, 0);





}


/* Deformation the reconstructed shifted srep's mean surface along the first three eigenmodes.
 * The training sample is stored into a 3*n by N matrix, where n is point number, N is srep number.
 * Here we have 3*n variables, N observations.
 * If set up_str, down_str, crest_str all to 0, the mean is the same as deformOriginalSurfaceMean().
*/
void shiftedsrep::deformShiftedSurfaceMean(string up_str, string down_str, string crest_str, int interpolationLevel){

    int pointNum = calPointsNum(interpolationLevel);
    cout<<"-----------------------------pointNum is: "<<pointNum<<endl;

    MatrixType X(3*pointNum, this->srepNum); // 3*n by 30 matrix, each row is feature, each column is a s-rep.
    X.fill(0);

    // Get the boundary points on each side
    surfacePointsEachSide(up_str, 0, this->varStand, X, interpolationLevel);
    surfacePointsEachSide(down_str, 1, this->varStand, X, interpolationLevel);
    surfacePointsEachSide(crest_str, 2, this->crestAtomNum, X, interpolationLevel);

    saveMatrix("/NIRAL/work/ltu/WorkSpace/May15/Pablo2/lib/vtksrep/Correspondence/deform_mean/X.txt", X);

    deformAlongEigMode(X, pointNum, this->srepNum);
}





/* Get the points on the boundary of each side. up or down side has 39 points, crest has 28. total 106 points for each srep's boundary.
 * boundaryPoints: storing all the 30 sreps' side's boundary points.
 * sspn: up or down side points number.
*/
void shiftedsrep::surfacePointsEachSide(string vars_str, int side, int varNums, MatrixType &X, int interpolationLevel){
    stringstream ss(vars_str);

    // Move and update each side of each s-rep
    int srepNum = this->quadFigList.size();
    int sspn = getSSPN(interpolationLevel);

    vector<double> vars;

    for(unsigned int q = 0; q < srepNum; q++){ //Loop each srep
        // Get up side points
        string token;
        int i = 0;
        while (i < varNums) {
            ss >> token;
            vars.push_back(atof(token.c_str()));
            //cout<<"---------------------------var "<< vars[i]<<endl;
            i++;
        }
        //cout<<"---------------------------var size: "<< vars.size()<<endl;

        // Shift the srep and get the shifted surface points
        vtkSmartPointer< vtkPoints > surfacePos = vtkSmartPointer< vtkPoints >::New();
        getSurfacePoints(this->quadFigList[q], side, vars, varNums, surfacePos, interpolationLevel);

        // Store the points into matrix
        int pNum = surfacePos->GetNumberOfPoints();

        if(side==0){
            for(unsigned int j = 0; j < pNum; j++){
                double p[3];
                surfacePos->GetPoint(j, p);
                for(unsigned int dim = 0; dim<3; dim++){
                    X[j*3+dim][q] = p[dim];
                }
            }
        }
        else if(side==1){
            for(unsigned int j = 0; j < pNum; j++){
                double p[3];
                surfacePos->GetPoint(j, p);
                for(unsigned int dim = 0; dim<3; dim++){
                    X[sspn*3+j*3+dim][q] = p[dim];
                }
            }
        }
        else if(side==2){
            for(unsigned int j = 0; j < pNum; j++){
                double p[3];
                surfacePos->GetPoint(j, p);
                for(unsigned int dim = 0; dim<3; dim++){
                    X[sspn*3*2+j*3+dim][q] = p[dim];
                }
            }
        }
        else{
            cout<<"Message from shiftedsrep::surfacePointsEachSide:: Invaild side type!!"<<endl;
            cout<<"---"<<side<<endl;
            exit(1);
        }

        // Clear for next
        vars.clear();
    }
}



/* Get surface point and save to the reference passing para surfacePoints.
 * This function is diff from getPoints() by getting the surface point under specific interpolation level.
*/
void shiftedsrep::getSurfacePoints(M3DQuadFigure* quadFig, int side, vector<double> vars, int varNums, vtkSmartPointer< vtkPoints > surfacePos,
                            int interpolationLevel){
    // Convert coeffs from vector to pointer
    double *coeff = new double [varNums]; // storing up or down spokes variables.
    for(unsigned int i = 0; i < varNums; i++){
        coeff[i] = vars[i];
    }

    M3DQuadFigure* tempQuadFig = dynamic_cast<M3DQuadFigure*>(quadFig->clone()); //Move all base on original srep.

    visualizesrepsurface vss(this->rX, this->rY, this->rZ, this->opacity);

    switch (side){
    case 0 : // up or down
    case 1 :
    {
        slidestandardspokes mStandardSpoke;
        // Shift spokes by given vars.
        mStandardSpoke.updateSpokeInfo(tempQuadFig, coeff, side);

        // Get spoke info
        vss.getUpOrDownSurfacePoints(tempQuadFig, surfacePos, interpolationLevel, side);
    }
        break;
    case 2 : // crest
    {
        slidecrestspokes mCrestSpoke;
        // Shift crest spokes
        visualizecrest vcrest;
        vtkSmartPointer<vtkSRep> srepFig = vcrest.getSrepFig(tempQuadFig);
        mCrestSpoke.updateCrestSpokeInfo(srepFig, tempQuadFig, coeff); //Move base on original srep.

        // Get spoke info of the moved crest
        vss.getCrestSurfacePoints(tempQuadFig, surfacePos, interpolationLevel);
    }
        break;
    }

    delete coeff;
}



/* col: is the shape index, which shape to be saved.
 * Each column in mShape is a srep.
*/
void shiftedsrep::saveToVTK(int pointNum, MatrixType mShape, string outputFile, int mode){
    // Save each new shape points to vtk file
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New(); // same as vtkPolyData* polyData;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for(unsigned int i = 0; i < pointNum; i++){
        int r = i*3;
        points->InsertNextPoint(mShape[r][mode], mShape[r+1][mode], mShape[r+2][mode]);
    }

    polyData->SetPoints(points);
    vtkSmartPointer<vtkPolyDataWriter> polyDataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    polyDataWriter->SetInput(polyData);

    polyDataWriter->SetFileName(outputFile.c_str());
    polyDataWriter->Update();
}




/* Read the surface point saved in vtk file, and plot these points.
*/
void shiftedsrep::drawSurfacePoint(const char* filename, vtkSmartPointer<vtkRenderer> renderer){
    // Read in the source points set
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New(); // same as vtkPolyData* polyData_source;

    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(filename);
    reader->Update();
    polyData = reader->GetOutput();
    cout<<"-----------------filename is: "<<filename<<endl;

    vtkSmartPointer< vtkPoints > surfacePoints = vtkSmartPointer< vtkPoints >::New();

    // Read in the source points set
    for(unsigned int i = 0; i < polyData->GetNumberOfPoints(); i++){
        double p[3];
        polyData->GetPoint(i,p);

        surfacePoints->InsertNextPoint(p[0], p[1], p[2]);
    }

    // Draw the points (surface).
    visualizesrepsurface vss(this->rX, this->rY, this->rZ, this->opacity);
    vss.drawPointsSetSurface_zcolor(surfacePoints, renderer);
}




/* Save matrix to .txt file.*/
void shiftedsrep::saveMatrix(const char* filename, MatrixType matrix){

    std::ofstream fout;
    fout.open(filename);

    if(fout)  {
        for(int i =0; i< matrix.rows();i++){
            for(int j=0; j<matrix.columns(); j++){
                fout<< matrix[i][j]<<" ";
            }
            fout<<endl;
        }

        cout<<"Successfully saved matrix to: "<<filename<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();
}

















