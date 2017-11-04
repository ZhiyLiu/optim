#include "slidecrestspokes.h"

slidecrestspokes::slidecrestspokes()
{
}





slidecrestspokes::slidecrestspokes(int interpolationLevel){
    this->interpolationLevel = interpolationLevel;
}




/* Move crest spokes by given vars, compute the moved crest regulartiy entropy and move srep's p r u values.
 * vars: a array holing all the changing deltau or deltav for this srep. (Each srep crest has varCrest variables.)
*/
double slidecrestspokes::moveCrestSpokes(M3DQuadFigure *quadFig, vtkSmartPointer<vtkSRep> srepfig, double *vars, DoubleVec subTVs){

    // Shift spokes. Update quadFig' crest spoke information using the moved spokes' p r u values.
    updateCrestSpokeInfo(srepfig, quadFig, vars);

    // Compute the shifted srep's regularity entropy using the changed quadFig and new computed srepfig.
    double entropy = computeRegularityEntropyOfCrest(quadFig, subTVs);

    return entropy;
}



/* GetInterpolatedSpoke function return a spoke direction, this direction is the vector of "boundary point - skeletal point",
 * its coordinate is smaller, such as: <0.00479999, -0.0391918, -0.0290593>. In pable and CPNS, we normalize each spoke to be a unit vector,
 * that is: divid x,y,z by sqrt(x^2+y^2+Z^2) seperately, so the size of this vector is 1, breifly name as be normalized to 1.
 * Compute the shifted spoke and update to its original index's primitive.
 * Set variable (dv or du) to primitive. Move this crest spoke[u][v] forward (var>0) or backward (var<0).
 * If var <0, we use its previous cellid to got the new spoke position.
 *
 * cellIndex: the index of the spoke hub along the crest curve.
 * var: the movement distance for this spoke.
*/
void slidecrestspokes::shiftCrestSpokeAndSetToPrimitive(vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> shiftCrestSpokes,
                                                        M3DQuadFigure *quadFig, double u, double v, int cellIndex, double var){
//    int rowNums = quadFig->getRowCount();
//    int colNums = quadFig->getColumnCount();
//    int crestSpokeNum = 2*colNums + (rowNums -2)*2;

    double slide_t = var; //the shift distance.
    int cellid = cellIndex; // crest atom index.

    if(var < 0){ //move -0.4, equals move its previous cell to 1-0.4=0.6.
        cellid = cellIndex - 1;
        slide_t = 1 + var;

        // The (0,0) not move, no need this.
//        // The first cellid, if backward, will use the last cellid, that's (crestSpokeNum - 1).
//        if(u==0 && v==0){
//            cellid = crestSpokeNum - 1;
//        }
    }

    if(slide_t<0.0000001) //sometimes 0 is represent in tiny double; if cellid=19 and slide_t = 1.32955e-316, the spokedir will be NAN, why??
        slide_t = 0;


    // Get the interpolated skeletal position at this point
    vtkSRep::VNLType spokeHub = shiftCrestSpokes->GetInterpolatedPoint(cellid, slide_t);
    // Get interpolated spoke direction at this point, without normalize.
    vtkSRep::VNLType upSpokeDir = shiftCrestSpokes->GetInterpolatedSpoke(cellid, slide_t, 0); //up spoke
    vtkSRep::VNLType crestSpokeDir = shiftCrestSpokes->GetInterpolatedSpoke(cellid, slide_t, 0.5); //crest spoke
    vtkSRep::VNLType downSpokeDir = shiftCrestSpokes->GetInterpolatedSpoke(cellid, slide_t, 1); //down spoke

    //set spoke hub position to new value.
    quadFig->getPrimitivePtr(u,v)->setX(spokeHub[0], spokeHub[1], spokeHub[2]);
    //cout<<"---- cellid is: "<<cellid<<", slide_t is: "<<slide_t<<" ---------------------------spokeHub is: "<<spokeHub<<endl;

    Vector3D unitSpokeDir;
    double spokeLength;

    //set up spoke length & dirction.
    spokeLength = sqrt(upSpokeDir[0]*upSpokeDir[0] + upSpokeDir[1]*upSpokeDir[1] + upSpokeDir[2]*upSpokeDir[2]);
    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setR0(spokeLength);
    //cout<<"---------------------------------the up spokeLength is: "<<spokeLength<<endl;
    // nomalize the spoke direction to be unit vector.
    unitSpokeDir.set(upSpokeDir[0]/spokeLength,upSpokeDir[1]/spokeLength,upSpokeDir[2]/spokeLength);
    //cout<<"---------------------------------the unit up spokeDir is: "<<unitSpokeDir<<endl;
    quadFig->getPrimitivePtr(u,v)->setU0(unitSpokeDir);


    //set down spoke length & dirction.
    spokeLength = sqrt(downSpokeDir[0]*downSpokeDir[0] + downSpokeDir[1]*downSpokeDir[1] + downSpokeDir[2]*downSpokeDir[2]);
    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setR1(spokeLength);
    //cout<<"---------------------------------the down spokeLength is: "<<spokeLength<<endl;
    // nomalize the spoke direction to be unit vector.
    unitSpokeDir.set(downSpokeDir[0]/spokeLength,downSpokeDir[1]/spokeLength,downSpokeDir[2]/spokeLength);
    //cout<<"---------------------------------the unit down spokeDir is: "<<unitSpokeDir<<endl;
    quadFig->getPrimitivePtr(u,v)->setU1(unitSpokeDir);


    //set crest spoke length & dirction.
    spokeLength = sqrt(crestSpokeDir[0]*crestSpokeDir[0] + crestSpokeDir[1]*crestSpokeDir[1] + crestSpokeDir[2]*crestSpokeDir[2]);
    dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(u,v))->setREnd(spokeLength);
    //cout<<"---------------------------------the crest spokeLength is: "<<spokeLength<<endl;
    // nomalize the spoke direction to be unit vector.
    unitSpokeDir.set(crestSpokeDir[0]/spokeLength,crestSpokeDir[1]/spokeLength,crestSpokeDir[2]/spokeLength);
    //cout<<"---------------------------------the unit crest spokeDir is: "<<unitSpokeDir<<endl;
    dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(u,v))->setUEnd(unitSpokeDir);

    // if use the following functions, is slow than the above code...
//    //set up spoke length & dirction.
//    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setR0(calculateSpokeLength(upSpokeDir));
//    quadFig->getPrimitivePtr(u,v)->setU0(calculateUnitSpokeDir(upSpokeDir));


//    //set down spoke length & dirction.
//    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setR1(calculateSpokeLength(downSpokeDir));
//    quadFig->getPrimitivePtr(u,v)->setU1(calculateUnitSpokeDir(downSpokeDir));


//    //set crest spoke length & dirction.
//    dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(u,v))->setREnd(calculateSpokeLength(crestSpokeDir));
//    dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(u,v))->setUEnd(calculateUnitSpokeDir(crestSpokeDir));


//    // Set deltau & deltav. (Not need for compute, but for look and check.)
//    if(u==0 || u ==rowNums-1) { //Only v change
//        dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setDeltaU2(0);
//        dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setDeltaV2(var);
//    }
//    if(v==0 || v ==colNums-1) { //Only u change
//        dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setDeltaU2(var);
//        dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setDeltaV2(0);
//    }
//    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setDeltaU1(0);
//    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setDeltaV1(0);
//    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setDeltaU0(0);
//    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setDeltaV0(0);
}


double slidecrestspokes::calculateSpokeLength(vtkSRep::VNLType newSpoke){
    double squareDis = newSpoke[0]*newSpoke[0] + newSpoke[1]*newSpoke[1] + newSpoke[2]*newSpoke[2];

    return sqrt(squareDis);
}


// nomalize the spoke direction to be unit vector.
Vector3D slidecrestspokes::calculateUnitSpokeDir(vtkSRep::VNLType newSpoke){
    double spokeLength = calculateSpokeLength(newSpoke);

    Vector3D unitSpokeDir;

    unitSpokeDir.set(newSpoke[0]/spokeLength,newSpoke[1]/spokeLength,newSpoke[2]/spokeLength);

    return unitSpokeDir;
}



/* Crest spoke can only move along the crest line. each spoke only v or u move, so each spoke only has one variable.
 * The corner four atom not move.
 * varArraySrep: size is 24.
 * The four corner spokes has variables always set to 0, fixed: (0, 0), (0, colNums-1), (rowNums-1, colNums-1), (rowNums-1, 0).
*/
void slidecrestspokes::updateCrestSpokeInfo(vtkSmartPointer<vtkSRep> srepfig, M3DQuadFigure *quadFig, double * varArraySrep){

    // New a crest interplator....
    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> shiftCrestSpokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
    shiftCrestSpokes->SetInput(srepfig);
    // Here can set to any level. doesn't matter. We only need to get the moveed crest atom's spokes info. Its position index is fixed.
    shiftCrestSpokes->SetInterpolationLevel(0); // test: set to 0, 1, 2, 3. output GetInterpolatedSpoke(cellid, slide_t, 0.5), same!!!
    shiftCrestSpokes->Update();

    int rowNums = quadFig->getRowCount();
    int colNums = quadFig->getColumnCount();

    // Output all the variables for all the sreps in this movement.
    /*cout<<"---------------------varArraySrep for this srep is : "<<endl;
    for(unsigned j = 0; j < crestSpokeNum; j++){
        cout<<varArraySrep[j]<<"  ";
    }
    cout<<endl;*/

    // Get the crest atoms one by one, along the fold curve in counter clockwise. Starting from (0,0).
    unsigned int u = 0;
    unsigned int v = 0;

    int varIndex = 0; //count the varArraySrep, same as the crest atom id

    // primitive (0,0), (0, colNums-1) are fixed, no need to be updated.

    // u=0, v=1,2,...,colNums-1. Bottom, only v change.
    for(v = 1; v < colNums-1; v++){

        // Get the spoke information at the new position. And update primitive[u][v] to this new values.
        shiftCrestSpokeAndSetToPrimitive(shiftCrestSpokes, quadFig, u, v, varIndex+1, varArraySrep[varIndex]);

        varIndex++;
    }

    // primitive (rowNums-1,colNums - 1) is fixed, no need to be updated.
    // u=1,2,...,rowNums-1, v=colNums - 1, only u change.
    v = colNums - 1;
    for(u = 1; u < rowNums-1; u++){

        // Get the spoke information at the new position. And update primitive[u][v] to this new values.
        shiftCrestSpokeAndSetToPrimitive(shiftCrestSpokes, quadFig, u, v, varIndex+2, varArraySrep[varIndex]);

        varIndex++;
    }

    // u=rowNums-1, v=1,2,...,colNums-2
    u = rowNums-1;
    for(v = colNums-2; v > 0; v--){

        // Get the spoke information at the new position. And update primitive[u][v] to this new values.
        shiftCrestSpokeAndSetToPrimitive(shiftCrestSpokes, quadFig, u, v, varIndex+3, varArraySrep[varIndex]);

        varIndex++;
    }

    // primitive (rowNums-1,0) is fixed, no need to be updated.
    // u=rowNums-1,...,2,1, v=0
    v = 0;
    for(u = rowNums-2; u > 0; u--){

        // Get the spoke information at the new position. And update primitive[u][v] to this new values.
        shiftCrestSpokeAndSetToPrimitive(shiftCrestSpokes, quadFig, u, v, varIndex+4, varArraySrep[varIndex]);

        varIndex++;
    }

    //cout<<"-----varIndex is (shoulde be 24): "<<varIndex<<endl;

//    // Set all the up and down spokes deltaU & deltaV to 0. No need, Just for m3d file easy look. If not set to 0, some values very strange
//    // 3.8103891997142292e-316; 2.7e21... not effect the result, because don't use that one, but why these values?
//    // Loop the interior atoms, which only have up and down spokes.
//    for(unsigned int i =1; i<rowNums-1; i++){
//        for(unsigned int j =1; j<colNums-1; j++){
//            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j))->setDeltaU1(0);
//            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j))->setDeltaV1(0);
//            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j))->setDeltaU0(0);
//            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j))->setDeltaV0(0);
//            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j))->setDeltaU2(0);
//            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j))->setDeltaV2(0);
//        }
//    }
}







/* Compute the crest spokes regularity entropy.
 * Regularity entropy was divided into three kinds of entropies:
 * angle entropy, horizonal edge length entropy, vertical edge length entropy.
 * The three type of entropies correspondence to 3 feature matrix:
 * angle feature matrix is 3-by-24. (24 is the quad number).
 * horizonal edge length feature matrix is 1-by-24.
 * vertical edge length feature matrix is 2-by-24.
*/
double slidecrestspokes::computeRegularityEntropyOfCrest(M3DQuadFigure *quadFig, DoubleVec subTVs){

    double regEntropy = 0;
    int crestAtomNums = (quadFig->getRowCount() + quadFig->getColumnCount()) * 2 - 4;

    // Horizontal edge feature & vertical edge feature
    MatrixType verEdgeFeatures(2, crestAtomNums);
    MatrixType horEdgeFeatures(1, crestAtomNums);

    vtkSmartPointer< vtkPoints > points_s = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkPoints > points_b = vtkSmartPointer< vtkPoints >::New();

    calcrestregularityfeatures regFeature(quadFig, subTVs, this->interpolationLevel, crestAtomNums);

    regFeature.getSubQuadsPositionOfCrestRegion(points_s, points_b);

    regularityentropy regE;

    //Calculate the edges
    regFeature.calculateCrestEdges_2(points_s, points_b, verEdgeFeatures, horEdgeFeatures);

    // Compute the horizonal edge length entropy
    regEntropy += regE.calculateEntropy(horEdgeFeatures, 0.01);

    // Compute the vertical edge length entropy
    regEntropy += regE.calculateEntropy(verEdgeFeatures, 0.01);


    MatrixType angleFeatures(3, crestAtomNums);

    // Angle feature
    regFeature.calculateCrestAngles_2(points_b, angleFeatures);


//    for(unsigned dim = 0; dim<angleFeatures.rows(); dim++){
//        for(unsigned int i=0; i<angleFeatures.columns(); i++){
//            cout<<angleFeatures[dim][i]<<"  ";
//        }
//        cout<<endl;
//    }
    // Compute the angle entropy
    regEntropy += regE.calculateEntropy(angleFeatures, 0.01);

    return regEntropy;
}






