/* Sliding the standard spokes (up spokes or down spokes).
 * The crest atom can only move along the skeletal sheet (only in u or v direction).
 * The interior atom can move anywhere with both u and v.
 * Currently the four corner atoms were fixed, no sliding.
 *
 * Liyun Tu
 * Apr 13, 2014
*/

#include "slidestandardspokes.h"



slidestandardspokes::slidestandardspokes()
{
}

slidestandardspokes::slidestandardspokes(int side, int interpolationLevel){
    this->side = side;
    this->interpolationLevel = interpolationLevel;
}


/* This function don't save any intermedia sreps, the moved srep's p r u values is stored in a matrix and passed to matlab for
 * geometry entropy compute.
 * Input a srep and its variables, output its moved entropy.
 * side; up(0), down(1).
 * vars: a array holing all the changing deltau and deltav for this srep. (Each srep has varStand*2 + varCrest variables.)
*/
double slidestandardspokes::moveStandardSpokes(M3DQuadFigure* quadFig, double *vars, DoubleVec subUVs){
    // Shift spokes. Update this->curQuadFig' up or down spoke information using the moved spokes' p r u values.
    updateSpokeInfo(quadFig, vars, this->side);

    // Compute the shifted srep's regularity entropy
    double entropy = computeRegularityEntropy(quadFig, subUVs); //compute the entropy directly, not passing matrix to matlab anymore.

    return entropy; // return this srep's up or down entropy
}




void slidestandardspokes::setMovedSpokeInfoToPrimitive(M3DQuadFigure *quadFig, double newU, double newV, double primitiveIndexU,
                                                       double primitiveIndexV){
    //cout<<"Set moved spoke info for primitive["<<primitiveIndexU<<","<<primitiveIndexV<<"]"<<endl;

    //get the correspondence delta u , v of each of the four points.
    M3DQuadInterpolater standsideintp;

    //point, a 3D vector, which store the x, y, z coordinate of this point.
    M3DSpoke newspoke = standsideintp.interpolateQuadSpoke2(quadFig, newU, newV, side);

    //interpolate and get new information for moved spokes.
    //cout<<"------------interpolating a spoke at ["<<newU<<", "<<newV<<"]"<<endl;
    //cout<<"------------The interpolated spoke tail is: "<<newspoke->getX().getX()<<"  "<< newspoke->getX().getY()<<"  "<<newspoke->getX().getZ()<<endl;

    //set spoke hub position to new value.
    quadFig->getPrimitivePtr(primitiveIndexU,primitiveIndexV)->setX(newspoke.getX().getX(), newspoke.getX().getY(), newspoke.getX().getZ());

    if(side==0){
        //set up spoke dirction.
        quadFig->getPrimitivePtr(primitiveIndexU,primitiveIndexV)->setU0(newspoke.getU());

        //set up spoke length.
        dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(primitiveIndexU,primitiveIndexV))->setR0(newspoke.getR());
    }
    else{
        //set down spoke dirction.
        quadFig->getPrimitivePtr(primitiveIndexU,primitiveIndexV)->setU1(newspoke.getU());

        //set down spoke length.
        dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(primitiveIndexU,primitiveIndexV))->setR1(newspoke.getR());
    }
}




/* Update the 46 spokes' deltau, deltav and its corresponding primitive's r, u, x, y, z coordinates.
 * We move a spoke with a step deltau or deltav, get a new spoke information through interpolateQuadSpoke method, and save this new spoke information
 * into the original primitive for simple.
 * 4 Coner atoms do not move: spoke[0][0], spoke[rowNum-1][0], spoke[0][colNum-1], spoke[rowNum-1][colNum-1]. No need to set value.
 * inputPathName: is the original srep, whiche stored in /NIRAL/work/ltu/WorkSpace/Pablo2_Dec12/lib/vtksrep/Correspondence/models/testSet_Original
 * outputPathName: is the update srep, only used to calculate the entropy. These files are temporately, which changing in every loop of optimization.
 **/
void slidestandardspokes::updateSpokeInfo(M3DQuadFigure *quadFig, double * varArraySrep, int side){
    /*cout<<"---------the variables for spokes are: "<<endl;
        for(unsigned j =0; j< 30; j++){
            cout<<varArraySrep[j]<<"  ";
        }
        cout<<endl;*/
    int rowNum = quadFig->getRowCount();
    int colNum = quadFig->getColumnCount();

    int varArraySrepIndex = 0; //count the varArraySrep.

    double newU, newV;

    if(side==0) {
        //Crest spokes only u move: spoke[1][0],...,spoke[rowNum-2][0], spoke[1][colNum-1],...,spoke[rowNum-2][colNum-1].
        //***********************************************************************************************
        for(int i = 1; i< rowNum-1;i++){
            //For Primitive[1][0] to spoke[rowNum-2][0].
            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,0))->setDeltaU0(varArraySrep[varArraySrepIndex]);
            newU = i + varArraySrep[varArraySrepIndex];
            newV = 0;

            setMovedSpokeInfoToPrimitive(quadFig, newU, newV, i, 0);

            varArraySrepIndex++;

            //For spoke[1][colNum-1],...,spoke[rowNum-2][colNum-1].
            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,colNum-1))->setDeltaU0(varArraySrep[varArraySrepIndex]);
            newU = i + varArraySrep[varArraySrepIndex];
            newV = colNum-1;

            setMovedSpokeInfoToPrimitive(quadFig, newU, newV, i, colNum-1);

            varArraySrepIndex++;
        }

        //Crest spoke only v move: spoke[0][1], spoke[0][2], spoke[0][3], …, spoke[0][11]; spoke[2][1], spoke[2][2], spoke[2][3], …, spoke[2][11].
        //***********************************************************************************************
        for(int i=1; i< colNum-1;i++){
            //For Primitive[0][1] to spoke[0][colNum-2].
            newU = 0;
            newV = i + varArraySrep[varArraySrepIndex];

            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(0,i))->setDeltaV0(varArraySrep[varArraySrepIndex]);

            setMovedSpokeInfoToPrimitive(quadFig, newU, newV, 0, i);
            varArraySrepIndex++;


            //For Primitive[rowNum-1][1] to spoke[rowNum-1][colNum-2].
            newU = rowNum-1;
            newV = i + varArraySrep[varArraySrepIndex];

            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(rowNum-1,i))->setDeltaV0(varArraySrep[varArraySrepIndex]);

            setMovedSpokeInfoToPrimitive(quadFig, newU, newV, rowNum-1, i);
            varArraySrepIndex++;
        }

        //Interior spokes both u and v move: spoke[1][1], spoke[1][2], spoke[1][3], ..., spoke[1][11].
        //***********************************************************************************************
        //varArraySrep[24] - varArraySrep[34] store the deltaU value for Primitive[1][1] - Primitive[1][11].
        //varArraySrep[35] - varArraySrep[45] store the deltaV value for Primitive[1][1] - Primitive[1][11].
        for(int i = 1; i < rowNum-1;i++){
            for(int j=1; j<colNum-1;j++){
                newU = i + varArraySrep[varArraySrepIndex];
                dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j))->setDeltaU0(varArraySrep[varArraySrepIndex]);
                varArraySrepIndex++;

                newV = j + varArraySrep[varArraySrepIndex];
                dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j))->setDeltaV0(varArraySrep[varArraySrepIndex]);
                varArraySrepIndex++;

                setMovedSpokeInfoToPrimitive(quadFig, newU, newV, i, j);
            }
        }
    }
    else {
        //Crest spoke only u move: spoke[1][0],...,spoke[rowNum-2][0], spoke[1][colNum-1],...,spoke[rowNum-2][colNum-1].
        //***********************************************************************************************
        for(int i = 1; i< rowNum-1;i++){
            //For Primitive[1][0] to spoke[rowNum-2][0].
            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,0))->setDeltaU1(varArraySrep[varArraySrepIndex]);
            newU = i + varArraySrep[varArraySrepIndex];
            newV = 0;

            setMovedSpokeInfoToPrimitive(quadFig, newU, newV, i, 0);

            varArraySrepIndex++;


            //For spoke[1][colNum-1],...,spoke[rowNum-2][colNum-1].
            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,colNum-1))->setDeltaU1(varArraySrep[varArraySrepIndex]);
            newU = i + varArraySrep[varArraySrepIndex];
            newV = colNum-1;

            setMovedSpokeInfoToPrimitive(quadFig, newU, newV, i, colNum-1);

            varArraySrepIndex++;
        }

        //Crest spoke only v move: spoke[0][1], spoke[0][2], spoke[0][3], …, spoke[0][11]; spoke[2][1], spoke[2][2], spoke[2][3], …, spoke[2][11].
        //***********************************************************************************************
        for(int i=1; i< colNum-1;i++){
            //For Primitive[0][1] to spoke[0][colNum-2].
            newU = 0;
            newV = i + varArraySrep[varArraySrepIndex];

            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(0,i))->setDeltaV1(varArraySrep[varArraySrepIndex]);

            setMovedSpokeInfoToPrimitive(quadFig, newU, newV, 0, i);
            varArraySrepIndex++;


            //For Primitive[rowNum-1][1] to spoke[rowNum-1][colNum-2].
            newU = rowNum-1;
            newV = i + varArraySrep[varArraySrepIndex];

            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(newU,i))->setDeltaV1(varArraySrep[varArraySrepIndex]);

            setMovedSpokeInfoToPrimitive(quadFig, newU, newV, rowNum-1, i);
            varArraySrepIndex++;
        }


        //Interior spokes both u and v move: spoke[1][1], spoke[1][2], spoke[1][3], ..., spoke[1][11].
        //***********************************************************************************************
        for(int i = 1; i < rowNum-1;i++){
            for(int j=1; j<colNum-1;j++){
                newU = i + varArraySrep[varArraySrepIndex];
                dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j))->setDeltaU1(varArraySrep[varArraySrepIndex]);
                varArraySrepIndex++;

                newV = j + varArraySrep[varArraySrepIndex];
                dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j))->setDeltaV1(varArraySrep[varArraySrepIndex]);
                varArraySrepIndex++;

                setMovedSpokeInfoToPrimitive(quadFig, newU, newV, i, j);
            }
        }
    }    
}





/* Compute the regularity entropy
 * Regularity entropy was divided into three kinds of entropies:
 * angles entropy, horizonal edge length entropy, vertical edge length entropy.
 * The three type of entropies correspondence to 3 feature matrix:
 * angle feature matrix is 4-by-24. (24 is the quad number).
 * horizonal edge length feature matrix is 2-by-24.
 * vertical edge length feature matrix is 2-by-24.
*/
double slidestandardspokes::computeRegularityEntropy(M3DQuadFigure* quadFig, DoubleVec subUVs){
    double regEntropy = 0.0;

    int colNum = quadFig->getColumnCount();
    int rowNum = quadFig->getRowCount();

    int horEdgeNum = (rowNum-1) * colNum; //26
    int verEdgeNum = (colNum-1) * rowNum; //36

    MatrixType horEdgeFeatureMatrix(2, horEdgeNum);
    MatrixType verEdgeFeatureMatrix(2, verEdgeNum);

    calregularityfeatures regFeature(quadFig, this->interpolationLevel, this->side, subUVs);

    vtkSmartPointer< vtkPoints > hubPosition_b = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkPoints > hubPosition_s = vtkSmartPointer< vtkPoints >::New();
    regFeature.getSubQuadVertexesPosition(hubPosition_b, hubPosition_s);


    regularityentropy regE;

    //For boundary edges
    regFeature.calculateEdges_method3(hubPosition_b, horEdgeFeatureMatrix, verEdgeFeatureMatrix, 0);

    //For skeletal edges
    regFeature.calculateEdges_method3(hubPosition_s, horEdgeFeatureMatrix, verEdgeFeatureMatrix, 1);

    // horizonal edge length entropy
    regEntropy += regE.calculateEntropy(horEdgeFeatureMatrix, 0.01);

    // vertical edge length entropy
    regEntropy += regE.calculateEntropy(verEdgeFeatureMatrix, 0.01);


    int quadNum = (rowNum-1) * (colNum-1);
    //Calculate quads angles.
    MatrixType angleFeatureMatrix(6, quadNum);

    // angles on boundary
    regFeature.calculateAnglesForStandardSide(hubPosition_b, angleFeatureMatrix, 0);

    // angles on skeletal
    regFeature.calculateAnglesForStandardSide(hubPosition_s, angleFeatureMatrix, 1);

    // the angle entropy
    regEntropy += regE.calculateEntropy(angleFeatureMatrix, 0.01);

    return regEntropy;
}
