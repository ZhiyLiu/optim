/* Use variables gotten from optimizer as deltaU, deltaV. by adding this to u v, through interpolate we can got new p, r and u.
 * Compute the regularity features using the new spoke info(stored in an given folder(moved_sreps) in a m3d format)
 * At the meanwhile, the m3d files will be the input to compute the geometry entropy.
*/


#include "movespokes.h"


//using namespace std;

movespokes::movespokes()
{
}

movespokes::~movespokes()
{
}


movespokes::movespokes(string rootDir, int spokeNum, int srepNum, int side, int interpolationLevel) {
    this->rootDir = rootDir;
    this->interpolationLevel = interpolationLevel;

    toolsfunc tools;
    // Delete files under temp_data folder
    tools.deleteFiles(tools.connectStr(rootDir, "/temp_data/"));

    this->side = side;

    this->pMatrix.set_size(spokeNum*3, srepNum);
    this->uMatrix.set_size(spokeNum*3, srepNum);
    this->rMatrix.set_size(spokeNum, srepNum);

}





/* Return the up or down or crest regularity entropy. (0: up, 1: down, 2: crest)
 * quadFigList: the input to the iteration, original model, don't change during iteration.
 * coeff: all the sreps' changing variables set to this move.
 * registryPointerList: holding the list of Registry of the input sreps, used to write the new model to m3d files. By doing this, can save a lot
 * of time during iteration no need to read m3d file again and again...
 * side: 0: up, 1: down, 2: crest. Should pass right coeff to each type!!
 * coeff_this: For crest entropy, each spoke has a variables, four corner spokes has 0 values.
 *             For up or down entropy, four corner is fixed, don't have variables.
*/
double movespokes::calculateRegEntropy(const double * coeff, M3DQuadFigure* shiftingQuadFig, DoubleVec subVs,
                                       int varsNum, vector<M3DQuadFigure *> quadFigList, vector<vtkSmartPointer<vtkSRep> > srepfigList){
    double *coeff_this = new double [varsNum];
    double regEntropy = 0.0; //sum of all the srep's regularity entropies.

    //For each srep, update sreps one by one.
    for(int q = 0; q < quadFigList.size(); q++){

        // update the current fig info to quadFig for shifting, because we shouldn't change the original srep list.
        copySrepFig(quadFigList[q], shiftingQuadFig);

        //Split coeff to get variables for srep q.
        for(unsigned int j=0; j<varsNum; j++){
            coeff_this[j] = coeff[q*varsNum + j];
        }

        double entropy = 0.0;

        if(side == 2) {  // Shifting crest spokes, compute crest regularity entropy.
            // the subdivided t v are stored in the first entry.
            slidecrestspokes mCrestSpoke(interpolationLevel);
            entropy = mCrestSpoke.moveCrestSpokes(shiftingQuadFig, srepfigList[q], coeff_this, subVs);
        }
        else { // Shifting up or down spokes, compute up or down regularity entropy.
            slidestandardspokes mStandardSpoke(side, interpolationLevel);
            entropy = mStandardSpoke.moveStandardSpokes(shiftingQuadFig, coeff_this, subVs);
        }

        // Sum the three kinds of entropy of each srep, sum regularity of all the input sreps.
        regEntropy += entropy;

        storeSpokesToMatrix(q, shiftingQuadFig);

        // Save the optimized results (moved m3d file)
        //saveShiftedSrep(shiftingQuadFig, inputSreps[q]);
    }

    toolsfunc tools;
    //composite pMatrix, rMatrix, uMatrix and save to a fixed file: tools.connectStr(rootDir,"/temp_data/GEO_input.txt").
    string geoFileName = tools.connectStr(rootDir,"/temp_data/GEO_input.txt");
    saveGEO_Input_Matrix(geoFileName);

    // Delete pointer, release memory.
    delete [] coeff_this;

    // Return sum of the three part of entropies.
    std::cout<<"----regEntropy: "<<regEntropy<<std::endl;

    return regEntropy;
}




/* Generate the p r u matrix. *
 * Input a value reference, after loop each srep, the pMatrix, rMatrix, uMatrix will have values.
 * pMatrix: the first spokeNum*3 rows are up spokes, the following spokeNum*3 rows are down spokes, the rest are crest spokes.*
 *
 *
 * When compute geometry entropy, matlab will load this matrix.txt*/
void movespokes::storeSpokesToMatrix(int srepIndex, M3DQuadFigure* quadFig){
    int index = 0;

    switch(this->side){
    case 0: // For up
    case 1: // or down entropy.
        // Loop each spoke
        for(unsigned int i=0; i< quadFig->getRowCount(); i++){
            for(unsigned int j=0; j< quadFig->getColumnCount(); j++){
                M3DQuadPrimitive* prim = dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j));

                if(!prim){ //If point to nothing
                    std::cout<<"Message from: movespokes::saveStandSpokesToMatrix:: the prim pointer point is NULL!! EXIT_FAILURE!!"<<std::endl;
                    exit(1);
                }

                int pIndex = index * 3;
                this->pMatrix(pIndex,srepIndex) = prim->getX().getX();
                this->pMatrix(pIndex+1,srepIndex) = prim->getX().getY();
                this->pMatrix(pIndex+2,srepIndex) = prim->getX().getZ();

                // For up spokes.
                if(this->side == 0) {
                    this->uMatrix(pIndex,srepIndex) = prim->getU0().getX();
                    this->uMatrix(pIndex+1,srepIndex) = prim->getU0().getY();
                    this->uMatrix(pIndex+2,srepIndex) = prim->getU0().getZ();

                    this->rMatrix(index,srepIndex) = prim->getR0();
                }

                // For down spokes.
                else {
                    this->uMatrix(pIndex,srepIndex) = prim->getU1().getX();
                    this->uMatrix(pIndex+1,srepIndex) = prim->getU1().getY();
                    this->uMatrix(pIndex+2,srepIndex) = prim->getU1().getZ();

                    this->rMatrix(index,srepIndex) = prim->getR1();
                }

                index++;
            }
        }
        break;
    case 2: //For crest entropy.
        for(unsigned int u=0; u< quadFig->getRowCount(); u++){
            for(unsigned int v=0; v< quadFig->getColumnCount(); v++){

                if(u==0 || u==(quadFig->getRowCount()-1) || v==0 || v==(quadFig->getColumnCount()-1)){

                    M3DQuadEndPrimitive* primEnd = dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(u,v));

                    if(!primEnd){ //If point to nothing
                        std::cout<<"Message from: movespokes::saveSpokesToMatrix:: the primEnd pointer point is NULL!! EXIT_FAILURE!!"<<std::endl;
                        exit(1);
                    }

                    int pIndex = index * 3;

                    this->pMatrix(pIndex,srepIndex) = primEnd->getX().getX();
                    this->pMatrix(pIndex+1,srepIndex) = primEnd->getX().getY();
                    this->pMatrix(pIndex+2,srepIndex) = primEnd->getX().getZ();

                    this->uMatrix(pIndex,srepIndex) = primEnd->getUEnd().getX();
                    this->uMatrix(pIndex+1,srepIndex) = primEnd->getUEnd().getY();
                    this->uMatrix(pIndex+2,srepIndex) = primEnd->getUEnd().getZ();

                    this->rMatrix(index,srepIndex) = primEnd->getREnd();

                    index++;
                }
            }
        }
        break;
    }
}



/* Composite pMatrix, rMatrix, uMatrix to GEO_input matrix. And save GEO_input to .txt file.*/
void movespokes::saveGEO_Input_Matrix(string geoFileName){
    std::ofstream fout;
    fout.open(geoFileName.c_str());
   //cout<<"-----------------------------------geoFileName is: "<<geoFileName<<endl;

    if(fout)  {
        for(int i =0; i< this->pMatrix.rows();i++){
            for(int j=0; j< this->pMatrix.columns(); j++){
                fout<< this->pMatrix[i][j]<<" ";
            }
            fout<<std::endl;
        }
        //fout<<"---------------------------------------------the following is crest spokes: "<<endl;
        for(int i =0; i< this->rMatrix.rows();i++){
            for(int j=0; j< this->rMatrix.columns(); j++){
                fout<< this->rMatrix[i][j]<<" ";
            }
            fout<<std::endl;
        }
        //fout<<"---------------------------------------------crest spokes end! "<<endl;
        for(int i =0; i< this->uMatrix.rows();i++){
            for(int j=0; j< this->uMatrix.columns(); j++){
                fout<< this->uMatrix[i][j]<<" ";
            }
            fout<<std::endl;
        }

        //cout<<"Successfully saved p r u values to: "<<filename<<endl;
    }
    else
        std::cerr<<"Write out failed, cannot open the file!"<<geoFileName<<std::endl;

    fout.close();
}



/* save the optimized results (moved m3d file)
*/
void movespokes::saveShiftedSrep(M3DQuadFigure* quadFig, string filepath, string rootDir){
    Registry *registry = new Registry();
    registry->readFromFile(filepath.c_str(), false);

    toolsfunc tools;
    // Get filename without extension and path from a full path.
    string bn = tools.getBaseFileName(filepath);

    //new sreps each time the optimizer moving the primitive. this folder changing timely in the iteration.
    string movedSrepPath = rootDir + string("/moved_sreps/") + bn + ".m3d";

    writeSpokeInfoToFile(movedSrepPath.c_str(), *registry, quadFig);
}



/* Copy spoke info from the current fig to the shifting fig.
 * do similar thing as: quadFig = dynamic_cast<M3DQuadFigure*>(quadFigList[q]->clone()); // memory leak!!
 * but the above one has memory leak. The following function do not leak.
 */
void movespokes::copySrepFig(M3DQuadFigure* quadFig, M3DQuadFigure* shiftingFig){   
    // get the current quadFig info and set to shiftingQuadFig
    for(unsigned int i = 0; i < quadFig->getRowCount(); i++){
        for(unsigned int j = 0; j < quadFig->getColumnCount(); j++){
            M3DQuadPrimitive* prim_curr = dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j));
            Vector3D hubpos = prim_curr->getX();

            // set skeletal point position
            dynamic_cast<M3DQuadPrimitive*>(shiftingFig->getPrimitivePtr(i,j))->setX(hubpos.getX(), hubpos.getY(), hubpos.getZ());

            // set spoke length
            dynamic_cast<M3DQuadPrimitive*>(shiftingFig->getPrimitivePtr(i,j))->setR0(prim_curr->getR0());
            dynamic_cast<M3DQuadPrimitive*>(shiftingFig->getPrimitivePtr(i,j))->setR1(prim_curr->getR1());

            // set spoke direction
            dynamic_cast<M3DQuadPrimitive*>(shiftingFig->getPrimitivePtr(i,j))->setU0(prim_curr->getU0());
            dynamic_cast<M3DQuadPrimitive*>(shiftingFig->getPrimitivePtr(i,j))->setU1(prim_curr->getU1());

            // for crest spokes
            if(i==0 || i==(quadFig->getRowCount()-1) || j==0 || j==(quadFig->getColumnCount()-1)){
                double rEnd = dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(i,j))->getREnd();
                Vector3D uEnd = dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(i,j))->getUEnd();

                dynamic_cast<M3DQuadEndPrimitive*>(shiftingFig->getPrimitivePtr(i,j))->setREnd(rEnd);
                dynamic_cast<M3DQuadEndPrimitive*>(shiftingFig->getPrimitivePtr(i,j))->setUEnd(uEnd);
            }
        }
    }
}






/* This function will move spokes and save the moved srep to a file. But, because save to a m3d file and let CPNS read from m3d file
 * is very painful especially for large population. So this is not a good implement.
 * So we use moveStandardSpokesFast() instead of this one. The difference between this function and moveStandardSpokesFast() is:
 * moveStandardSpokesFast do not save the model files during iteration.
 * This function is used to generate some test sreps or call to save some needed intermidiate models. Because we do not need to save every temp
 * sreps, that cost too much time. So for large population, we use moveSpokesFast() instead.
 *
 * But, sometimes we will need to save the temprately srep's...
 * This function is used after optimization, we got the minimum objective function and its corresponding variables, we use this function to
 * generate the srep's m3d file for visualization.
 *
*/
/* varArray store all the varNums params. varArray[i*varNumEachSrep+j] is the ith srep's jth parameter.
*/
void movespokes::saveMovedModel(string rootDir, int count,  Registry * qRegistry, string quadFileName, M3DQuadFigure* quadFig){

    //new sreps each time the optimizer moving the primitive. this folder changing timely in the iteration.
    string movedSrepPath = rootDir + string("/moved_sreps/") + quadFileName;

    writeSpokeInfoToFile(movedSrepPath.c_str(), *qRegistry, quadFig);


   //For test only: copy a model out for demo how it was moved.
   /*if(fileName =="3.m3d"){
       cout<<"-----try to copy: "<<fileName<<endl;
       string outputPath = tools.connectFilePath(rootDir, "/models/modelChanging/3_", count, ".m3d");
       string strTemp = string("cp ") + tools.connectStr(rootDir,"/models/moved_sreps/3.m3d ") + outputPath;

        int returnValue = system(strTemp.c_str());
        if (returnValue != 0 ) std::cout << "Message from moveSpokes: Cannot copy that changing model file!" << endl;
    }*/
}


/* Write the new spoke information into file.
 * The inputPath is the original srep folder;
 * The outputPath is the folder where the new srep will be stored.
 * srepIndex: the index of current srep.
 * The original srep info is first read into a registry object, and then we update the information of this curQuadFig, and save this registry.
*/
void movespokes::writeSpokeInfoToFile(const char * outputPath, Registry &registry, M3DQuadFigure* quadFig){
    std::cout<<"----write to: "<<outputPath<<std::endl;
    char figureStr[1024];
    //sprintf(figureStr, "model.figure[%d]", figNum);
    sprintf(figureStr, "model.figure[%d]", 0);//In m3d file, every figure is identified as figure[0], cannot change the index, why?
    //cout<<"Message from movespokes.cpp; model.figure[%d] is: "<<figureStr<<endl;
    quadFig->writeFigure(figureStr,registry); //write the figure into registry.

    //overwrite the input file to save the changes for next loop.
    registry.writeToFile(outputPath);
    std::cout<<"----finished! "<<std::endl;
}
