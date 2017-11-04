#include "srepprocrustes.h"

srepprocrustes::srepprocrustes()
{
}


srepprocrustes::srepprocrustes(const char* srepFolder)
{
    this->srepFolder = string(srepFolder);
}



/* srepFolder: directory where a list of s-reps which we will performance alignment on.
*/
int srepprocrustes::procrustesAlignment( ){
    cout<<"+++++++++++ Prepare to do Procrustes alignment to sreps under: "<<this->srepFolder<<" ++++++++++++"<<endl;
    vector<string> srepPathNames = tls.getModelname(this->srepFolder);
    int interpolationLevel = 2;
    //Align every mesh to the mean mesh.
    ProcrustesFilterType::Pointer procrustesFilter = ProcrustesFilterType::New();
    procrustesFilter->SetNumberOfInputs(srepPathNames.size());
    procrustesFilter->SetUseInitialAverageOn(); //set to on, use the actual mean.
    procrustesFilter->SetUseScalingOn();
    procrustesFilter->SetUseSingleIterationOff();
    procrustesFilter->SetUseNormalizationOn();//on by default.
    procrustesFilter->SetAlignRotationOn();

    if(!srepPathNames.empty()){
        cout<<"-------------there are: "<<srepPathNames.size()<<"sreps in this folder."<<endl;
        for(int q =0; q<srepPathNames.size();q++){
            // Create a alignsrep object
            alignsrep as(srepPathNames[q].c_str(), interpolationLevel);

            // Get points position datas and areas.
            as.weightAreaToBoundaryPoints();

            vector<Vector3D> bPoints = as.getBoundaryPoints();
            vector<Vector3D> sPoints = as.getPointsOnSkeletal(); // Get skeletal points.
            //vector<double> wAreas = as.getBoundaryPointsAreas();

            // Save these values to mesh.
            int mesh_b_index = 0;
            int mesh_s_index =0;

            //Per srep's boundary points compsite its own mesh.
            this->mesh_b = MeshType::New();
            //this->mesh_b->GetPoints()->Reserve(this->quadFig->getRowCount() * this->quadFig->getColumnCount());
            // Create points
            PointType point;
            this->mesh_s = MeshType::New();
            //this->mesh_s->GetPoints()->Reserve(this->quadFig->getRowCount() * this->quadFig->getColumnCount());

            //Loop each entry in bPoints, get the x,y,z coordinate of this srep's boundary points.
            for(unsigned i = 0; i< bPoints.size(); i++){
                point[0]= bPoints[i].getX();
                point[1]= bPoints[i].getY();
                point[2]= bPoints[i].getZ();
                this->mesh_b->SetPoint(mesh_b_index, point);
                mesh_b_index++;
             }

            //Loop each entry in sPoints, get the x,y,z coordinate of this srep's skeletal points.
            for(unsigned i = 0; i< sPoints.size(); i++){
                point[0]= sPoints[i].getX();
                point[1]= sPoints[i].getY();
                point[2]= sPoints[i].getZ();
                this->mesh_s->SetPoint(mesh_s_index, point);
                mesh_s_index++;
            }

            //Output the points for check.
            /*std::cout << "Mesh points of skeletal is: " << this->mesh_s->GetNumberOfPoints() << std::endl;
            PointsIterator pointIterator = this->mesh_s->GetPoints()->Begin();

            //Output these mesh points.
            PointsIterator end = this->mesh_s->GetPoints()->End();
            while( pointIterator != end )   {
              MeshType::PointType p = pointIterator.Value();  // access the point
              std::cout << p << std::endl;                    // print the point
              ++pointIterator;                                // advance to next point
            }*/

            // Get the center of mass for this srep
           /* PointType centerOfMass;
            for( int dim = 0; dim < 3; dim++ ) {
                centerOfMass[dim] =0;
            }
            PointsIterator pntItX;
            for( pntItX = this->mesh_s->GetPoints()->Begin(); pntItX != this->mesh_s->GetPoints()->End(); ++pntItX ) {
                for( int dim = 0; dim < 3; dim++ ) {
                    centerOfMass[dim] += pntItX.Value()[dim];
                }
            }
            cout<<"-------------------------this->mesh_s->GetPoints()->Size() is: "<<this->mesh_s->GetPoints()->Size()<<endl;
            for( int dim = 0; dim < 3; dim++ ) {
                centerOfMass[dim] /= this->mesh_s->GetPoints()->Size();
            }

            // Center this mesh to the origin
            cout<<"------------------------the centered input skeletal mesh is:"<<endl;
            for( pntItX = this->mesh_s->GetPoints()->Begin(); pntItX != this->mesh_s->GetPoints()->End(); ++pntItX ) {
                for( int dim = 0; dim < 3; dim++ ) {
                    pntItX.Value()[dim] -= centerOfMass[dim];
                }
                cout<<pntItX.Value()<<endl;
            }
*/

            // Add this srep's boundary mesh into procrustesFilter
            //procrustesFilter->SetInput(q, this->mesh_b);
            // Test only...
            procrustesFilter->SetInput(q, this->mesh_s);

            // Store this srep's skeletal points (mesh_s) into vector list.
            this->mesh_s_list.push_back(this->mesh_s);
        }

        // Execute Procrustes alignment.
        procrustesFilter->Update();

        /*cout<<"the 0 th skeletal after transform is: "<<endl;
        MeshType::Pointer _RegisteredMesh = procrustesFilter->GetOutput(0);
        PointsIterator pointIterator = _RegisteredMesh->GetPoints()->Begin();
        //Output these mesh points.
        PointsIterator end = _RegisteredMesh->GetPoints()->End();
        while( pointIterator != end )   {
          MeshType::PointType p = pointIterator.Value();  // access the point
          std::cout << p << std::endl;                    // print the point
          ++pointIterator;                                // advance to next point
        }*/


        // Get the transforms of each srep's boundary points, and apply to its own skeletal points.
        for(unsigned i =0; i<srepPathNames.size();i++){

            /*this->affineTransform = procrustesFilter->GetTransform(i);
            procrustesFilter->PrintTransform(this->affineTransform);

            MeshType::Pointer outputMesh = procrustesFilter->GetOutput(i);*/

            // Create a Filter
            FilterType::Pointer filter = FilterType::New();

            // Get affineTransform used to align mesh_b i.
            this->affineTransform = procrustesFilter->GetTransform(i);

            // Connect the inputs. (Apply ith affineTransform to mesh_s i.)
            filter->SetInput(this->mesh_s_list[i]);
            filter->SetTransform(this->affineTransform);

            // Execute the filter
            filter->Update();
            //print out the translation detail information, offset...
            std::cout << "Filter: " << filter;

            // Get the Smart Pointer to the Filter Output
            MeshType::Pointer outputMesh = filter->GetOutput();

            cout << "Output mesh_s has " << outputMesh->GetNumberOfPoints()<< " points." <<endl;

            // Get the the point container
            MeshType::PointsContainerPointer transformedPoints = outputMesh->GetPoints();
            PointsContainerType::ConstIterator it = transformedPoints->Begin();
            while( it != transformedPoints->End() ) {
              PointType p = it.Value();
              std::cout.width( 5 ); std::cout << p[0] << ", ";
              std::cout.width( 5 ); std::cout << p[1] << ", ";
              std::cout.width( 5 ); std::cout << p[2] << std::endl;
              ++it;
            }

            // Write the transformed x,y,z of the skeletal sheet into its m3d files.
            writeAlignToSkeletal(srepPathNames[i], outputMesh);
        }

        return EXIT_SUCCESS;
    }
    else{
        cout<<"You input a invalid srep file, please check the path or name of this srep!"<<endl;
        return EXIT_FAILURE;
    }
}




/* Save the aligned skeletal information into its srep file(.m3d file).
 *
*/
void srepprocrustes::writeAlignToSkeletal(string filePathName, MeshType::Pointer outputMesh){
    cout<<"----------write....."<<filePathName<<" to /models/alignment/"<<endl;//need change.

    Registry registry;
    registry.readFromFile(filePathName.c_str(), false);

    M3DQuadFigure* quadFig = tls.GetQuadFigure(filePathName.c_str());

    //Get the the point container
    MeshType::PointsContainerPointer transformedPoints = outputMesh->GetPoints();
    PointsContainerType::ConstIterator it = transformedPoints->Begin();

    //int count=0;
    //Loop each primitive, get the x,y,z coordinate of this srep's boundary.
    for(unsigned u =0; u<quadFig->getRowCount(); u++){
        for(unsigned v =0; v<quadFig->getColumnCount(); v++){
            if( it != transformedPoints->End() ) {
              PointType p = it.Value();
              // Set spoke hub position to new value.
              quadFig->getPrimitivePtr(u,v)->setX(p[0],p[1],p[2]);
              ++it;
              //count++;
              //cout<<"------------------loop time is: "<<count<<endl;
            }
        }
    }

    // Write the changed quadFig to m3d file.
    char figureStr[1024];
    //sprintf(figureStr, "model.figure[%d]", figNum);
    sprintf(figureStr, "model.figure[%d]", 0);//In m3d file, every figure is identified as figure[0], cannot change the index, why?
    //cout<<"Message from movespokes.cpp; model.figure[%d] is: "<<figureStr<<endl;
    quadFig->writeFigure(figureStr,registry); //write the figure into registry.

    string filename = tls.getPathName(filePathName) + string("/output/")+ tls.getFileNameWithExtension(filePathName);
    //cout<<"-------------filename is:---"<<filename<<endl;

    //overwrite the input file to save the changes for next loop.
    registry.writeToFile(filename.c_str());

}

