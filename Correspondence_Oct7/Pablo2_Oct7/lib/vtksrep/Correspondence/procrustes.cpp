#include "procrustes.h"



using namespace std;

procrustes::procrustes()
{
   // mesh_Xp = MeshType::New();
}



/* Doing ordinary procrustes alignment. Conform Y to X.
 * X: target (mean). 3-by-n.
 * Y: source. 3-by-n. The matrix to be transformed.
 */
procrustes::MeshType::Pointer procrustes::OPA(MeshType::Pointer MX, MeshType::Pointer MY){
    // Convert MeshType X, Y to MatrixType
    MatrixType X, Y;
    X.set_size(3, this->pointsNum);
    Y.set_size(3, this->pointsNum);
    PointsIterator pntItX;
    int i =0;
    PointsIterator pntItY = MY->GetPoints()->Begin();
    for( pntItX = MX->GetPoints()->Begin(); pntItX != MX->GetPoints()->End(); ++pntItX ) {
        for( int dim = 0; dim < 3; dim++ ) {
            X[dim][i] = pntItX.Value()[dim];
            Y[dim][i] = pntItY.Value()[dim];
        }
        i++;
        ++pntItY;
    }

    // Center of mass
    MatrixType Xmean, Ymean;
    Xmean = getCenterOfMass(X);
    Ymean = getCenterOfMass(Y); // (0,0,0) infact no need to center here, because already centered in GPA.
    // Output center of mass for check
    cout<<"------------the Center of mass of Yc (source) is: "<<endl;
    for(unsigned int i =0; i<Ymean.rows();i++){
        for(unsigned int j=0; j<Ymean.columns();j++){
            cout<<Ymean[i][j]<<"  ";
        }
        cout<<endl;
    }
    cout<<"------------the Center of mass of Xc (target) is: "<<endl;
    for(unsigned i =0; i<Xmean.rows();i++){
        for(unsigned j=0; j<Xmean.columns();j++){
            cout<<Xmean[i][j]<<"  ";
        }
        cout<<endl;
    }

    // Center the srep
    MatrixType Xc = centerThisSrep(X, Xmean);
    MatrixType Yc = centerThisSrep(Y, Ymean);

    // Output X and Y for check
    cout<<"------------the centered Yc (source) is: "<<endl;
    for(unsigned int i =0; i<Yc.rows();i++){
        for(unsigned int j=0; j<Yc.columns();j++){
            cout<<Yc[i][j]<<"  ";
        }
        cout<<endl;
    }
    cout<<"------------the centered Xc (target) is: "<<endl;
    for(unsigned i =0; i<Xc.rows();i++){
        for(unsigned j=0; j<Xc.columns();j++){
            cout<<Xc[i][j]<<"  ";
        }
        cout<<endl;
    }

    // SVD(Singular Value Decomposition)
    MatrixType x1 = Yc * Xc.transpose();
    cout<<"+++++++++++++++++++++++++++++++++++++++++++++the x1 is: "<<endl;
    for(int i=0; i<x1.rows();i++){
        for(int j=0; j<x1.columns();j++){
            cout<<x1[i][j]<<"  ";
        }
        cout<<endl;
    }
    vnl_svd<double> svd(x1);
    cout<<"+++++++++++++++++++++++++++++++++++++++++++++ svd.U() is: "<<endl;
    for(int i=0; i<svd.U().rows();i++){
        for(int j=0; j<svd.U().columns();j++){
            cout<<svd.U()[i][j]<<"  ";
        }
        cout<<endl;
    }
    cout<<"+++++++++++++++++++++++++++++++++++++++++++++ svd.V() is: "<<endl;
    for(int i=0; i<svd.V().rows();i++){
        for(int j=0; j<svd.V().columns();j++){
            cout<<svd.V()[i][j]<<"  ";
        }
        cout<<endl;
    }
    this->trans_R = svd.V() * svd.U().transpose();

    cout<<"+++++++++++++++++++++++++++++++++++++++++++++ this->trans_R is: "<<endl;
    for(int i=0; i<this->trans_R.rows();i++){
        for(int j=0; j<this->trans_R.columns();j++){
            cout<<this->trans_R[i][j]<<"  ";
        }
        cout<<endl;
    }

    MatrixType x2 = Xc.transpose() * this->trans_R * Yc;   //x2 is n-by-n
    /*cout<<"+++++++++++++++++++++++++++++++++++++++++++++the x2 is: "<<endl;
    for(int i=0; i<x2.rows();i++){
        for(int j=0; j<x2.columns();j++){
            cout<<x2[i][j]<<"  ";
        }
        cout<<endl;
    }*/
    double x2Trace = 0;
    for(unsigned int i = 0; i < x2.rows(); i++ ) {
        x2Trace += x2[i][i];
    }
    MatrixType x3 = Yc.transpose() * Yc;
    /*cout<<"+++++++++++++++++++++++++++++++++++++++++++++the x3 is: "<<endl;
    for(int i=0; i<x3.rows();i++){
        for(int j=0; j<x3.columns();j++){
            cout<<x3[i][j]<<"  ";
        }
        cout<<endl;
    }*/
    double x3Trace = 0;
    for( unsigned int i = 0; i < x3.rows(); i++ ) {
        x3Trace += x3[i][i];
    }
    this->trans_scale = x2Trace / x3Trace;
    cout<<"-----------------trans_scale is: "<<trans_scale<<endl;

    // Translation
    this->trans_T = Xmean - this->trans_scale * this->trans_R * Ymean; // should be 0, because all the data has been centered to origin.
    cout<<"------------Translation this->trans_T is: "<<endl;
    for(unsigned i =0; i<this->trans_T.rows();i++){
        for(unsigned j=0; j<this->trans_T.columns();j++){
            cout<<this->trans_T[i][j]<<"  ";
        }
        cout<<endl;
    }

    // repmat trans_T from 3-by-1 to 3-by-n
    int pointsNum = X.columns();
    MatrixType trans_T_n;
    trans_T_n.set_size(3, pointsNum);
    for(unsigned int i=0; i<3; i++){
        for(unsigned int j=0; j<pointsNum; j++){
            trans_T_n[i][j] = this->trans_T[i][0];
        }
    }

    // Transformed Y, 3-by-n.
    MatrixType Yp = this->trans_scale * trans_R * Y + trans_T_n;

    // Convert Yp to MeshType
    MeshType::Pointer trans_Y = MeshType::New();
    PointType point;
    cout<<"----the transformed srep is: "<<endl;
    for(unsigned int k = 0; k < Yp.columns(); k++){
        for(unsigned int dim = 0; dim < 3; dim++){
            point[dim] = Yp[dim][k];
            cout<<Yp[dim][k]<<"  ";
        }
        cout<<endl;
        trans_Y->SetPoint(k, point); // save point to mesh.
    }

    return trans_Y;
}




/* mesh_x_list holding a set of srep's boundary points, need to do procrustes matching.
 * mesh_x_list[i] storing one srep's boundary points. There are N sample sreps, i=1,2,...,N.
*/
/*void procrustes::GPA(vector<MeshType::Pointer> mesh_x_list){

    this->srepNums = mesh_x_list.size();
    this->pointsNum = mesh_x_list[0]->GetNumberOfPoints();

    // Initialize this->mesh_X. Center each srep in mesh_X.
    vector<MeshType::Pointer> mesh_Xp_temp = mesh_x_list;//centerSrepsToOrigin(mesh_x_list);

    // Make a copy of mesh_X. while loop for computing the mean, only mesh_Xp change.
    vector<MeshType::Pointer> source = mesh_Xp_temp;
    vector<MeshType::Pointer> mesh_Xpmean_input;
    cout<<"11111111111111111111111111111111111111111111"<<endl;
    double Gold = 1e16;
    double G = 1e15;
    double tol = 1e-6;
    int counter =0;
    do {
        Gold = G;

        // Loop each srep under the criterial. Minimise the sum of squared norms of pairwise differences.
        for(unsigned n=0; n< this->srepNums; n++){
            // Put the N-1 srep(not including the current srep mesh_x_list[n]) to compute this->mesh_Xpmean.
            for(unsigned i=0; i< this->srepNums; i++){
                if(i!=n)
                    mesh_Xpmean_input.push_back(mesh_Xp_temp[i]);
            }

            // Compute the mean srep over the N-1 srep. (Not including Xn. Because conform Xn to Xmean.). Storing in this->mesh_Xpmean.
            MeshType::Pointer mesh_Xpmean = calculateMean(mesh_Xpmean_input);
            mesh_Xpmean_input.clear(); // clear vector for next loop.

            // Conform nth srep to the mean srep. Replace the nth srep in this->mesh_Xp[n] with the transformed one.
            mesh_Xp_temp[n] = OPA(mesh_Xpmean, source[n]);
        }       

        counter++;
        cout<<"+++++++++++++++++++++counter is: "<<counter<<endl;
        G = getG(mesh_Xp_temp);
        //cout<<"----The pairwise difference is: "<<G<<endl;
    } while(counter<250);//(Gold-G > tol);//

    //cout<<"---------counter is: "<<counter<<", Gold-G=" << Gold-G<<endl;

    // this->mesh_Xp is the final transformed sreps..
    // write get and set method to output the transform information and transformed sreps(this->mesh_Xp).

    this->mesh_Xp = mesh_Xp_temp;
}*/



void procrustes::GPA(vector<MeshType::Pointer> mesh_x_list){

    this->srepNums = mesh_x_list.size();
    this->pointsNum = mesh_x_list[0]->GetNumberOfPoints();

    // Initialize this->mesh_X. Center each srep in mesh_X.
    vector<MeshType::Pointer> mesh_Xp = centerSrepsToOrigin(mesh_x_list);

    // Make a copy of mesh_X. while loop for computing the mean, only mesh_Xp change.
    vector<MeshType::Pointer> mesh_Xo = mesh_Xp;

    double Gold = 1e16;
    double G = 1e15;
    double tol = 1e-6;
int counter =0;
    do {
        Gold = G;

        // Compute the mean srep over the N-1 srep. (Not including Xn. Because conform Xn to Xmean.).
        MeshType::Pointer mesh_Xpmean = calculateMean(mesh_Xp);

        // Loop each srep under the criterial. Minimise the sum of squared norms of pairwise differences.
        for(unsigned int n=0; n< this->srepNums; n++){
            // Conform nth srep to the mean srep. Replace the nth srep in mesh_Xp with the transformed one.
            mesh_Xp[n] = OPA(mesh_Xpmean, mesh_Xo[n]);
        }

        counter++;

        G = getG(mesh_Xp);
        cout<<"----The pairwise difference is: "<<G<<endl;
    } while(counter<3);//(Gold-G > tol);

    // mesh_Xp is the final transformed sreps..
    // write get and set method to output the transform information and transformed sreps(mesh_Xp).

    this->mesh_Xp = mesh_Xp;

}



/* Center all the sample sreps by minus its own centroid.*/
vector<procrustes::MeshType::Pointer> procrustes::centerSrepsToOrigin(vector<MeshType::Pointer> mesh_x_list){

    //vector<MeshType::Pointer> mesh_Xp;

    PointsContainerPointer pos;
    PointsIterator pntIt;

    // Loop each srep.
    for(unsigned int i = 0; i< this->srepNums; i++){
        double x=0;
        double y=0;
        double z=0;
        // Loop each point on this srep
        pos = mesh_x_list[i]->GetPoints();
        for( pntIt = pos->Begin(); pntIt != pos->End(); ++pntIt ) {
            x += pntIt.Value()[0];
            y += pntIt.Value()[1];
            z += pntIt.Value()[2];
        }
        // Compute the centroid of each srep.
        x = x/this->pointsNum;
        y = y/this->pointsNum;
        z = z/this->pointsNum;        

        // Center this srep to the original. Storing the new xyz to mesh_x_list[i].
        for( pntIt = pos->Begin(); pntIt != pos->End(); ++pntIt ) {
            //cout<<"before center, pntIt.Value() is: "<<pntIt.Value()<<endl;
            pntIt.Value()[0] -= x;
            pntIt.Value()[1] -= y;
            pntIt.Value()[2] -= z;
            //cout<<"after center, pntIt.Value() is: "<<pntIt.Value()<<endl;
        }

        // Storing this centered srep into mesh_Xp
        //mesh_Xp.push_back(mesh_x_list[i]);
    }

    return mesh_x_list;
}






/* Get the sum of squared errors (pairwise differences betweent sreps in in mesh_x_list).
 * The difference between 2 srep is the sum of distance between the correspondence points on the 2 sreps.
 * Go through each srep, calculate its distance with each of the others srep.
 * 1 to 2,3,..,n; 2 to 3,4,...,n; 3 to 4,5,...,n; ...; n-1 to n. pairwise. no repeat. 1 to 2 is same as 2 to 1, only compute once.
*/
double procrustes::getG(vector<MeshType::Pointer> mesh_x_list){

    double G =0;
    for(unsigned int i=0; i<this->srepNums; i++){
        for(unsigned int j=i+1; j<this->srepNums; j++){
            G = G + getFroDistance(mesh_x_list[i], mesh_x_list[j]);
        }
    }

    G = G/this->srepNums;

    return G;
}


/* Compute the difference(distance) between mesh A and B. Which is the Frobenius Norm(sum of squared) of A-B.
*/
double procrustes::getFroDistance(MeshType::Pointer A, MeshType::Pointer B){

    PointsContainerPointer aPoints = A->GetPoints();
    PointsIterator pointItA; // A iterator

    PointsIterator pointItB = B->GetPoints()->Begin();

    double diff, squaredDiff; // difference between two points.

    for( pointItA = aPoints->Begin(); pointItA != aPoints->End(); ++pointItA ) {
        for( int dim = 0; dim < 3; dim++ ) {
            diff = pointItA.Value()[dim] - pointItB.Value()[dim];  // access the point
            squaredDiff += diff * diff;
        }
        ++pointItB;                               // advance to next point
    }
}



/* Get center of mass. return a 3-by-1 matrix holding the xyz position of this srep center.*/
procrustes::MatrixType procrustes::getCenterOfMass(MatrixType srepPos){
    MatrixType center;
    center.set_size(3,1);

    int pointsNum = srepPos.columns();

    for(unsigned int dim=0; dim<3; dim++){
        double rowSum = 0;
        for(unsigned int j=0; j<pointsNum; j++){
            rowSum += srepPos[dim][j];
        }
        center[dim][0] = rowSum/pointsNum;
    }

  return center;
}


/* Center the srep. Minus each point in X by COM.
 * X: srep
 * COM: X's center of mass, a point, a 3-by-1 matrix, holding the xyz coordinate of the center of mass.
*/
procrustes::MatrixType procrustes::centerThisSrep(MatrixType X, MatrixType COM){

    for(unsigned int i=0; i<3; i++){
        for(unsigned int j=0; j<X.columns(); j++){
            X[i][j] = X[i][j] - COM[i][0];
        }
    }

    return X;
}



/* Compute the mean srep by giving N-1 sreps.
 * this->mesh_Xmean_input holding N-1 srep's boundary points.
*/
procrustes::MeshType::Pointer procrustes::calculateMean(vector<MeshType::Pointer> mesh_Xp){
    unsigned int n = mesh_Xp.size(); // n = N-1
    cout<<"=============================The srep number used to compute mean is: "<<n<<endl;
    PointsContainerPointer inputSrep;
    PointsIterator pntIt;

    MeshType::Pointer Xpmean = MeshType::New();

    /*PointsContainerPointer meanPoints = this->mesh_Xpmean->GetPoints();//->Reserve(this->pointsNum);;
    PointsIterator meanIt; // mean iterator

    // Loop each of the N-1 srep. Sum the correspondence xyz of the n sreps, get a mesh storing 3*n double.
    for(unsigned i = 0; i< n; i++){
        inputSrep = mesh_Xp[i]->GetPoints();
        pntIt = inputSrep->Begin();


        for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt ) {
            for( int dim = 0; dim < 3; dim++ ) {
                meanIt.Value()[dim] += pntIt.Value()[dim];
            }
            ++pntIt;
        }
    }

    // Get the mean xyz.
    for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt ) {
        for( int dim = 0; dim < 3; dim++ )  {
            meanIt.Value()[dim] /= n;
        }
    }*/


    MatrixType mean(3,this->pointsNum);
    // Initial to 0.
    /*for(unsigned int dim =0; dim<3; dim++){
        for(unsigned int l =0; l<this->pointsNum; l++){
            mean[dim][l] =0;
        }
    }*/
    mean.fill(0);

    // Loop each of the N-1 srep. Sum the correspondence xyz of the n sreps, get a mesh storing 3*n double.
    for(unsigned int k = 0; k< n; k++){
        inputSrep = mesh_Xp[k]->GetPoints();        

        int i =0; // point index
        cout<<"------------------------the points used to compute mean is: "<<endl;
        for( pntIt = inputSrep->Begin(); pntIt != inputSrep->End(); ++pntIt ) {
            for( int dim = 0; dim < 3; dim++ ) {
                mean[dim][i] += pntIt.Value()[dim];                
            }
            cout<<pntIt.Value()<<endl;
            i++;
        }
    }

    // Get the mean xyz. And store into mesh_Xpmean
    PointType point;
    for(unsigned int j=0; j<this->pointsNum; j++) {
        for( int dim = 0; dim < 3; dim++ )  {
            point[dim] = mean[dim][j]/n;
        }
        Xpmean->SetPoint(j, point);
    }

    //Output the points for check.
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~the mean srep points is: " << Xpmean->GetNumberOfPoints() << std::endl;
    procrustes::PointsIterator pointIterator = Xpmean->GetPoints()->Begin();

    //Output these mesh points.
    procrustes::PointsIterator end = Xpmean->GetPoints()->End();
    while( pointIterator != end )   {
     procrustes::MeshType::PointType p = pointIterator.Value();  // access the point
      std::cout << p << std::endl;                    // print the point
      ++pointIterator;                                // advance to next point
    }

    return Xpmean;
}



vector<procrustes::MeshType::Pointer> procrustes::getTrans_X(){    

    // Because pablo use (0,1) coordinate, here we plus each point with the center (0.5,0.5,0.5)
    PointsIterator pntItX;
    PointsContainerPointer X;
    for(unsigned i=0; i<this->mesh_Xp.size(); i++){
        X = this->mesh_Xp[i]->GetPoints();
        for(pntItX = X->Begin(); pntItX != X->End(); ++pntItX){
            for( int dim = 0; dim < 3; dim++ ) {
                pntItX.Value()[dim] += 0.5;
            }
        }
    }

    return this->mesh_Xp;
}

procrustes::MatrixType procrustes::getTrans_Rot(){
    return this->trans_R;
}

double procrustes::getTrans_Sca(){
    return this->trans_scale;
}

procrustes::MatrixType procrustes::getTrans_Tra(){
    return this->trans_T;
}






/* Save matrix to .txt file.*/
void procrustes::saveMeshVector(const char* filename, vector<MeshType::Pointer> X){

    int srepNum = X.size();
    int pointsNum = X[0]->GetNumberOfPoints();
    // Save MeshType X to MatrixType
    MatrixType X_matrix;
    X_matrix.set_size(pointsNum*3, srepNum);

    PointsIterator pntItX;
    for(unsigned k=0; k< X.size();k++){ //number of sreps.
        int pointIndex =0; //loop each point.
        for( pntItX = X[k]->GetPoints()->Begin(); pntItX != X[k]->GetPoints()->End(); ++pntItX ) {
            for( int dim = 0; dim < 3; dim++ ) {
                X_matrix[pointIndex*3 + dim][k] = pntItX.Value()[dim];
            }
            pointIndex++;
        }
    }

    // Cout X_matrix to file.
    std::ofstream fout;
    fout.open(filename);

    if(fout)  {
        for(int i =0; i< X_matrix.rows();i++){
            for(int j=0; j< X_matrix.columns(); j++){ // Number of points.
                fout<< X_matrix[i][j]<<" ";
            }
            fout<<endl;
        }

        cout<<"Successfully saved matrix to: "<<filename<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();
}


/* Get the final transformed N srep's mean srep.
 * ProcrustesMean (also referred to as the Frechet mean)
*/
procrustes::MatrixType procrustes::getProcrustesMean(){

    PointsContainerPointer inputSrep;
    PointsIterator pntIt;

    MatrixType procrustesMean;
    procrustesMean.set_size(3,this->pointsNum);

    // Loop each of the N-1 srep. Sum the correspondence xyz of the n sreps, get a mesh storing 3*n double.
    for(unsigned k = 0; k< this->srepNums; k++){
        inputSrep = this->mesh_Xp[k]->GetPoints();

        int i =0; // point index
        for( pntIt = inputSrep->Begin(); pntIt != inputSrep->End(); ++pntIt ) {
            for( int dim = 0; dim < 3; dim++ ) {
                procrustesMean[dim][i] += pntIt.Value()[dim];
            }
            i++;
        }
    }

    for(unsigned j=0; j<this->pointsNum; j++) {
        for( int dim = 0; dim < 3; dim++ )  {
            procrustesMean[dim][j] /= this->srepNums;

            // Pablo's center is (0.5,0.5,0.5). We move the procrustesMean to Pablo's center.
            procrustesMean[dim][j] += 0.5;
        }
    }

    return procrustesMean;
}










