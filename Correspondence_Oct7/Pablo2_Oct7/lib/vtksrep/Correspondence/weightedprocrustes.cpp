/* Performance a generalized procrustes analysis to sreps.
 * X: 3*k*n, storing all the nth srep's points set. Each srep is a 3*k matrix, 3 is x y z dimesion, k is points number.
 * W: n*k, storing all the nth srep's area weight, each row is a srep, each column is a point. *
 * realWeight: If true, we add a real area weight correspondence to each point. we can doing procrustes focus on shape not the representation.
 * If false (weightedprocrustes pro(false)), means we use equal weight = 1 to each points.
 * We can use this class to do procrustes to any points set, first get tanslation, scaling and rotaion information from the input data set (X),
 * and then call function applyTransform(X_s) can apply this transform to any given dataset (X_s).
 *
 * Liyun Tu
 * Apr 1, 2014.
*/



#include "weightedprocrustes.h"



using namespace std;


weightedprocrustes::weightedprocrustes()
{
}


weightedprocrustes::weightedprocrustes(bool realWeight)
{
    this->realWeight = realWeight;
}


/* Generalized procrustes analysis for weighted srep points.
 * X being an m-by-k-by-n sample of boundary points for all the sreps.
 * W is a weight matrix holding the area weight for each boundary points on each srep. In W, each row is a srep.
 * m : dimension of space that landmarks lie
 * j : number of landmarks
 * i : number of samples
*/
vector<MatrixType> weightedprocrustes::weightedGPA(vector<MeshType::Pointer> X, VectorPoints W){

    vector<MatrixType> transformed_X; // holding the transformed sreps boundary points.
    // Step 1: Compute the weighted centroid of each srep
    // sum the area of each srep.
    vector<double> srepAreas;
    PointType point;
    for(unsigned i=0; i<W.size(); i++){
        double srepArea = 0;
        for(unsigned j=0; j<W[0].size(); j++){
            srepArea += W[i][j];
        }
        srepAreas.push_back(srepArea); //srepAreas[i] is the ith srep's area.
    }    

    PointsContainerPointer pnt;
    PointsIterator pntIt;
    // Multiply each point's xyz with its area weight, divide by its srep areas. Each srep's centroid is the sum of its weighted points.
    MeshType::Pointer Xbar = MeshType::New();
    for(unsigned i=0; i<X.size(); i++){ // srep index
        // Each iteration the point begin from 0
        for( int dim = 0; dim < 3; dim++ ) {
            point[dim] = 0;
        }
        pnt = X[i]->GetPoints();
        int j =0; // point index
        cout<<"----Import waring: Currently set all area weight to 1, and use 1/(n-1) instead of 1/n. If use real area weight,\n"
            <<" you should uncomment the 63th line W[i][j]/srepAreas[i] in weightedprocrustes::weightedGPA(): "<<endl;
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            for( int dim = 0; dim < 3; dim++ ) {
                if(this->realWeight) // If use real weight, you shoud uncomment this line.
                    point[dim] += pntIt.Value()[dim]*W[i][j]/srepAreas[i];
                else // Use area weight = 1, use 1/(n-1) instead of 1/n;
                    point[dim] += pntIt.Value()[dim]*W[i][j]/(srepAreas[i]-1);
            }            
            j++;            
        }
        cout<<endl;
        Xbar->SetPoint(i, point); // Xbar[i] storing ith srep's centroid.
    }
    this->centroid = Xbar; // save each srep's centroid. centroid[i] storing ith srep's centroid.


    // Step 2: Center each srep to the origin
    vector<MeshType::Pointer> Xp;
    PointsIterator XbarPosIt = Xbar->GetPoints()->Begin();
    for(unsigned i=0; i<X.size(); i++){ // srep index
        MeshType::Pointer centerPos = MeshType::New();
        pnt = X[i]->GetPoints();        
        int j =0; // point index
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            for( int dim = 0; dim < 3; dim++ ) {
                point[dim] = pntIt.Value()[dim] - XbarPosIt.Value()[dim];
            }
            centerPos->SetPoint(j, point);            
            j++;
        }
        Xp.push_back(centerPos);
        this->mesh_Xp.push_back(centerPos); // save Xp for display.
        XbarPosIt++;
    }
    // Output the centered srep
    /*for(unsigned i=0; i<Xp.size(); i++){ // srep index
        pnt = Xp[i]->GetPoints();
        cout<<"----------------The Xp["<<i<<"] is: "<<endl;
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            cout<<pntIt.Value()<<endl;
        }
    }*/


    // Step 3: Compute the scale factor of each s-rep using the centered srep.
    vector<double> theta;
    for(unsigned i = 0; i< Xp.size(); i++){ // srep index
        double theta_i = 0; // this srep's weighted PDM size.
        pnt = Xp[i]->GetPoints();
        int j =0; // point index
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            double p = 0;
            for( int dim = 0; dim < 3; dim++ ) {
                p += pntIt.Value()[dim] * pntIt.Value()[dim]; // Frobenius norm(2-norm)
            }
            if(this->realWeight)
                p = p * W[i][j] / srepAreas[i];
            else // All weight=1, use 1/(n-1) instead of 1/n;
                p = p * W[i][j] / (srepAreas[i]-1);

            theta_i += p; // theta_i is the ith srep's size.
            j++;
        }

        theta.push_back(sqrt(theta_i)); //theta[i] storing the ith srep's scalar.
    }
    this->scaleVector = theta; // save each srep's scale factor. scaleVector[i] is ith srep's scale factor.
    // Output theta
    /*cout<<"------------theta for each srep is: "<<endl;
    for(unsigned i=0; i<theta.size();i++){
        cout<<theta[i]<<"  ";
    }
    cout<<endl;*/


    // Step 4: Compute the scaled, centered srep
    vector<MeshType::Pointer> Q;
    for(unsigned i = 0; i< Xp.size(); i++){ // srep index
        pnt = Xp[i]->GetPoints();
        MeshType::Pointer scalePos = MeshType::New();
        int j =0; // point index
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            for( int dim = 0; dim < 3; dim++ ) {
                point[dim] = pntIt.Value()[dim]/theta[i];
            }
            scalePos->SetPoint(j, point);
            j++;
        }
        Q.push_back(scalePos); // sreps in Q have equal size(unit size?).
        this->mesh_Q.push_back(scalePos); // save Q for display.
    }
    // Output scaled srep
    /*cout<<"------------The scaled srep is: "<<endl;
    for(unsigned i=0; i<Q.size(); i++){ // srep index
        pnt = Q[i]->GetPoints();
        cout<<"----------------The Q["<<i<<"] is: "<<endl;
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            cout<<pntIt.Value()<<endl;
        }
    }*/

    // Prepare the de-scale factor for later. For pablo display. should de-scale the points back to (0,1) space.
    double avgTheta =0;
    for(unsigned i=0; i<theta.size();i++){
        avgTheta += log(theta[i]);
    }
    avgTheta /= theta.size();
    double thetaBar = exp(avgTheta);


    // Step 5: Rotation    
    // For each srep in Q
    for(unsigned j = 0; j< Q.size(); j++){ // srep index

        // Step 5.1A: Compute the m*m covariance matrix(C) for each s-rep.
        int pointsNum = Q[0]->GetNumberOfPoints(); //----------------------------------------------------------Q, Xp,
        // Copy this srep's mesh into matrix
        MatrixType Qj; // jth srep after centering, scaling.
        Qj.set_size(3, pointsNum);
        pnt = Q[j]->GetPoints(); //----------------------------------------------------------------------------Q, Xp,
        int k =0; // point index
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            for( int dim = 0; dim < 3; dim++ ) {
                Qj[dim][k] = pntIt.Value()[dim];
            }
            k++;
        }
        // Create weight matrix.
        MatrixType Wj; // diganal matrix holing the area weight of jth srep in the dignoal entries.
        Wj.set_size(pointsNum, pointsNum);
        for(unsigned r=0; r<pointsNum;r++){
            for(unsigned c=0; c<pointsNum;c++){
                Wj[r][c] = 0; //initial all entries to 0.
            }
        }
        for(unsigned i =0; i<pointsNum; i++){
            Wj[i][i] = W[j][i]; // the jth srep's ith point weight.
        }
        // Compute the covariance matrix for jth srep.
        MatrixType Cj = Qj * Wj * Qj.transpose();

        // Output Cj
        /*cout<<"----------the Cj is: "<<endl;
        for(unsigned r=0; r<Cj.rows();r++){
            for(unsigned c=0; c<Cj.columns();c++){
                cout<<Cj[r][c]<<"  ";
            }
            cout<<endl;
        }*/

        // Step 5.1B: Compute the eigenvalue and eigenvector of Cj, sort in correspondence descent order.
        vnl_symmetric_eigensystem<double> eig(Cj);
        // Out D and V
        /*cout<<"D is: "<<endl; //eigenValues, default in increasing order.
        for(unsigned r=0; r<eig.D.rows();r++){
            cout<<eig.D(r,r)<<endl;
        }
        cout<<"V is: "<<endl; //eigenVector, sequence is correspondence to its eigenValues.
        for(unsigned r=0; r<3;r++){
            for(unsigned c=0; c<3;c++){
                cout<<eig.V[r][c]<<"  ";
            }
            cout<<endl;
        }*/
        // Sort the eigenValues in decreasing order. And reorder its corrspondence eigenvector.
        vnl_diag_matrix<double> D(3,3); //eigenValues in decreasing order.
        D(0,0) = eig.D(2,2);
        D(1,1) = eig.D(1,1);
        D(2,2) = eig.D(0,0);
        // reorder its correspondence eigenvector by swap column 1 and 3.
        MatrixType V(3,3);
        for(unsigned int r=0; r<3;r++){
            V[r][0] = eig.V[r][2];
            V[r][1] = eig.V[r][1];
            V[r][2] = eig.V[r][0];
        }
        // If the sign of the first element in V was not positive, then multiply by -1.
        // This is used to avoid the reflection?? because some srep's V is not same as matlab compute.
        // determinent <1, reflection, first col is V1, second col is V2, third col is V3. sum of each col, is the determinant of V1,V2,V3??
        for(unsigned int c=0; c<3; c++){
            double detOfCol = 0;
            for(unsigned int r=0; r<3; r++){
                detOfCol += V[r][c];
            }
            if(detOfCol<0){
                for(unsigned int r=0; r<3; r++){
                    V[r][c] *= -1;
                }
            }
         }
        cout<<"D is: "<<endl; //eigenValues, default in increasing order.
        for(unsigned r=0; r<eig.D.rows();r++){
            cout<<D(r,r)<<endl;
        }
        cout<<"V is: "<<endl; //eigenVector, sequence is correspondence to its eigenValues.
        for(unsigned r=0; r<3;r++){
            for(unsigned c=0; c<3;c++){
                cout<<V[r][c]<<"  ";
            }
            cout<<endl;
        }
        // Check if the eigenVector are orthogonal. V=[v1,v2,v3], if orthogonal: dot(v1,v1)=1, dot(v1,v2)=0, dot(v1,v3)=0
        MatrixType VTV = V.transpose() * V;
        cout<<"------------the VTV is: "<<endl;
        for(unsigned int j=0; j<VTV.rows(); j++){
            for(unsigned int k=0; k<VTV.columns(); k++){
                cout<<VTV[j][k]<<"  ";
            }
            cout<<endl;
        }

        // Step 5.2: Apply the rotation to data set(here apply to skeletal points).
        MatrixType Rot = V.transpose(); // the rotation matrix is the transpose of the eigenvector.
        /*cout<<"------------the Rot is: "<<endl;
        for(unsigned int j=0; j<Rot.rows(); j++){
            for(unsigned int k=0; k<Rot.columns(); k++){
                cout<<Rot[j][k]<<"  ";
            }
            cout<<endl;
        }*/
        // Apply Rot to this srep's boundary points.
        MatrixType rotPoints = Rot * Qj;
        /*cout<<"------------the rotated boundary points is: "<<endl;
        for(unsigned int j=0; j<rotPoints.rows(); j++){
            int count=0;
            for(unsigned int k=0; k<rotPoints.columns(); k++){
                cout<<rotPoints[j][k]<<"  ";
                count++;
            }
            cout<<endl;
            cout<<"-------------There are "<<count<<"points."<<endl;
        }*/

        // For pablo display. should de-scale the points back to (0,1) space.
        for(unsigned int i=0; i<3; i++){
            for(unsigned int j=0; j<pointsNum; j++){
                rotPoints[i][j] *= thetaBar;

                // For pablo display, we move the origin to the center of pablo's window (0.5,0.5,0.5)
                rotPoints[i][j] += 0.5;
            }
        }

        // Store this transformed sreps' boundary to list
        transformed_X.push_back(rotPoints);

        // Store this srep's rotation matrix into list.
        this->rotationMatrix.push_back(Rot); // save rotation to list.
    }

    return transformed_X;
}





/* Apply the rotationMatrix, scaleVector and centroid to another points set (can be anydata sets). Here apply to skeletal points set and
 * descale and move to (0.5,0.5,0.5) for pablo display.
 * Here the the rotationMatrix, scaleVector and centroid is gotten from the srep's boundary points.
 * apply to each srep's own skeletal points.
 * corr_X[i] holding the ith srep's skeletal points.
*/
vector<MatrixType> weightedprocrustes::applyTransform(vector<MeshType::Pointer> corr_X){
    vector<MatrixType> transformedPoints;

    PointsContainerPointer pnt;
    PointsIterator pntIt;
    PointType point;

    PointsIterator centroidPntIt = this->centroid->GetPoints()->Begin();
    for(unsigned int i =0; i<corr_X.size(); i++){ // srep index
        int pointsNum = corr_X[0]->GetNumberOfPoints();

        // Center corr_X to origin
        pnt = corr_X[i]->GetPoints();
        MeshType::Pointer corr_Xp = MeshType::New();
        int j =0; // point index
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            for( int dim = 0; dim < 3; dim++ ) {
                point[dim] = pntIt.Value()[dim] - centroidPntIt.Value()[dim];
            }
            corr_Xp->SetPoint(j, point);
            j++;
        }
        centroidPntIt++; // exceed to next srep's centroid.

        // Output the centered srep
        /*pnt = corr_Xp->GetPoints();
        cout<<"----------------The Xp["<<i<<"] is: "<<endl;
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            cout<<pntIt.Value()<<endl;
        }*/


        // Scaling the centered srep.
        pnt = corr_Xp->GetPoints();
        MatrixType corr_Q(3, pointsNum);
        j = 0; // point index
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            for( int dim = 0; dim < 3; dim++ ) {
                corr_Q[dim][j] = pntIt.Value()[dim]/this->scaleVector[i];
            }
            j++;
        }
        // Output the centered, scaled srep
        /*for(unsigned int j=0; j<corr_Q.rows(); j++){
            for(unsigned int k=0; k<corr_Q.columns(); k++){
                cout<<corr_Q[j][k]<<"  ";
            }
        }*/

        // Rotate the centered, scaled srep.
        MatrixType rotPoints = this->rotationMatrix[i] * corr_Q;

        // For pablo display. should de-scale the points back to (0,1) space.
        double avgTheta =0;
        for(unsigned i=0; i<this->scaleVector.size();i++){
            avgTheta += log(this->scaleVector[i]);
        }
        avgTheta /= this->scaleVector.size();
        double thetaBar = exp(avgTheta);

        for(unsigned int i=0; i<3; i++){
            for(unsigned int j=0; j<pointsNum; j++){
                rotPoints[i][j] *= thetaBar;

                // For pablo display, we move the origin to the center of pablo's window (0.5,0.5,0.5)
                rotPoints[i][j] += 0.5;
            }
        }

        transformedPoints.push_back(rotPoints);
    }

    return transformedPoints;
}




/* Compute the procrustes distance over the aligned sreps samples.
 * alignedSreps: holding all the aligned points of each srep.
 * Get the sum of squared errors (pairwise differences betweent sreps in in srepList).
 * The difference between 2 srep is the sum of distance between the correspondence points on the 2 sreps.
 * Go through each srep, calculate its distance with each of the others srep.
 * 1 to 2,3,..,n; 2 to 3,4,...,n; 3 to 4,5,...,n; ...; n-1 to n. pairwise. no repeat. 1 to 2 is same as 2 to 1, only compute once.
*/
double weightedprocrustes::procrustesDistance(vector<MatrixType> srepList){

    double G =0;
    // Loop each pair of sreps(no repeat, AB and BA is same, only compute once.).
    for(unsigned int i = 0; i < srepList.size(); i++){ //srep A index
        for(unsigned int j = i+1; j < srepList.size(); j++){ //srep B index, only look at those follows A to keep no repeat.
            G = G + getFroDistance(srepList[i], srepList[j]);
        }
    }

    G = G/srepList.size();

    return G;
}




/* Compute the difference(distance) between srep A and srep B. Which is the Frobenius Norm(sum of squared) of A-B.
 * A and B have same size. Both with 3 rows, k columns.
*/
double weightedprocrustes::getFroDistance(MatrixType A, MatrixType B){

    double diff = 0;
    double squaredDiff =0; // difference between two points.

    // For each point
    for(unsigned int i =0; i< A.columns(); i++){
        for( int dim = 0; dim < 3; dim++ ) {
            diff = A[dim][i] - B[dim][i];  // access the point
            squaredDiff += diff * diff;
        }
    }

    return squaredDiff;
}



/* Compute the new spokes( r and u) between correspondence boundary point and skeletal point.
 * Return a matrix list, holding all the spokes' direction for each srep.
 * newSpokes[i] is the ith srep's spoke information. each spoke has xyz direction.
*/
vector<MatrixType> weightedprocrustes::getNewSpokesDirection(vector<MatrixType> bp, vector<MatrixType> sp){
    vector<MatrixType> newSpokesDirection;

    // For each srep
    for(unsigned int i = 0; i < bp.size(); i++){ //srep index
        // For each pair of correspondence points on boundary and skeletal
        MatrixType bPoints = bp[i]; // this srep's boundary points. 3*k.
        MatrixType sPoints = sp[i]; // this srep's skeletal points. 3*k.

        // This srep's spoke direction.
        MatrixType spokeVector = bPoints - sPoints;
        int pointNum = spokeVector.columns();

        // Normalize each spoke's direction to a unit vector
        for(unsigned int m=0; m< pointNum; m++){
            double spokeRadiu=0;
            for(unsigned int dim=0; dim<3; dim++){
                spokeRadiu += spokeVector[dim][m] * spokeVector[dim][m];
            }
            for(unsigned int dim=0; dim<3; dim++){
                spokeVector[dim][m] /= sqrt(spokeRadiu);
            }
        }

        // Save this srep's spokes.
        newSpokesDirection.push_back(spokeVector);
    }

    return newSpokesDirection;
}


/* For each pair of correspondence points on boundary and skeletal
*/
VectorPoints weightedprocrustes::getNewSpokesRadius(vector<MatrixType> bp, vector<MatrixType> sp){
    VectorPoints newSpokesRadius;

    // For each srep
    for(unsigned int i = 0; i < bp.size(); i++){ //srep index
        MatrixType bPoints = bp[i]; // this srep's boundary points. 3*k.
        MatrixType sPoints = sp[i]; // this srep's skeletal points. 3*k.
       /* cout<<"--------------------the bp is: "<<endl;
        for(unsigned int o=0; o<bPoints.rows(); o++){
            for(unsigned int m=0; m<bPoints.columns();m++){
                cout<<bPoints[o][m]<<"  ";
            }
            cout<<endl;
        }
        cout<<"--------------------the sp is: "<<endl;
        for(unsigned int o=0; o<sPoints.rows(); o++){
            for(unsigned int m=0; m<sPoints.columns();m++){
                cout<<sPoints[o][m]<<"  ";
            }
            cout<<endl;
        }*/
        // This srep's spoke direction.
        MatrixType spokeVector = bPoints - sPoints;
        /*cout<<"--------------------the spokeVector is: "<<endl;
        for(unsigned int o=0; o<spokeVector.rows(); o++){
            for(unsigned int m=0; m<spokeVector.columns();m++){
                cout<<spokeVector[o][m]<<"  ";
            }
            cout<<endl;
        }*/
        vector<double> spokeRadius;
        // Normalize each spoke's direction to a unit vector
        for(unsigned int m=0; m< spokeVector.columns(); m++){
            double spokeRadiu=0;
            for(unsigned int dim=0; dim<3; dim++){
                spokeRadiu += spokeVector[dim][m] * spokeVector[dim][m];
            }
            spokeRadius.push_back(sqrt(spokeRadiu));
        }

        // Save this srep's spokes.
        newSpokesRadius.push_back(spokeRadius);
    }

    return newSpokesRadius;
}


/* Return the centered points set.*/
vector<MatrixType> weightedprocrustes::getXp(){
    vector<MatrixType> matrix_Xp;

    for(unsigned int j = 0; j< this->mesh_Xp.size(); j++){ // srep index
        // Get points number
        int pointsNum = this->mesh_Xp[0]->GetNumberOfPoints();
        // Copy this srep's mesh into matrix
        MatrixType temp; // jth srep after centering, scaling.
        temp.set_size(3, pointsNum);
        PointsContainerPointer pnt;
        PointsIterator pntIt;
        pnt = this->mesh_Xp[j]->GetPoints();
        int k =0; // point index
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            for( int dim = 0; dim < 3; dim++ ) {
                temp[dim][k] = pntIt.Value()[dim];
            }
            k++;
        }
        matrix_Xp.push_back(temp);
    }

    return matrix_Xp;
}


/* Return the centered points set.*/
vector<MatrixType> weightedprocrustes::getQ(){
    vector<MatrixType> matrix_Q;

    for(unsigned int j = 0; j< this->mesh_Q.size(); j++){ // srep index
        // Get points number
        int pointsNum = this->mesh_Q[0]->GetNumberOfPoints();
        // Copy this srep's mesh into matrix
        MatrixType temp; // jth srep after centering, scaling.
        temp.set_size(3, pointsNum);
        PointsContainerPointer pnt;
        PointsIterator pntIt;
        pnt = this->mesh_Q[j]->GetPoints();
        int k =0; // point index
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            for( int dim = 0; dim < 3; dim++ ) {
                temp[dim][k] = pntIt.Value()[dim];
            }
            k++;
        }
        matrix_Q.push_back(temp);
    }

    return matrix_Q;
}


