#include "toolsfunc.h"

#include <sys/types.h>
#include <dirent.h>             /* struct DIR, struct dirent, opendir().. */

ControlParms * globalControl;	// Read the user's preferences file
int globalVerbosity;			// Current verbosity level of Pablo


toolsfunc::toolsfunc()
{
}



M3DQuadFigure* toolsfunc::GetQuadFigure(const char* figfilename){
    globalControl = new ControlParms(NULL, -1, false);	// Ignore user's preferences
    globalVerbosity = 0;
    globalControl->setDefault(OutputVerbosity, globalVerbosity);
    globalControl->setDefault(ReorderModels, 0);
    globalControl->setDefault(SmoothImages, false);
    globalControl->setDefault(ConvertImages, false);
    globalControl->setDefault(ByteOrder, 1);
    globalControl->setDefault(CompressImages, true);
    globalControl->setDefault(ShowLandmarks, true);

    P3DControl* p3d = new P3DControl(10);

    // m3d file was read in and stored in a M3DQuadFigure object
    p3d->read(figfilename, false);
    //std::cout<<"finish read m3d file"<<std::endl;
    M3DObject* m3dobject = p3d->getObjectPtr();
    M3DFigure * figure = m3dobject->getFigurePtr( 0 ) ;

    return dynamic_cast<M3DQuadFigure*>( figure );
}



string toolsfunc::connectFilePath(string rootDir, const char* specificDir, int changingIndex, const char* extension){
    //transform integer to string, and connect two string together.
    char num[20];
    sprintf(num, "%d", changingIndex);

    string path = rootDir + string(specificDir) + num + string(extension);

    return path;
}




string toolsfunc::conStr(string rootDir, const char* specificDir, string fileName, const char* extension) const{

    return rootDir + string(specificDir) + fileName + string(extension);
}




string toolsfunc::connectStr(string str1, const char* str2) const{
    return str1 + string(str2);
}



void toolsfunc::deleteFiles(string path){

    //string cmd = "exec rm -r " + path + "/*";
    string cmd = "find " + path + " -type f | xargs rm -rf";

    int returnValue1;

    try {
        returnValue1 = system(cmd.c_str());
    }
    catch(exception e){
            cout<<"The initial feature folder is null, no need to delete!"<<endl;
       }


    if (returnValue1 != 0 ){
        cout << "Error when delete regularity features files from folder: " << cmd << endl;
    }
}




/* Get all the filename under a folder into a vector. */
vector< std::string > toolsfunc::getModelname(string filefolder){
    vector< std::string > modelnamevector;
    DIR* dir = opendir(filefolder.c_str());
    if (!dir) {
        perror("opendir");
        exit(1);
    }

    //this structure is used for storing the name of each entry in turn.
    struct dirent* entry;

    //read the directory's contents, print out the name of each entry. scanning the directory using the readdir(), The first call returns the first entry of the directory.
    //Each successive call returns the next entry in the directory. When all entries have been read, NULL is returned.
    //printf("Directory contents:\n");
    while ( (entry = readdir(dir)) != NULL) {

        char * pattern1 = ".m3d";
        char * pattern2 = "~";
        //check if the pattern matchs. skip the "." and ".." entries, to avoid loops. See http://www.cpp-home.com/tutorials/107_6.htm
        /*if (entry->d_name && !strstr(entry->d_name, pattern) && strstr(entry->d_name, ".m3d")
                && (strcmp(entry->d_name, ".") != 0) && (strcmp(entry->d_name, "..") != 0)) {*/
            //printf("%s\n", "", entry->d_name); // d_name is a null-terminated character array, containing the name (NOT the path) of the entry.
        if (entry->d_name && !strstr(entry->d_name, pattern2) && strstr(entry->d_name, pattern1)) {            
            //Connect two string together.
            string path = filefolder + entry->d_name;
            modelnamevector.push_back(path);            
        }
    }

    //closedir() will return '0' on success, or '-1' if it failed.
    if (closedir(dir) == -1) {
        perror("closedir");
        exit(1);
    }

    return modelnamevector;
}




/* Save featurematrix to a txt file.
 * path: is the derectory and file name of the txt.
 * If the file not exist, it will create one.
 **/
void toolsfunc::savefeatures(string path, VectorTrainingSetFeaturesType featurematrix){

    const char* filename = const_cast<char*>(path.c_str());

    std::ofstream fout;
    fout.open(filename);

    if (fout)  {
        for (int i=0;i< featurematrix.size();i++){
            for (int j=0;j< featurematrix[i].size();j++)
                fout<< featurematrix[i][j]<<" ";
            fout<<endl;
        }
        //cout<<"Successfully saved to: "<<filename<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();
}



/* Combine 3 vectors to a matrix, save in .txt file. The 3 vectors has same columns.*/
void toolsfunc::saveMatrix3(string path, vector<double> A, vector<double> B, vector<double> C){

    const char* filename = const_cast<char*>(path.c_str());

    std::ofstream fout;
    fout.open(filename);

    if (fout)  {
        for (int i=0;i< A.size();i++){
            fout<< A[i]<<" ";
        }
        fout<<endl;
        for (int i=0;i< B.size();i++){
            fout<< B[i]<<" ";
        }
        fout<<endl;
        for (int i=0;i< C.size();i++){
            fout<< C[i]<<" ";
        }
        fout<<endl;
        //cout<<"Successfully saved feature to: "<<filename<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();
}



/* Combine 2 vectors to a matrix, save in .txt file. The 2 vectors has same columns.*/
void toolsfunc::saveMatrix2(string path, vector<double> A, vector<double> B){

    const char* filename = const_cast<char*>(path.c_str());

    std::ofstream fout;
    fout.open(filename);

    if (fout)  {
        if(A.size()!=B.size()){
            cout<<"Get wrong features!! Msg from toolsfunc: the two vector has different columns!"<<endl;
            return;
        }
        for (int i=0;i< A.size();i++){
            fout<< A[i]<<" ";
        }
        fout<<endl;
        for (int i=0;i< B.size();i++){
            fout<< B[i]<<" ";
        }
        fout<<endl;

        //cout<<"Successfully saved feature to: "<<filename<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();
}



/* Bilinear interpolate between two points:
 * the line between point1 and point2 is divided into 2^interpolationLevel pieces.
 * step: the pieces a line will be divided into;
 * the resultpoints vector size will be step+1, because step sub-line has step+1 points.
*/
void toolsfunc::splitLine(double point1, double point2, int step, vector<double> &resultpoints){
    // Clear vector
    resultpoints.clear();

    double stepdistance = (point2 - point1)/step;

    //resultpoints is an array, contain step+1 double value (step+1 points).
    for(unsigned i =0; i<step+1; i++){
        resultpoints.push_back(point1 + i*stepdistance);
    }
}




/* Split a quad into many sub quads. By bilinear interpolating to u and v seperately for each quad point.
 * Return a double vector that store all the u or v coordinate of the sub-quad.
 * p0, p1, p2, p3 is u or v coordinate of the four vertex of the quad.
 * first consider only the u coordinate, then do the same to v coordinate.
 * The bilinear interpolating works this way:
 * Firstly, divide the row line p0p3 and p1p2 each into n(n equals to 2^interpolationLevel) pieces;
 * Secondly, divide the n+1 column line into n pieces.
 * It return all the points contained in the quad in a column first order.
 */
vector<double> toolsfunc::splitQuad(double p0, double p1, double p2, double p3, int step){

    //do bilinear interpolate to the two rows.
    //for the first row, p0p3
    vector<double> firstrow, secondrow;
    splitLine(p0, p3, step, firstrow);

    //for the second row, p1p2
    splitLine(p1, p2, step, secondrow);

    //do bilinear interpolate to the step+1 coloum line, each line is divided into 2^interpolationLevel pieces.
    vector<double> columntemp;
    vector<double> quadpoints;  //store the final subquads vertexes.

    //loop each pair of points in firstrow[i] and secondrow[i].
    for(unsigned i =0; i<firstrow.size(); i++){
        //interpolate to each pair of column
        splitLine(firstrow[i], secondrow[i], step, columntemp);

        //store the column's step+1 points into vector quadpoints.
        for(unsigned j=0; j<columntemp.size(); j++){
            //cout<<"columntemp ["<<j<<"] is: "<<columntemp[j]<<endl;
            quadpoints.push_back(columntemp[j]);
        }

        //clear the vector for next loop.
        columntemp.clear();
    }

    firstrow.clear();
    secondrow.clear();

    return quadpoints;
}

/* No repeat points.
*/
vector<double> toolsfunc::splitQuad_2(double p0, double p1, double p2, double p3, int step){
    //do bilinear interpolate to the two rows.
    //for the first row, p0p3
    vector<double> firstrow, secondrow;
    splitLine(p0, p3, step, firstrow);

    //for the second row, p1p2
    splitLine(p1, p2, step, secondrow);

    //do bilinear interpolate to the step+1 coloum line, each line is divided into 2^interpolationLevel pieces.
    vector<double> columntemp;
    vector<double> quadpoints;  //store the final subquads vertexes.

    //loop each pair of points in firstrow[i] and secondrow[i].
    for(unsigned i =0; i<firstrow.size()-1; i++){
        //interpolate to each pair of column
        splitLine(firstrow[i], secondrow[i], step, columntemp);

        //store the column's step points into vector quadpoints.
        for(unsigned j=0; j< columntemp.size()-1; j++){
            //cout<<"columntemp ["<<j<<"] is: "<<columntemp[j]<<endl;
            quadpoints.push_back(columntemp[j]);
        }

        //clear the vector for next loop.
        columntemp.clear();
    }

    firstrow.clear();
    secondrow.clear();

    return quadpoints;
}


/*get filename with extension from a full path.*/
string toolsfunc::getFileNameWithExtension(string path){

    string fileName;

    int charIndex = path.find_last_of('/');
    if(charIndex != string::npos)
        fileName = path.substr(charIndex+1);

    if(fileName.empty()){
        cout<<"Message from toolsfunc::getFileNameWithExtension: Not a valid srep file.";
        EXIT_FAILURE;
    }
    else
        return fileName;
}



/*get path from a full path.*/
string toolsfunc::getPathName(string path){

    string pathName;

    int charIndex = path.find_last_of('/');
    if(charIndex != string::npos)
        pathName = path.substr(0, charIndex);

    if(pathName.empty()){
        cout<<"Message from toolsfunc::getFileNameWithExtension: Not a valid srep file.";
        EXIT_FAILURE;
    }
    else
        return pathName;
}



/*get file path without extension from a full path.*/
string toolsfunc::splitExtension(string fileName){

    int charIndex = fileName.find_last_of('.');

    fileName = fileName.substr(0, charIndex);

    if(fileName.empty()){
        cout<<"Message from toolsfunc::getFileNameWithExtension: Not a valid srep file.";
        EXIT_FAILURE;
    }
    else
        return fileName;
}



/* Get filename without extension and path from a full path.*/
string toolsfunc::getBaseFileName(string path){

    int charIndex = path.find_last_of('/');

    int lastDotIndex = path.find_last_of(".");
    int length = lastDotIndex - charIndex;

    string filename = path.substr(charIndex+1, length-1); // substr(0,4); //Get substr from index 0, length 4 char.

    if(filename.empty()){
        cout<<"Message from toolsfunc::getFileNameWithExtension: Not a valid srep file.";
        EXIT_FAILURE;
    }
    else
        return filename;
}


bool toolsfunc::quadAreaMin(Vector3D *point){
    //the first method: diagonal of p0p2
    double diagonalp0p2 = sqrt(pow(point[0].getX()-point[2].getX(), 2) + pow(point[0].getY()-point[2].getY(), 2) + pow(point[0].getZ()-point[2].getZ(), 2));
    //the second method: diagonal of p1p3
    double diagonalp1p3 = sqrt(pow(point[1].getX()-point[3].getX(), 2) + pow(point[1].getY()-point[3].getY(), 2) + pow(point[1].getZ()-point[3].getZ(), 2));

    //the four sides of each sub-quad.
    //left edge: p0p1
    double p0p1 = sqrt(pow(point[0].getX()-point[1].getX(), 2) + pow(point[0].getY()-point[1].getY(), 2) + pow(point[0].getZ()-point[1].getZ(), 2));
    //bottom edge: p1p2
    double p1p2 = sqrt(pow(point[1].getX()-point[2].getX(), 2) + pow(point[1].getY()-point[2].getY(), 2) + pow(point[1].getZ()-point[2].getZ(), 2));
    //right edge: p2p3
    double p2p3 = sqrt(pow(point[2].getX()-point[3].getX(), 2) + pow(point[2].getY()-point[3].getY(), 2) + pow(point[2].getZ()-point[3].getZ(), 2));
    //top edge: p0p3
    double p0p3 = sqrt(pow(point[0].getX()-point[3].getX(), 2) + pow(point[0].getY()-point[3].getY(), 2) + pow(point[0].getZ()-point[3].getZ(), 2));

    //calculate triangle area using Heron's formula
    double s1 = (diagonalp0p2 + p0p1 + p1p2)/2; //half of the sum of three edges.
    double s2 = (diagonalp0p2 + p0p3 + p2p3)/2;
    double s3 = (diagonalp1p3 + p0p1 + p0p3)/2;
    double s4 = (diagonalp1p3 + p1p2 + p2p3)/2;

    //two triangle's area using the first method, diagonalp0p2
    double size1 = sqrt(fabs(s1*(s1-p0p1)*(s1-p1p2)*(s1-diagonalp0p2)));
    double size2 = sqrt(fabs(s2*(s2-p0p3)*(s2-p2p3)*(s2-diagonalp0p2)));

    //two triangle's area using the second method, diagonalp1p3
    double size3 = sqrt(fabs(s3*(s3-p0p1)*(s3-p0p3)*(s3-diagonalp1p3)));
    double size4 = sqrt(fabs(s4*(s4-p1p2)*(s4-p2p3)*(s4-diagonalp1p3)));

    //the smaller one of the two method as the area for the sub-quad
    if(size1+size2 < size3+size4)
        return true;
    else
        return false;
}



/* Save matrix to .txt file.*/
void toolsfunc::save_vector_to_txt(const char* filename, vector<double> data){

    std::ofstream fout;
    fout.open(filename);

    if(fout)  {
        for(int i =0; i< data.size();i++){
            fout<< data[i]<<" ";
        }
        fout<<endl;
        cout<<"Successfully saved data to: "<<filename<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();
}


/* duration = (double)(finish - start) / CLOCKS_PER_SEC;
 *get seconds this way.  You know, total clocks divided by how many clocks in a second.
 **/
double toolsfunc::RunningTime(clock_t time1, clock_t time2){
    std::cout << "------------------CLOCKS_PER_SEC: " << CLOCKS_PER_SEC << " ms";
    double t=time1 - time2;
    double time = (t*1000)/CLOCKS_PER_SEC; // milliseconds
    return time;
}



/* Return all the u v coordinate, no repeat!!
 * interpolatedU: store the u coordinate of all the sub-quad of quad[q].
 * interpolatedV: store the v coordinate of all the sub-quad of quad[q].
*/
void toolsfunc::getUVCoordinate(vector<double> &interpolatedU, vector<double> &interpolatedV, int rowNum, int colNum, int step){
    vector<double> subquadpoint_v, subquadpoint_u;
    interpolatedU.clear();
    interpolatedV.clear();

    // Loop each quad, get the top-left subpositions
    for(unsigned int i = 0; i < rowNum -1; i++){ //the row number of the quads. its 3.
        for(unsigned int j = 0; j < colNum -1; j++){//coloums, its 13.
            //four points of a quad with delta u, v. counter clockwise.
            double quadpointsu[4], quadpointsv[4];
            //left-top point
            quadpointsu[0] = i;
            quadpointsv[0] = j;
            //left-bottom point
            quadpointsu[1] = i + 1;
            quadpointsv[1] = j;
            //right-bottom point
            quadpointsu[2] = i + 1;
            quadpointsv[2] = j + 1;
            //right-top point
            quadpointsu[3] = i;
            quadpointsv[3] = j + 1;

            //given four points of a quad (p0, p1, p2, p3), split each side of the quad into 2^interpolationLevel sub-line.
            //first, get the subquads's u coordinate, in column first order.
            subquadpoint_u = splitQuad_2(quadpointsu[0],quadpointsu[1],quadpointsu[2],quadpointsu[3], step);

            // Get the 16 subquads's v coordinate, in column first order.
            subquadpoint_v = splitQuad_2(quadpointsv[0],quadpointsv[1],quadpointsv[2],quadpointsv[3], step);

            for(unsigned int m = 0; m < subquadpoint_u.size(); m++){
                interpolatedU.push_back(subquadpoint_u[m]);
                interpolatedV.push_back(subquadpoint_v[m]);
            }

            subquadpoint_v.clear();
            subquadpoint_u.clear();
        }
    }

    // Get the bottom edge (rowNum -1) sub pieces.
    for(unsigned int j = 0; j < colNum -1; j++){//coloums, its 13.
        vector<double> subpieces;
        splitLine(j, j+1, step, subpieces);

        for(unsigned int m = 0; m < subpieces.size()-1; m++){
            interpolatedU.push_back(rowNum -1);
            interpolatedV.push_back(subpieces[m]);
        }
    }

    // Get the the right edge (colNum -1) sub pieces
    for(unsigned int i = 0; i < rowNum -1; i++){//coloums, its 13.
        vector<double> subpieces;
        splitLine(i, i+1, step, subpieces);

        for(unsigned int m = 0; m < subpieces.size()-1; m++){
            interpolatedU.push_back(subpieces[m]);
            interpolatedV.push_back(colNum -1);
        }
    }

    // Add the top-right corner.
    interpolatedU.push_back(rowNum -1);
    interpolatedV.push_back(colNum -1);
}




/* Return the spoke number of the up or down side's.
 * level = 0, there are 39 standard spoke;
 * level = 1, there are 125 standard spoke;
 * level = 2, there are 441 standard spoke;
*/
int toolsfunc::getSSSN(int interpolationLevel, int rowNum, int colNum){ // Get the standard side point number.

    int step = pow((double)2, (double)interpolationLevel);

    // Standard side
    int standQuads = (rowNum-1)*(colNum-1);
    int spokeEachQuad = step * step;
    int spokeStandSide = standQuads * spokeEachQuad + (colNum-1)*step + (rowNum-1)*step + 1;

    return spokeStandSide;
}


/* Return the spoke number of the crest side.
 * level = 0, there are 28 fold spoke;
 * level = 1, there are ;
 * level = 2, there are ;
*/
int toolsfunc::getCSSN(int interpolationLevel, int rowNum, int colNum){ // Get the standard side point number.
    // Crest atom number
    int crestQuads = rowNum*2 + colNum*2 - 4;

    // return the number of fold spokes.
    if(interpolationLevel==0)
        return crestQuads;
    else {
        int step = pow((double)2, (double)interpolationLevel);

        // Standard side
        int spokeEachQuad = (step+1)*step;
        int spokeCrestSide = crestQuads * spokeEachQuad;

        return spokeCrestSide;
    }
}



/* Save the srep points to a vtk file.
 * points: store all the points of a srep.
*/
void toolsfunc::saveVtkPointsToVTKFile(vtkSmartPointer<vtkPoints> points, string outputFile){

    // Save each new shape points to vtk file
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New(); // same as vtkPolyData* polyData;

    polyData->SetPoints(points);
    vtkSmartPointer<vtkPolyDataWriter> polyDataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    polyDataWriter->SetInput(polyData);

    polyDataWriter->SetFileName(outputFile.c_str());
    polyDataWriter->Update();
}


//void toolsfunc::write_text_to_log_file( const char* logFileName, std::string text ) {
//    std::ofstream log_file(logFileName, std::ios_base::out | std::ios_base::app );
//    log_file << text << std::end;
//}




/*compute the dot product as the cosin angle of p0p1 and p0p2*/
double toolsfunc::dotProductAngle(double * point_c, double * point_h, double * point_v){

    // vertical vector
    Vector3D a = getVector(point_c, point_v);

    // horizontal vector
    Vector3D b = getVector(point_c, point_h);

    // normalize to unit vector
    a.normalize();
    b.normalize();

    // dot product of a and b
    return a * b;
}


/*compute the dot product as the cosin angle of p0p1 and p0p2*/
double toolsfunc::dotProductAngle_2(Vector3D point_c, Vector3D point_h, Vector3D point_v){

    // vertical vector
    Vector3D a = point_h - point_c;

    // horizontal vector
    Vector3D b = point_v - point_c;

    // normalize to unit vector
    a.normalize();
    b.normalize();

    // dot product of a and b
    return a * b;
}



/*compute the dot product as the cosin angle of p0p1 and p0p2*/
Vector3D toolsfunc::trangularNormal_vec3D(Vector3D point_c, Vector3D point_h, Vector3D point_v){

    // vertical vector
    Vector3D a = point_h - point_c;

    // horizontal vector
    Vector3D b = point_v - point_c;

    return crossProduct(a, b);
}

/*compute the dot product as the cosin angle of p0p1 and p0p2*/
Vector3D toolsfunc::trangularNormal_doublePot(double * point_c, double * point_h, double * point_v){

    // vertical vector
    Vector3D a = getVector(point_c, point_v);

    // horizontal vector
    Vector3D b = getVector(point_c, point_h);

    return crossProduct(a, b);
}



Vector3D toolsfunc::crossProduct(Vector3D a, Vector3D b){
    Vector3D normalEndPoint;

    normalEndPoint.setX(a.getY()*b.getZ() - a.getZ()*b.getY());
    normalEndPoint.setY(a.getZ()*b.getX() - a.getX()*b.getZ());
    normalEndPoint.setZ(a.getX()*b.getY() - a.getY()*b.getX());

    // dot product of a and b
    return normalEndPoint;
}


//Vector3D toolsfunc::unitVector(Vector3D vecPara){

//    double length = sqrt(pow(vecPara.getX(), 2) + pow(vecPara.getY(), 2) + pow(vecPara.getZ(), 2));

//    return vecPara/length;
//}



Vector3D toolsfunc::getVector(double * point1, double * point2){
    Vector3D vec;

    vec.set(point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]);

    return vec;
}


//calculate the distance between two points.
double toolsfunc::lengthofedges(Vector3D point[2]){
    double x0 = point[0].getX();
    double y0 = point[0].getY();
    double z0 = point[0].getZ();
    //cout<<"point[0] is:"<<x0<<"; "<<y0<<"; "<<z0<<endl;

    double x1 = point[1].getX();
    double y1 = point[1].getY();
    double z1 = point[1].getZ();
    //cout<<"point[1] is:"<<x1<<"; "<<y1<<"; "<<z1<<endl;

    double distance = sqrt(pow(x1-x0, 2) + pow(y1-y0, 2) + pow(z1-z0, 2));

    //cout<<"distance between two edges is: "<<distance<<endl;

    return distance;
}
