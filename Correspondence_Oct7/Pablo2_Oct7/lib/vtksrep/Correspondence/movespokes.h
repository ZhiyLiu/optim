#ifndef MOVESPOKES_H
#define MOVESPOKES_H


#include "toolsfunc.h"
#include "M3DQuadInterpolater.h"
#include "calregularityfeatures.h"

//#include <itkMatrix.h>
#include <itkVector.h>

#include "regularityentropy.h"
#include "slidestandardspokes.h"
#include "slidecrestspokes.h"


/* Given deltau and deltav, move the (up or down) spoke u v coordinate by this step. The moved spoke information was gotten from
 * interpolation method. After got the new hub position, radius and spoke dir, we save this new values to its primitive. Here the primitive
 * index keep consistence with the old primitive's index.
*/

//using namespace std;


// For ellipsoids
/*const int srepNum = 50;     //input srep number

const int atomNums = 27;    //the total atoms on each srep.
const int varCrest = 16; // crest spoke number, each spoke has one variable.

const int varStand = 30; // variables for standard spokes (up or down spokes).
const int varNumEachSrep = varStand *2 + varCrest; // variables for current srep.
const int varNums = srepNum*varNumEachSrep;
const int crestSpokeNum = 20;*/


// For lateral ventricels
/*const int srepNum = 100;     //input srep number

const int atomNums = 39;    //the total atoms on each srep.
const int varCrest = 28; // crest spoke number, each spoke has one variable. The corner four cannot move, variables for them is set to 0.

const int varStand = 46; // variables for standard spokes (up or down spokes).
const int varNumEachSrep = varStand *2 + varCrest; // variables for current srep.
const int varNums = srepNum*varNumEachSrep;*/
//const int crestSpokeNum = 28;




typedef vnl_matrix<double> MatrixType;


class movespokes
{
public:

    typedef itk::Point< double, 3 > PointType;
    typedef vector< PointType > VectorPointType;
    typedef vector< VectorPointType > SrepVectorPointType;

    typedef vector< double > doubleVector;
    typedef vector< doubleVector > SrepDoubleVectorType;

    movespokes();
    ~movespokes();

    movespokes(string rootDir, int spokeNum, int srepNum, int side, int interpolationLevel);
    //constructor with parameters not convinience??.
    double calculateRegEntropy(const double * coeff, M3DQuadFigure* quadFig, DoubleVec subVs,
                                           int varsNum, vector<M3DQuadFigure *> quadFigList, vector<vtkSmartPointer<vtkSRep> > srepfigList);

    void saveMovedModel(string rootDir, int count,  Registry * qRegistry, string fileName, M3DQuadFigure* quadFig);
    void writeSpokeInfoToFile(const char * outputPath, Registry &registry, M3DQuadFigure* quadFig);

    void copySrepFig(M3DQuadFigure* quadFig, M3DQuadFigure* shiftingFig);

    void saveShiftedSrep(M3DQuadFigure* quadFig, string filepath, string rootDir);

    void saveGEO_Input_Matrix(string geoFileName);
    void storeSpokesToMatrix(int srepIndex, M3DQuadFigure* quadFig);

private:

    string rootDir; //The folder where input sreps stored.

    int side;
    int interpolationLevel;

    // Prepare data for geometry entropy. Save crest spokes' p r u values into pMatrix, uMatrix, rMatrix.
    MatrixType pMatrix;// storing all the spokes hub position, in sequence: up spoke, down spoke, crest spoke.
    MatrixType uMatrix;
    MatrixType rMatrix;

};

#endif // MOVESPOKES_H
