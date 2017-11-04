#include "objectivefunction.h"

objectivefunction::objectivefunction()
{
}


objectivefunction::objectivefunction(const char* rootDir, int side, int varsNum, int spokeNum, M3DQuadFigure* shiftingQuadFig,
                                     int interpolationLevel, vector< std::string > inputSreps){
    this->varsNum = varsNum;
    this->spokeNum = spokeNum;
    this->rootDir = rootDir;    
    this->interpolationLevel = interpolationLevel;
    this->side = side;    
    this->shiftingQuadFig = shiftingQuadFig;

    toolsfunc tools;
    //Initialized the srep quadFig list
    for(unsigned int i = 0; i < inputSreps.size(); i++){
        //Initialize quadFigList, which holding the original sreps, the changing during each iteration will based on this.
        std::cout<<"---------Reading:"<<inputSreps[i]<<std::endl;
        M3DQuadFigure* qFig = tools.GetQuadFigure(inputSreps[i].c_str());
        this->quadFigList.push_back(qFig);

        // Initialize srepfig.
        vtkSmartPointer<vtkReadSRep> readsrep = vtkSmartPointer<vtkReadSRep>::New();
        readsrep->SetFileName(inputSreps[i].c_str());
        readsrep->Update();
        vtkSmartPointer<vtkSRep> srepfig = readsrep->GetOutput();
        this->srepfigList.push_back(srepfig);
    }

    uvmap mp;
    // holding the 0~1 subdivision of the line.
    this->subVs = mp.getUVSubdivision(interpolationLevel);
}




/* Compute the total entropy. */
double objectivefunction::objectiveFunction(const double * wholeCoeff, double w1, double w2) const{
    //Step 1: Move spokes, compute the regularity entropy and save geometry matrix to .txt.
    movespokes mSpokes(rootDir, spokeNum, this->quadFigList.size(), side, this->interpolationLevel);
    double regEntropy = mSpokes.calculateRegEntropy(wholeCoeff,this->shiftingQuadFig, subVs, varsNum, quadFigList, srepfigList);

    toolsfunc tools;

    //Step 2: Call matlab script to compute geometry entropy.
    //matlab -nodesktop -nodisplay -nojvm -nosplash -r "dataDir='/work/ltu/WorkSpace/Sep_3_newFeature/downSide/temp_data/';geoEntropy"
    string s = string("matlab -nodesktop -nodisplay -nojvm -nosplash -singleCompThread -r ") + '"' + string("dataDir='") + tools.connectStr(rootDir, "/temp_data/");

    string script_str;
    if(this->side == 0)
        script_str = s + string("';upSide") + '"';
    else if(this->side == 1)
        script_str = s + string("';downSide") + '"';
    else
        script_str = s + string("';crestSide") + '"';

    int returnValue1 = system(script_str.c_str());
    if (returnValue1 != 0 ) std::cout << "something has happened in Matlab, following value returned: " << returnValue1 << std::endl;


    //Step 3: Read the geo entropy from .txt
    string geoEntropyResultPath = tools.connectStr(tools.connectStr(rootDir, "/temp_data/"), "geoEntropyResult.txt");
    ifstream myfile(geoEntropyResultPath.c_str());
    string line;
    double geoEntropy = 0;
    if (myfile.is_open()) {
        while ( getline (myfile,line) )
        {
            geoEntropy = atof(line.c_str());
            std::cout <<"----geoEntropy: "<<geoEntropy<<std::endl;
        }
        myfile.close();
    }

    else std::cout << "Unable to open file: " << geoEntropyResultPath << std::endl;

    double objectFuc = w1*geoEntropy - w2*regEntropy;
    std::cout <<"----objective function: "<<objectFuc<<std::endl;

    return objectFuc;
}


int objectivefunction::getIterationCount(int iteCounter) const{

    return iteCounter+1;
}


