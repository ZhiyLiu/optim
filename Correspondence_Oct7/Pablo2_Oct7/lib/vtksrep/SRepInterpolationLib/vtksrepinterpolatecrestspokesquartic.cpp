#include "vtksrepinterpolatecrestspokesquartic.h"


#include "vtkSmartPointer.h"


#include "vtkCellArray.h"
#include "vtkLine.h"
#include "vtkQuad.h"

#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vnl/vnl_least_squares_function.h"
#include "vtkAppendPolyData.h"

#include "minimizecurvaturefunction.h"

#include "vnl/algo/vnl_levenberg_marquardt.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkSRepInterpolateCrestSpokesQuartic);



vtkSRepInterpolateCrestSpokesQuartic::vtkSRepInterpolateCrestSpokesQuartic()
{
    m_InterpolationLevel = 0;
    m_SRepOutput = 0;

    m_AtomId = -1;
    m_SpokeType = -1;
    m_Gamma_t = 1;
    m_Gamma_theta = 1;
    m_UseAllSpokes = false;
    m_CyclicSpokes = false;
    m_CyclicCurve = true;
    m_DoubleEndCrestAtoms = false;
    m_Inplanenterpolate = 0;
    m_P0 = VNLType(3);
    m_P0.fill(0);    
    m_SphericalCaps = false;
}

vtkSRepInterpolateCrestSpokesQuartic::~vtkSRepInterpolateCrestSpokesQuartic()
{    
    m_Spokesinterpolator.clear();
}


// Superclass method to update the pipeline
int vtkSRepInterpolateCrestSpokesQuartic::RequestData(vtkInformation* request,
                        vtkInformationVector** inputVector,
                        vtkInformationVector* outputVector){



    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkSRep *input = dynamic_cast<vtkSRep*>(vtkSRep::SafeDownCast(inInfo->Get(vtkSRep::DATA_OBJECT())));
    m_Input = input;

    m_SRepOutput = vtkSmartPointer<vtkSRep>::New();


    vtkSRep::VectorIdsType crestids;

    /*if(m_DoubleEndCrestAtoms){
        for(unsigned i = 0; i < 4; i++){
            vtkSRep::VectorIdsType crestidsside = input->GetCrestMedialAtomsIds(i);
            //crestids.push_back(crestidsside[0]);
            //crestids.push_back(crestidsside[0]);
            crestids.insert(crestids.end(), crestidsside.begin(), crestidsside.end());
        }
        //crestids.insert(crestids.begin(), crestids[0]);

        //cout<<crestids.size()<<endl;

    }else{
        crestids = input->GetCrestMedialAtomsIds();
    }*/

    crestids = input->GetCrestMedialAtomsIds();


    if(crestids.size() == 0){
        cout<<"ERROR: No crest atoms found!!!"<<endl;
        return 0;
    }

    vtkSRep::VectorVNLType crestpositions = input->GetCrestMedialAtoms(crestids);
    vtkSRep::VectorVNLType crestderivatives = input->GetCrestMedialAtomsDerivatives(crestpositions, false, m_CyclicCurve);
    vtkSRep::VectorVNLType crestnormals = input->GetMedialSheetNormals(crestids);
    vtkSRep::VectorVNLType crestderivativesnormals;

    for(unsigned i = 0; i < crestderivatives.size(); i++){


        //VNLType vect = vnl_cross_3d(crestnormals[i], crestderivatives[i]);
        //vect = vnl_cross_3d(vect, crestnormals[i]);

        VNLType vect = crestderivatives[i] - dot_product(crestderivatives[i], crestnormals[i]) * crestnormals[i];

        crestderivativesnormals.push_back(vect);
    }



    int atomPos = -1;
    if(m_AtomId!=-1){
        for(unsigned i = 0; i < crestids.size(); i++){
            if(m_AtomId == crestids[i]){
                atomPos = i;
            }
        }
        if(atomPos == -1){
            if(this->GetDebug()){
                cout<<"ERROR: Crest interpolation called for an atom that doesn't belong to the crest or end atom!!!"<<endl;
                return 0;
            }else{
                return 1;
            }
        }
    }

    m_MedialCrestCurveInterpolator = vtkSmartPointer< vtkInterpolateCurve >::New();
    m_MedialCrestCurveInterpolator->SetInterpolationLevel(m_InterpolationLevel);
    m_MedialCrestCurveInterpolator->SetCurvePoints(crestpositions);
    m_MedialCrestCurveInterpolator->SetCurveDerivatives(crestderivativesnormals);
    m_MedialCrestCurveInterpolator->SetCyclic(m_CyclicCurve);
    m_MedialCrestCurveInterpolator->Update();

    if(this->GetDebug())
        cout<<"medialcrest="<<m_MedialCrestCurveInterpolator->GetNumberOfPoints()<<endl;



    vector< vtkSRep::SPOKES_TYPE > spokestype;

    if(m_CyclicSpokes){        
        for(unsigned i = 0; i < input->GetNumberOfSpokes(); i++){
            spokestype.push_back((vtkSRep::SPOKES_TYPE)i);
        }
    }else{        
        spokestype.push_back(vtkSRep::TOP_SPOKE);
        spokestype.push_back(vtkSRep::CREST_SPOKE);
        spokestype.push_back(vtkSRep::BOTTOM_SPOKE);
    }

    m_Spokesinterpolator.clear();

    for(unsigned i = 0; i < spokestype.size(); i++){

        vtkSRep::VectorVNLType spokes = input->GetSpokes(crestids, spokestype[i], false);
        vtkSRep::VectorVNLType spokesderivatives = input->GetCrestMedialAtomsDerivatives(spokes, false, m_CyclicCurve);

        vtkSmartPointer< vtkInterpolateCurve > spokeinterpolate = vtkSmartPointer< vtkInterpolateCurve >::New();
        spokeinterpolate->SetInterpolationLevel(m_InterpolationLevel);
        spokeinterpolate->SetCurvePoints(spokes);
        spokeinterpolate->SetCurveDerivatives(spokesderivatives);
        spokeinterpolate->SetCyclic(m_CyclicCurve);
        spokeinterpolate->Update();        

        m_Spokesinterpolator.push_back(spokeinterpolate);

    }


    int crestpos = -1;
    if(m_AtomId != -1){
        for(unsigned i = 0; i < crestids.size() && crestpos == -1; i++){
            if(m_AtomId == crestids[i]){
                crestpos = i;
            }
        }
    }

    int spokepos = -1;
    if(m_SpokeType != -1){
        for(unsigned i = 0; i < spokestype.size() && spokepos == -1; i++){
            if(m_SpokeType == spokestype[i]){
                spokepos = i;
            }
        }
    }

    InterpolateCrest(spokestype, m_MedialCrestCurveInterpolator, m_Spokesinterpolator, crestpos, spokepos);

    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkSmartPointer<vtkCellArray> interpolatedcellarray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> interpolatedpoints = vtkSmartPointer<vtkPoints>::New();

    vtkSRep::VectorSRepIdsType pointsIds;

    for(unsigned i = 0; i < m_SRepOutput->GetNumberOfPoints(); i++){
        pointsIds.push_back(vtkSRep::VectorIdsType());
        double point[3];
        m_SRepOutput->GetPoint(i, point);
        vtkSRep::VNLType p0(point, 3);
        vtkSRep::VectorVNLType spokes = m_SRepOutput->GetSpokes(i);

        for(unsigned j = 0; j < spokes.size(); j++){
            vtkSRep::VNLType p;
            p = p0 + spokes[j];
            pointsIds[i].push_back(interpolatedpoints->InsertNextPoint(p[0], p[1], p[2]));
        }

    }

    for(unsigned i = 0; i < pointsIds.size() - 1; i++){
         for(unsigned j = 0; j < pointsIds[i].size() - 1; j++){

             vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
             quad->GetPointIds()->SetId(0, pointsIds[i][j]);
             quad->GetPointIds()->SetId(1, pointsIds[i+1][j]);
             quad->GetPointIds()->SetId(2, pointsIds[i+1][j+1]);
             quad->GetPointIds()->SetId(3, pointsIds[i][j+1]);

             //quad->Print(cout);

             interpolatedcellarray->InsertNextCell(quad);

         }
     }


    if(m_SphericalCaps){



    }


    //output->DeepCopy(appendpoly->GetOutput());
    output->SetPoints(interpolatedpoints);
    //output->SetLines(interpolatedcellarray);
    output->SetPolys(interpolatedcellarray);



    return 1;

}


void vtkSRepInterpolateCrestSpokesQuartic::InterpolateCrest(vector<vtkSRep::SPOKES_TYPE > spokestype,
                                                            vtkSmartPointer< vtkInterpolateCurve > medialcrestcurveinterpolator,
                                                            vector< vtkSmartPointer< vtkInterpolateCurve > > spokesinterpolator,
                                                            int crestpos, int spokepos){

    vtkIdType srepid0 = -1;
    vtkIdType srepidprev = -1;

    vtkSRep::VectorVNLType inplanespokes(spokestype.size());
    vtkSRep::VectorVNLType inplanederivativesspokes(spokestype.size());


    vtkSmartPointer<vtkCellArray> srepoutinterpolatedcellarray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> srepoutinterpolatedpoints = vtkSmartPointer<vtkPoints>::New();
    vtkSRep::RadiusVectorType srepoutallradius;
    vtkSRep::SpokesVectorType srepoutallspokes;

    if(m_Gamma_t > 1){
        m_Gamma_t = 1;
    }
    if(m_Gamma_theta > 1){
        m_Gamma_theta = 1;
    }

    double step = pow((double)2, (double)m_InterpolationLevel);
    double stepsize = 1/step;

    vtkSRep::MapUVCoordToIdType mapuvtoid;

    int istart = 0;
    int iend = medialcrestcurveinterpolator->GetNumberOfPoints();

    double t0 = 0;
    double tend = 1;

    bool modifystart = false;

    if(crestpos != -1){
        modifystart = true;
        if(crestpos == 0){
            istart = 0;
            iend = 1;
            t0 = 0;
            tend = m_Gamma_t;
            modifystart = false;
        }else if(crestpos == medialcrestcurveinterpolator->GetNumberOfPoints()){
            istart = crestpos - 1;
            iend = crestpos;
            t0 = m_Gamma_t;
            if(t0 < stepsize){
                t0 = stepsize;
            }
            tend = 1;
            modifystart = false;
        }else{
            istart = crestpos - 1;
            iend = crestpos + 1;
        }

        if(m_CyclicCurve && crestpos == medialcrestcurveinterpolator->GetNumberOfPoints()){
            istart = crestpos - 1;
            iend = crestpos;
            t0 = m_Gamma_t;
            if(t0 < stepsize){
                t0 = stepsize;
            }
            tend = 1;
            modifystart = false;
        }
    }

    for(unsigned i = istart; i < iend && i < medialcrestcurveinterpolator->GetNumberOfPoints(); i++){

        if(modifystart){
            if(i == istart){
                t0 = 1.0 - m_Gamma_t;
                if(t0 < stepsize){
                    t0 = stepsize;
                }
                tend = 1;
            }else if(i == iend - 1){
                t0 = 0;
                tend = m_Gamma_t;
                if(tend < stepsize){
                    tend = 0;
                }
            }
        }

        for(double t = t0; t <= tend; t+=stepsize){

            VNLType p0 = medialcrestcurveinterpolator->EvaluateFunction(i, t);
            //vtkIdType id0 = interpolatedpoints->InsertNextPoint(p0[0], p0[1], p0[2]);

            pair<double, double> uvcoord((double)i + t, 0);

            if(mapuvtoid.find(uvcoord) == mapuvtoid.end()){

                if(srepid0 != -1){
                    srepidprev = srepid0;
                }
                srepid0 = srepoutinterpolatedpoints->InsertNextPoint(p0[0], p0[1], p0[2]);
                mapuvtoid[uvcoord] = srepid0;

                if(srepidprev != -1){
                    vtkSmartPointer<vtkLine> srepline = vtkSmartPointer<vtkLine>::New();
                    srepline->GetPointIds()->SetId(0, srepidprev);
                    srepline->GetPointIds()->SetId(1, srepid0);
                    srepoutinterpolatedcellarray->InsertNextCell(srepline);
                }

                for(unsigned j = 0; j < spokesinterpolator.size(); j++){

                    VNLType spoke = spokesinterpolator[j]->EvaluateFunction(i, t);
                    //VNLType radius = radiusinterpolator[j]->EvaluateFunction(i, t);

                    VNLType p1 = spoke;//*radius[0];

                    inplanespokes[j] = p1;

                    /*vtkIdType id1 = interpolatedpoints->InsertNextPoint(p1[0], p1[1], p1[2]);

                    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                    line->GetPointIds()->SetId(0, id0);
                    line->GetPointIds()->SetId(1, id1);
                    interpolatedcellarray->InsertNextCell(line);*/
                }
                //Calculate the derivatives in "plane" of the three spokes
                inplanederivativesspokes = m_Input->GetCrestMedialAtomsDerivatives(inplanespokes, false, m_CyclicSpokes);

                //fit the crest curve for the plane
                vtkSmartPointer< vtkInterpolateCurve > inplaneinterpolate = vtkSmartPointer< vtkInterpolateCurve >::New();
                inplaneinterpolate->SetInterpolationLevel(m_InterpolationLevel);
                inplaneinterpolate->SetCurvePoints(inplanespokes);
                inplaneinterpolate->SetCurveDerivatives(inplanederivativesspokes);
                inplaneinterpolate->SetCyclic(m_CyclicSpokes);
                inplaneinterpolate->Update();


                vtkSRep::VectorVNLType srepoutspokes;
                vtkSRep::VectorDoubleType srepoutradius;

                vtkSRep::MapUVCoordToIdType mapuvtoidspokes;

                int jstart = 0;
                int jend = inplaneinterpolate->GetNumberOfPoints();

                double theta0 = 0;
                double thetaend = 1;
                bool modifystarttheta = false;

                if(spokepos != -1){

                    modifystarttheta = true;

                    if(spokepos == 0){
                        jstart = 0;
                        jend = 1;
                        theta0 = 0;
                        thetaend = m_Gamma_theta;
                        modifystarttheta = false;
                    }else if(spokepos == inplaneinterpolate->GetNumberOfPoints()){
                        jstart = spokepos - 1;
                        jend = spokepos;
                        theta0 = m_Gamma_theta;
                        if(theta0 < stepsize){
                            theta0 = stepsize;
                        }
                        thetaend = 1;
                        modifystarttheta = false;
                    }else{
                        jstart = spokepos - 1;
                        jend = spokepos + 1;
                    }

                    if(m_CyclicSpokes && spokepos == inplaneinterpolate->GetNumberOfPoints()){
                        jstart = spokepos - 1;
                        jend = spokepos;
                        theta0 = m_Gamma_theta;
                        if(theta0 < stepsize){
                            theta0 = stepsize;
                        }
                        thetaend = 1;
                        modifystart = false;
                    }
                }

                for(int j = jstart; j < jend; j++){

                    if(modifystarttheta){
                        if(j == jstart){
                            theta0 = 1.0 - m_Gamma_theta;
                            if(theta0 < stepsize){
                                theta0 = stepsize;
                            }
                            thetaend = 1;
                        }else if(j == jend - 1){
                            theta0 = 0;
                            thetaend = m_Gamma_theta;
                            if(thetaend < stepsize){
                                thetaend = 0;
                            }
                        }
                    }

                    for(double theta = theta0; theta <= thetaend; theta+=stepsize){

                        pair<double, double> uvcoordspoke((double)j + theta, 0);

                        if(mapuvtoidspokes.find(uvcoordspoke) == mapuvtoidspokes.end()){

                            VNLType p1 = inplaneinterpolate->EvaluateFunction(j, theta);

                            //vtkIdType id1 = interpolatedpoints->InsertNextPoint(p1[0], p1[1], p1[2]);
                            //vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                            //line->GetPointIds()->SetId(0, id0);
                            //line->GetPointIds()->SetId(1, id1);
                            //interpolatedcellarray->InsertNextCell(line);

                            //if(theta != 0 || j != 0){
                                vtkSRep::VNLType srepoutspoke = p1;// - p0;
                                srepoutradius.push_back(srepoutspoke.magnitude());
                                srepoutspoke = srepoutspoke.normalize();

                                srepoutspokes.push_back(srepoutspoke);

                                mapuvtoidspokes[uvcoordspoke] = srepoutspokes.size() - 1;
                            //}
                        }
                    }
                }

                srepoutallspokes.push_back(srepoutspokes);
                srepoutallradius.push_back(srepoutradius);

                srepoutspokes.clear();
                srepoutradius.clear();
            }
        }
    }

    m_SRepOutput->SetPoints(srepoutinterpolatedpoints);
    m_SRepOutput->SetLines(srepoutinterpolatedcellarray);
    m_SRepOutput->SetAllRadius(srepoutallradius);
    m_SRepOutput->SetAllSpokes(srepoutallspokes);
}

vtkSRep::VNLType vtkSRepInterpolateCrestSpokesQuartic::GetInterpolatedSpoke(vtkIdType cellid, double t, double v){

    vtkSRep::VectorVNLType inplanespokes;
    vtkSRep::VectorVNLType inplanederivativesspokes;

    VNLType p0 = GetInterpolatedPoint(cellid, t);
    if(m_P0 == p0){
        return GetInterpolatedSpoke(v);
    }
    m_P0 = p0;

    for(unsigned j = 0; j < m_Spokesinterpolator.size(); j++){

        VNLType spoke = m_Spokesinterpolator[j]->EvaluateFunction(cellid, t);
        VNLType p1 = spoke;//*radius[0];// + m_P0;

        inplanespokes.push_back(p1);

    }


    //Calculate the derivatives in "plane" of the three spokes
    for(unsigned j = 0; j < inplanespokes.size(); j++){
        if( j == 0){
            VNLType der = inplanespokes[j];
            if(m_CyclicSpokes){
                der = (inplanespokes[j+1] - inplanespokes[inplanespokes.size()-1])/2.0;
            }else{

                if(inplanespokes.size() == 3){

                    VNLType N = inplanespokes[j] - inplanespokes[j+2];
                    N.normalize();

                    der = inplanespokes[j+1] - inplanespokes[j];
                    der = der - dot_product(der, N) * N;
                }else{
                    der.fill(0);
                }
            }
            inplanederivativesspokes.push_back(der);
        }else if(j == inplanespokes.size() - 1){
            VNLType der = inplanespokes[j];
            if(m_CyclicSpokes){
                der = (inplanespokes[0] - inplanespokes[j-1])/2.0;
            }else{
                if(inplanespokes.size() == 3){
                    VNLType N = inplanespokes[j] - inplanespokes[j-2];
                    N.normalize();

                    der = inplanespokes[j] - inplanespokes[j-1];
                    der = der - dot_product(der, N) * N;
                }else{
                    der.fill(0);
                }
            }
            inplanederivativesspokes.push_back(der);
        }else{
            inplanederivativesspokes.push_back((inplanespokes[j+1] - inplanespokes[j-1])/2.0);
        }
    }

    //fit the crest curve for the plane
    m_Inplanenterpolate = vtkSmartPointer< vtkInterpolateCurve >::New();
    m_Inplanenterpolate->SetCurvePoints(inplanespokes);
    m_Inplanenterpolate->SetCurveDerivatives(inplanederivativesspokes);
    m_Inplanenterpolate->SetCyclic(m_CyclicSpokes);
    m_Inplanenterpolate->Update();

    return GetInterpolatedSpoke(v);
}



vtkSRep::VNLType vtkSRepInterpolateCrestSpokesQuartic::GetInterpolatedSpoke(double v){

    return m_Inplanenterpolate->EvaluateFunction(v);// - m_P0;

}

vtkSRep::VNLType vtkSRepInterpolateCrestSpokesQuartic::GetInterpolatedSpoke(double t, double v){
    //update cell with t
    double tc = t*m_MedialCrestCurveInterpolator->GetNumberOfPoints();
    double tcfl = floor(tc);
    vtkIdType coef = tcfl;
    double tcoef = tc - tcfl;

    if(coef >= m_MedialCrestCurveInterpolator->GetNumberOfPoints()){
        coef = m_MedialCrestCurveInterpolator->GetNumberOfPoints()-1;
        tcoef = 1;
    }

    return GetInterpolatedSpoke(coef, tcoef, v);
}

vtkSRep::VNLType vtkSRepInterpolateCrestSpokesQuartic::GetInterpolatedPoint(vtkIdType cellid, double t){
    return m_MedialCrestCurveInterpolator->EvaluateFunction(cellid, t);
}
vtkSRep::VNLType vtkSRepInterpolateCrestSpokesQuartic::GetInterpolatedPoint(double t){
    return m_MedialCrestCurveInterpolator->EvaluateFunction(t);
}
