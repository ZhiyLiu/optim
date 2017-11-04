#include "calquantum.h"

#define PI 3.1415926535

//------------------------------------------------------------------------------------------------------------
//code added by Xiao Yang for computing the uniformity
//
//
//
//what we need is 14 eigenvalues and from that form the 14 sigmas. the 14 eigenvalues
//forms from 7 pairs:
//1. area of medial quad and upper surface;
//2. area of medial quad and lower surface;
//3. length of vertical length of medial and upper;
//4. length of vertical length of medial and lower;
//5. length of horizental length of medial and upper;
//6. length of horizental length of medial and lower;
//7. distance between points on the edge of the medial sheet and on the crest

//function used to calculate the eigenvalue of the covariance matrix
void calEigen(double matrix[][2], vector<double> data1, vector<double> data2, int datanum, double &eigen1, double &eigen2)
{
        matrix[0][0] = matrix[0][1] = matrix[1][0] = matrix[1][1] = 0;
        for(int i = 0; i < datanum; i++)
        {
                matrix[0][0] += data1[i] * data1[i] / (datanum-1);
                matrix[0][1] += data1[i] * data2[i] / (datanum-1);
                matrix[1][0] = matrix[0][1];
                matrix[1][1] += data2[i] * data2[i] / (datanum-1);
        }

        //get the eigenvalue simply by solving the equation:
        //lambda^2 - (m11+m00)*lambda + m00m11 - m01^2 = 0;
        eigen1 = ((matrix[0][0] + matrix[1][1]) + sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        eigen2 = ((matrix[0][0] + matrix[1][1]) - sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
}



double calquantum(vtkSmartPointer<vtkSRep> srepfig, vtkSmartPointer<vtkSRepInterpolateMedialSheet> medialsheetinterpolator, vtkSmartPointer<vtkSRep> srepcrest, vtkSmartPointer<vtkSRepInterpolateMedialSpokesHermite> medialspokesinterpolator)
{
//the function to calculate the quantum of s-rep features. the inputs are:
//srepfig: stores the original s-rep
//medialsheetinterpolator: stores the information for the medial sheet interpolation
//srepcrest: stores the information about the interpolated crest points.
//medialspokesinterpolator: stores the information of the spokes of the interpolated medial sheet points


        unsigned CellNum = srepfig->GetNumberOfCells();//the number of the quads(cells) in the medial sheet
        cout<<"number of the quads(cells) in the medial sheet is: "<<CellNum<<endl; // 39 atoms, has 24 quads.
        int interlevel = medialsheetinterpolator->GetInterpolationLevel();//the interpolation level
        cout<<"the interpolation level: "<<interlevel<<endl;
        vector<double> areas, upperareas, lowerareas;//calculate the quads for medial sheet and upper and lower surface

        vector<double> medialhori, medialver, upperhori, upperver, lowerhori, lowerver;//these are for recording the vertical and the horizonal side of all the subquads

        vector<double> medialsubquad, uppersubquad, lowersubquad; //these are for storing each subquads

        bool usesubquad = true; // the parameter of whether use subquad instead of quad in calculating feature 1 and 2

        double step = pow((double)2, (double)interlevel);// if interlevel==2, step==4
       // cout<<"the step is: "<<step<<endl;

        //loop for each of the 24 quads(cells)
        for(unsigned i = 0; i < CellNum; i++)
        {
                areas.push_back(0);
                upperareas.push_back(0);
                lowerareas.push_back(0);
                //calculate the area of each quad(cell)
                for(double u = 0; u <= step-1; u++)
                {
                        for(double v = 0; v <= step-1; v++)
                        {
                                vnl_vector<double> point[4];
                                //get the for points for the sub-quad
                                point[0] = medialsheetinterpolator->GetInterpolatedPoint((vtkIdType)i, u/step, v/step);
                                point[1] = medialsheetinterpolator->GetInterpolatedPoint((vtkIdType)i, (u+1)/step, v/step);
                                point[2] = medialsheetinterpolator->GetInterpolatedPoint((vtkIdType)i, (u+1)/step, (v+1)/step);
                                point[3] = medialsheetinterpolator->GetInterpolatedPoint((vtkIdType)i, u/step, (v+1)/step);

                                //get the two method for segmenting the sub-quad to two triangles
                                double midline1 = sqrt(pow(point[0][0]-point[2][0], 2) + pow(point[0][1]-point[2][1], 2) + pow(point[0][2]-point[2][2], 2));
                                double midline2 = sqrt(pow(point[1][0]-point[3][0], 2) + pow(point[1][1]-point[3][1], 2) + pow(point[1][2]-point[3][2], 2));

                                //calculate the four sides for each sub-quad. Note that when the sides are get,
                                //they will also be recorded for calculation of feature 3-6
                                double line11= sqrt(pow(point[0][0]-point[1][0], 2) + pow(point[0][1]-point[1][1], 2) + pow(point[0][2]-point[1][2], 2));
                                medialhori.push_back(line11);
                                double line12= sqrt(pow(point[0][0]-point[3][0], 2) + pow(point[0][1]-point[3][1], 2) + pow(point[0][2]-point[3][2], 2));
                                medialver.push_back(line12);
                                double line21= sqrt(pow(point[2][0]-point[1][0], 2) + pow(point[2][1]-point[1][1], 2) + pow(point[2][2]-point[1][2], 2));
                                medialver.push_back(line21);
                                double line22= sqrt(pow(point[2][0]-point[3][0], 2) + pow(point[2][1]-point[3][1], 2) + pow(point[2][2]-point[3][2], 2));
                                medialhori.push_back(line22);


                                //compute size using Heron's formula
                                double p1 = (midline1+line11+line21)/2;
                                double p2 = (midline1+line12+line22)/2;
                                double p3 = (midline2+line11+line12)/2;
                                double p4 = (midline2+line21+line22)/2;

                                //segmentation method 1, using midline 1
                                double size1 = sqrt(abs(p1*(p1-midline1)*(p1-line11)*(p1-line21)));
                                double size2 = sqrt(abs(p2*(p2-midline1)*(p2-line12)*(p2-line22)));
                                //segmentation method 2, using midline 2
                                double size3 = sqrt(abs(p3*(p3-midline2)*(p3-line11)*(p3-line12)));
                                double size4 = sqrt(abs(p4*(p4-midline2)*(p4-line21)*(p4-line22)));

                                //get the smaller one of two as the area for the sub-quad
                                areas[i] += min(size1+size2, size3+size4);
                                medialsubquad.push_back(min(size1+size2, size3+size4));

                                //compute upper side of the surface, basiclly the same procedure as calculating the medial sheet sub-quads
                                vnl_vector<double> upperpoint[4];
                                upperpoint[0] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 0, u/step, v/step);
                                upperpoint[1] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 0, (u+1)/step, v/step);
                                upperpoint[2] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 0, (u+1)/step, (v+1)/step);
                                upperpoint[3] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 0, u/step, (v+1)/step);
                                for(int pointnum = 0; pointnum < 4; pointnum++)
                                {
                                        for(int pointdim = 0; pointdim < 3; pointdim++)
                                        {
                                                //calculate the point on the upper surface by adding the
                                                //point on the medial sheet with the spoke vector
                                                upperpoint[pointnum][pointdim] += point[pointnum][pointdim];
                                        }
                                }


                                midline1 = sqrt(pow(upperpoint[0][0]-upperpoint[2][0], 2) + pow(upperpoint[0][1]-upperpoint[2][1], 2) + pow(upperpoint[0][2]-upperpoint[2][2], 2));
                                midline2 = sqrt(pow(upperpoint[1][0]-upperpoint[3][0], 2) + pow(upperpoint[1][1]-upperpoint[3][1], 2) + pow(upperpoint[1][2]-upperpoint[3][2], 2));
                                line11= sqrt(pow(upperpoint[0][0]-upperpoint[1][0], 2) + pow(upperpoint[0][1]-upperpoint[1][1], 2) + pow(upperpoint[0][2]-upperpoint[1][2], 2));
                                upperhori.push_back(line11);
                                line12= sqrt(pow(upperpoint[0][0]-upperpoint[3][0], 2) + pow(upperpoint[0][1]-upperpoint[3][1], 2) + pow(upperpoint[0][2]-upperpoint[3][2], 2));
                                upperver.push_back(line12);
                                line21= sqrt(pow(upperpoint[2][0]-upperpoint[1][0], 2) + pow(upperpoint[2][1]-upperpoint[1][1], 2) + pow(upperpoint[2][2]-upperpoint[1][2], 2));
                                upperhori.push_back(line21);
                                line22= sqrt(pow(upperpoint[2][0]-upperpoint[3][0], 2) + pow(upperpoint[2][1]-upperpoint[3][1], 2) + pow(upperpoint[2][2]-upperpoint[3][2], 2));
                                upperver.push_back(line22);

                                p1 = (midline1+line11+line21)/2;
                                p2 = (midline1+line12+line22)/2;
                                p3 = (midline2+line11+line12)/2;
                                p4 = (midline2+line21+line22)/2;

                                size1 = sqrt(abs(p1*(p1-midline1)*(p1-line11)*(p1-line21)));
                                size2 = sqrt(abs(p2*(p2-midline1)*(p2-line12)*(p2-line22)));
                                size3 = sqrt(abs(p3*(p3-midline2)*(p3-line11)*(p3-line12)));
                                size4 = sqrt(abs(p4*(p4-midline2)*(p4-line21)*(p4-line22)));

                                upperareas[i] += min(size1+size2, size3+size4);
                                uppersubquad.push_back(min(size1+size2, size3+size4));

                                //compute lower side of the surface
                                vnl_vector<double> lowerpoint[4];
                                lowerpoint[0] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 1, u/step, v/step);
                                lowerpoint[1] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 1, (u+1)/step, v/step);
                                lowerpoint[2] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 1, (u+1)/step, (v+1)/step);
                                lowerpoint[3] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 1, u/step, (v+1)/step);
                                for(int pointnum = 0; pointnum < 4; pointnum++)
                                {
                                        for(int pointdim = 0; pointdim < 3; pointdim++)
                                        {
                                                lowerpoint[pointnum][pointdim] += point[pointnum][pointdim];
                                        }
                                }
                                midline1 = sqrt(pow(lowerpoint[0][0]-lowerpoint[2][0], 2) + pow(lowerpoint[0][1]-lowerpoint[2][1], 2) + pow(lowerpoint[0][2]-lowerpoint[2][2], 2));
                                midline2 = sqrt(pow(lowerpoint[1][0]-lowerpoint[3][0], 2) + pow(lowerpoint[1][1]-lowerpoint[3][1], 2) + pow(lowerpoint[1][2]-lowerpoint[3][2], 2));
                                line11= sqrt(pow(lowerpoint[0][0]-lowerpoint[1][0], 2) + pow(lowerpoint[0][1]-lowerpoint[1][1], 2) + pow(lowerpoint[0][2]-lowerpoint[1][2], 2));
                                lowerhori.push_back(line11);
                                line12= sqrt(pow(lowerpoint[0][0]-lowerpoint[3][0], 2) + pow(lowerpoint[0][1]-lowerpoint[3][1], 2) + pow(lowerpoint[0][2]-lowerpoint[3][2], 2));
                                lowerver.push_back(line12);
                                line21= sqrt(pow(lowerpoint[2][0]-lowerpoint[1][0], 2) + pow(lowerpoint[2][1]-lowerpoint[1][1], 2) + pow(lowerpoint[2][2]-lowerpoint[1][2], 2));
                                lowerhori.push_back(line21);
                                line22= sqrt(pow(lowerpoint[2][0]-lowerpoint[3][0], 2) + pow(lowerpoint[2][1]-lowerpoint[3][1], 2) + pow(lowerpoint[2][2]-lowerpoint[3][2], 2));
                                lowerhori.push_back(line22);

                                p1 = (midline1+line11+line21)/2;
                                p2 = (midline1+line12+line22)/2;
                                p3 = (midline2+line11+line12)/2;
                                p4 = (midline2+line21+line22)/2;

                                size1 = sqrt(abs(p1*(p1-midline1)*(p1-line11)*(p1-line21)));
                                size2 = sqrt(abs(p2*(p2-midline1)*(p2-line12)*(p2-line22)));
                                size3 = sqrt(abs(p3*(p3-midline2)*(p3-line11)*(p3-line12)));
                                size4 = sqrt(abs(p4*(p4-midline2)*(p4-line21)*(p4-line22)));



                                lowerareas[i] += min(size1+size2, size3+size4);
                                lowersubquad.push_back(min(size1+size2, size3+size4));
                        }
                }
        }

        //compute the data of the crest points
        unsigned crestpointnum = srepcrest->GetNumberOfPoints();
        //cout<<"crest points number is: "<<crestpointnum<<endl; // 140

        vector<double> pointdistances;
        vector<double> crestpointdistances;

        //compute the deviation of distance between interpolated points on the spoke
        for(unsigned i = 0; i < crestpointnum; i++)
        {
                unsigned next = (i+1)%crestpointnum;
                double point[3], nextpoint[3];
                srepcrest->GetPoint(i, point);
                srepcrest->GetPoint(next, nextpoint);
                pointdistances.push_back(sqrt(pow(nextpoint[0]-point[0], 2) + pow(nextpoint[1]-point[1], 2) + pow(nextpoint[2] - point[2], 2)));
                //here we get the position of the points on the crest from three things:
                //1. position of the point on the edge of the medial sheet;
                //2. direction of the crest spoke;
                //3. length of the crest spoke
                vtkSRep::VectorVNLType currentspokes = srepcrest->GetSpokes(i);
                vtkSRep::VectorVNLType nextspokes = srepcrest->GetSpokes(next);
                vtkSRep::VectorDoubleType radius = srepcrest->GetSpokesRadius(i);
                unsigned spokenum = currentspokes.size()/2;
                //cout<<"spokenum is: "<<spokenum<<endl;//the output is 4
                vtkSRep::VNLType p = currentspokes[spokenum] * radius[spokenum];
                for(int num = 0; num < 3; num++)
                {
                        point[num] += p[num];
                }

                radius = srepcrest->GetSpokesRadius(next);
                p = nextspokes[spokenum] * radius[spokenum];
                for(int num = 0; num < 3; num++)
                {
                        nextpoint[num] += p[num];
                }
                crestpointdistances.push_back(sqrt(pow(nextpoint[0]-point[0], 2) + pow(nextpoint[1]-point[1], 2) + pow(nextpoint[2] - point[2], 2)));
        }




        //substract the mean from the data
        double medialsheetmean = 0;
        double uppermean = 0;
        double lowermean = 0;
        for(int i = 0; i < CellNum; i++)
        {
                medialsheetmean += areas[i]/CellNum;
                uppermean += upperareas[i]/CellNum;
                lowermean += lowerareas[i]/CellNum;
        }
        for(int i = 0; i < CellNum; i++)
        {
                areas[i] -=medialsheetmean;
                upperareas[i] -= uppermean;
                lowerareas[i] -= lowermean;
        }

        double medialsubquadmean = 0;
        double uppersubquadmean = 0;
        double lowersubquadmean = 0;
        int subquadnum = medialsubquad.size();
        for(int i = 0; i < subquadnum; i++)
        {
                medialsubquadmean += medialsubquad[i]/subquadnum;
                uppersubquadmean += uppersubquad[i]/subquadnum;
                lowersubquadmean += lowersubquad[i]/subquadnum;
        }
        for(int i = 0; i < subquadnum; i++)
        {
                medialsubquad[i] -=medialsubquadmean;
                uppersubquad[i] -= uppersubquadmean;
                lowersubquad[i] -= lowersubquadmean;
        }


        int linesize = medialver.size();
        int medialvermean = 0, medialhorimean = 0, uppervermean = 0, upperhorimean = 0, lowervermean = 0, lowerhorimean = 0;
        for(int i = 0; i < linesize; i++)
        {
                medialvermean += medialver[i];
                medialhorimean += medialhori[i];
                uppervermean += upperver[i];
                upperhorimean += upperhori[i];
                lowervermean += lowerver[i];
                lowerhorimean += lowerhori[i];
        }

        medialvermean /= linesize;
        medialhorimean /= linesize;
        uppervermean /= linesize;
        upperhorimean /= linesize;
        lowervermean /= linesize;
        lowerhorimean /= linesize;

        for(int i = 0; i < linesize; i++)
        {
                medialver[i] -= medialvermean;
                medialhori[i] -= medialhorimean;
                upperver[i] -= uppervermean;
                upperhori[i] -= upperhorimean;
                lowerver[i] -= lowervermean;
                lowerhori[i] -= lowerhorimean;
        }



        double pointsmean = 0;
        double crestpointsmean= 0;

        for(unsigned i = 0; i < crestpointnum; i++)
        {
                pointsmean += pointdistances[i];
                crestpointsmean += crestpointdistances[i];
        }
        pointsmean /= crestpointnum;
        crestpointsmean /= crestpointnum;

        for(unsigned i = 0; i < crestpointnum; i++)
        {
                pointdistances[i] -= pointsmean;
                crestpointdistances[i] -= crestpointsmean;
        }



        //form the matrix adn calculate the eigenvalues
        double sigma[14];
        double matrix[2][2];
        if(usesubquad)
        {
                calEigen(matrix, medialsubquad, uppersubquad, subquadnum, sigma[0], sigma[1]);
                calEigen(matrix, medialsubquad, lowersubquad, subquadnum, sigma[2], sigma[3]);
        }
        else
        {
                calEigen(matrix, areas, upperareas, CellNum, sigma[0], sigma[1]);
                calEigen(matrix, areas, lowerareas, CellNum, sigma[2], sigma[3]);
        }



        //calculate feature 3
        calEigen(matrix, medialver, upperver, linesize, sigma[4], sigma[5]);

        //calculate feature 4
        calEigen(matrix, medialver, lowerver, linesize, sigma[6], sigma[7]);


        //calculate feature 5
        calEigen(matrix, medialhori, upperhori, linesize, sigma[8], sigma[9]);

        //calculate feature 6
        calEigen(matrix, medialhori, lowerhori, linesize, sigma[10], sigma[11]);


        //calculate feature 7
        calEigen(matrix, pointdistances, crestpointdistances, crestpointnum, sigma[12], sigma[13]);

        //sum up the result to get the entropy
        double feature_entropy = 0;

        for(int i = 0; i < 14; i++)
        {
            cout<<"The eigen values is: "<<sigma[i]<<endl;
                feature_entropy += log(sigma[i]);
        }
        feature_entropy += 14 * (0.5 + 0.5 * log(PI));

/* print each entrophy. not necessary
 *        for(int i = 0; i < 7; i++)
 *        {
 *                cout << 2*0.5 + 2*0.5*log(PI) + log(sigma[i*2]) + log(sigma[i*2+1]) << endl;
 *        }
 *
 */
        return feature_entropy;
}
