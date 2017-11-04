/* This class is the costfunction for itkOnePlusOneEvolutionaryOptimizer.
 * Simple example reference to : http://www.na-mic.org/svn/Slicer3-lib-mirrors/trunk/Insight/Testing/Code/Numerics/itkPowellOptimizerTest.cxx
 * We only need to write our own cost function in GetValue( const ParametersType & parameters ), here the parameters is the initial variables
 * we passed in.
 *
 * Liyun Tu
 * Apr 26, 2014.
*/


#ifndef ONEPLUSONECOSTFUNCTION_H
#define ONEPLUSONECOSTFUNCTION_H

#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include "itkNormalVariateGenerator.h"
#include <vnl/vnl_math.h>

#include "objectivefunction.h"


//using namespace std;



class oneplusonecostfunction : public itk::SingleValuedCostFunction
{
public:

  typedef oneplusonecostfunction       Self;
  typedef itk::SingleValuedCostFunction   Superclass;
  typedef itk::SmartPointer<Self>         Pointer;
  typedef itk::SmartPointer<const Self>   ConstPointer;
  itkNewMacro( Self );
  itkTypeMacro( oneplusonecostfunction, SingleValuedCostFunction );

  //enum { SpaceDimension=200 };

  typedef Superclass::ParametersType      ParametersType;
  typedef Superclass::DerivativeType      DerivativeType;
  typedef Superclass::MeasureType         MeasureType ;

  oneplusonecostfunction();
  oneplusonecostfunction(const char* rootDir, vector< std::string > inputSreps, int interpolationLevel, int regEntType);


  void GetDerivative( const ParametersType & , DerivativeType &  ) const { }


  // Return Current cost function Value, we write our cost function here.
  MeasureType  GetValue( const ParametersType & parameters ) const;

  unsigned int GetNumberOfParameters(void) const     {
      return this->totalVars;
  }

  MeasureType doIteration();

//  int callOnePlusOneOptimizer(const char* rootDir, vector< std::string > inputSreps, int interpolationLevel, int regEntType, bool initialOpt,
//                              const char* logFileName);


private:
    int totalVars;
    const char* rootDir;

//    //holding the input sreps. Reading in before iteration.
//    vector<M3DQuadFigure *> quadFigList;
//    // hoding srepfig, used for crest spoke interpolate.
//    vector<vtkSmartPointer<vtkSRep> > srepfigList;

//    int varsNum;
//    int spokeNum;
//    vector< std::string > inputSreps;
//    int interpolationLevel;
//    int side;
//    DoubleVecVec subUs;
//    DoubleVecVec subVs;

//    M3DQuadFigure* shiftingQuadFig; // as a template, storing a s-rep for shifting.


    objectivefunction objFunc;
    int srepNum;
        int varsNum;


};



#endif // ONEPLUSONECOSTFUNCTION_H




