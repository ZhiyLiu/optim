/* vim: set syn=cpp :
 *
 * Swig wrapper for the optima library
 * - Rohit Saboo
 */
%{
#include "ControlKeys.h"
#include "control_parms_defaults.h"
#include "ControlParms.h"

// paul_code includes
#include "libs.h"
#include "support.h"
#include "optima.h"
#include "BrentLinearMethod.h"
#include "ConjugateGradientMethod.h"
#include "EvolutionaryStrategy.h"
#include "SimplexMethod.h"
#include "LogManager.h"
#include "Registry.h"
#include "matrix.h"
#include "problems.h"
%}

%include swig-helpers.i

%include libs.h
%include support.h
%include optima.h
%include BrentLinearMethod.h
%include ConjugateGradientMethod.h
%include EvolutionaryStrategy.h
%include SimplexMethod.h
%include LogManager.h
%include Registry.h
%include matrix.h
%include problems.h
