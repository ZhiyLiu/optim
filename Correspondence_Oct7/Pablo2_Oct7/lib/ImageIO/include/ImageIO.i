/* vim: set syn=cpp :
 *
 * Swig wrapper for the ImageIO library
 * - Rohit Saboo
 */
%{
#include "ControlKeys.h"
#include "control_parms_defaults.h"
#include "ControlParms.h"

// ImageIO includes
#include "AllImageIO.h"
%}

%include swig-helpers.i

%include AllImageIO.h
