#ifndef CONTROL_PARMS_H
#define CONTROL_PARMS_H



// This include file defines enum control_keys_t, the application
// name, the control names array and its length.  These are in a
// separate file for easy replacement by other applications built
// on the Pablo libraries.
#include "ControlKeys.h"

// The class from which ControlParms is derived.
#include "ControlParmsAccess.h"


/*
    ControlParms is a light-weight, inline class provided for using
    class ControlParmsAccess in a particular application.  It
    eliminates the need to pass the application name and array
    of control names to ControlParmsAccess.  Also it provides
    a shorter type name in the application's code.
*/
class ControlParms : public ControlParmsAccess
{

    public:

        ControlParms(const char * ctrl_dir = NULL, int verbosity = 0,
            bool load = true)
            : ControlParmsAccess(APPLICATION_NAME, NUM_CONTROL_ITEMS, CONTROL_NAMES,
            ctrl_dir, verbosity, load)
                { }
        ~ControlParms() { }

    private:

        ControlParms(const ControlParms &);	// Not implemented

};



// The global ControlParms object
extern ControlParms * globalControl;	// See pablo.cpp

// The current verbosity level of Pablo
extern int globalVerbosity;		// See pablo.cpp





#endif	/* CONTROL_PARMS_H */

