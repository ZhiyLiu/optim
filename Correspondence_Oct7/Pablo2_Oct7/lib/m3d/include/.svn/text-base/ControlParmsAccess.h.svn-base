#ifndef CONTROL_PARMS_ACCESS_H
#define CONTROL_PARMS_ACCESS_H


/*
    This class provides a facility for managing user preferences,
    which are separately stored using class ControlFile.

    To use this class, the application should initialize a single global
    ControlParmsAccess object before main() is entered.  It will obtain
    any values for any parameters set in a control file, with a name
    made for the application and depending on the operating system.  See
    ControlFile.h for a description of how this is done.

    Later, the ControlParmsAccess object will be queried for settings
    used to decide among various display options in the execution of the
    program.  Calls to the various read() functions will always return
    the user-specified preferences, if they exist.  If not, the default
    value is returned.  Note that before a value can be read, the
    default should have been set.

    If any value is changed, by calling one of the write() functions,
    then the control file will be rewritten when the program exits.

    The final argument of the constructor should usually be ignored.
    If it is false, no control parameters will be set.  This is used
    for running programs from the command line and without showing
    any windows.

    A requirement of this class is that enum control_keys_t be separately
    defined.  For every control parameter name passed to the constructor,
    there must be a corresponding entry in this enum.  The ctrl_names
    array provided with to the constructor must have a 1:1 correspondence
    to the enum.

    See ControlParms.h for a class derived from this class.
*/


class ControlParmsAccess
{

public:

    // Constructor.  Application_name is the name of the program.
    // Num_ctrl_names is the length of ctrl_names, which is the list
    // of parameter names to be stored in the control file.  Ctrl_dir
    // is an alternate path where the control file may be stored.  The
    // verbosity argument controls the amount of printed output.
	ControlParmsAccess(const char * application_name, int num_ctrl_names,
            const char ** ctrl_names, const char * ctrl_dir = NULL,
            int verbosity = 0, bool load = true);
	~ControlParmsAccess();

	void setVerbosity(int verbosity);   // 0, 1, or 2

	const char * path() const;	// Name of the control file

	// Functions to set the parameter defaults.
	void setDefault(control_keys_t key, const char * deflt);
	void setDefault(control_keys_t key, bool deflt);
	void setDefault(control_keys_t key, int deflt);
	void setDefault(control_keys_t key, double deflt);

	// Functions to set the parameter defaults and simultaneously read
	// the value (returning the default, if there is no value).
	char * read(control_keys_t key, const char * deflt);
	bool read(control_keys_t key, bool deflt);
	int read(control_keys_t key, int deflt);
	double read(control_keys_t key, double deflt);

	// Functions to read the values.
	char * readString(control_keys_t key);
	bool readBool(control_keys_t key);
	int readInt(control_keys_t key);
	double readDouble(control_keys_t key);

	// Functions to change the values.
	void write(control_keys_t key, const char * val);
	void write(control_keys_t key, bool val);
	void write(control_keys_t key, int val);
	void write(control_keys_t key, double val);

	void remove(control_keys_t key);


private:

	union value_t {
		char * s;
		bool b;
		int i;
		double d;
	};

	void * control_file;	// ControlFile object

	int * changed;
	bool * loaded;
	value_t * values;
	bool * allocated;
        int num_items;
        const char ** control_names;

	ControlParmsAccess(const ControlParmsAccess &);		// Not implemented

};




#endif	/* CONTROL_PARMS_ACCESS_H */

