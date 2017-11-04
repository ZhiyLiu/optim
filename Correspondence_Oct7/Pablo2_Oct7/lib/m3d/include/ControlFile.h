#ifndef CONTROL_FILE_H
#define CONTROL_FILE_H


/*
    This class provides a facility for storing user preferences in a
    "control file".  For application X, this file will be named ".X"
    on UNIX or LINUX systems and "X.txt" on Windows 95/98/NT/2000.
    Under UNIX or LINUX, this file should be in the user's home
    directory.  Under Windows NT/2000, it should be in the user's
    Profile directory.  Under Windows 95/98, a file common to all
    users may be designated by changing the definition of
    WINDOWS_9x_HOME in ../src/ControlFile.cpp.  Otherwise, only the
    current directory will be searched.

    In all cases, if the constructor is called with an argument
    containing the path of a directory, then the control file will
    be expected to be in that directory.

    A global ControlFile object will be initialized before main() is
    entered (see ControlParmsAccess.h).  It will read the control
    file, if present.  Later, this object will be queried for settings
    used to decide among various display options in the execution of
    the program.

    Calls to the various read() functions will always return the user-
    specified preferences, if they exist.  If not, the default value is
    returned.

    The final argument of the constructor should usually be ignored.
    If it is false, the control file will not be loaded.  This is
    used for running a program from the command line and without
    showing any windows.

*/

class ControlFile
{

public:

	ControlFile(const char * application_name,
            const char * ctrl_dir = NULL, int verbosity = 0,
            bool load = true);
	~ControlFile();

	void setVerbosity(int verbosity);

	const char * controlFilename() { return path.c_str(); }

	const char * read(const char * key, const char * deflt);
	bool read(const char * key, bool deflt);
	int read(const char * key, int deflt);
	double read(const char * key, double deflt);

	bool write(const char * key, const char * val);
	bool write(const char * key, bool val);
	bool write(const char * key, int val);
	bool write(const char * key, double val);

	bool remove(const char * key);


private:

	void * preferences;
	std::string path;
	bool modified;
	bool ignore;
	int loudness;
	std::string control_file_name;

	ControlFile(const ControlFile &);	// Not implemented

};




#endif	/* CONTROL_FILE_H */

