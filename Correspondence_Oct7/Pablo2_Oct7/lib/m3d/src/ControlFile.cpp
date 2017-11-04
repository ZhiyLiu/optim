
#include <stdlib.h>
#include <iostream>
#include <string>
#include "Registry.h"
#include "ControlFile.h"

#ifdef _WIN32
#include <windows.h>
#endif

/*	Implementation Note:  This file works for UNIX, LINUX, Windows NT and
	Windows 95/98/ME.  Under 95/98/ME, it will only look for the control file
	in the directory from which the program is run.  It would be possible to
	change this to always get the file from the directory containing the
	executable, but then the globalControl object in pablo.cpp would have
	to be allocated from inside main(), so the path of the executable
	could be passed to this class.  To accomplish the same thing as the
	program is now written, the value of the WINDOWS_9x_HOME definition
	below must be changed to the specific path desired, before this file
	is compiled.
*/


// Path to global directory on Windows 95/98
#define WINDOWS_9x_HOME	"."		/* May be changed at installation time. */

// Environment variable name giving path to user's directory
#ifdef _WIN32
#define HOME_DIR	"USERPROFILE"
#else
#define HOME_DIR	"HOME"
#endif


using namespace std;




ControlFile::ControlFile(const char * application_name,
    const char * ctrl_dir, int verbosity, bool load)
{
    const char * home;

	modified = false;
	loudness = verbosity;
	ignore = ! load;

#ifdef _WIN32
    control_file_name = string(application_name);
    control_file_name.append(".txt");
#else
    control_file_name = string(".");
    control_file_name.append(application_name);
#endif

#ifdef _WIN32
    OSVERSIONINFO osvi;

    osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
    GetVersionEx(&osvi);

    if(osvi.dwPlatformId == VER_PLATFORM_WIN32_WINDOWS)
        home = WINDOWS_9x_HOME;
    else
        home = getenv(HOME_DIR);
#else
	home = getenv(HOME_DIR);
#endif

    preferences = NULL;
    if(home == NULL)
        return;

    if (ctrl_dir == NULL)
		path = home;
	else
		path = ctrl_dir;

	if (path.length() == 0) {
		preferences = NULL;
		return;
	}

	path.append("/");
    path.append(control_file_name);

	Registry * reg;
	try {
		reg = new Registry;
		if (load) {
			reg->readFromFile(path.c_str());
			if (loudness == 1)
				cout << "Preferences file: " << path << ".\n";
		}
	}
	catch (RException excp) {
		string error = excp.message();
		string open_msg("Not a valid file.");
		if (error == open_msg) {
			delete reg;
			reg = NULL;
			if (loudness >= 0)
				cout << "No preferences file was found.\n";
		}
		else {
			excp.print(cout);
		}
	}

	if (reg == NULL) {
		// No preferences file was found - create one, if possible
		try {
			reg = new Registry;
			if (! reg ->writeToFile(path.c_str()))
				cout << "Warning: preferences file may be corrupted" << endl;
		}
		catch (RException excp) {
			if (loudness >= 0)
				cout << excp.message() << '\n';
			cout << "Could not create preferences file - preferences will not be saved."
				<< endl;
		}
	}

	preferences = (void *) reg;
}


void ControlFile::setVerbosity(int verbosity)
{
	loudness = verbosity;
}


ControlFile::~ControlFile() {
	if (preferences != NULL) {
		Registry * reg = (Registry *) preferences;

		if (modified && ! ignore) {
			try {
				reg ->writeToFile(path.c_str());
				if (loudness == 1)
					cout << "Saved preferences file.\n";
			}
			catch (RException excp) {
				excp.print(cout);
			}
		}
		delete reg;
	}
}


const char * ControlFile::read(const char * key, const char * deflt) {
	Registry * reg = (Registry *) preferences;
	if (reg != NULL)
		return reg->getStringValue((char *) key, deflt);
	else return deflt;
}


bool ControlFile::read(const char * key, bool deflt) {
	Registry * reg = (Registry *) preferences;
	if (reg != NULL)
		return reg->getBooleanValue((char *) key, deflt);
	else return deflt;
}


int ControlFile::read(const char * key, int deflt) {
	Registry * reg = (Registry *) preferences;
	if (reg != NULL)
		return reg->getIntValue((char *) key, deflt);
	else return deflt;
}


double ControlFile::read(const char * key, double deflt) {
	Registry * reg = (Registry *) preferences;
	if (reg != NULL)
		return reg->getDoubleValue((char *) key, deflt);
	else return deflt;
}


bool ControlFile::write(const char * key, const char * val) {
	Registry * reg = (Registry *) preferences;
	if (reg == NULL)
		return false;
	else {
		reg->setStringValue((char *) key, val);
		modified = true;
		return true;
	}
}


bool ControlFile::write(const char * key, bool val) {
	Registry * reg = (Registry *) preferences;
	if (reg == NULL)
		return false;
	else {
		reg->setBooleanValue((char *) key, val);
		modified = true;
		return true;
	}
}


bool ControlFile::write(const char * key, int val) {
	Registry * reg = (Registry *) preferences;
	if (reg == NULL)
		return false;
	else {
		reg->setIntValue((char *) key, val);
		modified = true;
		return true;
	}
}


bool ControlFile::write(const char * key, double val) {
	Registry * reg = (Registry *) preferences;
	if (reg == NULL)
		return false;
	else {
		reg->setDoubleValue((char *) key, val);
		modified = true;
		return true;
	}
}


bool ControlFile::remove(const char * key)
{
	Registry * reg = (Registry *) preferences;
	if (reg == NULL)
		return false;
	else {
		reg->deleteKey((char *) key);
		modified = true;
		return true;
	}
}

