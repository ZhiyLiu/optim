// Dispatch to an application's main function,
//  which parses argv/argc. and print resultant return code.

#ifdef _WIN32
#include <process.h>
#else
#include <sys/types.h>
#include <unistd.h>
#endif
#include <iostream>
#include <fstream>
#include "main.h"

using namespace std;


extern bool getarg(char * shortOption, char * longOption, bool takesArg,
	string & opt, int & argNum, bool & error, bool verbose, int argc,
	char **argv);

bool saveInfo(const char * path, int info)
{
    std::ofstream outFile;
	outFile.open(path, ios::out & ios::trunc);
	if (outFile.is_open()) {
	    outFile << info << endl;
	    outFile.close();
		return true;
	}
	else {
        cout << "Error: unable to open file " << path << endl;
		return false;
	}
}

// Print return code from app.
// The app may call exit() only for critical errors, not for normal
//  failure returns, so that the string "exit status" can be used
//  to determine if the app crashed or not.
// Other implementations of main() are possible which never call
//  exit(), such as a pablo server or an alternate GUI.
int main(int argc, char *argv[])
{
    bool error;
    bool verbose;

	error = false;
	verbose = false;
	for (int i = 1; i < argc; i++) {
		string opt = argv[i];

#ifdef _WIN32
		if (opt[0] != '-' && opt[0] != '/')
#else
		if (opt[0] != '-')
#endif
			continue;

		if (getarg("lv", "verbose", false, opt, i, error, false, argc, argv))
		    verbose = true;

		// Undocumented option for use by ConStruct or external scripts.
		// This causes the process id to be written into the designated
		// file when the program starts up.  This can be used later to
		// kill the process.
		if (getarg("tp", "pid", true, opt, i, error, false, argc, argv)) {
			if (error) {
				cout << "Option -tp (-pid) requires a filename argument\n";
				return -1;
			}
#ifdef _WIN32
			if (saveInfo(argv[i], _getpid()))
#else
			if (saveInfo(argv[i], getpid()))
#endif
			{
                if (verbose)
				    cout << "Pid file: " << argv[i] << endl;
				break;
			}
			else
				return -1;	// No pid file was created
		}
	}

	int exitCode = pablo(argc, argv);

	// Collect all possible exits; guarantee stdout is flushed
	if (exitCode != 0)
		cout << "Exit status: " << exitCode << endl;
#ifndef AE2_BUILD
	else
		cout << "Pablo finished." << endl;
#endif




	if (touchFile != NULL)
		(void) saveInfo(touchFile, exitCode);

	// dibyendu
	
	//cin.get() ;


	return exitCode;
}

