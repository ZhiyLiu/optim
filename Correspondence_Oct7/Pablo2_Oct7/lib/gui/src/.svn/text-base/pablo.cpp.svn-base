#include <string>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#ifdef WIN32
#include <new.h>
#endif
#include <FL/gl.h>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Tooltip.H>
#include "make_windows.h"
#include "P3DUserInterface.h"
#include "ControlParms.h"
#include "Registry.h"
#include "control_parms_defaults.h"
#include "pablo_version.h"
#include "globalBuildParms.h"
#include "Tuning.h"
#ifdef OPTIMIZATION_VISUALIZER
#include "OptVisualizerUI.h"
#endif

#include "BrentLinearMethod.h"


//#define DEBUG
//#define CHECK_DOMAIN	/* Uncomment to use user-domain checking in Binary Pablo */

#ifdef BINARY
// Whether or not a license key is required to run Binary Pablo
const bool keyRequired = false;    // Change to true for all distributed copies
#endif

// This is the official release number of the last tagged release (last major
// release).  It is recorded in file releases, at the top of the Pablo2 SVN tree.
const char * release = "Pablo-1-4-0";	// The latest number from file 'releases'

// The copyright must also be changed in initAboutPabloWindow().
const char * copyright = "(c) 2002-2009 MIDAG, UNC-CH.\n\
 Licensed for research and educational use only.";

#define DELAY	0.001
#define LONG_DELAY	0.1


using namespace std;


#ifdef BINARY
#ifdef CHECK_DOMAIN

// List of domains where Binary Pablo may run without a key.
// To have a name permanently added to this list, send email
// to gash@cs.unc.edu.
char * domains[] = {
	"UNC-COMPSCI",
	"RAD-ONC-PC",
	NULL	// Required terminator
};
#endif

// License key management
bool require_key;		// Is a key required to run Binary Pablo?
string thisLicVersion("app=Binary-Pablo-License,ver=1.03");
bool pk;

// Filename containing key info when require_key is true
#ifdef _WIN32
string keyFilename = "c:\\bpkey.dat";
#else
string keyFilename = "bpkey.dat";
#endif

// For generating a key from components
unsigned long magicNumber = 8675309;
unsigned long magicOffset = 5;

#endif	/* BINARY */

// To have Binary Pablo load a script automatically add these:
// project->settings->[GUI project]->debug->program arguments
//  -po -cs BPscript.txt
// project->settings->[GUI project]->debug->working dir
//  [dir containing BPscript.txt]

// New parameter and tuning mechanism "tuneVals":
// 1) tuneVals is a global class containing an array of weighting and
//    other configuration DOUBLEs. It's initialized to default values,
//    also kept in the class, which are used only if no other info
//    is supplied by the user, so you can depend on the weights being
//    set at any time (they're set before main() is invoked).
//    Each value is named to associate it with named values in files.
// 2) command-line argv is parsed:
//    2a) the '-cs' switch spec's the 'control script', a registry-style
//        text file containing values that are deposited into tuneVals
//        as well as image, model, etc that are kept outside tuneVals.
//    2b) the other argv switches are values that are likely to be
//        change at run-time, eg, model, image, etc. These will override
//        values set in step (2a).
// 3) most lower-level routines have been rewritten to read tuneVals
//    directly instead of getting their values through a Set* accessor or
//    constructor, and to write their match values into tuneVals directly.
//    A few routines still pass values.
// 4) print tuneVals match values.
// 5) routines exist to write tuneVals out as a registry-style file,
//    so that it may be used in a future run. This includes the command-line
//    changes applied to a tuneVals. You can only get to this by recompiling,
//    for now.
// 6) the old landmark routines were disabled but are in the code. They should
//    have been cleaned out to prevent confusion, but are intimately integrated
//    into the GUI and so this was saved for LATER.

// PABLO VERSION
// The previous "versioning capability" identifies the date of the last cvs checkout,
// but we are not currently checking binary pablo into cvs.
// Instead, the "born-on" time stamp is represented by a float: YYYYMMDD.HH which tells
// the feature set as accepted by the tuning param's.
// Since binary pablo is to be released infrequently, the rule for matching
// a registry's version against the program's version is that
// the date has to match (YYYYMMDD). The hour portion (HH) is for emergency
// purposes, not new features, so it's the same version.

// FEATURES DOCUMENTED ONLY IN PROG GUIDE
// the following are left undocumented for internal UNC use:
// -id  distance map image: -po will read this if present, write it if absent.
//      it contains the unnormalized 3-D distance map to the nearest image edge.
// -pt  prints tab-delimited control script key usage for direct insertion into
//      the user manual
// -pk  generates a keyFile. see the current magicNumber, which tells the program to
//      write over the keyfile instead of read it. this is only used to
//      generate key files at a site, for EMERGENCY OR INTERNAL USE ONLY.
//      UNIMPLEMENTED: user can specify the keyfile's name.


P3DUserInterface * ui;
#ifdef OPTIMIZATION_VISUALIZER
OptVisualizerUI * ovUI;
#endif


bool all_iconified = false;
bool windowing = false;
bool startUpWait = true;
Fl_Window ** icons = NULL;
int iconCount = 0;
int errorCount = 10;
char * headerString;
extern char * longChangesString;
bool constructExec;
char * touchFile;

extern void deleteUserInterface(P3DUserInterface * userInterface);


#ifndef PRODUCTION_VERSION
#ifdef _WIN32
int _matherr(_exception * excp)
#else
int _matherr(struct ::exception * excp)
#endif
{
	if (--errorCount > 0) {
		cout << "Math error: trapped in pablo.cpp" << endl;
		if (errorCount == 1)
			cout << "Suppressing further occurrences of this message" << endl;
	}
	return 0;	// continue execution
}
#endif	/* PRODUCTION_VERSION */


#ifdef _MSC_VER
int newFailureTrap(size_t size)
{
	cout << "Error: New() failed to allocate the requested space of " << size << " bytes\n";
	abort();
	return 0;
}
#else	/*  _MSC_VER */
void newFailureTrap()
{
	cout << "Error: New() failed to allocate the requested space; ";
	cout << "terminating Pablo\n";
	abort();
}
#endif	/*  _MSC_VER */


int globalHandler(int event)
{
	switch (event)
	{
		case FL_SHORTCUT:
			return ui->mainWindow->handle(event);
			break;

		default:
			break;
	};

	return 0;
}


void iconifyAllWindows()
{
	Fl_Window * window;
	int i;

	if (all_iconified)
		return;

	// Get list of windows to be iconified.  The extra looping is
	// necessary to avoid multiple calls to iconize() per window.
	window = Fl::first_window();
	while (window != NULL) {
		iconCount++;
		window = Fl::next_window(window);
	}
	icons = new Fl_Window *[iconCount];
	window = Fl::first_window();
	for (i = 0; i < iconCount; i++) {
		icons[i] = window;
		window = Fl::next_window(window);
	}

	// Iconify the list
	for (i = 0; i < iconCount; i++) {
		window = icons[i];
		window->iconize();
	}

	all_iconified = true;
}


void deiconifyAllWindows()
{
	if (! all_iconified)
		return;

	for (int i = 0; i < iconCount; i++)
		icons[i]->show();

	delete [] icons;
	icons = NULL;
	iconCount = 0;
	all_iconified = false;
	return;
}


void check_deiconify(); 	// Forward declaration

void check_iconify()
{
	Fl_Window ** windows;
	bool have_icon, restore;
	int n;

	windowing = true;
	ui->mainWindow->setHideCallback(NULL);
	restore = false;

	switch (globalControl->readInt(IconifyMode))
	{
		case 0: 
			restore = true;
			break;
		case 1: 
			if (! ui->mainWindow->visible())	// The Main Window is always shown
				iconifyAllWindows();
			else
				restore = true;
			break;
		case 2: 
			windows = ui->mainWindow->madeWindows(n);
			have_icon = false;
			for (int i = 0; i < n; i++) {
				if (windows[i]->shown() && ! windows[i]->visible()) {
					have_icon = true;
					break;
				}
			}
			if (have_icon)
				iconifyAllWindows();
			else
				restore = true;
			break;
	}
	if (restore)
		ui->mainWindow->setHideCallback(check_iconify);
	else
		ui->mainWindow->setShowCallback(check_deiconify);
	windowing = false;
}


void check_deiconify()
{
	bool restore;
	int i;

	windowing = true;
	ui->mainWindow->setShowCallback(NULL);
	restore = false;

	switch (globalControl->readInt(IconifyMode))
	{
		case 0: 
			restore = true;
			break;
		case 1: 
			if (ui->mainWindow->visible())	// The Main Window is always shown
				deiconifyAllWindows();
			else {
				for (i = 0; i < iconCount; i++) {
					if (! (icons[i]->shown() && ! icons[i]->visible())) {
						icons[i] = icons[iconCount - 1];
						iconCount--;
						if (iconCount == 0) {
							delete [] icons;
							icons = NULL;
							all_iconified = false;
						}
						else
							restore = true;
						break;
					}
				}
			}
			break;
		case 2: 
			deiconifyAllWindows();
			break;
	}
	if (restore)
		ui->mainWindow->setShowCallback(check_deiconify);
	else
		ui->mainWindow->setHideCallback(check_iconify);
	windowing = false;
}


void timeoutCallback(void *)
{
	if (ui->callback->exiting())
		return;

	if (windowing) {
		// Just a precaution
		Fl::repeat_timeout(DELAY, timeoutCallback);
		return;
	}

	if (all_iconified)
		Fl::repeat_timeout(LONG_DELAY, timeoutCallback);
	else
		Fl::repeat_timeout(DELAY, timeoutCallback);

	if (ui->callback->imagePlaneSliderFinishedMoving()) {
//        cout << "+\n" << flush;
		ui->modelView->updateImagePlanes(true);
		return;
	}

	if (ui->modelView->rockingMode()) {
		ui->modelView->redraw();
		return;
	}

	if (ui->modelView->screenDumpScheduled > 0 ) {
		if (--ui->modelView->screenDumpScheduled == 0)
			ui->modelView->screenDump();
		else
			Fl::repeat_timeout(LONG_DELAY, timeoutCallback);
	}

	if (startUpWait) {
		if (! ui->modelWindow->shown()) {
			Fl::repeat_timeout(DELAY, timeoutCallback);
			return;
		}

		// Changing the width forces the resize to propagate to the
		// groups the window contains.  Merely calling redraw() has
		// no effect.
		ui->modelWindow->resize(ui->modelWindow->x(), ui->modelWindow->y(),
			ui->modelWindow->w() - 1, ui->modelWindow->h());
		ui->modelWindow->resize(ui->modelWindow->x(), ui->modelWindow->y(),
			ui->modelWindow->w() + 1, ui->modelWindow->h());

		startUpWait = false;
	}
}


void triggerModelWindowRedraw()
{	
	Fl::repeat_timeout(DELAY, timeoutCallback);
	startUpWait = true;
}


static void printInputsAsTDTable()		// For those who don't know, TD apparently stands for "tag - data"
{
#ifdef BINARY
	cout << "image = ;" << "\t" << "-ii" << "\t"
		<< "[raw3 FILE] binary image: contains intensities 0 and 1 in cubic voxels" << endl;
	cout << "model = ;" << "\t" << "-im" << "\t"
		<< "[m3d FILE] shape model to be fit to binary image" << endl;
	cout << "landmarkModel = ;" << "\t" << "-il" << "\t"
		<< "[m3d FILE] landmark model corresponding to binary image" << endl;
	cout << "outModel = ;" << "\t" << "-om" << "\t"
		<< "[m3d FILE] fitted model" << endl;
	cout << "outTile = ;" << "\t" << "-ot" << "\t"
		<< "[BYU FILE] surface of fitted model. See tileSurfaceLevel and tileQuads" << endl;
	cout << "outImageTile = ;" << "\t" << "-oit" << "\t"
		<< "[BYU FILE] simple marching cubes tiling of binary image.  Always uses Quads" << endl;
	cout << "outModelImage = ;" << "\t" << "-omi" << "\t"
		<< "[RAW3 FILE] scan-convert model to image" << endl;
	cout << "outSimTransModel = ;" << "\t" << "-ostm" << "\t"
		<< "[m3d FILE] similarity transform model after init or figure stage(s).  See doMethoOfMoments and figureIterations" << endl;

	// undocumented
	//	cout << "distanceMap = ;" << "\t" << "-id" << "\t"
	//		<< "[raw3 FILE] distance map image. Read if exists; written if not" << endl;
#else
	cout << "image = ;" << "\t" << "-ii" << "\t"
		<< "[raw3 FILE] gray/binary image: optimize without distance map" << endl;
	cout << "model = ;" << "\t" << "-im" << "\t"
		<< "[m3d FILE] shape model to be fit to image" << endl;
	cout << "landmarkModel = ;" << "\t" << "-il" << "\t"
		<< "[m3d FILE] landmark model corresponding to image positions and tolerances of such" << endl;
	cout << "outModel = ;" << "\t" << "-om" << "\t"
		<< "[m3d FILE] fitted model" << endl;
	cout << "outTile = ;" << "\t" << "-ot" << "\t"
		<< "[BYU FILE] surface of fitted model.  See tileSurfaceLevel and tileQuads" << endl;
	cout << "outModelImage = ;" << "\t" << "-omi" << "\t"
		<< "[RAW3 FILE] scan-convert model to image" << endl;
#endif
}

#ifdef BINARY

// license and key management
// - findValidKey is from immap::main.cxx
// - changed keyfile name from ivkey.dat to bpkey.dat

// ret true if token is found on seps-delimited string searchee
bool findStringDelimited(const char * searchee, const char * token,
						 const char * seps, bool debug = false)
{
	if (! searchee || !* searchee || ! token || ! *token || ! seps || ! *seps)
		return false;

	// strtok alters the string, so work on a duplicate
	char * searcheeCopy = strdup(searchee);

	if (debug)
		cout << " List is \"" << searchee << "\"" << endl;
	char * thisToken = strtok(searcheeCopy, seps);
	while (thisToken)
	{
		if (debug)
			cout << " Checking (\"" << thisToken << "\" == \"" << token << "\")" << endl;
		if (! strcmp(thisToken, token))
		{
			if (debug)
				cout << " Found" << endl;
			free(searcheeCopy);
			return true;
		}
		thisToken = strtok(NULL, seps);
	}
	if (debug)
		cout << " Not found " << endl;
	if (searcheeCopy)
		free(searcheeCopy);
	return false;
}


#ifdef _WIN32
#include <winsock.h>
#else
#include <unistd.h>
#endif
#include <iostream>
#include <fstream>
#include <ctime>

// acquires and/or produces a valid key for the computer,
//  to allow program to run
//
// RETURN
// 0: key approved (pass)
// 1: key denied, reason printed to cerr (fail)
// else key denied, returns windows status code from WSAStartup
//
// METHOD
// if (keyFile not found)
//  if (secret -pk switch used)
//   allow user to type keyFile info (key becomes a password)
//  else fail
// else
//  read version from keyFile
//  fail if version is not THIS version
//  read key, hostList, expireDate from keyFile
//  calc if host is on hostList (but don't act on it yet)
// if key is magicNumber, write file and pass
// calc productKey from hostList, expireDate
// STANDARD TESTS:
// fail if host is not on hostlist
// fail if (productKey != key || expireDate < today)
// warn if (expireDate is within a month of today)
// pass
//

int findValidKey(bool pk, const char * cd_arg, bool debug)
{
#ifdef _WIN32	/* Required gethostname() to work on windows */
	WORD wVersionRequested;
	WSADATA wsaData;
	int err;

	wVersionRequested = MAKEWORD(2, 2);

	err = (int) WSAStartup(wVersionRequested, &wsaData);
	if (err == SOCKET_ERROR) 
	{
		cerr << "Error: winsock failed; see status code for error number" << endl;
		return err;
	}
#else
	const char * home = getenv("HOME");
#endif

	// Get host name
	char nameBuffer[10000];
	gethostname(nameBuffer, sizeof(nameBuffer));
	string hostName(nameBuffer);
	char hostList[10000];	// list of |-delimited hostnames
	bool onHostList;		// hostName found in hostList

	// Get current date as a 8-digit integer YYYYMMDD
	int todaysDate;
	int expirationDate;
	time_t ltime;				// local time, in seconds from 1970
	time(&ltime);
	struct tm *today = localtime(&ltime);	// convert seconds to year,month,day

	todaysDate = 10000 * (1900 + today->tm_year)	// [0..115] -> YYYY0000
		+ 100 * (1 + today->tm_mon)					// [0..11]  ->     MM00
		+ today->tm_mday;							// [1..31]  ->       DD

	// The client's site name and location, no newlines
	char site[1000];		// arbitrary string

	// Key files have a version that tells us what to expect in the key file
	char licVersionBuffer[1024];
	string licVersion(licVersionBuffer);

	double key;		// From the key file

	// Try to load key file
	char *del = "|";	// delimiter between hostnames in a host list
#ifdef _WIN32
	ifstream input(keyFilename.c_str());
#else
	string path;
	ifstream input;

	// Unix/Linux processing.
	// First look in the cd_arg directory
	if (cd_arg != NULL) {
		path = cd_arg;
		path += "/";
		path += keyFilename;
		input.open(path.c_str());
	}
	if (! input.is_open()) {
		// Next look in the current directory
		input.clear();
		path = keyFilename;
		input.open(path.c_str());
		if (! input.is_open()) {
			// Next look in the user's home directory
			input.clear();
			path = home;
			path += "/";
			path += keyFilename;
			input.open(path.c_str());
			if (! input.is_open()) {
				// Next look in the root directory
				input.clear();
				path = "/";
				path += keyFilename;
				input.open(path.c_str());
			}
		}
	}
#endif
	if (! input.is_open())
    {
		cerr << "Warning: no license found for host " << hostName << " in " << keyFilename << endl;

		if (pk)
		{
			licVersion = thisLicVersion;

			// Get key from user
			cerr << "--> enter key: ";
			cin >> key;

			cerr << "--> enter host list -- delimit with " << del << " ";
			cin >> hostList;

			cerr << "--> enter expiration date (YYYYMMDD): ";
			cin >> expirationDate;

			cerr << "--> enter client and site name, no spaces (eg: Joe_Smith,Schmo_Hospital): ";
			cin >> site;		// no spaces!
		}
		else
			return 1;	// Fail
    }
	else {
		// Get key from file
		input >> licVersion;

		if (licVersion != thisLicVersion) {
			cout << "Warning: license has wrong version:" << endl;
			cout << " Need " << thisLicVersion << " -- Have " << licVersion << endl;
			return 1;	// Fail
		}

		input >> key;
		input >> hostList;
		input >> expirationDate;
		input >> site;
		input.close();

		onHostList = findStringDelimited(hostList, hostName.data(), del, debug);
    }

	// Generate key
	//	key = checksum of all components
	char checksumBuffer[1000] = "";
	strcat(checksumBuffer, hostList);
	char ed[100];
	sprintf(ed, "%ld", expirationDate);
	strcat(checksumBuffer, ed);
	strcat(checksumBuffer, site);

	int len = strlen(checksumBuffer);
	unsigned long crcSum = 0;
	for (unsigned int i = 0; i < len; i++)
		crcSum += i * static_cast<unsigned long>(checksumBuffer[i]);

	// The all-important key calculation equation; this is
	// the ONLY code in the entire app that calc's a key
	// and is tied to the license keyFile version number.
	// Check for overflowing an unsigned long, just to
	// tell when it happens (for debugging)
	double dproduct = double(crcSum) * int(magicNumber/1000) + double(magicOffset);
	unsigned long product = int(dproduct) % (ULONG_MAX-1);
	if ((unsigned long)dproduct >= ULONG_MAX && debug)
		cerr << "Warning: key overflow (" << ULONG_MAX << ")" << endl;

	if (debug) {
		// Print all digits (not exponential format) and nothing past the decimal place
		ios_base::fmtflags f = cerr.flags();	// record default
		cerr.setf (ios_base::fixed, ios_base::floatfield);
		cerr.precision(0);

		cerr << " keyFile=" << key
			<< " as=" << crcSum
			<< " ed=" << expirationDate
			<< " mo=" << magicOffset
			<< " mn=" << magicNumber 
		    << " dp=" << dproduct
			<< " lp=" << product
			<< endl;

		cerr.flags(f);		// reset to default
	}

	if (key == magicNumber) {
		// Write correct key to file
		std::ofstream output(keyFilename.data());
		if (output.fail())
		{
			cerr << "Error opening key file (" << keyFilename << ") for writing." << endl;
			return 0;	// pass
		}

		cerr << "Writing key file: " << keyFilename << endl;
		output << licVersion << endl;

		// Print all digits (not exponential format) and nothing past the decimal place
		ios_base::fmtflags f = output.flags();	// record default
		output.setf (ios_base::fixed, ios_base::floatfield);
		output.precision(0);

		output << dproduct << endl;
		output.flags(f);		// Reset to default

		output << hostList << endl;
		output << expirationDate << endl;
		output << site << endl;
		output.close();

		if (output.fail())
			cerr << "Error: could not write key file. Pass." << endl;

		return 0;		// pass
	}

	string contact = "See http://titan.radonc.unc.edu/colab ";
	if (! onHostList)
	{
		cerr << "Error: hostname " << hostName << " is not licensed. "
			<< contact << "to extend your license." << endl;
		return 1;	// fail
	}

	if (dproduct != key)
	{
		cerr << "Error: Invalid key for hostname " << hostName << ". "
			<< contact << "for a license." << endl;
		return 1;	// fail
	}

	if (expirationDate < todaysDate)
	{
		cerr << "Error: your license has expired for hostname " << hostName << ". "
			<< contact << "for an updated license." << endl;
		cerr << " Today is " << todaysDate << endl;
		cerr << " Expiration date is " << expirationDate << endl;

		return 1;		// fail
	}

	if ((expirationDate - 100) < todaysDate)
	{
		cerr << "Warning: Your license expires in a month for hostname " << hostName << ". "
			<< contact << "for an updated license." << endl;
		cerr << " Today is " << todaysDate << endl;
		cerr << " Expiration date is " << expirationDate << endl;

		return 0;	// pass, with a warning
	}

	cout << "Licensed to " << site << ".   Expires " << expirationDate << endl;


	return 0;		// pass
}

#endif	/* BINARY */


bool getarg(char * shortOption, char * longOption, bool takesArg, string & opt,
    int & argNum, bool & error, bool verbose, int argc, char **argv)
{
#ifdef _WIN32
	if (opt[0] == '-' || opt[0] == '/')
#else
	if (opt[0] == '-')
#endif
	{
		string str = string(opt.data() + 1);

		if (str == shortOption || (longOption != NULL && str == longOption)) {
			if (takesArg) {
#ifdef _WIN32
#ifdef __CYGWIN32__
				if (++argNum >= argc || argv[argNum][0] == '-')	// Allow paths beginning with '/'
#else
				if (++argNum >= argc || argv[argNum][0] == '-' || argv[argNum][0] == '/')
#endif
#else
				if (++argNum >= argc || argv[argNum][0] == '-')
#endif
				{
					cout << "Error: Option " << opt << " requires an argument" << endl;
					error = true;
					return false;	// Necessary argument missing
				}
				else
					if (verbose)
						cout << opt << " " << argv[argNum] << endl;
			}
			return true;
		}
		else
			return false;	// No str match
	}
	return false;	// No leading - or /
}


// These two global variables are declared in ControlParms.h
ControlParms * globalControl;	// Read the user's preferences file
int globalVerbosity;			// Current verbosity level of Pablo


int pablo(int argc, char *argv[])
{
    int i;
    bool halt, error;
    bool left_handed_mouse;
    bool convertFormat;
    bool convertAltFormat = false;
    int compressImage;
    bool align, optimize;
    bool runThinPlateSpline;
    bool interactive;
    char * scriptFile;
    char * alignPtsFile;
    bool ignore;
    bool debug;		// Unpublished internal option: super-verbose, for debugging
    bool verbose;	
    bool silent;
    bool versionShown;
    bool allow_aboutPablo;
    int newVersion;
    string cd_arg;
    P3DControl * control;
#ifdef TASKING
    string taskFileName;
#endif
#if defined(PRODUCTION_VERSION) && ! defined(BINARY)
    time_t ltime;
#endif
#ifdef BINARY
    bool regError;
    bool resampleImage;
    bool nonInteractive;
    char * landmarkModelFile;
#else
    bool regularize;
    double reg_steps;
    int reg_iterations;
#endif


#ifdef DEBUG
    cout << "Executing " << argv[0] << '\n';
#endif
    // Set up the version-dependent header string
    headerString = new char[18 + strlen(release) + strlen(copyright) + 1];
#ifdef BINARY
    strcpy(headerString, "Binary ");
#else
    headerString[0] = '\0';
#endif
    strcat(headerString, release);
    strcat(headerString, " ");
    strcat(headerString, copyright);
    // Clear the 1st dash in the release from headerString
    char * p = strstr(headerString, "Pablo");
    if (p) {
        p += 5;
        if (*p == '-')
            *p = ' ';
    }

    debug = false;
    verbose = false;
    silent = false;
    versionShown = false;
    constructExec = false;
    touchFile = NULL;

    halt = false;
    error = false;
    compressImage = 0;		// 0:skip  1:compress  2:uncompress
    convertFormat = false;
    align = false;
    optimize = false;
    runThinPlateSpline = false;
    scriptFile = NULL;
    ignore = false;
#ifdef BINARY
    regError = false;
    resampleImage = false;
    nonInteractive = false;
    landmarkModelFile = NULL;
#else
    regularize = false;
    // Regularization defaults
    reg_steps = 0.01;
    reg_iterations = 40;
#endif
    interactive = false;
    left_handed_mouse = false;
    allow_aboutPablo = true;

    string imageFile;
    string outputImageFile;
    string modelFile;
    string outputModelFile;
    string xformFile;
    string outXformFile;
    string maskFile;

    Registry scriptParms(25, 10);
    control = NULL;

#ifdef BINARY

    pk = false;

#ifdef CHECK_DOMAIN
    // This allows users at UNC or other designated locations to
    // run Binary Pablo without a key
    i = 0;
    require_key = true;
    while (domains[i] != NULL) {
        if (0 != strcmp(getenv("USERDOMAIN"), domains[i]))
            i++;
        else {
            require_key = false;
            break;
        }
    }
#else
    require_key = keyRequired;
#endif

#endif	/* BINARY */

    // -----------------------------  Start of options parsing -----------------------------

    for (i = 1; i < argc; i++) {
        string opt = argv[i];

#ifdef _WIN32
        if (opt == "-h" || opt == "/h" || opt == "-help" || opt == "/help")
#else
        if (opt == "-h" || opt == "-help")
#endif
        {
            cout << headerString << "\n\n";
#ifdef AE2_BUILD
            return -1;	// Printing Help is not allowed in the AE2 runs
#endif

            // Note: when options are added to the command line, they also should be
            // added to bpTuning or gpTuning, or the spell checking when processing
            // option -cs may fail.

            // The following options are intentionally not documented in the help
            // output:
            //		Binary Pablo:  -pt, -tf, -id, -pt. -pk
            //		Regular Pablo:	-pt, -tf, -construct, -qv

#ifdef BINARY

            cout << "Syntax: 'bin_pablo [option] ...', with options as described below.\n";
            cout << "\nInput options:\n";            
            cout << "  -ii or -image FILE          Load image file.\n";
			cout << "  -id or -distanceMap FILE    Load distance map file.\n";
            cout << "  -il or -landmarkModel FILE  Load landmark model file.\n";
            cout << "  -im or -model FILE          Load model file.  Instead of using '-im' or\n";
            cout << "                                  '-model', the model file name may simply be\n";
            cout << "                                  given at the end of the command.\n";
            cout << "  -ix or -transform FILE      Load the similarity transform from the specified\n";
            cout << "                                  file.  Any transformation in the model file,\n";
            cout << "                                  will be ignored.\n";
            cout << "  -ic or -contour FILE        Load contour file.\n";
            cout << "\nInput control options:\n";
            cout << "  -cd directory               Specifies location of the pablo.txt control file.\n";
            cout << "  -caf or -altFormat FILE     Convert image file type based on extension and halt.\n";
            cout << "  -cafc or -altFormatc FILE   Like -caf but compresses the resulting file.\n";
            cout << "  -cf or -format N            Convert image file version to N (use -image) and\n";
            cout << "                                  halt.\n";
            cout << "  -cs FILE                    Load and execute script file.\n";
            cout << "  -ci or -ignore              Ignore any unknown parameters found in a script.\n";
            cout << "                                  Must be used before -cs to have effect.\n";
#ifdef TASKING
            cout << "  -t or task FILE             Load the indicated text file into the task window.\n";
#endif
            cout << "\nOutput options:\n";
            cout << "  -om or -outModel FILE       Save fitted model; omit for match value probe.\n";
            cout << "  -omi or -outModelImage FILE Save fitted model's binary image; omit for match\n";
            cout << "                                  value probe.\n";			
			cout << "                                   This option helps to make a video of the optimization process\n";
            cout << "  -ot or -outTile FILE        Save BYU tiles of fitted model.\n";
            cout << "  -oit or -outImageTile FILE  Save BYU tiles of binary surfaces within image.\n";
            cout << "  -ostm or\n";
            cout << "    -outSimTransModel FILE    Save similarity transform after init/fig stage.\n";
            cout << "  -ox or -outTransform FILE   Save any transformation produced.  Used only for alignment.\n";
			cout << "  -vf or -videoFolder DIRECTORY Save models in intermediate optimization phases\n";
            cout << "\nProcessing Options (ignored when running interactively):\n";			
            cout << "  -pc or -compress FILE       Compress the image file, then quit.\n";
            cout << "  -ph FILE                    Print header of image file; then load image and\n";
            cout << "                                  continue.\n";
            cout << "  -pa or -align FILE      Align using points from the file, then quit; open no windows.\n";
            cout << "  -po or -optimize        Optimize, then quit; open no windows.\n";
            cout << "  -ps                         Run landmark based diffeomorphism helper, then quit; open no windows.\n";
            cout << "  -pr                         Print control script defaults, then quit; open no\n";
            cout << "                                  windows.\n";
            cout << "  -pd                         Print control script descriptions, then quit;\n";
            cout << "                                  open no windows.\n";
            cout << "  -pq or -quit                Verify arg syntax, then quit; open no windows.\n";
            cout << "  -pQ                         Quit after performing all -p requests; no windows\n";
            cout << "                                  will be displayed.\n";
            cout << "  -pi or -resample FILE       Isotropically resample the image file, then quit;\n";
            cout << "                                  open no windows.\n";
            cout << "  -pm or -extent FILES        Print max extent of all image files, then quit;\n";
            cout << "                                  open no windows.\n";
            cout << "                              Note: -pm consumes the rest of the command line.\n";
            cout << "  -pu or -uncompress FILE     Uncompress the image file, then quit.\n";
#ifndef _WIN32
            cout << "  -l or -left                 Mouse is left-handed.\n";
#endif
            cout << "\nLogging options:\n";
            cout << "  -h or -help                 Print help, then quit.\n";
            cout << "  -hc                         Print change history, then quit.\n";
            //cout << "  -lf FILE                    Print to logFile [stdout]\n";
            cout << "  -lV or -version             Print the version number (does not terminate).\n";
            cout << "  -lv or -verbose             Verbose: print more than the usual info.  After\n";
            cout << "                                  -lv is seen, the rest of the command-line\n";
            cout << "                                  is printed.\n";
            cout << "  -lq or -quiet               Quiet: print less than the usual info, mostly\n";
            cout << "                                  only errors.\n";
            cout << "  -I or -interactive          Show the scripted optimization running in the\n";
            cout << "                                  GUI.  This may only be used with -cs and -po.\n";
            cout << "\nAny tuning value can be used as a switch.  See -pr.\n";
            cout << "The last value seen overrides previous values.\n\n";
            cout << "Pathnames accepted: Windows absolute, Windows relative (to bin dir).\n\n";
            cout << "For details, see the Binary Pablo Manual.\n";

#else	/* BINARY */

            cout << "Syntax: 'pablo [option] ... [[-model] file]', with options as described below.\n";
            cout << "\nOptions:\n";
            cout << "  -h or -help             Produces this output.  Pablo will then terminate.\n";
            cout << "  -hc                     Print change history, then quit.\n";
            cout << "  -cd directory           Specifies location of the pablo.txt control file.\n";
            cout << "  -cf or -format N        Convert image file version to N (use -image) and halt.\n";
            cout << "  -cs FILE                Load and execute script file.\n";
            cout << "  -ci or -ignore          Ignore any unknown parameters found in a script.\n";
            cout << "                              Must be used before -cs to have effect.\n";
            cout << "  -ii or -image FILE      Load image file.\n";
            cout << "  -ik or -mask FILE       Load simple mask model file.\n";
            cout << "  -im or -model FILE      Load model file.  Instead of using '-im' or '-model',\n";
            cout << "                              the model file name may simply be given at the\n";
            cout << "                              end of the command.\n";
            cout << "  -ip or -pga FILE        Use the PGA data from the specified model file,\n";
            cout << "                              discarding any PGA data loaded by -im or -model.\n";
            cout << "  -ipr or -PGAModel FILE  Load a residue model with either figure or atom PGA\n";
            cout << "                              statistics.\n";
            cout << "  -ix or -transform FILE  Load the similarity transform from the specified\n";
            cout << "                              file.  Any transformation in the model file, will\n";
            cout << "                              be ignored.\n";
            cout << "  -ih or -hist FILE       Load the surface patch ensemble histogram from the\n";
            cout << "                              specified file.\n";
            cout << "  -ihm or -histMean FILE  Load an SPE file to override the Mean Vector loaded\n";
            cout << "                              by -hist.  The -hist and -histMean files should\n";
            cout << "                              have identical vector lengths.\n";
            cout << "  -ic or -contour FILE    Load contour file.\n";
            cout << "  -iqf or -QFData FILE    Load additional image features.\n";
            cout << "  -oh or -histOut FILE    Save the generated surface patch ensemble histogram\n";
            cout << "                              from training\n";
            cout << "  -histPatchImageFIG_NUM FILE\n";
            cout << "                          Save file of histogram patches, as an image file,\n";
            cout << "                              generated during histogram training for figure\n";
            cout << "                              FIG_NUM.\n";
            cout << "  -histBinaryImageFIG_NUM FILE\n";
            cout << "                          Load a binary image for histogram training of the\n";
            cout << "                              specified figure.\n";
            cout << "  -mpBinaryImageFIG_NUM FILE\n";
            cout << "                          Load a binary image for MultiPatch training of the\n";
            cout << "                              specified figure.\n";
            cout << "  -I or -interactive      Show the scripted optimization running in the GUI.\n";
            cout << "                              This may only be used with -po.\n";
            cout << "  -lV or -version         Print the version number.\n";
            cout << "  -lv or -verbose         Verbose: print more than the usual info.  After -lv\n";
            cout << "                              is seen, the rest of the command-line is printed.\n";
            cout << "  -lq or -quiet           Quiet: print less than usual, mostly only errors.\n";
#ifndef _WIN32
            cout << "  -l or -left             Mouse is left-handed.\n";
#endif
            cout << "  -om or -outModel FILE   Name of fitted model to be written after optimization.\n";
            cout << "  -ot or -outTile FILE    Generate a BYU tile set of the fitted model.\n";
            cout << "  -ox or -outTransform FILE   Save any transformation produced.  Used only for alignment.\n";
            cout << "  -pc or -compress FILE   Compress the image file (use -ii), then quit.\n";
            cout << "  -pd                     Print control script descriptions, then quit;\n";
            cout << "                              open no windows.\n";
            cout << "  -pa or -align FILE      Align using points from the file, then quit; open no windows.\n";
            cout << "  -po or -optimize        Optimize, then quit; open no windows.\n";
            cout << "  -pr                     Print control script defaults, then quit; open no\n";
            cout << "                              windows.\n";
            cout << "  -pq or -quit            Stop after verifying options' syntax; open no windows.\n";
            cout << "  -pu or -uncompress FILE Uncompress the image file (use -ii), then quit.\n";
            cout << "  -r or -reg [N S]        Regularize the model N steps of size S and halt.\n";
            cout << "                              Default values are 40 and 0.01, respectively.\n";
#ifdef TASKING
            cout << "  -t or -task FILE        Load the indicated text file into the task window.\n";
#endif

#endif	/* BINARY */
            cout << endl;
            halt = true;
            continue;
        }

#ifdef AE2_BUILD
        // Undocumented option required for running with unflipped (PLUNC) inputs.
        // This avoids accidentally running the AE2 version of Pablo on flipped
        // (Comp Sci) data.
        if (getarg("construct", NULL, false, opt, i, error, verbose, argc, argv)) {
            constructExec = true;
            continue;
        }
#endif
        // Undocumented option processed in main()
        if (getarg("tp", "pid", true, opt, i, error, false, argc, argv)) {
            // Ignore; this option is processed in main()
            continue;
        }

        if (getarg("cd", NULL, true, opt, i, error, verbose, argc, argv)) {
            cd_arg = argv[i];
            continue;
        }

        if (getarg("caf", "altFormat", true, opt, i, error, verbose, argc, argv)) {
            convertAltFormat = true;
            compressImage = 2;
            outputImageFile = argv[i];
            continue;
        }

        if (getarg("cafc", "altFormatc", true, opt, i, error, verbose, argc, argv)) {
            convertAltFormat = true;
            compressImage = 1;
            outputImageFile = argv[i];
            continue;
        }

        if (getarg("cf", "format", true, opt, i, error, verbose, argc, argv)) {
            convertFormat = true;
            newVersion = atoi(argv[i]);
            continue;
        }

        if (getarg("ci", "ignore", false, opt, i, error, verbose, argc, argv)) {
            ignore = true;
            continue;
        }

        if (opt == "-hc") {
            cout << longChangesString << endl;
            halt = true;
            continue;
        }

        if (getarg("ii", "image", true, opt, i, error, verbose, argc, argv)) {
            imageFile = argv[i];
            scriptParms.setStringValue("Image", argv[i]);
            continue;
        }

        if (getarg("im", "model", true, opt, i, error, verbose, argc, argv)) {
            modelFile = argv[i];
            scriptParms.setStringValue("Model", argv[i]);
            continue;
        }

        if (getarg("ipr", "PGAModel", true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("PgaResidueModel", argv[i]);
            continue;
        }

        if (getarg("ip", "pga", true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("PgaModel", argv[i]);
            continue;
        }

        if (getarg("ix", "transform", true, opt, i, error, verbose, argc, argv)) {
            xformFile = argv[i];
            scriptParms.setStringValue("Transform", argv[i]);
            continue;
        }

#ifndef BINARY
        if (getarg("ih", "hist", true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("Hist", argv[i]);
            continue;
        }

        if (getarg("ihm", "histMean", true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("HistMean", argv[i]);
            continue;
        }


        //if (getarg("dqf", NULL, true, opt, i, error, verbose, argc, argv)) {
        //	scriptParms.setStringValue("Dqf", argv[i]);
        //	continue;
        //}

        if (getarg("mpConfig", NULL, true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("MpConfig", argv[i]);
            continue;
        }

        if (getarg("mpfDebug", NULL, true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("MPFDebug", argv[i]);
            continue;
        }

        if (getarg("mpOut", NULL, true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("MpOut", argv[i]);
            continue;
        } 

        if (getarg("dqfConfig", NULL, true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("DqfConfig", argv[i]);
            continue;
        }

        if (getarg("dqfOut", NULL, true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("DqfOut", argv[i]);
            continue;
        }

        if (getarg("oh", "histOut", true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("HistOut", argv[i]);
            continue;
        }

        if (getarg("iqf", "QFData", true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("QFData", argv[i]);
            continue;
        }

        // This will work for up to 19 objects
        int f;
        for (f = 0; f < 20; f++) {
            char key[256];

            sprintf(key, "histPatchImage%d", f);
            if (getarg(key, NULL, true, opt, i, error, verbose, argc, argv)) {
                sprintf(key, "HistPatchImage%d", f);
                scriptParms.setStringValue(key, argv[i]);
                break;
            }

            sprintf(key, "histBinaryImage%d", f);
            if (getarg(key, NULL, true, opt, i, error, verbose, argc, argv)) {
                sprintf(key, "HistBinaryImage%d", f);
                scriptParms.setStringValue(key, argv[i]);
                break;
            }

            sprintf(key, "mpBinaryImage%d", f);
            if (getarg(key, NULL, true, opt, i, error, verbose, argc, argv)) {
                sprintf(key, "MPBinaryImage%d", f);
                scriptParms.setStringValue(key, argv[i]);
                break;
            }

            sprintf(key, "dqf%d", f);
            if (getarg(key, NULL, true, opt, i, error, verbose, argc, argv)) {
                sprintf(key, "Dqf%d", f);
                scriptParms.setStringValue(key, argv[i]);
                break;
            }

            sprintf(key, "dqfImage%d", f);
            if (getarg(key, NULL, true, opt, i, error, verbose, argc, argv)) {
                sprintf(key, "DqfImage%d", f);
                scriptParms.setStringValue(key, argv[i]);
                break;
            }

        }
        if (f < 20)
            continue;

#endif

#ifdef AE2_BUILD
        // Undocumented option for use by ConStruct only
        if (getarg("qv", NULL, false, opt, i, error, verbose, argc, argv)) {
            cout << release << ' ' << revision << '\n';
            versionShown = true;
            continue;
        }
#endif

        // Undocumented option for use by ConStruct or external scripts.
        // This causes the execution status to be written into the designated
        // file when the program exits.  See main.cpp.
        if (getarg("tf", "touch", true, opt, i, error, false, argc, argv)) {
            touchFile = argv[i];
            if (verbose)
                cout << "Touchfile: " << touchFile << endl;
            continue;
        }

        if (getarg("lV", "version", false, opt, i, error, verbose, argc, argv)) {
            cout << headerString << "\n\n";
            cout << "Pablo version: " << release << ' ' << revision << "\n\n";
#ifdef BINARY
            cout.precision(2);
            cout.setf(ios_base::fixed, ios_base::floatfield);
            cout << "Tuning param version: " << tuningWt(ScriptVersion) << '\n';
            cout.setf(~ios_base::fixed, ios_base::floatfield);
            cout << "Accepts license version: " << thisLicVersion << "\n\n";
#endif
            versionShown = true;
            continue;
        }

        // The option is unpublished, for use at UNC only
        if (getarg("lD", "debug", false, opt, i, error, verbose, argc, argv)) {
            if (silent)
                cout << "Ignoring " << opt << " because of -lq or -quiet" << endl;
            else if (verbose)
                cout << "Ignoring " << opt << " because of -lv or -verbose" << endl;
            else
                debug = true;
            continue;
        }

        if (getarg("lv", "verbose", false, opt, i, error, verbose, argc, argv)) {
            if (silent)
                cout << "Ignoring " << opt << " because of -lq or -quiet" << endl;
            else if (debug)
                cout << "Ignoring " << opt << " because of -lD or -debug" << endl;
            else
                verbose = true;
            continue;
        }

        if (getarg("lq", "quiet", false, opt, i, error, verbose, argc, argv)) {
            if (verbose)
                cout << "Ignoring " << opt << " because of -lv or -verbose" << endl;
            else if (debug)
                cout << "Ignoring " << opt << " because of -lD or -debug" << endl;
            else
                silent = true;
            continue;
        }

        if (getarg("I", "interactive", false, opt, i, error, verbose, argc, argv)) {
            interactive = true;
            continue;
        }

#ifndef _WIN32
        if (getarg("l", "left", false, opt, i, error, verbose, argc, argv)) {
            left_handed_mouse = true;
            continue;
        }
#endif

        if (getarg("pc", "compress", true, opt, i, error, verbose, argc, argv)) {
            outputImageFile = argv[i];
            compressImage = 1;
            continue;
        }

        if (getarg("pq", "quit", false, opt, i, error, verbose, argc, argv)) {
            halt = true;
            continue;
        }

        if (getarg("cs", NULL, true, opt, i, error, verbose, argc, argv)) {
            // Load the script into a registry, so later command-line options
            // may change the values
            scriptFile = argv[i];
            try {
                scriptParms.readFromFile(scriptFile, true);
            }
            catch (RException excp) {
                excp.print(std::cout);
                std::cout << "Error: invalid script file: " << scriptFile << '\n';
                return -1;
            }

            // Make sure the script doesn't contain spelling errors.
            // To do this, lookup every script entry in tuneVals.
            // Note that this will only consider command-line arguments
            // that are processed before -cs is.
            int nkeys = scriptParms.getKeyArraySize() - 1;
            char ** keys = new char *[nkeys + 1];
            if (scriptParms.getValueKeys(keys) != nkeys) {
                // This should never occur
                cout << "Error: failed to locate all keys; aborting script" << endl;
                return -1;
            }
            bool unknown = false;
            for (int k = 0; k < nkeys; k++) {
                if (tuneVals.indexByLongName(keys[k]) < 0) {
                    if (! unknown)
                        cout << '\n';
                    cout << "Warning: unknown script specification: " << keys[k] << '\n';
                    unknown = true;
                }
                keys[k] = NULL; 
            }
            delete [] keys;
            if (unknown) {
                cout << endl;
                if (! ignore)
                    return -1;
            }

            if (! silent)
                cout << "Loaded control script " << scriptFile << '\n';
            if (verbose)
                cout << "  Script specifications:\n";

            // Override current tuning values with the script values.
            // First verify that the version number is compatible.
            // See comment on PABLO VERSION, above.
            double proposedVersion = scriptParms.getDoubleValue(tuningLongName(ScriptVersion),
                tuningWt(ScriptVersion));
            if ((int) proposedVersion > (int) tuningWt(ScriptVersion)) {
                cout << "  Warning: script specifies " << tuningLongName(ScriptVersion) << " = ";
                printVersion(proposedVersion);
                cout << " -- expecting ";
                printVersion(tuningWt(ScriptVersion), 0);
                cout << ".xx ... proceed with caution\n";
            }

            // Next loop over other script entries to update the tuning parameters
            for (int w = 1; w < TUNE_COUNT; w++) {
                if (! tuningLocked(w)) {	// Silently ignore changes to locked vals
                    const char * key = tuningLongName(w);
                    if (tuneVals.isString(w)) {
                        // Strings have no default value in tuneVals
                        const char * newStr = scriptParms.getStringValue(key, "");
                        if (strlen(newStr) > 0) {
                            (void) tuneVals.setString(w, newStr);
                            if (verbose)
                                cout << "    " << key << " = " << newStr << '\n';
                        }
                    }
                    else if (tuneVals.isList(w)) {
                        int len;

                        // This makes the array unreadable a second time
                        int * list = scriptParms.getIntArray(key, &len);
                        // Lists have no default value in tuneVals
                        if (len > 0) {
                            // Copy the array into class Tuning
                            if (tuneVals.setList(w, len, list) == false) {
                                cout << "Error processing list item in script file: parameter name is "
                                    << tuningLongName(w) << endl;
                                return -1;
                            }
                            if (verbose) {
                                cout << "    " << key << " = " << list[0];
                                for (int i = 1; i < len; i++)
                                    cout << ", " << list[i];
                                cout << '\n';
                            }
                        }
                    }
                    else {
                        double lastWt = tuningWt(w);
                        double newWt = scriptParms.getDoubleValue(key, lastWt);
                        if (lastWt != newWt) {
                            (void) tuneVals.setWeight(w, newWt);
                            if (verbose)
                                cout << "    " << key << " = " << newWt << '\n';
                        }
                    }
                }
            }

#ifndef BINARY
            // Convert parameter TileFiguresList into TileFiguresMask
            int len;
            unsigned long mask = 0;
            const int * tileFigs = tuningList(TileFiguresList, len);
            if (len != 0) {
                for (int i = 0; i < len; i++)
                    mask |= (1 << tileFigs[i]);

                if (tuningWt(TileFiguresMask) != 0.0) {
                    cout << "Error: TileFiguresMask and TileFiguresList are mutually exclusive"
                        << endl;
                    return -1;
                }

                (void) tuneVals.setWeight(TileFiguresMask, mask);
            }
#endif

            continue;
        }

#ifdef BINARY

        if (getarg("il", "landmarkModel", true, opt, i, error, verbose, argc, argv)) {
            landmarkModelFile = argv[i];
            scriptParms.setStringValue("LandmarkModel", argv[i]);
            continue;
        }

		// dibyendu - this should be used now - 07/28/2011
        if (getarg("id", "distanceMap", true, opt, i, error, verbose, argc, argv)) {		// undocumented
            scriptParms.setStringValue("DistanceMap", argv[i]);
            continue;
        }

#else

        if (getarg("ik", "mask", true, opt, i, error, verbose, argc, argv)) {
            maskFile = argv[i];
            scriptParms.setStringValue("SimpleMaskFile", argv[i]);
            continue;
        }

#endif

        if (getarg("ic", "contour", true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("Contour", argv[i]);
            continue;
        }

        if (getarg("om", "outModel", true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("OutModel", argv[i]);
            outputModelFile = argv[i];
            continue;
        }

        if (getarg("ot", "outTile", true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("OutTile", argv[i]);
            continue;
        }

        if (getarg("ox", "outTransform", true, opt, i, error, verbose, argc, argv)) {
            //			scriptParms.setStringValue("outTransform", argv[i]);
            outXformFile = argv[i];
            continue;
        }

        if (getarg("vf", "videoFolder", true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("VideoFolder", argv[i]);
            continue;
        }


#ifdef BINARY

        if (getarg("ostm", "outSimTransModel", true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("OutSimTransModel", argv[i]);
            continue;
        }

        if (getarg("oit", "outImageTile", true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("OutImageTile", argv[i]);
            continue;
        }

        if (getarg("omi", "outModelImage", true, opt, i, error, verbose, argc, argv)) {
            scriptParms.setStringValue("OutModelImage", argv[i]);
            continue;
        }

        if (getarg("ph", NULL, true, opt, i, error, verbose, argc, argv)) {
            globalControl = new ControlParms(NULL, -1, false);
            setGlobalControlDefaults();
            globalVerbosity = (verbose ? 1 : 0);
            if (silent)
                globalVerbosity = -1;
            else if (debug)
                globalVerbosity = 2;
            globalControl->setDefault(OutputVerbosity, globalVerbosity);
            globalControl->setDefault(ReorderModels, 0);
            globalControl->setDefault(SmoothImages, false);
            globalControl->setDefault(ConvertImages, false);
            globalControl->setDefault(ConvertImageFormat, false);
            globalControl->setDefault(ByteOrder, 1);	// Native
            globalControl->setDefault(CompressImages, true);

            // undoLen not significant for image manipulations, so set it to 1
            if (control == NULL)
                control = new P3DControl(0);

            imageFile = argv[i];
            // For speed, only load header
            if (! control->loadImage(imageFile.data(), false, false, true))
                return 0;

            // Print like print_image_info() but format like a registry
            Image3D * image = control->getImagePtr();
            cout << "image = " << imageFile << ";" << endl;
            cout << "# voxel counts" << endl;
            cout << "xDim = " << image->getXDim() << ";" << endl;
            cout << "yDim = " << image->getYDim() << ";" << endl;
            cout << "zDim = " << image->getZDim() << ";" << endl;

            cout << "# voxel center to voxel center distance (cm)" << endl;
            cout << "xSpacing = " << image->getXSpacing() << ";" << endl;
            cout << "ySpacing = " << image->getYSpacing() << ";" << endl;
            cout << "zSpacing = " << image->getZSpacing() << ";" << endl;

            cout << "# center of first voxel (cm)" << endl;
            cout << "xOrigin = " << image->getXWorldOrigin() << ";" << endl;
            cout << "yOrigin = " << image->getYWorldOrigin() << ";" << endl;
            cout << "zOrigin = " << image->getZWorldOrigin() << ";" << endl;

            cout << "# center of last voxel (cm)" << endl;
            cout << "xOpposite = " << image->getXImageBound() << ";" << endl;
            cout << "yOpposite = " << image->getYImageBound() << ";" << endl;
            cout << "zOpposite = " << image->getZImageBound() << ";" << endl;

            cout << "# distance from first voxel to last voxel (cm)" << endl;
            cout << "xExtent = " << image->getXExtent() << ";" << endl;
            cout << "yExtent = " << image->getYExtent() << ";" << endl;
            cout << "zExtent = " << image->getZExtent() << ";" << endl;

            cout << "# max distance in any dimension (cm)" << endl;
            cout << "maxExtent = " << image->maxExtent() << ";" << endl;

            continue;
        }

        if (getarg("pi", "resample", true, opt, i, error, verbose, argc, argv)) {
            imageFile = argv[i];
            resampleImage = true;
            continue;
        }

        if (getarg("pk", NULL, false, opt, i, error, verbose, argc, argv)) {
            // Undocumented; for MIDAG internal use only
            // generating keyFile, so make sure one doesn't already exist
            pk = true;
            ifstream input(keyFilename.data());
            if (! input.fail()) {
                cout << "Error: cowardly refusing to overwrite keyFile ("
                    << keyFilename << "). Delete it manually." << endl;
                return -1;
            }
            input.close();

            halt = true;
            continue;
        }

        if (getarg("pm", "extent", true, opt, i, error, verbose, argc, argv)) {
            float maxExtent = 0, extent;
            allow_aboutPablo = false;
            globalControl = new ControlParms(NULL, -1, false);
            globalVerbosity = (verbose ? 1 : 0);
            if (silent)
                globalVerbosity = -1;
            else if (debug)
                globalVerbosity = 2;
            globalControl->setDefault(OutputVerbosity, globalVerbosity);
            globalControl->setDefault(ReorderModels, 0);
            globalControl->setDefault(SmoothImages, false);
            globalControl->setDefault(ConvertImages, false);
            globalControl->setDefault(ConvertImageFormat, false);
            globalControl->setDefault(ByteOrder, 1);	// Native
            globalControl->setDefault(CompressImages, true);

            // undoLen not significant for image manipulations, so set it to 1
            if (control == NULL)
                control = new P3DControl(0);

            // Could consume up to the next switch...
            int iCount = 0;		// image count
            while (i < argc) {
                imageFile = argv[i];
                extent = control->calcMaxExtent(imageFile.c_str());
                if (extent > maxExtent)
                    maxExtent = extent;
                if (globalVerbosity >= 1)
                    cout << imageFile << " has extent " << extent << endl;
                i++;
                iCount++;
            }
            if (globalVerbosity > -1)
                cout << "Maximum extent over " << iCount << " images: "
                << maxExtent << " (world units)" << endl;
            return 0;
        }

        if (getarg("pQ", NULL, false, opt, i, error, verbose, argc, argv)) {
            nonInteractive = true;
            continue;
        }

        if (getarg("ps", "spline", false, opt, i, error, verbose, argc, argv)) {
            optimize = true;
            runThinPlateSpline = true;
            continue;
        }

#endif	/* BINARY */

        if (getarg("pd", "printDescriptions", false, opt, i, error, verbose, argc, argv)) {
            tuneVals.printDesc();
            halt = true;
            continue;
        }

        if (getarg("pa", "align", true, opt, i, error, verbose, argc, argv)) {
            alignPtsFile = argv[i];
            align = true;
            continue;
        }

        if (getarg("po", "  optimize", false, opt, i, error, verbose, argc, argv)) {
            optimize = true;
            continue;
        }

        if (getarg("pr", "printDefaults", false, opt, i, error, verbose, argc, argv)) {
            tuneVals.printAsRegistry();
            halt = true;
            continue;
        }

        if (getarg("pt", "printTable", false, opt, i, error, verbose, argc, argv)) {
            // Undocumented; for MIDAG internal use only
            // make sure that this header line and printInputsAsTDTable
            // and printAsTDTable all print the same number of columns
            cout << "Registry Name and Default (case independent)"
                << "\tSwitch Alias\tDescription" << endl;
            printInputsAsTDTable();
            tuneVals.printAsTDTable();
            halt = true;
            continue;
        }

        if (getarg("pu", "uncompress", true, opt, i, error, verbose, argc, argv)) {
            outputImageFile = argv[i];
            compressImage = 2;
            continue;
        }

#ifndef BINARY

        if (getarg("r", "reg", false, opt, i, error, verbose, argc, argv)) {
            // Use of arguments is optional, so 'false' was passed to getarg()
            if (i + 2 < argc && argv[i + 1][0] != '-' && argv[i + 2][0] != '-')
            {
                // Located 2 arguments for -reg option
                reg_iterations = atoi(argv[++i]);
                reg_steps = atof(argv[++i]);
            }
            regularize = true;
            continue;
        }

#endif

#ifdef TASKING
        if (getarg("t", "task", true, opt, i, error, verbose, argc, argv)) {
            taskFileName = argv[i];
            continue;
        }
#endif	/* TASKING */

        // Any tuning value can be a switch (case independent).
        // Command-line tuning parameters override script values, if they
        // are placed after the -cs option on the command.  If placed
        // before -cs, then a script value overrides the command-line value.
        // The ScriptVersion parameter is a special case; it cannot be changed,
        // but using it verifies its value.
#ifdef _WIN32
        if (opt[0] == '-' || opt[0] == '/')
#else
        if (opt[0] == '-')
#endif
        {
            int nparm;

            nparm = (int) tuneVals.indexByLongName(opt.data() + 1);
            if (nparm < 0) {
                cout << "Error: Unknown option: " << opt << "\nTry using -help" << endl;
                error = true;
                continue;
            }
#ifdef _WIN32
#ifdef __CYGWIN32__
            if (++i >= argc || argv[i][0] == '-')	// Allow paths beginning with '/'
#else
            if (++i >= argc || argv[i][0] == '-' || argv[i][0] == '/')
#endif
#else
            if (++i >= argc || argv[i][0] == '-')
#endif
            {
                cout << "Error: Option " << opt << " requires an argument" << endl;
                error = true;
                return false;
            }
            else
                if (verbose)
                    cout << opt << " " << argv[i] << endl;

            const char * tag = opt.data() + 1;	// Tuning parameter's name
            // Silently ignore changes to locked values
            if (! tuningLocked(nparm)) {
                // Set the tuning value
                if (tuneVals.isList(nparm)) {
                    // List values are not supported
                    cout << "Error: Option " << opt
                        << " cannot be set on the command line" << endl;
                    error = true;
                    return false;
                }
                if (tuneVals.isString(nparm)) {
                    //char * last = new char[1 + sizeof(tuningStr(nparm))];
                    char * last = new char[1 + strlen(tuningStr(nparm))];
                    strcpy(last, tuningStr(nparm));
                    tuneVals.setString(nparm, argv[i]);
                    // Replace any script value
                    scriptParms.setStringValue(tag, argv[i]);
                    // Print when non-default values are specified
                    if (verbose && 0 == strcmp(last, tuningStr(nparm)))
                        cout << "  switch changes " << tuningLongName(nparm)
                        << " = " << tuningWt(nparm) << endl;
                    delete [] last;
                }
                else {
                    double last = tuningWt(nparm);
                    tuneVals.setWeight(nparm, atof(argv[i]));
                    // Replace any script value
                    scriptParms.setDoubleValue(tag, atof(argv[i]));
                    // Print when non-default values are specified
                    if (verbose && last != tuningWt(nparm))
                        cout << "  switch changes " << tuningLongName(nparm)
                        << " = " << tuningWt(nparm) << endl;
                }
            }
            if (nparm == ScriptVersion) {
                // Verify version number is compatible
                double proposedVersion = atof(argv[i]);
                if ((int) proposedVersion != (int) tuningWt(nparm)) {
                    cout << "Warning: using " << tuningLongName(nparm) << " = " 
                        << (int) tuningWt(nparm)
                        << ", which differs from command-line value of "
                        << (int) proposedVersion << endl;
                }
            }
            continue;
        }

        // The model file name may be placed after all options 
#ifdef _WIN32
#ifdef __CYGWIN32__
        if (i != argc - 1 || opt[0] == '-')    // Allow paths beginning with '/'
#else
        if (i != argc - 1 || opt[0] == '-' || opt[0] == '/')
#endif
#else
        if (i != argc - 1 || opt[0] == '-')
#endif
        {
            cout << "Error: Unknown option: " << opt << "\nTry using -help or -pd"
                << endl;
            error = true;
            continue;
        }
        else {
            modelFile = opt;
            scriptParms.setStringValue("Model", opt.data());
        }
    }

    // -----------------------------  End of options parsing -----------------------------

#ifdef BINARY

    // Here's where the normal processing begins, so introduce me
    if (! silent && ! versionShown)
        cout << headerString << endl;

    // Check for key file
    if (require_key && findValidKey(pk, cd_arg.c_str(), verbose)) {
        cerr << "Press control-c or close the window to exit" << endl;
        while (1) {};		// Infinite loop for user to view message
    }

    if (optimize && ! scriptFile) {
        cout << "Error: -po requires -cs to be specified" << endl;
        error = true;
    }

#endif	/* BINARY */

    if (compressImage && imageFile.length() <= 0) {
        cout << "Options -pc and -pu require using option -image" << endl;
        error = true;
    }

    if (convertFormat && imageFile.length() <= 0) {
        cout << "Option -format requires using option -image" << endl;
        error = true;
    }

    if (interactive && ! optimize) {
        cout << "Option -I requires using option -po" << endl;
        error = true;
    }

    if (interactive && align) {
        cout << "Option -I cannot be used with option -pa" << endl;
        error = true;
    }

    if (align && optimize) {
        cout << "Option -po cannot be used with option -pa" << endl;
        error = true;
    }

#ifndef BINARY
    if (regularize && modelFile.length() <= 0) {
        cout << "Cannot regularize without a model" << endl;
        halt = true;
    }

    if (regularize && scriptFile == NULL) {
        cout << "Cannot regularize and use a script" << endl;
        halt = true;
    }
#endif

    if (halt || error) {
        if (halt)
            return 0;
        else {
#ifdef BINARY
            cout << "Use -h for help" << endl;
#endif
            return -1;
        }
    }

    // Install out-of-memory handler
#ifdef _MSC_VER
    _set_new_handler(newFailureTrap);
#else
    set_new_handler(newFailureTrap);
#endif

    // Seed random number generator
    // Note: P3DControl::runBinaryPablo() resets this later for BP 
#if defined(PRODUCTION_VERSION) && ! defined(BINARY)
    time(&ltime);
    srand(ltime);
    rand();
    rand();
#else
    srand(12345);
#endif

    pauly::BrentLinearMethod::setDefaultBrentLinearSearchBoundFactor(
        tuningWt(BrentLinearSearchBoundFactor) );

    // Set the global controls, including the default for OutputVerbosity
    if (optimize && ! interactive)	// See if this is a scripted optimization run
    {
        bool ret;

        globalControl = new ControlParms(NULL, -1, false);	// AGG: Why not load preferences?
        globalVerbosity = (verbose ? 1 : 0);
        if (silent)
            globalVerbosity = -1;
        else if (debug)
            globalVerbosity = 2;

        setGlobalControlDefaults();

        // AGG: Should all of these be here?
        globalControl->setDefault(OutputVerbosity, globalVerbosity);
        globalControl->setDefault(ReorderModels, 0);
        globalControl->setDefault(SmoothImages, false);
        globalControl->setDefault(ConvertImages, false);
        globalControl->setDefault(ConvertImageFormat, false);
        globalControl->setDefault(ByteOrder, 1);	// Native
        globalControl->setDefault(CompressImages, true);
        globalControl->setDefault(ScaleImages, true);
        globalControl->setDefault(SimTransformSeparate, false);
        globalControl->setDefault(SimTransformMatrix, false);

        // undoLen not significant for scripting, so set it to 1
        if (control == NULL)
            control = new P3DControl(0);
#ifdef AE2_BUILD
#ifdef BINARY
        if (! constructExec || ! control->verifyApplication(AE2Binary)) {
            cout << "Error: Bad configuration; you need to recompile Binary Pablo using";
            cout << " the correct\n    configuration." << endl;
            return 1;
        }
#else
        if (! constructExec || ! control->verifyApplication(AE2)) {
            cout << "Error: Bad configuration; you need to recompile Pablo using";
            cout << " the correct\n    configuration." << endl;
            return 1;
        }
#endif
#endif

#ifdef BINARY
        if (runThinPlateSpline)
            ret = (0 == control->runThinPlateSpline(scriptParms));
        else
		{
			cout << "Running Binary Pablo..." << endl;
            ret = control->runBinaryPablo(scriptParms);
		}
#else
        ret = control->runPablo(scriptParms);
#endif
        delete control;
        delete globalControl;
        return (ret ? 0 : 1);
    }

    // Set the global controls, including the default for OutputVerbosity
    if (align)	// Perform an alignment run
    {
        bool ret;

        globalControl = new ControlParms(NULL, -1, false);	// AGG: Why not load preferences?
        globalVerbosity = (verbose ? 1 : 0);
        if (silent)
            globalVerbosity = -1;
        else if (debug)
            globalVerbosity = 2;

        setGlobalControlDefaults();

        // AGG: Should all of these be here?
        globalControl->setDefault(OutputVerbosity, globalVerbosity);
        globalControl->setDefault(ReorderModels, 0);
        globalControl->setDefault(SmoothImages, false);
        globalControl->setDefault(ConvertImages, false);
        globalControl->setDefault(ConvertImageFormat, false);
        globalControl->setDefault(ByteOrder, 1);	// Native
        globalControl->setDefault(CompressImages, true);
        globalControl->setDefault(ScaleImages, true);
        globalControl->setDefault(SimTransformSeparate, false);
        globalControl->setDefault(SimTransformMatrix, false);

        // undoLen not significant for scripting, so set it to 1
        if (control == NULL)
            control = new P3DControl(0);
#ifdef AE2_BUILD
#ifdef BINARY
        if (! constructExec || ! control->verifyApplication(AE2Binary)) {
            cout << "Error: Bad configuration; you need to recompile Binary Pablo using";
            cout << " the correct\n    configuration." << endl;
            return 1;
        }
#else
        if (! constructExec || ! control->verifyApplication(AE2)) {
            cout << "Error: Bad configuration; you need to recompile Pablo using";
            cout << " the correct\n    configuration." << endl;
            return 1;
        }
#endif
#endif

        ret = control->align(alignPtsFile, modelFile.c_str(), imageFile.c_str(),
            0.04, false, false, false, outXformFile.c_str());
        if (outputModelFile.size() > 0) {
            control->write(outputModelFile.c_str(), NULL);
            cout << "Saved model " << outputModelFile << endl;
        }

        delete control;
        delete globalControl;
        return (ret ? 0 : 1);
    }

    // Non-scripted execution follows

#ifdef AE2_BUILD
    // For this version, interactive execution is prohibited.
    // No error message is printed, since this is not an
    // authorized execution.
    cout << "Non-interactive version\n";
    return 1;
#endif

    int loudness = (silent ? -1 : 0);
    if (verbose)
        loudness = 1;

    // Set the global controls.
    // Note that the default for OutputVerbosity is set here.
#ifdef BINARY
    if (compressImage || convertFormat || 
        resampleImage || convertAltFormat)
#else
    if (compressImage || convertFormat || regularize)
#endif
    {
        allow_aboutPablo = false;
        globalControl = new ControlParms(NULL, -1, false);	// Ignore user's preferences
        globalVerbosity = (verbose ? 1 : 0);
        if (silent)
            globalVerbosity = -1;
        else if (debug)
            globalVerbosity = 2;
        globalControl->setDefault(OutputVerbosity, globalVerbosity);
        globalControl->setDefault(ReorderModels, 0);

#ifdef BINARY

        globalControl->setDefault(SmoothImages, false);
        //		globalControl->setDefault(ConvertImageFormat, false);
        globalControl->setDefault(ConvertImages, false);
        globalControl->setDefault(ByteOrder, 1);	// Native
        globalControl->setDefault(CompressImages, true);
        globalControl->setDefault(ShowLandmarks, true);

        // The undoLen is insignificant for image manipulations, so set it to 1
        if (control == NULL)
            control = new P3DControl(0);

        if (resampleImage) {
            // Load, re-sample, save image
            if (globalVerbosity > 0)
                cout << "Loaded image " << imageFile << endl;
            bool ok = control->loadImage(imageFile.data());
            if (! ok) {
                cout << "Error: did not load image " << imageFile << endl;
                return -1;	// fatal
            }

            if (tuningWt(BpCoordCubeImage) > 0.0) {		// padding enabled?
                control->resampleImageUnitCubic(tuningWt(BpCoordMaxExtent));
                if (globalVerbosity > -1)
                    cout << "Re-sampled and padded image " << imageFile << endl;
            }
            else {
                control->resampleImage();
                if (globalVerbosity > -1)
                    cout << "Re-sampled image " << imageFile << endl;
            }

            if (globalVerbosity > 0)
                cout << "Saving image " << imageFile << endl;
            ok = control->saveImage(imageFile.data());
            if (! ok) {
                cout << "Error: did not save image " << imageFile << endl;
                return -1;	// fatal
            }

            return 0;
        }
#endif	/* BINARY */
    } else {
        // Interactive use
        if (cd_arg.length() > 0)
            globalControl = new ControlParms(cd_arg.data(), loudness);
        else
            globalControl = new ControlParms(NULL, loudness);

        globalVerbosity = globalControl->read(OutputVerbosity, 0);
        // Command-line options override the control file setting
        if (silent)
            globalVerbosity = -1;
        else if (verbose)
            globalVerbosity = 1;
        else if (debug)
            globalVerbosity = 2;

        // Special handling is required for some ControlFile messages
        if (globalVerbosity >= 0 && loudness < 1)
            cout << "Preferences file: " << globalControl->path() << ".\n";
        globalControl->setVerbosity(globalVerbosity);	// For "Saved preferences file"
    }

    int windowPropsPosn[WIN_PROPS_LEN][2];
    bool windowPropsOpen[WIN_PROPS_LEN];

    // Read window locations and set their defaults
    windowPropsPosn[MAIN_WIN][0] = globalControl->read(MainWindow_X, 638);
    windowPropsPosn[MAIN_WIN][1] = globalControl->read(MainWindow_Y, 0);
    windowPropsPosn[MODEL_WIN][0] = globalControl->read(ModelWindow_X, 638);
    windowPropsPosn[MODEL_WIN][1] = globalControl->read(ModelWindow_Y, 88);
    windowPropsPosn[DISPLAY_CONTROL_WIN][0] = globalControl->read(DisplayControlWindow_X, 0);
    windowPropsPosn[DISPLAY_CONTROL_WIN][1] = globalControl->read(DisplayControlWindow_Y, 202);
    windowPropsPosn[VISIBILITY_CONTROL_WIN][0] = globalControl->read(VisibilityControlWindow_X, 0);
    windowPropsPosn[VISIBILITY_CONTROL_WIN][1] = globalControl->read(VisibilityControlWindow_Y, 437);
    windowPropsPosn[CONSTRAINTS_WIN][0] = globalControl->read(ConstraintsWindow_X, 20);
    windowPropsPosn[CONSTRAINTS_WIN][1] = globalControl->read(ConstraintsWindow_Y, 100);
    windowPropsPosn[ATOM_EDITOR_WIN][0] = globalControl->read(AtomEditorWindow_X, 190);
    windowPropsPosn[ATOM_EDITOR_WIN][1] = globalControl->read(AtomEditorWindow_Y, 470);
    windowPropsPosn[CUT_PLANES_CONTROL_WIN][0] = globalControl->read(CutPlanesControlWindow_X, 390);
    windowPropsPosn[CUT_PLANES_CONTROL_WIN][1] = globalControl->read(CutPlanesControlWindow_Y, 330);
    windowPropsPosn[OPTIMIZER_CONTROL_WIN][0] = globalControl->read(OptimizerControlWindow_X, 336);
    windowPropsPosn[OPTIMIZER_CONTROL_WIN][1] = globalControl->read(OptimizerControlWindow_Y, 0);
    windowPropsPosn[ATOM_PLANE_WIN][0] = globalControl->read(AtomPlaneWindow_X, 10);
    windowPropsPosn[ATOM_PLANE_WIN][1] = globalControl->read(AtomPlaneWindow_Y, 574);
    windowPropsPosn[CREST_PLANE_WIN][0] = globalControl->read(CrestPlaneWindow_X, 300);
    windowPropsPosn[CREST_PLANE_WIN][1] = globalControl->read(CrestPlaneWindow_Y, 574);
    windowPropsPosn[BPERP_N_PLANE_WIN][0] = globalControl->read(BPerpNPlaneWindow_X, 200);
    windowPropsPosn[BPERP_N_PLANE_WIN][1] = globalControl->read(BPerpNPlaneWindow_Y, 559);
    windowPropsPosn[PORT_PLANE_WIN][0] = globalControl->read(PortPlaneWindow_X, 880);
    windowPropsPosn[PORT_PLANE_WIN][1] = globalControl->read(PortPlaneWindow_Y, 574);
    windowPropsPosn[STARBOARD_PLANE_WIN][0] = globalControl->read(StarboardPlaneWindow_X, 590);
    windowPropsPosn[STARBOARD_PLANE_WIN][1] = globalControl->read(StarboardPlaneWindow_Y, 574);
    windowPropsPosn[INVOLUTES_PLANE_WIN][0] = globalControl->read(InvolutesPlaneWindow_X, 490);
    windowPropsPosn[INVOLUTES_PLANE_WIN][1] = globalControl->read(InvolutesPlaneWindow_Y, 559);
    windowPropsPosn[ADD_QUAD_FIGURE_WIN][0] = globalControl->read(AddQuadFigureWindow_X, 648);
    windowPropsPosn[ADD_QUAD_FIGURE_WIN][1] = globalControl->read(AddQuadFigureWindow_Y, 52);
    windowPropsPosn[ATTACH_SUBFIGURE_WIN][0] = globalControl->read(AttachSubfigureWindow_X, 390);
    windowPropsPosn[ATTACH_SUBFIGURE_WIN][1] = globalControl->read(AttachSubfigureWindow_Y, 70);
    windowPropsPosn[PREFERENCES_EDITOR_WIN][0] = globalControl->read(PreferencesEditorWindow_X, 280);
    windowPropsPosn[PREFERENCES_EDITOR_WIN][1] = globalControl->read(PreferencesEditorWindow_Y, 500);
    windowPropsPosn[IMAGE_PREFS_EDITOR_WIN][0] = globalControl->read(ImagePrefsSubEditorWindow_X, 340);
    windowPropsPosn[IMAGE_PREFS_EDITOR_WIN][1] = globalControl->read(ImagePrefsSubEditorWindow_Y, 577);
    windowPropsPosn[EDIT_MODEL_PROPERTIES_WIN][0] = globalControl->read(EditModelPropsWindow_X, 452);
    windowPropsPosn[EDIT_MODEL_PROPERTIES_WIN][1] = globalControl->read(EditModelPropsWindow_Y, 470);
    windowPropsPosn[EDIT_LANDMARKS_WIN][0] = globalControl->read(EditLandmarksWindow_X, 80);
    windowPropsPosn[EDIT_LANDMARKS_WIN][1] = globalControl->read(EditLandmarksWindow_Y, 490);
    windowPropsPosn[ELONGATION_WIN][0] = globalControl->read(ElongationWindow_X, 680);
    windowPropsPosn[ELONGATION_WIN][1] = globalControl->read(ElongationWindow_Y, 670);
    windowPropsPosn[ABOUT_PABLO_WIN][0] = globalControl->read(AboutPablo_X, 308);
    windowPropsPosn[ABOUT_PABLO_WIN][1] = globalControl->read(AboutPablo_Y, 0);
    windowPropsPosn[OPTIMIZER_SETTINGS_WIN][0] = globalControl->read(PenaltyWeightsWindow_X, 336);
    windowPropsPosn[OPTIMIZER_SETTINGS_WIN][1] = globalControl->read(PenaltyWeightsWindow_Y, 640);
    windowPropsPosn[MODEL_SLIDE_SHOW_WIN][0] = globalControl->read(SlideShowWindow_X, 100);
    windowPropsPosn[MODEL_SLIDE_SHOW_WIN][1] = globalControl->read(SlideShowWindow_Y, 100);
    windowPropsPosn[CPNS_DEFORMATION_WIN][0] = globalControl->read(CPNSDeformWindow_X, 50);
    windowPropsPosn[CPNS_DEFORMATION_WIN][1] = globalControl->read(CPNSDeformWindow_Y, 120);
    windowPropsPosn[PGA_DEFORMATION_WIN][0] = globalControl->read(PGADeformWindow_X, 50);
    windowPropsPosn[PGA_DEFORMATION_WIN][1] = globalControl->read(PGADeformWindow_Y, 120);
#ifdef PCA_VISUALIZER
    windowPropsPosn[PCA_DEFORMATION_WIN][0] = globalControl->read(PCADeformWindow_X, 60);
    windowPropsPosn[PCA_DEFORMATION_WIN][1] = globalControl->read(PCADeformWindow_Y, 130);
#endif
#ifdef OPTIMIZATION_VISUALIZER
    windowPropsPosn[OPT_VISUALIZER_WIN][0] = globalControl->read(OptVisualizerWindow_X, 0);
    windowPropsPosn[OPT_VISUALIZER_WIN][1] = globalControl->read(OptVisualizerWindow_Y, 100);
#endif

    // Test menu windows
#ifndef PRODUCTION_VERSION
    windowPropsPosn[INTERPOLATED_PRIMITIVE_WIN][0] =
        globalControl->read(InterpolatedPrimitiveWindow_X, 190);
    windowPropsPosn[INTERPOLATED_PRIMITIVE_WIN][1] =
        globalControl->read(InterpolatedPrimitiveWindow_Y, 630);
    windowPropsPosn[REGULARIZER_WIN][0] =
        globalControl->read(RegularizerWindow_X, 100);
    windowPropsPosn[REGULARIZER_WIN][1] =
        globalControl->read(RegularizerWindow_Y, 230);
    windowPropsPosn[MATCH_SURFACES_WIN][0] =
        globalControl->read(MatchSurfaces_X, 120);
    windowPropsPosn[MATCH_SURFACES_WIN][1] =
        globalControl->read(MatchSurfaces_Y, 250);
#endif

    // Set most remaining parameter defaults (see below).
    // Note that all defaults must be set before the interface appears,
    // because the user may immediately try to edit the preferences.  Many
    // of these also have to be set before the P3DUserInterface object is
    // constructed.
    setGlobalControlDefaults();

    switch (compressImage)
    {
    case 2:
        globalControl->write(CompressImages, false);
        // Fall through
    case 1:
        globalControl->write(ConvertImages, true);
        break;
    }

    if (convertFormat) {
        globalControl->write(ConvertImageFormat, true);
        globalControl->write(ImageFormat, newVersion);
    }
    else
        globalControl->write(ConvertImageFormat, false);

    // The main and model windows are always open.  The constraints window cannot be open
    // initially, because a primitive must be marked for this.  The Edit Model Properties
    // window is modal and can only be open when in use.
    windowPropsOpen[DISPLAY_CONTROL_WIN] = globalControl->read(DisplayControlWindow_Open, false);
    windowPropsOpen[VISIBILITY_CONTROL_WIN] = globalControl->read(VisibilityControlWindow_Open, false);
    windowPropsOpen[CONSTRAINTS_WIN] = globalControl->read(ConstraintsWindow_Open, false);
    windowPropsOpen[ATOM_EDITOR_WIN] = globalControl->read(AtomEditorWindow_Open, false);
    windowPropsOpen[CUT_PLANES_CONTROL_WIN] = globalControl->read(CutPlanesControlWindow_Open, false);
    if (interactive)
        windowPropsOpen[OPTIMIZER_CONTROL_WIN] = true;
    else
        windowPropsOpen[OPTIMIZER_CONTROL_WIN] = globalControl->read(OptimizerControlWindow_Open, false);
    windowPropsOpen[CREST_PLANE_WIN] = globalControl->read(CrestPlaneWindow_Open, false);
    windowPropsOpen[ATOM_PLANE_WIN] = globalControl->read(AtomPlaneWindow_Open, false);
    windowPropsOpen[BPERP_N_PLANE_WIN] = globalControl->read(BPerpNPlaneWindow_Open, false);
    windowPropsOpen[PORT_PLANE_WIN] = globalControl->read(PortPlaneWindow_Open, false);
    windowPropsOpen[STARBOARD_PLANE_WIN] = globalControl->read(StarboardPlaneWindow_Open, false);
    windowPropsOpen[INVOLUTES_PLANE_WIN] = globalControl->read(InvolutesPlaneWindow_Open, false);
    windowPropsOpen[ADD_QUAD_FIGURE_WIN] = globalControl->read(AddQuadFigureWindow_Open, false);
    windowPropsOpen[ATTACH_SUBFIGURE_WIN] = globalControl->read(AttachSubfigureWindow_Open, false);
    windowPropsOpen[PREFERENCES_EDITOR_WIN] = globalControl->read(PreferencesEditorWindow_Open, false);
    windowPropsOpen[ELONGATION_WIN] = globalControl->read(ElongationWindow_Open, false);
    windowPropsOpen[ABOUT_PABLO_WIN] = globalControl->read(AboutPablo_Open, allow_aboutPablo);
    windowPropsOpen[OPTIMIZER_SETTINGS_WIN] = globalControl->read(PenaltyWeightsWindow_Open, false);
    windowPropsOpen[MODEL_SLIDE_SHOW_WIN] = globalControl->read(SlideShowWindow_Open, false);
	windowPropsOpen[CPNS_DEFORMATION_WIN] = globalControl->read(CPNSDeformWindow_Open, false);
    windowPropsOpen[PGA_DEFORMATION_WIN] = globalControl->read(PGADeformWindow_Open, false);
#ifdef OPTIMIZATION_VISUALIZER
    windowPropsOpen[OPT_VISUALIZER_WIN] = globalControl->read(OptVisualizerWindow_Open, false);
#endif

    // Test menu windows
#ifndef PRODUCTION_VERSION
    windowPropsOpen[INTERPOLATED_PRIMITIVE_WIN] =
        globalControl->read(InterpolatedPrimitiveWindow_Open, false);
    windowPropsOpen[REGULARIZER_WIN] = globalControl->read(RegularizerWindow_Open, false);
    windowPropsOpen[MATCH_SURFACES_WIN] = globalControl->read(MatchSurfaces_Open, false);
#endif

    // Set the default smoothness parameter for figures
    M3DFigure::setDefaultSurfaceTolerance(globalControl->readInt(SurfaceSmoothnessDefault));

    // Constructing the interface causes P3DView::(*modelview) to be constructed,
    // which results in the setting of defaults for parameters SurfaceStyle,
    // SurfaceLevel, PartialSurfaceLevel, RockingAngle, RockingIncrement,
    // BackgroundColor_R, BackgroundColor_G, BackgroundColor_B, LandmarksColor_R,
    // LandmarksColor_G, LandmarksColor_B, LandmarkNarrowWidth, and LandmarkWideWidth.
    // The P3DUserInterfaceCallback::init() function sets the PartialSurfaceRendering
    // parameter default.
    ui = new P3DUserInterface(windowPropsPosn, windowPropsOpen);

#ifdef OPTIMIZATION_VISUALIZER
    ovUI = new OptVisualizerUI(windowPropsPosn[OPT_VISUALIZER_WIN][0], windowPropsPosn[OPT_VISUALIZER_WIN][1]);
    ovUI->setP3DUserInterface(ui);
    ui->callback->setOptVisualizerUI(ovUI);
    ovUI->init();
    ui->windowsOptVisWin->show();	// Show visualizer on Windows Menu
#endif

#ifdef BINARY

    if (convertAltFormat) {
        // Load and save with a different name, inferring file type
        bool ok = control->loadImage(imageFile.c_str());
        if (! ok) {
            cout << "Error: did not load image " << imageFile << endl;
            return -1;	// fatal
        }
        ok = control->saveImage(outputImageFile.c_str());
        if (! ok) {
            cout << "Error: did not save image " << outputImageFile << endl;
            return -1;	// fatal
        }
        // if (globalVerbosity > 0)
        std::cout << "Converted " << imageFile 
                  << " to " << outputImageFile << std::endl;
        // Automatically noninteractive
        return 0;
    }

    if (compressImage || convertFormat) {
        if (imageFile.length() > 0)
            ui->callback->loadImage(imageFile.data(), true);	// compress/decompress image
        delete globalControl;
        return 0;
    }

    // COMPUTE AND SAVE THE MODEL TILES, THEN EXIT			
    // AGG: This appears to be non-interactive code that should
    // precede allocation of ui
    bool loadedImage = false;
    const char * filename = scriptParms.getStringValue("OutImageTile", NULL);
    if (filename && *filename)		// to gen tiles, read image voxels
    {
        if (imageFile.length() == 0)
        {
            cout << "Error: -outImageTile requires -image" << endl;
            return 1;
        }

        // TODO: add bitStack mask here
        ui->callback->loadImage(imageFile.data(), false); // do not compress/decompress image		
        loadedImage = true;
        ui->callback->generateTileSet();
        if (control != NULL)
            delete control;
        control = ui->callback->getControl();
        //		ui->callback->saveTileSet(filename);	AGG: original version in BP - Do I need a saveTileSet fn in callback?
        control->saveTileSet(filename, true);
        if (nonInteractive) {
            if (globalVerbosity > 0)
                cout << "Option -pQ inhibits interactive mode" << endl;
            return 0;		// OK return
        }
    }

    if (nonInteractive)
        return 1;

    // Load the image
    // TODO: add bitStack mask here
    if (imageFile.length() > 0 && ! loadedImage && scriptFile == NULL)
        ui->callback->loadImage(imageFile.data(), false);	// do not compress/decompress image

    if (modelFile.length() > 0 && scriptFile == NULL)
        if (! ui->callback->loadModel(modelFile.data()))
            return 1;

    if (landmarkModelFile && landmarkModelFile != NULL) {
        // seems to be a display bug -- warn user
        cout << "Warning: -il has display bugs; use the File Menu instead" << endl;
        ui->callback->loadLandmarkModel(landmarkModelFile);
    }

#else	/* BINARY */

    if (imageFile.length() > 0  && scriptFile == NULL)
        ui->callback->loadImage(imageFile.data(), false, false, compressImage || convertFormat);

    if (regularize) {
        if (! ui->callback->loadModel(modelFile.data()))
            return 1;
        ui->callback->regularizeModel(reg_iterations, reg_steps);
        ui->callback->saveModel();
    }

    if (compressImage || convertFormat || regularize) {
        delete globalControl;
        return 0;
    }

    if (modelFile.length() > 0 && scriptFile == NULL)
        if (! ui->callback->loadModel(modelFile.data()))
            return 1;

#endif	/* BINARY */

    if (xformFile.length() > 0 && scriptFile == NULL)
        if (! ui->callback->loadSimTransform(xformFile.data()))
            return 1;

#ifndef _WIN32
    if (left_handed_mouse)
        ui->callback->set_mouse_button_order(true);
#endif

    Fl::gl_visual(FL_DOUBLE | FL_RGB | FL_ALPHA | FL_DEPTH);
    Fl::visual(FL_ALPHA | FL_RGB | FL_DOUBLE);
    Fl::add_handler(globalHandler);

    ui->show();
#ifdef OPTIMIZATION_VISUALIZER
    if (windowPropsOpen[OPT_VISUALIZER_WIN])
        ui->callback->openOptimizationVisualizer();
#endif

    // Pablo should run without previewing .m3d files
    Fl_File_Chooser fc(".", "*.m3d", Fl_File_Chooser::SINGLE, NULL);
    int previewStatus = fc.preview();	// Save user's default preview status
    fc.preview(0);

    Fl_Tooltip::delay(0.0f);	// Delay before tooltip becomes visible

    ui->mainWindow->setHideCallback(check_iconify);

#ifdef TASKING
    Fl_Text_Buffer * textBuffer = new Fl_Text_Buffer;
    if (taskFileName.length() > 0) {
        textBuffer->appendfile(taskFileName.data());
        ui->taskDisplayText->buffer(textBuffer);
    }
#endif

    Fl::add_timeout(DELAY, timeoutCallback);

    if (interactive) {
        int matchMethod;
#ifdef BINARY
        matchMethod = 7;	// Binary: See OptimizerMatchMethod in gpTuning.cpp
#else
        matchMethod = (int) tuningWt(GpOptimizerMatchMethod);
#endif
        ui->callback->startInteractive(matchMethod);
    }

    bool ovFlag;
#ifdef OPTIMIZATION_VISUALIZER
    // If the optimization visualizer is turned on and a script file is
    // specified, then optimize.  Note: -po should not be used.
    ovFlag = true;
#else
    ovFlag = false;
#endif  /* OPTIMIZATION_VISUALIZER */

    if (interactive || ovFlag) {
#ifdef BINARY
        if (! ui->callback->getControl()->runBinaryPablo(scriptParms))
#else
        if (! ui->callback->getControl()->runPablo(scriptParms))
#endif  /* BINARY */
            cout << "Error: Script failed\n";
        if (interactive)
            ui->callback->stopInteractive();
    }

    int retn = Fl::run();	
#ifdef DEBUG
    cout << "run complete" << endl;
#endif

    fc.preview(previewStatus);	// Restore user's default preview status

    Fl::remove_timeout(timeoutCallback);

    deleteUserInterface(ui);
    delete ui;
#ifdef OPTIMIZATION_VISUALIZER
    delete ovUI;
#endif
    delete globalControl;

    //extern void testImage3D();
    //testImage3D();

    return retn;
}

