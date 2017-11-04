/*
	This is a program to determine the latest file modification date
	and time for a source tree, by scanning the Entries files in the
	CVS or .svn directories found in the tree.

	It was written for use in building Pablo, and for that must be run
	before Pablo is compiled to create the version file loaded by Pablo.
	Whenever Pablo is updated, it should be rerun.

	The first input consists of the name of the output file.  The
	remaining inputs comprise a list of the Entries file paths.  The
	output is a printed string giving the date and time of the last file
	checked in.

    There is a single option:
	-v variable     The name of the variable to be defined in the
	                output file.  If ommited, 'revision' is used.
	                For building ConStruct, constructRevision is used.

	If the program fails, it exits with a return code of 1.


    Maintenance Status:

	This program has evolved from earlier versions, starting with one
	for CVS and then going to XML-based SVN.  Originally all code was
	intermingled.  To simplify maintenance, it is separated in this
	version.  However, the OLD_SVN code has not been tested.  The
	SVN version has only been tested with subversion 1.6.6.  The CVS
	version was also tested.
*/


#include <iostream>
#include <fstream>
#include <strstream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

// -------------------------------
// Exactly one of the following 3 definitions should be uncommented

// Uncomment to build this program for CVS files
//#define CVS	

// Uncomment to build this program for old-style (XML) SVN files
//#define OLD_SVN
	
// Uncomment to build this program for newer SVN files
#define SVN	
// -------------------------------


//#define DEBUG


#define MAX_WORD_SIZE	64		// Length of a word
#define MAX_REC_SIZE	1024	// Bytes per full record


using namespace std;


unsigned long max_date = 0;
unsigned long max_time = ULONG_MAX;
#ifndef CVS
int revision;
#endif

void next_line(ifstream & infile) {
	char c;

	if (infile.eof())
		return;

	infile.get(c);
	while (! (c == '\n' || c == '\r') && ! infile.eof())
		infile.get(c);
}

char * findValidExt(const char * name)
{
	char * ext;

	// Locate the source-code extension
	ext = (char *) strstr(name, ".cpp");	// Solaris CC requires these casts
	if (ext == NULL) {
		ext = (char *) strstr(name, ".h");
		if (ext == NULL) {
			ext = (char *) strstr(name, ".c");
			if (ext != NULL) {
				// Have .c
				if (ext[2] != '\0')
					ext = NULL;
			}
		}
		else {
			// Have .h
			if (ext[2] != '\0')
				ext = NULL;
		}
	}
	else {
		// Have .cpp
		if (ext[4] != '\0')
			ext = NULL;
	}
	return ext;
}


// --------------------------- Start CVS SECTION ---------------------------

#ifdef CVS	

bool line(ifstream & infile) {
	char c;
	char name[MAX_WORD_SIZE];
	int i;
	char * ext;

	infile >> c;

	// Consider only lines beginning with a slash
	if (c == '/') {

		// Get the file name
		infile >> name[0];
		i = 0;
		while (name[i] != '/')
			infile >> name[++i];
		name[i] = '\0';

		// Check the extension
		ext = findValidExt(name);

		// Skip files not containing code
		if (ext == NULL)
			next_line(infile);
		else {
			// Found a record to be examined
			return true;
		}
	}
	else
		next_line(infile);

	return false;
}

void cvs_get_date_time(ifstream & infile, char * month, char * day, char * year, char * time) {
	char c;
	int i;
	char str[8];

	// A sample CVS date/time entry:  Fri Jun  9 20:39:48 2006

	// Skip revision number
	infile >> c;
	while (c != '/')
		infile >> c;

	// Skip day of the week, but test it for "Result of merge"
	i = 0;
	infile.get(c);
	if (c == '*') {
		// Locally added file
		month[i] = '\0';
		day[i] = '\0';
		time[i] = '\0';
		year[i] = '\0';
		return;
	}
	while (c != ' ') {
		str[i++] = c;
		infile.get(c);
	}
	str[i] = '\0';
	if (0 == strcmp(str, "Initial")) {
		// Locally added file
		month[i] = '\0';
		day[i] = '\0';
		time[i] = '\0';
		year[i] = '\0';
		return;
	}
	if (0 == strcmp(str, "Result")) {
		// Skip to the day of the week
		while (c != '+') {
			infile.get(c);
			if (c == '\n') {
				month[i] = '\0';
				day[i] = '\0';
				time[i] = '\0';
				year[i] = '\0';
				return;
			}
		}
		// Skip the day of the week
		while (c != ' ')
			infile.get(c);
	}

	infile >> month[0];
	i = 0;
	while (month[i] != ' ')
		infile.get(month[++i]);
	month[i] = '\0';

	infile >> day[0];
	i = 0;
	while (day[i] != ' ')
		infile.get(day[++i]);
	day[i] = '\0';

	infile >> time[0];
	i = 0;
	while (time[i] != ' ')
		infile.get(time[++i]);
	time[i] = '\0';

	infile >> year[0];
	i = 0;
	while (year[i] != '/')
		infile >> year[++i];
	year[i] = '\0';
}

unsigned long cvs_reformat(char * month, char * day, char * year, char * timeStr,
	unsigned long & time)
{
	char * months[] = {
		"Jan", "Feb", "Mar", "Apr", "May", "Jun",
		"Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
	};
	char str[3];
	int m, d, y;
	unsigned long date;

	// Compute an integer representation of the date
	for (m = 0; m < 12; m++)
		if (0 == strcmp(month, months[m]))
			break;
	m++;
	d = atoi(day);
	y = atoi(year);
	date = y*10000 + m*100 + d;

	// Compute an integer representation of the time.
	str[0] = timeStr[0];
	str[1] = timeStr[1];
	str[2] = '\0';
	time = 10000*atoi(str);
	str[0] = timeStr[3];
	str[1] = timeStr[4];
	time += 100*atoi(str);
	str[0] = timeStr[6];
	str[1] = timeStr[7];
	time += atoi(str);

#ifdef DEBUG
    cout << year << '-' << month << '-' << day << 'T' << timeStr
        << " --> " << date << ' ' << time << endl;
#endif
	return date;
}

bool process_cvs(int argc, char ** argv,
	char * month, char * day, char * time, char * year)
{
	char c;
	unsigned long d, t;
	int i;

	for (i = 0; i < argc; i++) {

		ifstream infile(argv[i]);
		if (! infile) {
			cout << "Cannot open file " << argv[i];
			return false;
		}

#ifdef DEBUG
		cout << "Command-line argument " << i << ": " << argv[i] << '\n';
#endif
		while (1) {
			infile >> c;
			if (infile.eof())   // This will skip blank lines
				break;
			else infile.putback(c);

			if (line(infile)) {
				cvs_get_date_time(infile, month, day, year, time);
//				cout << "\tmonth=" << month << " day=" << day << " year=" << year << " time=" << time << '\n';
				d = cvs_reformat(month, day, year, time, t);
#ifdef DEBUG
				cout << "date = " << d << " time = " << t << "   max_date = "
					<< max_date << " max_time = " << max_time << '\n';
#endif
				if (d > max_date) {
					max_date = d;
					max_time = t;
				}
				else {
					if (d == max_date) {
						if (t > max_time) {
							max_time = t;
						}
					}
				}
				next_line(infile);
			}
		}
	}
#ifdef DEBUG
	cout << "Final max_date = " << max_date << '\n';
	cout << "Final max_time = " << max_time << '\n';
#endif

	return true;
}

#endif	/* CVS */

// --------------------------- End CVS SECTION ---------------------------


// --------------------------- Start OLD_SVN SECTION ---------------------------

#ifdef OLD_SVN

char rec[MAX_REC_SIZE];
char com_date[MAX_WORD_SIZE];


// Verify that a concatenated XML record is for a source file
bool extract()
{
	char * s1;
	char * s2;
	char * ext;
	char c;
	int n;

	ext = NULL;
	if (0 != strncmp(rec, "entry", 5))
		return false;
	s1 = rec + 5;
	while (*s1 == ' ')
		s1++;

	n = 0;
	while (*s1) {
		if (0 == strncmp(s1, "committed-rev=\"", 15)) {
			s1 += 15;
			s2 = s1;
			while (*s2 != '\"' && *s2 != '\0')
				s2++;
			c = *s2;
			*s2 = '\0';
			int r = atoi(s1);
			if (r > revision)
				revision = r;
			*s2 = c;
			s1 = s2;
			n++;
		}
		else if (0 == strncmp(s1, "name=\"", 6)) {
			s1 += 6;
			s2 = s1;
			while (*s2 != '\"' && *s2 != '\0')
				s2++;
			c = *s2;
			*s2 = '\0';
			ext = findValidExt(s1);
			*s2 = c;
			s1 = s2;
			if (ext != NULL)
				n++;
		}
		else if (0 == strncmp(s1, "committed-date=\"", 16)) {
			s1 += 16;
			s2 = s1;
			while (*s2 != '\"' && *s2 != '\0')
				s2++;
			c = *s2;
			*s2 = '\0';
			strcpy(com_date, s1);
			*s2 = c;
			s1 = s2;
			n++;
		}
		// The following XML fields are skipped for speed
		else if (0 == strncmp(s1, "url=\"", 5)) {
			s1 += 5;
			while (*s1 != '\"' && *s1 != '\0')
				s1++;
		}
		else if (0 == strncmp(s1, "last-author=\"", 13)) {
			s1 += 5;
			while (*s1 != '\"' && *s1 != '\0')
				s1++;
		}
		else if (0 == strncmp(s1, "text-time=\"", 11)) {
			s1 += 5;
			while (*s1 != '\"' && *s1 != '\0')
				s1++;
		}
		else if (0 == strncmp(s1, "prop-time=\"", 11)) {
			s1 += 5;
			while (*s1 != '\"' && *s1 != '\0')
				s1++;
		}
		else if (0 == strncmp(s1, "kind=\"", 6)) {
			s1 += 5;
			while (*s1 != '\"' && *s1 != '\0')
				s1++;
		}
		else if (0 == strncmp(s1, "checksum=\"", 10)) {
			s1 += 5;
			while (*s1 != '\"' && *s1 != '\0')
				s1++;
		}
		if (n == 3)
			break;
		s1++;
		while (*s1 == ' ')
			s1++;
	}	
	return (ext == NULL ? false : true);
}

bool xml_line(ifstream & infile) {
	char c;
	int i;

	infile >> c;	// Skip whitespace and get first char

	// Start with lines beginning with a less than symbol
	if (c == '<') {

		// Concatenate all the lines of the record
		infile >> rec[0];
		i = 0;
		while (rec[i] != '<') {
			if (rec[i] != '\n' && rec[i] != '\r')
				i++;
			if (i == MAX_REC_SIZE - 1) {
				cout << "Error in version.exe: buffer overflow" << endl;
				break;
			}
			// The interior lines are assumed to begin with a space
			infile.get(rec[i]);
			if (infile.eof())
				break;
		}
		infile.putback(rec[i]);
		rec[i] = '\0';

		// Verify that the XML record is for a source file
		return extract();
	}
	else
		next_line(infile);

	return false;
}

void xml_get_date_time(ifstream & infile, char * month, char * day, char * year,
	char * time)
{
	char * s1;
	char * s2;

	// The SVN (XML) date/time format is yyyy-mm-ddThh:mm:ss.ssssssZ
	s1 = strchr(com_date, 'T');
	if (s1 == NULL) {
		cout << "Error in version: invalid date/time" << endl;
		time[0] = '\0';
		month[0] = '\0';
		day[0] = '\0';
		year[0] = '\0';
		return;
	}
	*s1 = '\0';

	strncpy(year, com_date, 4);
	year[4] = '\0';

	strncpy(month, com_date + 5, 2);
	month[2] = '\0';

	strcpy(day, com_date + 8);

	s1++;
	// Ignore the fractional seconds and time zone of SVN time fields
	s2 = s1 + 1;
	while (*s2 != '.')
		s2++;
	*s2 = '\0';
	strcpy(time, s1);
}

unsigned long old_svn_reformat(char * month, char * day, char * year, char * timeStr,
	unsigned long & time)
{
	char * months[] = {
		"Jan", "Feb", "Mar", "Apr", "May", "Jun",
		"Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
	};
	char str[3];
	int m, d, y;
	unsigned long date;

	// Compute an integer representation of the date
	for (m = 0; m < 12; m++)
		if (0 == strcmp(month, months[m]))
			break;
	m++;
	d = atoi(day);
	y = atoi(year);
	date = y*10000 + m*100 + d;

	// Compute an integer representation of the time.
	str[0] = timeStr[0];
	str[1] = timeStr[1];
	str[2] = '\0';
	time = 10000*atoi(str);
	str[0] = timeStr[3];
	str[1] = timeStr[4];
	time += 100*atoi(str);
	str[0] = timeStr[6];
	str[1] = timeStr[7];
	time += atoi(str);

	// Increment to cover the lost fractional seconds of SVN time fields
	time++;

#ifdef DEBUG
    cout << year << '-' << month << '-' << day << 'T' << timeStr
        << " --> " << date << ' ' << time << endl;
#endif
	return date;
}

bool process_old_svn(int argc, char ** argv,
	char * month, char * day, char * time, char * year)
{
	char c;
	unsigned long d, t;
	int i;

	for (i = 0; i < argc; i++) {

		ifstream infile(argv[i]);
		if (! infile) {
			cout << "Cannot open file " << argv[i];
			return false;
		}

#ifdef DEBUG
		cout << "Command-line argument " << i << ": " << argv[i] << '\n';
#endif
		while (1) {
			infile >> c;
			if (infile.eof())   // This will skip blank lines
				break;
			else infile.putback(c);

			if (c == '<') {
				cout << "File " << argv[i] << " does not appear to be an old-style SVN file\n";
				return false;
			}

			if (xml_line(infile)) {
				xml_get_date_time(infile, month, day, year, time);
//				cout << "\tmonth=" << month << " day=" << day << " year=" << year << " time=" << time << '\n';
				d = old_svn_reformat(month, day, year, time, t);
#ifdef DEBUG
				cout << "date = " << d << " time = " << t << "   max_date = "
					<< max_date << " max_time = " << max_time << '\n';
#endif
				if (d > max_date) {
					max_date = d;
					max_time = t;
				}
				else {
					if (d == max_date) {
						if (t > max_time) {
							max_time = t;
						}
					}
				}
			}
		}
	}
#ifdef DEBUG
	cout << "Final max_date = " << max_date << '\n';
	cout << "Final max_time = " << max_time << '\n';
#endif

	return true;
}

#endif	/* OLD_SVN */

// --------------------------- End OLD_SVN SECTION ---------------------------


// --------------------------- Start SVN SECTION ---------------------------

#ifdef SVN

char rec[MAX_REC_SIZE];
char com_date[MAX_WORD_SIZE];


bool checkForDate()
{
    int len = strlen(rec);
    if (len < 5) {
        for (int i = 0; i < len; i++)
            if (! isdigit(rec[i]))
                return false;
    }

    if (rec[4] != '-')
        return false;
    if (rec[len - 1] != 'Z')    // Requiring time zone Z
        return false;
    
    if (len >= MAX_WORD_SIZE)
        return false;

    strcpy(com_date, rec);  // For get_date_time()
    return true;    // Located a date line
}

// This function is used for the newer version of SVN entries files
bool line(ifstream & infile) {
	int i;

	infile >> rec[0];	// Skip whitespace and get first char

	// Look for lines that start with the year
	if (isdigit(rec[0])) {

		// Read the entire record
		i = 0;
		while (rec[i] != '\n' && rec[i] != '\r') {
			i++;
			if (i == MAX_REC_SIZE - 1) {
				cout << "Error in version.exe: buffer overflow" << endl;
				break;
			}
			infile.get(rec[i]); // Do not ignore newlines
			if (infile.eof())
				break;
		}
		infile.putback(rec[i]);
		rec[i] = '\0';

        return checkForDate();
	}
	else
		next_line(infile);

	return false;
}

void get_date_time(ifstream & infile, char * month, char * day, char * year,
	char * time)
{
	char * s1;
	char * s2;

	// The SVN (XML-like) date/time format is yyyy-mm-ddThh:mm:ss.ssssssZ
	s1 = strchr(com_date, 'T');
	if (s1 == NULL) {
		cout << "Error in version: invalid date/time" << endl;
		time[0] = '\0';
		month[0] = '\0';
		day[0] = '\0';
		year[0] = '\0';
		return;
	}
	*s1 = '\0';

	strncpy(year, com_date, 4);
	year[4] = '\0';

	strncpy(month, com_date + 5, 2);
	month[2] = '\0';

	strcpy(day, com_date + 8);

	s1++;
	// Ignore the fractional seconds and time zone of SVN time fields
	s2 = s1 + 1;
	while (*s2 != '.')
		s2++;
	*s2 = '\0';
	strcpy(time, s1);
}

unsigned long reformat(char * month, char * day, char * year, char * timeStr,
	unsigned long & time)
{
	char str[3];
	int m, d, y;
	unsigned long date;

	m = atoi(month);
	d = atoi(day);
	y = atoi(year);
	date = y*10000 + m*100 + d;

	// Compute an integer representation of the time.
	str[0] = timeStr[0];
	str[1] = timeStr[1];
	str[2] = '\0';
	time = 10000*atoi(str);
	str[0] = timeStr[3];
	str[1] = timeStr[4];
	time += 100*atoi(str);
	str[0] = timeStr[6];
	str[1] = timeStr[7];
	time += atoi(str);

#ifdef DEBUG
        cout << year << '-' << month << '-' << day << 'T' << timeStr
            << " --> " << date << ' ' << time << endl;
#endif
	return date;
}

bool process_svn(int argc, char ** argv,
	char * month, char * day, char * time, char * year)
{
	char c;
	unsigned long d, t;
	string str;
	bool ok, done;
	int i;

	for (i = 0; i < argc; i++) {

		ifstream infile(argv[i]);
		if (! infile) {
			cout << "Cannot open file " << argv[i];
			return false;
		}

#ifdef DEBUG
		cout << "Command-line argument " << i << ": " << argv[i] << '\n';
#endif
		infile >> c;
		if (infile.eof())   // This will skip blank lines
			break;
		else infile.putback(c);

		if (! isdigit(c)) {
			cout << "File " << argv[i] << " does not appear to be an SVN file\n";
			return 1;
		}
		next_line(infile);

		ok = true;
		done = false;
		str.erase();
		while (! done) {
			if (str == "dir")
				done = true;
			if (! getline(infile, str)) {
				ok = false;
				cout << "File " << argv[i] << " is too short\n";
				break;
			}
		}
		if (! ok)
			continue;

		int rev = atoi(str.c_str());
//		cout << "File " << argv[i] << " has revision " << rev << '\n';
		// Normally all files will have the same revision number, but someone
		// may have copied an old one in
		if (rev > revision)
				revision = rev;

		while (! infile.eof()) {
			if (line(infile)) {
				get_date_time(infile, month, day, year, time);
//				cout << "\tmonth=" << month << " day=" << day << " year=" << year << " time=" << time << '\n';
				d = reformat(month, day, year, time, t);
#ifdef DEBUG
				cout << "date = " << d << " time = " << t << "   max_date = "
					<< max_date << " max_time = " << max_time << '\n';
#endif
				if (d > max_date) {
					max_date = d;
					max_time = t;
				}
				else {
					if (d == max_date) {
						if (t > max_time) {
							max_time = t;
						}
					}
				}
			}
		}
	}
#ifdef DEBUG
	cout << "Final max_date = " << max_date << '\n';
	cout << "Final max_time = " << max_time << '\n';
#endif

	return true;
}

#endif	/* SVN */

// --------------------------- End SVN SECTION ---------------------------


int main(int argc, char ** argv)
{
	char month[4];
	char day[3];
	char time[9];
	char year[5];
	char * var;
	int i;
	char version[17];

	if (argc < 3) {
		cout << "Syntax: version [-v variable] output_file list_of_input_files\n";
		return 1;
	}

    if (argv[1][0] == '-') {
        if (0 != strcmp(argv[1], "-v")) {
            cout << "Unknown option: " << argv[1] << '\n';
            return 1;
        }
        var = argv[2];
        i = 3;
    }
    else {
        var = NULL;
        i = 1;
    }
	ofstream outfile(argv[i]);
	if (! outfile) {
		cout << "Cannot open file " << argv[i] << '\n';
		return 1;
	}
    i++;
	argc -= i;

#ifdef CVS
	if (! process_cvs(argc, argv + i, month, day, time, year))
		return 1;
#endif

#ifdef OLD_SVN
	revision = 0;
	if (! process_old_svn(argc, argv + i, month, day, time, year))
		return 1;
#endif

#ifdef SVN
	revision = 0;
	if (! process_svn(argc, argv + i, month, day, time, year))
		return 1;
#endif

	if (max_date == 0 || max_time == ULONG_MAX)
		return 1;

	// Produce the version as "yyyy/mm/dd hh:mm:ss"
	ostrstream vfile(version, 16);
	vfile << max_date;
    vfile.width(6);
    vfile << max_time;

	outfile << "// This file is generated and should not be edited\n\n";
    if (var != NULL)
		outfile << "static const char * " << var << " = \"";
    else
	    outfile << "static const char * revision = \"";
#ifndef CVS
	outfile << revision << ' ';
#endif
	outfile << version[0] <<  version[1] << version[2] <<  version[3] << '/' << version[4] <<  version[5] << '/'
		<< version[6] <<  version[7] << ' '
		<< version[8] <<  version[9] << ':' << version[10] <<  version[11] << ':'
		<< version[12] <<  version[13];
	outfile << "\";\n\n";

#ifdef DEBUG
	cout << "\nWrote version file:\n\n";
	cout << "// This file is generated and should not be edited\n\n";
    if (var != NULL)
		cout << "static const char * " << var << " = \"";
    else
	    cout << "static const char * revision = \"";
#ifndef CVS
	cout << revision << ' ';
#endif
	cout << version[0] <<  version[1] << version[2] <<  version[3] << '/' << version[4] <<  version[5] << '/'
		<< version[6] <<  version[7] << ' '
		<< version[8] <<  version[9] << ':' << version[10] <<  version[11] << ':'
		<< version[12] <<  version[13];
	cout << "\";\n\n";
#endif

	return 0;
}

