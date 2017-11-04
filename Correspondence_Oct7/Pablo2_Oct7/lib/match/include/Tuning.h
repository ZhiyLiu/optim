#ifndef _TUNING_H
#define _TUNING_H

/*	Tuning, a base class for classes that supply tuning values,
	primarily for optimization, to various configurations of Pablo.
	This capability was originally developed for Binary Pablo by
	Gregg Tracton and Derek Merck.  It was later extended for more
	general use by Graham Gash.

	The optimization tuning parameters for Pablo are settings used
	to compute and control optimizations.  Each parameter has names
	for code and script file references.

	There is only one instance of this class in Pablo, which is
	defined in file P3DControl.cpp as a global variable, so that any
	code may refer to it.

	The Tuning class stores an array of TuneVal_t values, each of
	which contains one of the following items in element 'weight':
		1. A double precision weight for a term to be used somewhere
			in the program (usually in optimization).  Although a
			double is used, this may serve to contain an int or
			even boolean value.
		2. A character string.  These are used mainly to store
			file names or paths.
		3. A list (1D array) of integers.
	A lock flag is provided, so that the calling code knows which
	values may be set by the user.  The lock flag is always checked
	before changing an entry, however attempts to violate it are
	usually ignored without printing an error.

	Regular weights contain a default value, while strings and
	lists do not.  The default for those is NULL.

	In most places in Pablo, this file should be included instead
	of the subclass files, bpTuning.h and gpTuning.h, since they
	are automatically included below.
*/

class TuneVal_t
{
	public:

		char * shortName;		// A 2 or 3 byte name, for compact code
		char * longName;		// A 20 to 30 byte name, for script files
		double defaultWeight;	// Default weight for weights that are doubles;
								// the default for strings and lists is NULL.
		char * desc;			// A short description for users
		union {
			double d;
			int * lst;
			char * str;
		} weight;
		int listLen;			// -1 when weight is a double; -2 if weight is a
								//     string; the length of the list otherwise.
		bool locked;			// A user may not change the weight when locked

		TuneVal_t() {
			shortName = NULL;
			longName = NULL;
			desc = NULL;
			listLen = -1;
		}
		~TuneVal_t() {
			// ShortName, longName and desc do not belong to this class
			if (listLen != -1)
				delete [] weight.str;	// Deletes string or list entries
		}
};


/*	Non-class utility function to print the version number
	to cout as "VVVVVVVV.XX", when the default precision is
	used.  The only other reasonable precision is 0, used in
	pablo.cpp.  The precision and I/O flags are left correct
	by this function.
*/
void printVersion(double weight, int precision = 2);


// The following definition limits the length of stored lists
#define MAX_TUNING_LIST_LEN		16


class Tuning
{
public:

	Tuning();
	~Tuning();

	// Associative search: return index number of a name
	int indexByLongName(const char * name) const;

	// Functions for setting individual tuning parameters
	bool setString(int index, const char * s);
	bool setList(int index, int len, int * newList);
	bool setWeight(int index, double weight);

	// Locked entries cannot be changed
	bool locked(int index) { return v[index].locked; }

	// Return values from a particular tuning entry
	const char * longName(int index) { return v[index].longName; }
	const char * shortName(int index) { return v[index].shortName; }
	const char * descr(int index) { return v[index].desc; }
	double weight(int index) { return v[index].weight.d; }
	const char * string(int index) { return v[index].weight.str; }
	const int * list(int index, int & listLen) {
		listLen = v[index].listLen;
		return v[index].weight.lst;
	}
	const TuneVal_t getEntry(int index) const { return v[index]; }

	// These report the type of an entry
	bool isList(int index) { return v[index].listLen >= 0; }
	bool isString(int index) { return v[index].listLen == -2; }
	bool isWeight(int index) { return v[index].listLen == -1; }

	// Reset to initial (default) values; strings and lists
	//		will be deleted
	void resetWeight(int index);
	void resetAllWeights() {
		for (int index = 0; index < numEntries; index++)
			resetWeight(index);
	}

	// Print entries to stdout
	void printAsRegistry();
	void printDesc();
	void printAsTDTable();	// A tab-delimited table; inserts easily
							//     into a table using Word's "text to
							//		table" command.
	void print(int index, bool full = false);	// For debugging only
	void printWeight(int i);		// for users to read

protected:

	TuneVal_t * v;	// Array of the parameter values
	int numEntries;
	int index;

	virtual bool initialize() = 0;

	// Functions for creating individual tuning parameters
	bool setNext(char * shortName, char * longName, double defaultWeight,
		char * desc, bool lock = false);
	bool setNextList(char * shortName, char * longName, char * desc,
		bool lock = false);
	bool setNextString(char * shortName, char * longName, char * desc,
		bool lock = false);
};

// An instance of one of the classes derived from this one is
//		instantiated in P3DControl at the global scope.



// Automatically include either bpTuning.h or gpTuning.h
#ifdef BINARY
#include "bpTuning.h"
#else
#include "gpTuning.h"
#endif



#endif

