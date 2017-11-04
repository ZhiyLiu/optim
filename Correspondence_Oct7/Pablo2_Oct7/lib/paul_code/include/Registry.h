/******************************************************************
 * REGISTRY Library                                               *
 ******************************************************************
 * Author:						Paul Yushkevich
 *
 * Date:							1998
 *
 * Description					Text I/O with a key/value/folder organization
 *
 *
 * Dependencies:				None
 ******************************************************************
 * Registry.h
 *	-----------
 *
 ******************************************************************/
#ifndef _Registry_H_
#define _Registry_H_

#ifdef _MSC_VER
# pragma warning(disable:4786)	/* '255' characters in the debug information */
#endif	/* _MSC_VER */

#include <iostream>
#include <stdio.h>
#include "libs.h"
#include <string>

#include <list>


// Begin namespace
NAMESPACE_PAULY_START

struct Hashentry {
	char * key;
	void * value;	// Either a string or a Registry
};


/*
    Hash table class associated with registry processing
*/
class Hashtable {
	private:

		int capacity;
		int increment;
		int size;
		int size_w_deletions;

		Hashentry ** data;

		unsigned int hash(const char * s) const;
		int findSpot(const char * key) const;
		void rehash();
		void getAllKeys(char ** returnArray, int & idx);	

	public:	

		Hashtable(int initialCapacity, int capacityIncrement);
		~Hashtable();

		void put(const char * key, const void * value);
		void * get(const char * key) const;
		Hashentry * entry(const char * key) const;
		void remove(const char * key);
		void removeAll();

		void printMainKeys();
		void printAllKeys();
        void dump(int indent = 2) const;     // For debugging only

		int countKeys() { return size; }
		void getKeys(char ** returnArray);

		int countAllKeys();
		void getAllKeys(char ** returnArray);
};


/*
    Exception associated with registry processing
*/
class RException {
		std::string m_error;

	public:
		RException(std::string error) : m_error(error) {}
		RException() {
            std::cout << "Except." << std::endl;
        }
		~RException() {}

		void print(std::ostream & s) const;
		std::string message() const { return m_error; }
};

inline std::ostream & operator << (std::ostream & s, const RException & exc)
{
   exc.print(s);
   return s;
}


/*
    Class Registry

	This class allows key value pairs to be stored in an ASCII file, called a
	registry file, which may be later read by this class.  The class instance
	may then be queried to obtain the values for specific keys.

	The keys should be printable ASCII strings.  The first letter of key
	names must be upper case, and if a lower case example is encountered in a
	registry file, it will be converted to upper case on input before other
	processing.  Registry files are intended to be hand-editable.

	There are three kinds of keys, formatted in registry files as follows:
		Regular value key:	keyName = value;
		Array key:			keyName = { ArrayLength ArrayType ArrayValues }
		Folder key:			keyName { ... }
	Usually the closing brace of folder keys is on a later line.  The
	ellipsis consists of numerous key-value pairs, including other folder
	keys.

	Registry files may also contain single-line comments indicated by the
	first non-whitespace character being the pound sign, '#'.  They may
	also contain statements for including other files, of the form:
		#include file-path
	See below and at the top of Registry.cpp for more information on include
	files and array keys.

	Functions for storing and retreiving keys from an instance of this
	class are provided.  Deletion of regular value keys is supported.
    Functions to read from or write to a registry file, as well as for
    performing a number of utility operations are also provided.

	A value of regular key may be boolean, integer, float, double or a
	string.  Almost any legal strings may be passed as values to the
	functions of this class.  In particular, strings should not contain
	numeric escape sequences (for example \012 or \0A to represent \n).
	Any non-printable ASCII characters, such as NEWLINE, will be encoded
	as 2-digit escape sequences in output registry files.  On input, the
	original ASCII string will then be recovered.  Thus 3-digit escape
	sequences will cause errors in the parsing of registry files and
	2-digit escape sequences will not be treated as one might expect.
*/

class ofstream;

#ifdef SWIG
%ignore	Registry(const Registry &);
#endif

class Registry {
	friend class Hashtable;

private:

	Hashtable table;

	char * includePathList;
	bool writable;
	const char ** orderingList;
	int arrayWidth;

	// Hashtable table;
	static char * tmpKey;

	bool readFile(const char * pathname, std::string & oss,
		bool caseIndep);
	void setValue(const char * key, char * value);
	char * getValue(const char * key);
	bool write(std::ostream * f, int indent);
	void writeEncodedString(std::ostream * f, char * s);
	int hexDigit(char c);
	void decodeString(char * s);
	char * read(char * file, std::string & oss, bool caseIndep);
	Registry * getFolder(const char * key);
    bool hasKeyPrivate(char * key);
	void adjustOrdering(char ** keys, int numKeys);
	// Returns a list of keys in the registry, optionally starting with the
	// specified prefix.  If simple is true the format is "key value" or
	// "key {array_value"; if false, it is "key=value" or "key={array_value}".
	// If ignore is true, the ordering list will be ignored.
	void collectKeys(char ** lst, int & lstlen, std::string prefix,
		bool simple, bool ignore);

	Registry(const Registry &);	// Not implemented

public:

	// For construction, length and width are merely hints.  Length is an
	// estimate of the total number of keys; width is an estimate of the
	// number of nested levels to be used in the names of the keys.  These
	// are used to optimize allocation of the hash table.
	Registry(int length = 100, int width = 5);
	Registry(const char * fname) throw(RException);		// Uses default length & width
    virtual ~Registry();

	void print() { table.printAllKeys(); }
	int size() { return table.countAllKeys(); }
	int mainSize() { return table.countKeys(); }
	void printMajorKeys() { table.printMainKeys(); }

	const char * getStringValue(const char * key, const char * defaultValue, ...);
	bool getBooleanValue(const char * key, bool defaultValue, ...);
	int getIntValue(const char * key, int defaultValue, ...);
	float getFloatValue(const char * key, float defaultValue, ...);
	double getDoubleValue(const char * key, double defaultValue, ...);

	// Arrays may only be read once and become property of the caller.
	// The length argument returns the number of items in the array.
	int * getIntArray(const char * key, int * length, ...);
	float * getFloatArray(const char * key, int * length, ...);
	double * getDoubleArray(const char * key, int * length, ...);

	bool cmpStringValue(char * key,char * cmpValue, ...);

	void setStringValue(const char * key, const char * value, ...);
	void setBooleanValue(const char * key, bool value, ...);
	void setIntValue(const char * key, int value, ...);
	void setFloatValue(const char * key, float value, ...);
	void setDoubleValue(const char * key, double value, ...);

	//	The array arguments in these functions should be created
	//	using new[] and will become property of the registry.
	void setIntArray(const char * key, int length, int * array, ...);
	void setFloatArray(const char * key, int length, float * array, ...);
	void setDoubleArray(const char * key, int length, double * array, ...);

    // Removes a key from the registry.  The key cannot be an array
    // and should always be fully qualified (i.e. not a folder).
	void deleteKey(const char * key);

	/*  The next function may be used to specify the number of elements
		written per line in arrays.  A default of 10 is used until
		changed.  Warning: the value specified will apply to all
		arrays written to an output file by this class.  In general,
		this should never be changed.
	*/
	void setArrayWidth(int width) { arrayWidth = width; }

	Registry & getSubFolder(const char * key, ...);
	bool hasKey(const char * key, ...);

   /* These methods allow us to list the keys in the registry.  First
	  two methods both take a string array that must be allocated to
	  hold getKeyArraySize() strings.  (The caller allocates the array).

	  The first method places all value-bearing keys in the folder into
	  the array and returns the number of keys.  The array is
	  NULL-terminated and getKeyArraySize() - 1 keys will be returned.
	  The second method returns the array and number of subfolder names
	  in that folder.  Keys returned are pointers into the registry data,
	  so they should not be modified or deleted.
    */
	int getValueKeys(char * keys[]);
	int getFolderKeys(char * keys[]);
	int getKeyArraySize();

	/* Constrain the sorting order of highest-level key names.  This
	   function affects the action of collectAllKeys() and writeToFile().
	   The order array should contain a NULL terminated set of zero-ended
	   strings specifying the major ordering of keys.  Specifically, the
	   list should contain the first group of keys or folders, in the order
	   that they should appear in an output file or list of key-value pairs,
	   produced by writeToFile() and collectAllKeys(), respectively.  Any
	   other major names will appear alphabetically (case-insensitive) after
	   the names provided.  The major key names must begin with a letter for
	   this ordering to work.  If no order list is provided, the keys are
	   returned in case-insensitive alphabetic order.  If no ordering is
	   specified, any existing ordering is cleared.
	*/
	void ordering(const char ** order = NULL);

	/* Collect all keys into a list.
       This function  returns an array of key-value strings.  The key is
       the fully qualified name of a key, and value is the key's value.
       The output will be in one of two possible formats, depending on
       the value of the second argument:
			simple				    True                False
			Single-valued keys		key value			key=value
			Array-valued keys		key {array_value	key={array_value}
       In all cases, terminal semicolons will be omitted.

       For keys referring to arrays that have already been read, the
       value will be of the form "{length type ...}".  If the array has not
       been read, the ellipsis will be replaced by the list of data values.
       Length returns the number of items in the returned array and type is
       an indication of the data type: 0 for integer, 1 for float and 2 for
       double.  Function ordering(), above, may be used to influence the
       order of primary key names in the list returned.  The array returned
	   should be deleted after use.
	*/
	char ** collectAllKeys(int & length, bool simple = false,
		bool ignore = false);

	/* Specify a list of paths to search for registry files
	   required by #include directives in the file being loaded.
	   The path list must have the form "dir1|dir2|dir3".
	   This path list may be superceded by an INCLUDEPATH key
	   in the registry file being loaded.  If no path is specified
	   by either method, the arguments of any #include directives
	   must be absolute paths.
	*/
	void setIncludePath(const char * pathlist);

	// Registry's containing #include cannot be written.
	bool isWritable() const { return writable; }

	/* Produce or replace a registry file.  The file produced will not
	   contain any #include directives, although they may be added
	   manually.  Registry objects loaded from files containing a
	   #include will not be written (see above).  Function ordering()
	   may be used to influence the order of primary key names in
	   the output file.  If non-fatal errors occur, they will be
	   printed to cout and false will be returned.  There are only
	   two of these, both resulting from corrupted or empty arrays.
	*/
	bool writeToFile(const char * pathname) throw(RException);

	// This is like writeToFile(), except that the output goes to cout.
	// Error messages that writeToFile() produces will not be printed,
	// but the returned value will be false.
	bool printAsRegistry();

	/* Add the contents of the specified file to the current Registry
	   object.  This may be used multiple times to merge the contents
	   of several files into the same Registry.  If the second argument
	   is true, the first letter of all key names in the loaded file
	   will be converted to upper case before they are stored in the
	   Registry object.
	*/
	void readFromFile(const char * pathname, bool caseIndep = false)
		throw(RException);

	// Discard the contents of the current Registry object.
	void clear();

	/* Utility function to convert a text representation of a registry
	   array to an actual array.  The array string should not have a
	   leading brace, '{'.  A pointer to the array will be returned,
	   and it length and type will be set.  Used by program flat.
	*/
	char * parseArray(char * arrayString, int & length, int & type);

    // For debugging: prints the hash table tree
    void dump() const { table.dump(); }
};


// Begin namespace
NAMESPACE_PAULY_END

#endif

