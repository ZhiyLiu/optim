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
 * registry.cpp
 *	-----------
 *
 ******************************************************************/
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Registry.h"


//#define DEBUG


using namespace std;


#ifndef S_ISREG
#define S_ISREG(m)	(((S_IFREG & m) != 0) ? true : false)
#endif

const int MAX_KEY_LEN = 64;		// Maximum length of a hash key

// Begin namespace
NAMESPACE_PAULY_START


void RException::print(std::ostream & s) const
{
	s << m_error << std::endl;
}


/*
	This hash table supports the storing of simple data values 
	specified using the syntax "key = value;".  The key and value
	are character strings.  White-space characters between these
	elements of this format will be surpressed.  The key cannot
	contain white space.  Neither the key nor value can contain
	semicolons.  Non-white-space characters following the
	semicolon will be considered to belong to the next key, so
	it is not necessary to put each key-value pair on a separate
	line.

	The hash table also supports the storing of arrays.  Special
	encoding of the hash key and value are needed in the later
	case.  Specifically, the final name of the key will begin
	with '@'.  The value stored with that key is a string
	consisting of three parts:
		length			The number of items in the array.
		type			The integer casting of an ArrayType_t value
							indicating the type of the array data.
		array address	The hex address of the actual array in
							memory.  The leading "0x" is omitted.
							Once an array is retrieved by the user
							of this class, the array address in
							the hash table record is set to 0.

	Deletion of regular keys is supported.  Deletion of array keys
	has not been implemented.  Deletion can only be done for fullly
    qualified keys.  Thus, if a.b.c is a key, it is not possible to
    delete a.b or a.
    
    Note that deletion has to be implemented carefully because of
    the collision resolution mechanism of the hash table.
    Initially, deleted keys are just turned into keys with empty
    strings in the key field.  Only for leaf keys can the value be
    deleted.  Whenever a rehash is done to enlarge the table, then
    they are actually removed, but only for the level rehashed.

	A simple mechanism is also provided for causing one registry
	file to include another registry file.  This may be done by
	inserting a line of the form "#include filename" in the first
	registry file, where "filename" is the name (or path) of the
	registry file to be included.  There cannot be any whitespace
	between the "#" and "include"; at least one space or tab must
	immediately follow "include".  The special key INCLUDEPATH
	may be used to specify a set of directories to be searched for
	the file to be included.  This key must preceed the "#include"
	directive and should have the syntax:
		"INCLUDEPATH = dir1|dir2|...|dirN",
	where the list of directories are absolute paths.  The
	setIncludePath() member function also may be used to specify
	the same list of directories.  If no directory list is
	provided, the filename must be an absolute pathname.

	Registry files may also contain comments, which take the form
	of "#" followed by arbitrary text to the end of the line,
	except that comments cannot begin like "#include" statements.
	Note that semicolons cannot be used to terminate comments;
	any text between a semicolon and the newline will be treated
	as part of the comment. If the registry file is read and
	then written back out, any comments will be lost.
*/

unsigned int Hashtable::hash(const char * s) const
{
	// Some good hashing function
	int hash = 0;
	while (*s) {
		hash ^= (*s);
		s++;
	}

	return hash;
}

int Hashtable::findSpot(const char *key) const
{
	int position = hash(key) % capacity;
	int i = 1;

	while (data[position] && strcmp(key,data[position]->key) != 0) {
		position = (position + i*i) % capacity;
		i++;
	}

	return position;
}

void Hashtable::rehash()
{
	Hashentry ** oldData = data;
	int oldCap = capacity;
    int i;

	// Reallocate the data
	capacity += increment;
	data = new Hashentry *[capacity];
	memset(data, 0, capacity*sizeof(Hashentry *));

	// Put each entry in its place in the new data array
	for (i = 0; i < oldCap; i++) {
		if (oldData[i]) {
			if (oldData[i]->key[0] != '\0') // Move valid entries
				data[findSpot(oldData[i]->key)] = oldData[i];
			else    // Do not move the marked entries
				size_w_deletions--;
		}
	}

	delete [] oldData;
}

Hashtable::Hashtable(int initialCapacity, int capacityIncrement)
{
	capacity = initialCapacity;
	increment = capacityIncrement;

	data = new Hashentry *[capacity];
	memset(data, 0, capacity*sizeof(Hashentry *));
	size = 0;
	size_w_deletions = 0;
}

Hashtable::~Hashtable() {
	removeAll();
	// dibyendu
	if( data != NULL )
		delete [] data;	
	data = NULL ;
}

void Hashtable::put(const char * key, const void * value)
{
	int index = findSpot(key);
	Hashentry * entry = data[index];

	if (entry == NULL) {
		entry = new Hashentry;
		data[index] = entry;
		entry->key = new char[strlen(key) + 1];
		strcpy(entry->key, key);

		size++;
		size_w_deletions++;
        if (size_w_deletions > capacity/2)
			rehash();
	}
	entry->value = (void *) value;
}

void * Hashtable::get(const char * key) const
{
	int index = findSpot(key);
	Hashentry * entry = data[index];
	return (! entry) ? entry : entry->value;
}

Hashentry * Hashtable::entry(const char * key) const
{
	int index = findSpot(key);
	Hashentry * entry = data[index];
	return entry;
}

void Hashtable::remove(const char * key)
{
	int i = findSpot(key);
	if (! data[i])
        return;    // Ignore missing keys
	if ((data[i]->key)[0] != '\0') {
		if ((data[i]->key)[0] == '=') {
			delete [] (char *) data[i]->value;
		    data[i]->value = NULL;
		}
		else if ((data[i]->key)[0] == '{' && (data[i]->key)[1] == '}') {
            // Do nothing to the value pointer, which is an address of the
            // folder at the next level.  It will be removed when rehash()
            // is next called.
		}
        else {
			cout << "Error deleting " << key << endl;   // Handles arrays
            return;
        }
		delete [] data[i]->key;
		data[i]->key = new char[1]; // Mark this as a removed key
		data[i]->key[0] = '\0';
		size--;
	}
	else if ((data[i]->key)[0] == '{' && (data[i]->key)[1] == '}')
		cout << "Error deleting " << key << endl;    // Handles folder keys
}

enum ArrayType_t { IntArray, FloatArray, DoubleArray };

void Hashtable::removeAll()
{
#ifdef DEBUG
	// This can be used to tune the arguments of the Registry constructor
	cout << "capacity = " << capacity << " \tsize = " << countKeys()
		<< " \tratio = " << capacity/float(1 + size) << " \tincrement = "
		<< increment << endl;
#endif
	for (int i = 0; i < capacity; i++) {
		if (data[i]) {
			if ((data[i]->key)[0] == '=')
				delete [] (char *) data[i]->value;
			else if ((data[i]->key)[0] == '@') {
				int length, type;
				void * array;

				// Delete array, if it is unread
				char * arrayStr = (char *) data[i]->value;
				sscanf(arrayStr, "%d %d %x", &length, &type, &array);
				if (array != NULL) {
                	switch(type) {
                      case DoubleArray: { delete [] ((double *) array); break; }
                      case FloatArray: { delete [] ((float *) array); break; }
                      case IntArray: { delete [] ((int *) array); break;}
                      default: { cout << "Deleting array of unknown type" << endl; }
                    }
 				}
				delete [] (char *) data[i]->value;
			}
			else if ((data[i]->key)[0] == '{' && (data[i]->key)[1] == '}') {
				Registry * reg = (Registry *) data[i]->value;
				delete reg;
			}
			else if ((data[i]->key)[0] != '\0')
				cout << "Error deleting" << endl;
			delete [] data[i]->key;
			
			delete data[i];
			data[i] = NULL;
		}
	}
	size = 0;
	size_w_deletions = 0;
}

void Hashtable::getKeys(char ** returnArray)
{
	int idx = 0;
	for (int i = 0; i < capacity; i++) {
		if (data[i] && data[i]->key[0] != '\0')
			returnArray[idx++] = data[i]->key;
	}
}

// This function must initially be called with idx = 0.  The size
// of returnArray must be sufficiently large to hold the keys.
void Hashtable::getAllKeys(char ** returnArray, int & idx)
{
	for (int i = 0; i < capacity; i++) {
		if (data[i]) {
			if (data[i]->key[0] != '\0') {
				returnArray[idx] = data[i]->key;
				idx++;
				if ((data[i]->key)[0] == '{' && (data[i]->key)[1] == '}') {
					Registry * reg = (Registry *) data[i]->value;
					reg->table.getAllKeys(returnArray, idx);
				}
			}
		}
	}
}

void Hashtable::getAllKeys(char ** returnArray)
{
	int idx = 0;
	getAllKeys(returnArray, idx);
}

int Hashtable::countAllKeys()
{
	int sz = 0;
	for (int i = 0; i < capacity; i++) {
		if (data[i] && (data[i]->key)[0] != '\0') {
			sz++;
			if ((data[i]->key)[0] == '{' && (data[i]->key)[1] == '}') {
				Registry * reg = (Registry *) data[i]->value;
				sz += reg->table.countAllKeys();
			}
		}
	}
	return sz;
}

void Hashtable::printMainKeys()
{
	char ** keys = new char *[countKeys()];
	getKeys(keys);
	for (int i = 0; i < countKeys(); i++)
		cout << keys[i] << '\n';
	delete [] keys;
}

void Hashtable::printAllKeys()
{
	int sz = countAllKeys();
	char ** keys = new char *[sz];
	getAllKeys(keys);
	for (int i = 0; i < sz; i++)
		cout << keys[i] << '\n';
	delete [] keys;
}

void Hashtable::dump(int indent) const
{
	char * ind = new char[indent + 1];
	memset(ind, ' ', indent);
	ind[indent] = '\0';

	for (int i = 0; i < capacity; i++) {
		if (data[i]) {
			if (data[i]->key[0] == '\0')
                cout << ind << "Empty node\n";
			else if ((data[i]->key)[0] == '{' && (data[i]->key)[1] == '}') {
                cout << ind << data[i]->key << '\n';
				Registry * reg = (Registry *) data[i]->value;
                reg->table.dump(indent + 2);
            }
            else if ((data[i]->key)[0] == '=') {
                cout << ind << data[i]->key << " = ";
                cout << (char *) data[i]->value << '\n';
            }
            else if ((data[i]->key)[0] == '@')
                cout << ind << data[i]->key << " = array\n";
            else
                cout << "Unknown entry detected\n";   // Should never print
		}
	}
}


// -----------------------------  Registry Class Follows  -----------------------------


int keySortFunc(const void * p1, const void * p2)
{
	char *s1 = *((char **) p1);
	char *s2 = *((char **) p2);

	// Compare the two strings up to the brackets - first folders and keys are compared
	// though
	while(1) {
		if(*s1 < *s2)
			return -1;
		if(*s2 < *s1)
			return 1;
		if(*s1 == 0)
			return 0;

		// Special case - compare past brackets - this routine will even sort things like
		// array[12][14][1]
		if(*s1 == '[') {
			char * b1 = strchr(s1, ']');
			char * b2 = strchr(s2, ']');
			if(b1 == NULL || b2 == NULL)
				return strcmp(s1, s2);
			*b1 = *b2 = 0;
			int a1 = atoi(s1 + 1);
			int a2 = atoi(s2 + 1);
			*b1 = *b2 = ']';

			if(a1 < a2)
				return -1;
			if(a1 > a2)
				return 1;

			s1 = b1;
			s2 = b2;
		}

		s1++; s2++;
	}
}

void Registry::setValue(const char * key, char * value)
{
	// See if there is a dot
	char * dot = (char *) strchr(key, '.');   // Solaris CC requires this cast
	if (dot) {
		// See if there is a folder for pre-dot stuff.
		*dot = 0;
		Registry * f = getFolder(key);

		// Create a folder if one is not present
		if (! f) {
			f = new Registry;
			char * fkey = new char[strlen(key) + 3];
			strcpy(fkey, "{}");
			strcat(fkey, key);
			table.put(fkey, f);
			delete[] fkey;
		}

		*dot = '.';
		f->setStringValue(dot + 1, value);
		delete [] value;
	}
	else {
		char *newKey = new char[strlen(key) + 2];
		if (key[0] == '@')
			newKey[0] = '\0';
		else
			strcpy(newKey, "=");
		strcat(newKey, key);
		Hashentry * entry = (Hashentry *) table.entry(newKey);
		if (entry) {
			// Replacing an old entry
			char * oldvalue = (char *) entry->value;
			delete [] oldvalue;
			entry->value = value;
		}
		else	// Adding a new entry
			table.put(newKey, value);
		delete [] newKey;
	}
}

char * Registry::getValue(const char * key)
{
	// If folder has dot, delegate to next folder.
	char * dot = (char *) strchr(key, '.');   // Solaris CC requires this cast
	
	if (dot) {
		*dot = 0;
		Registry * f = getFolder(key);
		*dot = '.';
		return f ? f->getValue(dot + 1) : NULL;
	}
	else {
		char * trueKey = new char[strlen(key) + 2];
		if (key[0] != '@')
			strcpy(trueKey, "=");
		else
			trueKey[0] = '\0';
		strcat(trueKey, key);
		char * s = (char *) table.get(trueKey);
		delete [] trueKey;
		return s;
	}
}

void Registry::deleteKey(const char * key)
{
	// If folder has dot, delegate to next folder.
	char * dot = (char *) strchr(key, '.');   // Solaris CC requires this cast
	if (dot) {
		*dot = 0;
		Registry * f = getFolder(key);
		if (f) {
			f->deleteKey(dot + 1);
			if (f->mainSize() == 0) {
		        char * trueKey = new char[strlen(key) + 3];
			    strcpy(trueKey, "{}");
			    strcat(trueKey, key);
			    table.remove(trueKey);
		        delete [] trueKey;
			}
		}
		*dot = '.';
	}
	else {
		char * trueKey = new char[strlen(key) + 2];
		if (key[0] != '@')
			strcpy(trueKey, "=");
		else
			trueKey[0] = '\0';
		strcat(trueKey, key);
		table.remove(trueKey);
		delete [] trueKey;
	}
}

int Registry::getValueKeys(char * keys[])
{
   int keyArrayIdx = 0;

   table.getKeys(keys);
   qsort(keys, table.countKeys(), sizeof(char *), keySortFunc);

	for (int i = 0; i < table.countKeys(); i++) {
		if (keys[i][0] == '=' || keys[i][0] == '@') {
			keys[keyArrayIdx] = keys[i] + 1;
         keyArrayIdx++;
      }
   }

   keys[keyArrayIdx] = NULL;
   return keyArrayIdx;
}

int Registry::getFolderKeys(char * keys[])
{
   int keyArrayIdx = 0;

   table.getKeys(keys);
   qsort(keys, table.countKeys(), sizeof(char *), keySortFunc);

	for (int i = 0; i < table.countKeys(); i++) {
		if (keys[i][0] == '{') {
			keys[keyArrayIdx] = keys[i] + 2;
         keyArrayIdx++;
      }
   }

   keys[keyArrayIdx] = NULL;

   return keyArrayIdx;
}

int Registry::getKeyArraySize()
{
   int nKeys = table.countKeys();
   return nKeys + 1;
}

// Write this folder to output stream, at given indent, sorted by key
bool Registry::write(std::ostream * f, int indent)
{
	//cout << "Registry::write started !" << endl ;

	int numKeys;
	bool ret;

	char * ind = new char[indent + 1];
	memset(ind, ' ', indent);
	ind[indent] = '\0';

	numKeys = table.countKeys();
	char ** keys = new char *[numKeys];
	table.getKeys(keys);

	// Sort keys using a simple sorting algorithm
	if (numKeys > 1)
	    qsort(keys, numKeys, sizeof(char *), keySortFunc);

	// Apply special ordering to major names.
	if (orderingList != NULL)
		adjustOrdering(keys, numKeys);

	// For recursive calls, ignore the ordering list
	const char ** ordering = orderingList;
	orderingList = NULL;

	ret = true;

	// Write each key out
	for (int i = 0; i < numKeys; i++) {
		char * key = keys[i];

		// Determine the kind of key to be written
		if (key[0] == '{') {	// Key is a folder
			//Move up to actual key
			key += 2;

			// Write folder name, bracket
			*f << ind << key << " {\n";

			// Write out folder contents
			Registry * fld = getFolder(key);
			fld->setArrayWidth(arrayWidth);
			fld->write(f, indent + 3);

			// Write the closing bracket
			*f << ind << "}\n";
		}
		else if (key[0] == '@') {	// Key is an array
			int length, i;
			char * array;
			int * ia;
			float * fa;
			double * da;
			int type;

			char * arrayStr = getValue(key);
			// Must ignore the first character
			if (3 != sscanf(arrayStr, "%d %d %x", &length, &type, &array)) {
				if (f != &cout)
					cout << "Warning: could not write key " << key + 1 << '\n';
				ret = false;
				continue;
			}
			if (array == NULL) {
				if (f != &cout)
					cout << "Warning: empty array for key " << key + 1 << '\n';
				ret = false;
				continue;
			}
			key++;
			*f << ind << key << " = { " << length << ' ' << type << '\n';

			switch ((ArrayType_t) type) {
				case IntArray:
					ia = (int *) array;
					*f << ' ' << ia[0];
					for (i = 1; i < length; i++) {
						if (i%arrayWidth == 0)
							*f << '\n';
						*f << ' ' << ia[i];
					}
					break;
				case FloatArray:
					fa = (float *) array;
					*f << ' ' << fa[0];
					for (i = 1; i < length; i++) {
						if (i%arrayWidth == 0)
							*f << '\n';
						*f << ' ' << fa[i];
					}
					break;
				case DoubleArray:
					da = (double *) array;
					*f << ' ' << da[0];
					for (i = 1; i < length; i++) {
						if (i%arrayWidth == 0)
							*f << '\n';
						*f << ' ' << da[i];
					}
					break;
			}
			*f << " };\n";
		}
		else if (key[0] == '=') {	// Key has a single value
			// Write a key - value pair
			key++;

			// Write encoded value
			*f << ind << key << " = ";
			writeEncodedString(f, getValue(key));
			*f << ";\n";
		}
	}
	f->flush();

	orderingList = ordering;	// Restore ordering list

	delete [] keys;
	delete [] ind;
	return ret;
}

/*	Apply special ordering to major names.
	This makes the first part of a key list identical to
	the ordering list, while retaining the sorted ordering
	of any remaining keys.  For example, if the sorted
	keys are a.1, b.2, c.3, and d.4, and the ordering list
	contains "b d", then the adjusted list of keys will
	become "b.2 d.4 a.1 c.3".
*/
void Registry::adjustOrdering(char ** keys, int numKeys)
{
	char * tmp;
	char * str;
	int pos, size;

	pos = 0;
	size = 0;
	while (orderingList[size] != NULL)
		size++;

	for (int o = 0; o < size; o++) {
		for (int k = 0; k < numKeys; k++) {
			str = keys[k];

			while (! isalpha(*str)) str++;
			if (0 == strcmp(str, orderingList[o])) {
				if (k != pos) {
					tmp = keys[k];
					for (int j = k; j > pos; j--)
						keys[j] = keys[j - 1];
					keys[pos] = tmp;
				}
				pos++;
				break;
			}
		}
	}
}

void Registry::collectKeys(char ** lst, int & lstlen, std::string prefix, bool simple, bool ignore)
{
	int numKeys;
	string folder;

	numKeys = table.countKeys();
	if (numKeys == 0) {
		cout << "Collected no keys!\n";	// Should not happen
		return;
	}
	char ** keys = new char*[numKeys];
	table.getKeys(keys);
	folder = prefix;

	// Sort keys using a simple sorting algorithm
	qsort(keys, numKeys, sizeof(char *), keySortFunc);

	// Apply special ordering to major names.
	if (! ignore && orderingList != NULL) 
		adjustOrdering(keys, numKeys);

	// Expand each key
	for (int i = 0; i < numKeys; i++) {
		char *key = keys[i];

		// Expand the folder
		if (key[0] == '{') {
			//Move up to actual key
			key += 2;

			// Include folder name in the prefix
			prefix += key;
			prefix += '.';

			// Write out folder contents
			Registry *fld = getFolder(key);
			// For recursive calls, ignore the ordering list
			fld->collectKeys(lst, lstlen, prefix, simple, true);
			prefix = folder;
		}
		else if (key[0] == '@') {
			int length, i;
			char * array;
			int * ia;
			float * fa;
			double * da;
			int type;
			char ltstr[32];
			int ltstrlen;
			string str;
			ostringstream ss;

			char * arrayStr = getValue(key);
			// Must ignore the first character
			if (3 != sscanf(arrayStr, "%d %d %x", &length, &type, &array)) {
				cout << "Warning: could not process key " << key + 1 << '\n';
				continue;
			}
			key++;

			// Copy the length and type
			ltstrlen = 0;
			while (arrayStr[ltstrlen] != ' ')
				ltstr[ltstrlen] = arrayStr[ltstrlen++];
			ltstr[ltstrlen] = arrayStr[ltstrlen++];
			while (arrayStr[ltstrlen] != ' ')
				ltstr[ltstrlen] = arrayStr[ltstrlen++];
			ltstr[ltstrlen] = '\0';

			// Typically, the array will already have been read, which means
			// it no longer belongs to this class.
			if (array != NULL) {
				switch ((ArrayType_t) type) {
					case IntArray:
						ia = (int *) array;
						for (i = 0; i < length; i++)
							ss << ' ' << ia[i];
						break;
					case FloatArray:
						fa = (float *) array;
						for (i = 0; i < length; i++)
							ss << ' ' << fa[i];
						break;
					case DoubleArray:
						da = (double *) array;
						for (i = 0; i < length; i++)
							ss << ' ' << da[i];
						break;
				}
				lst[lstlen] = new char[prefix.size() + strlen(key) + ltstrlen + ss.str().size() + 5];
			}
			else
				lst[lstlen] = new char[prefix.size() + strlen(key) + ltstrlen + 9];

			strcpy(lst[lstlen], prefix.c_str());
			strcat(lst[lstlen], key);
			if (simple)
				strcat(lst[lstlen], " {");
			else
				strcat(lst[lstlen], "={");
			strcat(lst[lstlen], ltstr);
			if (array != NULL)
				strcat(lst[lstlen], ss.str().c_str());
			else
				strcat(lst[lstlen], " ...");
			if (! simple)
				strcat(lst[lstlen], "}");
			lstlen++;
		}
		else if (key[0] == '=') {
			// Store the full key and encoded value
			key++;
			char * val = getValue(key);

			lst[lstlen] = new char[prefix.size() + strlen(key) + strlen(val) + 3];
			strcpy(lst[lstlen], prefix.c_str());
			strcat(lst[lstlen], key);
			if (simple)
				strcat(lst[lstlen], " ");
			else
				strcat(lst[lstlen], "=");
			strcat(lst[lstlen], val);
			lstlen++;
		}
	}

	delete [] keys;
}

char ** Registry::collectAllKeys(int & length, bool simple, bool ignore)
{
	string prefix("");
	length = table.countAllKeys();	// Maximum possible length
	char ** lst = new char *[length];
	length = 0;
	collectKeys(lst, length, prefix, simple, ignore);
	return lst;
}

/*	String encoder.

	This does not encode spaces, so can only be used on values.
	Characters are assumed to be in the range [0x00 - 0x7F].
	DEL encodes to "\7F", and characters less than SPACE encode
	to "\XX", where XX is their hexadecimal value.  Semicolon
	also must be encoded as an escape sequence, because of it's
	use to terminate key value pairs in registry files.

	The backspace itself is passed through unchanged, so that
	values of Registry entries may contain character escapes,
	such as \n, \t, etc.  However, an escape sequence of the
	numeric kind, such as \7F should never be put in a value,
	because this function inserts them in the output stream
	and decodeString() below will remove them.
*/
void Registry::writeEncodedString(std::ostream * f, char * s)
{
	for(int i = 0; i < strlen(s); i++) {
		char c = s[i];
		// AGG: Why did Paul encode tilda?
		if (c >= ' ' /* && c <= '~' */ && c != ';' /* && c != '\\'*/) {
			*f << c;
		}
		else {
			// Encode semicolons and non-printing ASCII characters
			*f << '\\';
			f->width(2);
			f->fill('0');
			*f << hex << (int) c << dec;
		}
	}
}	

int Registry::hexDigit(char c)
{
	if(c >= 'a')
		return 10+c-'a';
	if(c >= 'A')
		return 10+c-'A';
	return c-'0';
}

void Registry::setIncludePath(const char * pathlist)
{
	if (includePathList != NULL)
		delete [] includePathList;
	includePathList = new char[1 + strlen(pathlist)];
	strcpy(includePathList, pathlist);
}

void Registry::decodeString(char *s)
{
	char *dst = s;

	while (*s) {
		if (*s == '\\') {
			if (isdigit(*(s + 1)) && isdigit(*(s + 2))) {
				// Encoded escape sequence of the form \XX found (see above)
				*dst = hexDigit(*(s + 1))*16 + hexDigit(*(s + 2));
				s += 3;
			}
			else {
				// The backslash was escaping the next character, so was
				// not inserted by function writeEncodedString().  This
				// is mainly used to escape spaces in Windows pathnames.
				*dst = *s;
				s++;
			}
		}
		else {
			*dst = *s;
			s++;
		}
		dst++;
	}

	*dst = 0;
}

void Registry::clear()
{
	table.removeAll();
	if (includePathList != NULL) {
		delete [] includePathList;
		includePathList = NULL;
	}
}

char * Registry::parseArray(char * arrayString, int & length, int & type)
{
	// Convert the array string
	int i;
	char * str;
	bool fail;
	int * ia;
	float * fa;
	double * da;
	char * array;

	if (2 != sscanf(arrayString, "%d %d", &length, &type))
		return NULL;

	str = strtok(arrayString, " \n\t\r");	// length
	str = strtok(NULL, " \n\t\r");			// type
	fail = false;
	switch ((ArrayType_t) type) {
		case IntArray:
			ia = new int[length];
			array = (char *) ia;
			for (i = 0; i < length; i++) {
				str = strtok(NULL, " \n\t\r");
				if (1 != sscanf(str, "%d", ia + i)) {
					fail = true;
					break;
				}
			}
			break;

		case FloatArray:
			fa = new float[length];
			array = (char *) fa;
			for (i = 0; i < length; i++) {
				str = strtok(NULL, " \n\t\r");
				if (1 != sscanf(str, "%f", fa + i)) {
					fail = true;
					break;
				}
			}
			break;

		case DoubleArray:
			da = new double[length];
			array = (char *) da;
			for (i = 0; i < length; i++) {
				str = strtok(NULL, " \n\t\r");
				if (1 != sscanf(str, "%lf", da + i)) {
					fail = true;
					break;
				}
			}
			break;
	}
	if (fail) {
		delete [] array;
		return NULL;
	}
	else
		return array;
}

// Read this folder in
char * Registry::read(char * file, string & oss, bool caseIndep)
{
	// cout << "Registry::read started" << endl ;

	char * keyStart;
	char * keyEnd;
	char * key;

	while(1) {
		// We are expecting a key - find its beginnning and end
		keyStart = file + strspn(file, " \t\n\r");

		// If key starts with closing brace, we are done
		if (*keyStart == '}' || *keyStart == '\0')
			return keyStart;

		if (*keyStart == '#')
			keyEnd = keyStart + 1;
		else
			keyEnd = keyStart + strcspn(keyStart, " \t\n\r={");

		// Create a key string
		key = new char[keyEnd - keyStart + 1];
		strncpy(key, keyStart, keyEnd - keyStart);
		key[keyEnd - keyStart] = '\0';

		// Check the next character after the key
		char * nextChar = keyEnd + strspn(keyEnd, " \t\n\r");

		// See if this is a name-value assignment
		if (*nextChar == '=') {
			// We got us a key-value pair
			char * valueStart = nextChar + 1 + strspn(nextChar + 1, " \t\n\r");
			char * valueEnd = strchr(valueStart, ';');
			/*  Bug: if a semicolon is missing on the current value, the first
				one from a later key-value pair will be used.  Therefore, while
				the current key will be processed correctly, at least the next
				one will simply be skipped and the following error will not be
				reported.  The setValue() function cannot be changed to detect
				the problem.  The only solution is apparently to enforce their
				being one key-value pair per line of the file.
			*/
			if (valueEnd == NULL) {
				oss = "Key ";
				oss += key;
				oss += " not followed by semicolon terminated value.\n";
				delete [] key;
				return NULL;
			}
			*valueEnd = '\0';

			if (*valueStart != '{') {
				// Decode the value string
				decodeString(valueStart);

				// Set value
				char * value = new char[strlen(valueStart) + 1];
				strcpy(value, valueStart);

				// Set the value
				setValue(key, value);
			}
			else {
				// Convert the array string
				int length, type;
				char * array = parseArray(valueStart + 1, length, type);
				if (array == NULL) {
					oss = "Key ";
					oss += key;
					oss += " is an invalid array.\n";
					delete [] key;
					return NULL;
				}
				valueStart = strtok(NULL, " \n\t\r");
				if (! valueStart || *valueStart != '}') {
					delete [] array;
					oss = "Key ";
					oss += key;
					oss += " not terminated by closing brace.\n";
					delete [] key;
					return NULL;
				}

				// Set the value
				char * value = new char[32];
				sprintf(value, "%d %d %x", length, type, array);

				// Set the value
				char * fullKey = new char[strlen(key) + 2];
				strcpy(fullKey, "@");
				strcat(fullKey, key);
				setValue(fullKey, value);
				delete [] fullKey;
			}

			delete [] key;

			// Update the string pointer
			file = valueEnd + 1;
		}
		else if (*nextChar == '{') {
			// Create a folder
			Registry *sub = new Registry;
			char *folderEnd = sub->read(nextChar + 1, oss, caseIndep);

			// Folder end should be a closing brace
			if (folderEnd == NULL) {
				delete [] key;
				return NULL;
			}

			if (*folderEnd != '}') {
				oss = "Sub-folder ";
				oss += key;
				oss += " is terminated by '";
				oss += *folderEnd;
				oss += "', not a '}'\n";
				delete [] key;
				return NULL;
			}

			// Need to prepend key with "{}"
			char * newKey = new char[strlen(key) + 3];
			strcpy(newKey, "{}");
			strcat(newKey, key);
			delete [] key;

			// Add this folder
			table.put(newKey, sub);
			delete [] newKey;

			file = folderEnd + 1;
		}
		else if (0 == strcmp(key, "#include")) {
			// Note: trailing whitespace has been stripped from key
			char * path;
			char * value;
			char * directory;
			int filelen;

			delete [] key;
			writable = false;	// Prevent saving changes to this file

			char * end = strpbrk(nextChar, "\n\r \t");
			if (end == NULL) {
				oss = "#include not followed by newline.\n";
				return NULL;
			}
			*end = '\0';
			filelen = end - nextChar + 1;
			// Must copy the include path, because strtok will modify it
			value = getValue("INCLUDEPATH");
			if (value == NULL)
				value = includePathList;
			if (value == NULL) {
				path = NULL;
				directory = NULL;
			}
			else {
				path = new char[1 + strlen(value)];
				strcpy(path, value);
				directory = strtok(path, "|");
			}

			bool done = false;
			while (! done) {
				struct stat buf;
				int dirlen;

				if (directory == NULL)
					dirlen = 0;
				else
					dirlen = strlen(directory) + 1;
				char * filename = new char[dirlen + filelen];
				if (directory) {
					strcpy(filename, directory);
					strcat(filename, "/");
					directory = strtok(NULL, "|");
				}
				else {
					filename[0] = '\0';
					done = true;
				}
				strcat(filename, nextChar);
				if (stat(filename, &buf) != 0 || S_ISREG(buf.st_mode) == false) {
					delete [] filename;
					if (done) {
						oss = "Include file ";
						oss += nextChar;
						oss += " could not be found\n";
						delete [] path;
						return NULL;
					}
					continue;
				}

				bool ok = readFile(filename, oss, caseIndep);
				delete [] filename;
				if (! ok) {
					// An exception was thrown
					delete [] path;
					return NULL;
				}
				done = true;
			}
			delete [] path;

			// Update the string pointer
			file = end + 1;
		}
		else if (*keyStart == '#') {
			// If line starts with a comment character, which may be indented,
			// and that is not part of #include, skip to the next line.  Note
			// that comments are discarded, so if the file is rewritten, they
			// are lost.
			nextChar = keyEnd + strcspn(keyEnd, "\n\r");
			file = nextChar + strspn(nextChar, "\n\r");		// Skip the line
			delete [] key;
		}
		else {
			oss = "Key ";
			oss += key;
			oss += " not followed by an equal sign\n";
			delete [] key;
			return NULL;
		}
	}
}

/**
* Get a named registry folder.
*/
Registry *Registry::getFolder(const char *key)
{
	// If folder has dot, delegate to next folder.
	char * dot = (char *) strchr(key, '.');   // Solaris CC requires this cast
	if (dot) {
		*dot = 0;
		Registry * f = getFolder(key);
		*dot = '.';
		return f ? f->getFolder(dot + 1) : NULL;
	}
	else {
		char *trueKey = new char[strlen(key) + 3];
		strcpy(trueKey,"{}");
		strcat(trueKey,key);
		Registry * f = (Registry *)table.get(trueKey);
//		if(f==NULL) {
//			f=new Registry;
//			table.put(trueKey,f);
//		}
		//delete [] trueKey;
		return f;
	}
}

/**
* See if anything with this key exists
*/
bool Registry::hasKeyPrivate(char * key)
{
	// If folder has dot, delegate to next folder.
	char * dot = strchr(key, '.');
	if(dot) {
		*dot = 0;
		bool hasKey = hasKeyPrivate(key);
		if (hasKey) {
			hasKey = getFolder(key)->hasKeyPrivate(dot + 1);
		}
		*dot = '.';
		return hasKey;
	}
	else {
		char * trueKey = new char[strlen(key) + 3];
		strcpy(trueKey, "{}");
		strcat(trueKey, key);
		if (table.get(trueKey)) {
			delete [] trueKey;
			return true;
		}

		strcpy(trueKey, "=");
		strcat(trueKey, key);
		if (table.get(trueKey)) {
			delete [] trueKey;
			return true;
		}

		return false;
	}
}

Registry::Registry(int length, int width) : table(width, length/width)
{
#ifdef DEBUG
	cout << "New Registry at " << this << "; table = " << &table << '\n';
#endif
	includePathList = NULL;
	writable = true;
	orderingList = NULL;
	arrayWidth = 10;
}

Registry::Registry(const char * fname) throw(RException) : table(100, 40)
{
#ifdef DEBUG
	cout << "New Registry at " << this << "; table = " << &table << '\n';
#endif
	includePathList = NULL;
	writable = true;
	orderingList = NULL;
	arrayWidth = 10;
	(void) readFromFile(fname);
}

Registry::~Registry()
{
#ifdef DEBUG
	cout << "Discarded Registry at " << this << "; table = " << &table << '\n';
#endif
	clear();
}

/**
* Return value is NULL if value not found
*/
const char * Registry::getStringValue(const char * key, const char * defaultValue, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, defaultValue);
	vsprintf(tmpKey, key, val);
	va_end(val);

	char * value = getValue(tmpKey);
	if (value) 
		return value;

	return defaultValue;
}

bool Registry::cmpStringValue(char * key, char * cmpValue, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, cmpValue);
	vsprintf(tmpKey, key, val);
	va_end(val);

	char * value = getValue(tmpKey);

	if (value == NULL) {
		return (cmpValue == NULL);
	}
	if (cmpValue == NULL) {
		return false;
	}
	return (0 == strcmp(value, cmpValue));
}

/**
* Return value is -1 if value not found, 0 if found.
*/
int Registry::getIntValue(const char * key, int defaultValue, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, defaultValue);
	vsprintf(tmpKey, key, val);
	va_end(val);

	char * value = getValue(tmpKey);
	if (value)
		return atoi(value);

	return defaultValue;
}

/**
* Return value is -1 if value not found, 0 if found.
*/
float Registry::getFloatValue(const char * key, float defaultValue, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, defaultValue);
	vsprintf(tmpKey, key, val);
	va_end(val);

	char * value = getValue(tmpKey);
	if (value) 
		return (float) atof(value);

	return defaultValue;
}

/**
* Return value is -1 if value not found, 0 if found.
*/
double Registry::getDoubleValue(const char * key, double defaultValue, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, defaultValue);
	vsprintf(tmpKey, key, val);
	va_end(val);

	char * value = getValue(tmpKey);
	if (value) 
		return atof(value);

	return defaultValue;
}

/**
* Get a bool value
*/
bool Registry::getBooleanValue(const char * key, bool defaultValue, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, defaultValue);
	vsprintf(tmpKey, key, val);
	va_end(val);

	char * value = getValue(tmpKey);
	if (value) 
		return (atoi(value) > 0);

	return defaultValue;
}

/**
* Return value is NULL if array not found
*/
int * Registry::getIntArray(const char * key, int * length, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, length);
	vsprintf(tmpKey, key, val);
	va_end(val);

	char * fullKey = new char[strlen(tmpKey) + 2];
	strcpy(fullKey, tmpKey);
	char * p = strrchr(fullKey, '.');
	if (p == NULL)
		p = fullKey;
	else
		p++;
	strcpy(p, "@");
	p++;
	strcpy(p, tmpKey + (p - fullKey) - 1);

	char * arrayStr = getValue(fullKey);
	delete [] fullKey;
	if (arrayStr == NULL) {
		*length = 0;
		return NULL;
	}

	int type;
	int * array;
	int n = sscanf(arrayStr, "%d %d %x", length, &type, &array);
	if (n != 3 && type != (int) IntArray) {
		*length = 0;
		return NULL;
	}

	char * addr = strrchr(arrayStr, ' ');
	addr++;
	addr[0] = '0';	// Mark the array as undeletable
	addr[1] = '\0';

	return array;
}

/**
* Return value is NULL if array not found
*/
float * Registry::getFloatArray(const char * key, int * length, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, length);
	vsprintf(tmpKey, key, val);
	va_end(val);

	char * fullKey = new char[strlen(tmpKey) + 2];
	strcpy(fullKey, tmpKey);
	char * p = strrchr(fullKey, '.');
	p++;
	strcpy(p, "@");
	p++;
	strcpy(p, tmpKey + (p - fullKey) - 1);

	char * arrayStr = getValue(fullKey);
	delete [] fullKey;
	if (arrayStr == NULL) {
		*length = 0;
		return NULL;
	}

	int type;
	float * array;
	int n = sscanf(arrayStr, "%d %d %x", length, &type, &array);
	if (n != 3 && type != (int) FloatArray) {
		*length = 0;
		return NULL;
	}

	char * addr = strrchr(arrayStr, ' ');
	addr++;
	addr[0] = '0';	// Mark the array as undeletable
	addr[1] = '\0';

	return array;
}

/**
* Return value is NULL if array not found
*/
double * Registry::getDoubleArray(const char * key, int * length, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, length);
	vsprintf(tmpKey, key, val);
	va_end(val);

	char * fullKey = new char[strlen(tmpKey) + 2];
	strcpy(fullKey, tmpKey);
	char * p = strrchr(fullKey, '.');
	p++;
	strcpy(p, "@");
	p++;
	strcpy(p, tmpKey + (p - fullKey) - 1);

	char * arrayStr = getValue(fullKey); // here lies the problem
	
	
	delete [] fullKey;
	if (arrayStr == NULL) {
		*length = 0;
		return NULL;
	}
	
	

	int type;
	double * array;
	int n = sscanf(arrayStr, "%d %d %x", length, &type, &array);
	if (n != 3 && type != (int) DoubleArray) {
		*length = 0;
		return NULL;
	}

	char * addr = strrchr(arrayStr, ' ');
	addr++;
	addr[0] = '0';	// Mark the array as undeletable
	addr[1] = '\0';

	
	return array;
	
}

/**
* Place a copy of a string into the registry
*/
void Registry::setStringValue(const char * key, const char * value, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, value);
	vsprintf(tmpKey, key, val);
	va_end(val);

	// Create new value string 
	char * newValue = new char[strlen(value)+1];
	strcpy(newValue, value);
	setValue(tmpKey, newValue);
}

/**
* Set integer value
*/
void Registry::setIntValue(const char * key, int value, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, value);
	vsprintf(tmpKey, key, val);
	va_end(val);

	// Create new value string 
	char * str = new char[12];
	sprintf(str, "%d", value);
	setValue(tmpKey, str);
}

/**
* Set double value
*/
void Registry::setFloatValue(const char * key, float value, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, value);
	vsprintf(tmpKey, key, val);
	va_end(val);

	// Create new value string 
	char * str = new char[64];
	sprintf(str, "%.8g", value);
	setValue(tmpKey, str);
}

/**
* Set double value
*/
void Registry::setDoubleValue(const char * key, double value, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, value);
	vsprintf(tmpKey, key, val);
	va_end(val);

	// Create new value string 
	char * str = new char[64];
	sprintf(str, "%.17lg", value);
	setValue(tmpKey, str);
}

/**
* Set boolean value
*/
void Registry::setBooleanValue(const char * key, bool value, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, value);
	vsprintf(tmpKey, key, val);
	va_end(val);

	// Create new value string 
	char * str = new char[24];
	sprintf(str, "%d", (value) ? 1 : 0);
	setValue(tmpKey, str);
}

/**
* Set int array
*/
void Registry::setIntArray(const char * key, int length, int * array, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, array);
	vsprintf(tmpKey, key, val);
	va_end(val);

	// Create new value string 
	char * arrayStr = new char[32];
	sprintf(arrayStr, "%d %d %x", length, (int) IntArray, array);

	// Create new key string 
	char * fullKey = new char[strlen(tmpKey) + 2];
	strcpy(fullKey, tmpKey);
	char * p = strrchr(fullKey, '.');
	if (p == NULL)
		p = fullKey;
	else
		p++;
	strcpy(p, "@");
	p++;
	strcpy(p, tmpKey + (p - fullKey) - 1);
	setValue(fullKey, arrayStr);
	delete [] fullKey;
}

/**
* Set float array
*/
void Registry::setFloatArray(const char * key, int length, float * array, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, array);
	vsprintf(tmpKey, key, val);
	va_end(val);

	// Create new value string 
	char * arrayStr = new char[32];
	sprintf(arrayStr, "%d %d %x", length, (int) FloatArray, array);

	// Create new key string 
	char * fullKey = new char[strlen(tmpKey) + 2];
	strcpy(fullKey, tmpKey);
	char * p = strrchr(fullKey, '.');
	if (p == NULL)
		p = fullKey;
	else
		p++;
	strcpy(p, "@");
	p++;
	strcpy(p, tmpKey + (p - fullKey) - 1);
	setValue(fullKey, arrayStr);
	delete [] fullKey;
}

/**
* Set double array
*/
void Registry::setDoubleArray(const char * key, int length, double * array, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, array);
	vsprintf(tmpKey, key, val);
	va_end(val);

	// Create new value string 
	char * arrayStr = new char[32];
	sprintf(arrayStr, "%d %d %x", length, (int) DoubleArray, array);

	// Create new key string 
	char * fullKey = new char[strlen(tmpKey) + 2];
	strcpy(fullKey, tmpKey);
	char * p = strrchr(fullKey, '.');
	if (p == NULL)
		p = fullKey;
	else
		p++;
	strcpy(p, "@");
	p++;
	strcpy(p, tmpKey + (p - fullKey) - 1);
	setValue(fullKey, arrayStr);
	delete [] fullKey;
}

Registry & Registry::getSubFolder(const char * key, ...) {
	// Initialize arg list
	va_list val;
	va_start(val, key);
	vsprintf(tmpKey, key, val);
	va_end(val);

	return *getFolder(tmpKey);
}

bool Registry::hasKey(const char * key, ...) {
   // Initialize arg list
	va_list val;
	va_start(val, key);
	vsprintf(tmpKey, key, val);
	va_end(val);

	return hasKeyPrivate(tmpKey);
}

void Registry::ordering(const char ** order)
{
	orderingList = order;
	if (orderingList != NULL)
		if (*orderingList == NULL)
			orderingList = NULL;
}

/**
	Print the registry like writeToFile(), but to cout
	Errors will not be reported.
*/
bool Registry::printAsRegistry()
{
	return write(&cout, 0);
}

/**
* Write the folder to a disk
*/
bool Registry::writeToFile(const char * pathname) throw(RException)
{

	if (! writable) {
		string s("File ");
		s += string(pathname);
		s += " contains #include:\n    changes will not be recorded.";
		throw RException(s);
	}

	std::ofstream f;
	f.open(pathname);
	if (!f) {
		string s("Unable to open ");
		s += string(pathname);
		s += " for writing.";
		throw RException(s);
	}

	bool ret = write(&f, 0);	

	f.close();
	
	return ret;	

}

/**
* Read folder from file
*/
void Registry::readFromFile(const char * pathname, bool caseIndep)
	throw(RException)
{
	string oss;

	// cout << "Registry::readFromFile(" << pathname << ')' << endl;

#ifdef DEBUG
	cout << "Registry::readFromFile(" << pathname << ')' << endl;
#endif
	if (pathname == NULL) {
		string s("No input file.");
		throw RException(s);
		return;
	}
	if (readFile(pathname, oss, caseIndep) == false)
		throw RException(oss);
}

void upper(char * file)
{
	while(1) {
		// We are expecting a key - find its beginnning and end
		char * keyStart = file + strspn(file, " \t\n\r");

		// If key starts with closing brace, we are done
		if (*keyStart == '}' || *keyStart == '\0')
			return;

		char * keyEnd = keyStart + strcspn(keyStart, " \t\n\r={");

		// Note: both comments and include files are ignored here
		if (*keyStart == '#') {
			// If line starts with a comment character, which may be indented,
			// and that is not part of #include, skip to the next line.  Note
			// that comments are discarded, so if the file is rewritten, they
			// are lost.
			char * nextChar = keyEnd + strcspn(keyEnd, "\n\r");
			file = nextChar + strspn(nextChar, "\n\r");		// Skip the line
		}
		else {
			// Map the first character of the key, if necessary
			*keyStart = toupper(*keyStart);

			// Check the next character after the key
			char * nextChar = keyEnd + strspn(keyEnd, " \t\n\r");

			// See if this is a name-value assignment
			if (*nextChar == '=') {
				// We got us a key-value pair
				char * valueStart = nextChar + 1 + strspn(nextChar + 1, " \t\n\r");
				char * valueEnd = strchr(valueStart, ';');

				// Update the string pointer
				file = valueEnd+1;
			}
			else
				return;
		}
	}
}

/**
* Read folder from file.
  This returns false, if it fails.
*/
bool Registry::readFile(const char * pathname, string & oss, bool caseIndep)
{
    struct stat buf;
   // std::cout<<"Registry::readFile is: "<<pathname<<std::endl;
    if (stat(pathname, &buf) != 0 || S_ISREG(buf.st_mode) == false) {
		char msg[1000];
		sprintf(msg, "File does not exist (Registry:readFile): %s", pathname);
		oss = msg;
		return false;
	}

	FILE * f = fopen(pathname, "rb");
	if (! f) {
		oss = "Unable to open ";
		oss += string(pathname);
		oss += " for reading.";
		return false;
	}
	fseek(f, 0, SEEK_END);
	int size = ftell(f);
	fseek(f, 0, SEEK_SET);

	char * buffer = new char[size+1];
	fread(buffer, 1, size, f);
	buffer[size] = 0;
	fclose(f);

	if (caseIndep)
		upper(buffer);

	char * rc = read(buffer, oss, false);	// Parse the file's contents

	// cout << "Registry::read " << pathname << " completed successfully" << endl ;	

	delete [] buffer;

	return rc == NULL ? false : true;	
}

char *Registry::tmpKey = new char[1024];


// Begin namespace
NAMESPACE_PAULY_END


