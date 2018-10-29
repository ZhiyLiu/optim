#include <string.h>
#include <iostream>
#include "Tuning.h"
#include "Registry.h"

using namespace std;

#include <stdlib.h>

Tuning::Tuning()
{
	v = NULL;
	numEntries = 0;
	index = 0;
}

Tuning::~Tuning()
{
	if (v != NULL)
		delete [] v;
}

bool Tuning::setNext(char * shortName, char * longName, double defaultWeight, char * desc,
	bool lock)
{
	if (index >= numEntries)
		return false;
	v[index].shortName = shortName;
	v[index].longName = longName;
	v[index].defaultWeight = defaultWeight;
	v[index].desc = desc;
	v[index].locked = lock;
	index++;
	return true;
}

bool Tuning::setNextList(char * shortName, char * longName, char * desc,
	bool lock)
{
	if (index >= numEntries)
		return false;
	v[index].shortName = shortName;
	v[index].longName = longName;
	v[index].defaultWeight = 0.0;	// Not really needed
	v[index].weight.lst = NULL;
	v[index].listLen = 0;
	v[index].desc = desc;
	v[index].locked = lock;
	index++;
	return true;
}

bool Tuning::setNextString(char * shortName, char * longName, char * desc,
	bool lock)
{
	if (index >= numEntries)
		return false;
	v[index].shortName = shortName;
	v[index].longName = longName;
	v[index].defaultWeight = 0.0;	// Not really needed
	v[index].weight.str = NULL;
	v[index].listLen = -2;
	v[index].desc = desc;
	v[index].locked = lock;
	index++;
	return true;
}

bool Tuning::setWeight(int index, double weight) {
	if (! isWeight(index)) {
		cout << "Tuning parameter " << index << " is not a weight\n";
		return false;
	}
	if (! v[index].locked) {
		v[index].weight.d = weight;
		return true;
	}
	return false;
}

bool Tuning::setString(int index, const char * s)
{
	if (! isString(index)) {
		cout << "Tuning parameter " << index << " is not a string\n";
		return false;
	}
	if (! v[index].locked) {
		if (v[index].weight.str != NULL)
			delete [] v[index].weight.str;
		int l = strlen(s);
		v[index].weight.str = new char[l + 1];
		strcpy(v[index].weight.str, s);
		return true;
	}
	return false;
}

bool Tuning::setList(int index, int len, int * newList)
{
	if (! isList(index)) {
		cout << "Tuning parameter " << index << " is not a list\n";
		return false;
	}
	if (! v[index].locked) {
		if (len >= MAX_TUNING_LIST_LEN) {
			cout << "Maximum list length for tuning parameters exceeded" << endl;
			return false;
		}
		if (v[index].weight.lst == NULL)
			v[index].weight.lst = new int[MAX_TUNING_LIST_LEN];
		for (int i = 0; i < len; i++)
			v[index].weight.lst[i] = newList[i];
		v[index].listLen = len;
		return true;
	}
	return false;
}

// This ignores locking because it is a reinitialization
void Tuning::resetWeight(int index)
{
	if (v[index].listLen == -1)
		v[index].weight.d = v[index].defaultWeight;
	else {
		// Entry is a list or string; either way, a pointer
		delete [] v[index].weight.str;
		v[index].weight.str = NULL;
		if (v[index].listLen >= 0)
			v[index].listLen = 0;
	}
}

void Tuning::printAsRegistry()
{
	int i;

	// Print the Pablo version number first
	cout << v[0].longName << " = ";
	printVersion(v[0].weight.d);
	cout << '\n';

	Registry r(numEntries);
	for (i = 1; i < numEntries; i++) {
		if (isWeight(i))
			r.setDoubleValue(v[i].longName, v[i].weight.d);
		else if (isString(i)) {
			if (v[i].weight.str != NULL)
				r.setStringValue(v[i].longName, v[i].weight.str);
		}
		else {
			if (v[i].weight.lst != NULL)
				r.setIntArray(v[i].longName, v[i].listLen, v[i].weight.lst);
		}
	}
	(void) r.printAsRegistry();	// This does a flush
}

void Tuning::printWeight(int i)
{
	if (isString(i)) {
		char *str = v[i].weight.str;
		cout << ((str && *str) ? str : "NULL");
	} else if (isList(i)) {
		cout << "(";
		for (int ele=0; ele < v[i].listLen; ele++)
			cout << v[i].weight.str[i] << " ";
		cout << ")";
	} else if (isWeight(i))
		cout << v[i].weight.d;
}

// print long lines of text folded into lines of maxlen chars,
// broken at white space and indented:
// - maxlen=0 disables all line breaks: good for grep'ing results.
// - exclude newline from msg -- printFold adds them as needed.
void printFold(const char* msg, int maxlen)
{
	if (!maxlen) {
		cout << msg << endl;
		return;
	}

	// copy string to replace spaces with \0 via strtok
	char *str = strdup(msg);

	int lineLen = 0;	// chars seen on this line
	char *token = strtok(str, " ");
	while (token)
	{
		int tokenLen = strlen(token);

		// if this word exceeds this line's capacity, output line delimiter
		if ( (lineLen+tokenLen)>maxlen )
		{
			cout << endl << "    ";
			lineLen = 4;	// 4 spaces indentation
		}
		cout << token << " ";
		lineLen += tokenLen + 1;	// include space
		token = strtok(NULL, " ");
	}
	cout << endl;
	free(str);
}

void Tuning::printDesc()
{
	const int string_max_len = (int) v[LineLength].weight.d;
	
	for (int i = 0; i < numEntries; i++) {
		char msg[1000];
		sprintf(msg, "%s: %s", v[i].longName, v[i].desc);
		printFold(msg, string_max_len);
	}
}

void Tuning::printAsTDTable()
{
	int i;

	// Print the Pablo version number first
	cout << v[0].longName << " = ";
	printVersion(v[0].weight.d);
	cout << '\n';

	for (i = 1; i < numEntries; i++)
	{
		cout << v[i].longName << " = ";
		printWeight(i);
		cout << ";\t\t" << v[i].desc << '\n';
	}
	cout << flush;
}

int Tuning::indexByLongName(const char * name) const
{
	for (int i = 0; i < numEntries; i++) {
#ifdef _WIN32
		if (_stricmp(name, v[i].longName) == 0)
#else
//		if (strcasecmp(name, v[i].longName) == 0)
#endif
			return i;
	}
	return -1;
}

void Tuning::print(int index, bool full)			// AGG: This does not handle array entries properly
{
		cout << "Tuning entry " << index << " = {\n";
		cout << "    shortName = " << v[index].shortName << '\n';
		if (full)
			cout << "    longName = " << v[index].longName << '\n';
		cout << "    defaultWeight = " << v[index].defaultWeight << '\n';
		if (full)
			cout << "    desc = " << v[index].desc << '\n';
		cout << "    weight = " << v[index].weight.d << '\n';
		if (full)
			cout << "    locked = " << v[index].locked << '\n';
		cout << "}\n";
}

// Print a Pablo version number from the value stored in class Tuning
void printVersion(double weight, int precision)
{
	cout.precision(precision);
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout << weight;
	cout.setf(~ios_base::fixed, ios_base::floatfield);
	cout.precision(6);
}


