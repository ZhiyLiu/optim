
#include <stdlib.h>
#include <iostream>
#include <string>
#include "ControlParms.h"
#include "ControlFile.h"


using namespace std;

#define UNCHANGED	0
#define STRING		1
#define BOOLEAN		2
#define INTEGER		3
#define DOUBLE		4

#include <string.h>

ControlParmsAccess::ControlParmsAccess(const char * application_name,
    int num_ctrl_names, const char ** ctrl_names, const char * ctrl_dir,
    int verbosity, bool load)
{
    num_items = num_ctrl_names;
    control_names = ctrl_names;
    changed = new int [num_items];
    loaded = new bool [num_items];
    values = new value_t [num_items];
    allocated = new bool [num_items];

    control_file = (void *) new ControlFile(application_name, ctrl_dir,
        verbosity, load);
    for (int i = 0; i < num_items; i++) {
	loaded[i] = false;
	changed[i] = UNCHANGED;
	allocated[i] = false;
    }
}


void ControlParmsAccess::setVerbosity(int verbosity)
{
    ControlFile * cf = (ControlFile *) control_file;
    cf->setVerbosity(verbosity);
}


ControlParmsAccess::~ControlParmsAccess() {
    int key;

    ControlFile * cf = (ControlFile *) control_file;
    for (key = 0; key < num_items; key++) {
	    switch (changed[key]) {
	        case UNCHANGED:	break;

	        case STRING:	cf->write(control_names[key], values[key].s);
				    break;

	        case BOOLEAN:	cf->write(control_names[key], values[key].b);
				    break;

	        case INTEGER:	cf->write(control_names[key], values[key].i);
				    break;

	        case DOUBLE:	cf->write(control_names[key], values[key].d);
				    break;
	    }
    }

    delete (ControlFile *) control_file;
    for (key = 0; key < num_items; key++) {
	if (allocated[key])
	    delete [] values[key].s;
    }

    delete [] changed;
    delete [] loaded;
    delete [] values;
    delete [] allocated;
}

const char * ControlParmsAccess::path() const {
    ControlFile * cf = (ControlFile *) control_file;
    return cf->controlFilename();
}

void ControlParmsAccess::setDefault(control_keys_t key, const char * deflt) {
    ControlFile * cf = (ControlFile *) control_file;
    const char * str = cf->read(control_names[key], (char *) deflt);
    values[key].s = new char[1 + strlen(str)];
    strcpy(values[key].s, str);
    loaded[key] = true;
    allocated[key] = true;
}


void ControlParmsAccess::setDefault(control_keys_t key, bool deflt) {
    ControlFile * cf = (ControlFile *) control_file;
    values[key].b = cf->read(control_names[key], deflt);
    loaded[key] = true;
}


void ControlParmsAccess::setDefault(control_keys_t key, int deflt) {
    ControlFile * cf = (ControlFile *) control_file;
    values[key].i = cf->read(control_names[key], deflt);
    loaded[key] = true;
}


void ControlParmsAccess::setDefault(control_keys_t key, double deflt) {
    ControlFile * cf = (ControlFile *) control_file;
    values[key].d = cf->read(control_names[key], deflt);
    loaded[key] = true;
}


char * ControlParmsAccess::readString(control_keys_t key) {
    ControlFile * cf = (ControlFile *) control_file;
    if (! loaded[key]) {
	cout << "Key " << control_names[key]
	    << " was read before its default was specified" << endl;
	const char * str = cf->read(control_names[key], "");	// Arbitrary default
	values[key].s = new char[1 + strlen(str)];
	strcpy(values[key].s, str);
	loaded[key] = true;
	allocated[key] = true;
    }
    return values[key].s;
}


bool ControlParmsAccess::readBool(control_keys_t key) {
    ControlFile * cf = (ControlFile *) control_file;
    if (! loaded[key]) {
	cout << "Key " << control_names[key]
	    << " was read before its default was specified" << endl;
	values[key].b = cf->read(control_names[key], false);	// Arbitrary default
	loaded[key] = true;
    }
    return values[key].b;
}


int ControlParmsAccess::readInt(control_keys_t key) {
    ControlFile * cf = (ControlFile *) control_file;
    if (! loaded[key]) {
	cout << "Key " << control_names[key]
	    << " was read before its default was specified" << endl;
	values[key].i = cf->read(control_names[key], 0);	// Arbitrary default
	loaded[key] = true;
    }
    return values[key].i;
}


double ControlParmsAccess::readDouble(control_keys_t key) {
    ControlFile * cf = (ControlFile *) control_file;
    if (! loaded[key]) {
	cout << "Key " << control_names[key]
	    << " was read before its default was specified" << endl;
	values[key].d = cf->read(control_names[key], 0.0);	// Arbitrary default
	loaded[key] = true;
    }
    return values[key].d;
}


char * ControlParmsAccess::read(control_keys_t key, const char * deflt) {
    ControlFile * cf = (ControlFile *) control_file;
    if (! loaded[key]) {
	const char * str = cf->read(control_names[key], (char *) deflt);
	values[key].s = new char[1 + strlen(str)];
	strcpy(values[key].s, str);
	loaded[key] = true;
	allocated[key] = true;
    }
    return values[key].s;
}


bool ControlParmsAccess::read(control_keys_t key, bool deflt) {
    ControlFile * cf = (ControlFile *) control_file;
    if (! loaded[key]) {
	values[key].b = cf->read(control_names[key], deflt);
	loaded[key] = true;
    }
    return values[key].b;
}


int ControlParmsAccess::read(control_keys_t key, int deflt) {
    ControlFile * cf = (ControlFile *) control_file;
    if (! loaded[key]) {
	values[key].i = cf->read(control_names[key], deflt);
	loaded[key] = true;
    }
    return values[key].i;
}


double ControlParmsAccess::read(control_keys_t key, double deflt) {
    ControlFile * cf = (ControlFile *) control_file;
    if (! loaded[key]) {
	values[key].d = cf->read(control_names[key], deflt);
	loaded[key] = true;
    }
    return values[key].d;
}


void ControlParmsAccess::write(control_keys_t key, const char * val) {
    if (allocated[key])
	delete [] values[key].s;
    values[key].s = new char[strlen(val) + 1];
    strcpy(values[key].s, val);
    loaded[key] = true;
    allocated[key] = true;
    changed[key] = STRING;
}


void ControlParmsAccess::write(control_keys_t key, bool val) {
    values[key].b = val;
    loaded[key] = true;
    changed[key] = BOOLEAN;
}


void ControlParmsAccess::write(control_keys_t key, int val) {
    values[key].i = val;
    loaded[key] = true;
    changed[key] = INTEGER;
}


void ControlParmsAccess::write(control_keys_t key, double val) {
    values[key].d = val;
    loaded[key] = true;
    changed[key] = DOUBLE;
}


void ControlParmsAccess::remove(control_keys_t key)
{
	if (allocated[key]) {
	    delete [] values[key].s;
		allocated[key] = false;
	}
    loaded[key] = false;
	ControlFile * cf = (ControlFile *) control_file;
	cf->remove(control_names[key]);
	changed[key] = UNCHANGED;
}

