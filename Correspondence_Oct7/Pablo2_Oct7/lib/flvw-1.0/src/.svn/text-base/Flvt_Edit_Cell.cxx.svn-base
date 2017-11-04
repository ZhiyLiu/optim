
#include <string.h>
#include <FL/Flvt_Edit_Cell.H>
#include <stdio.h>							// AGG: TEMPORARY


int Flvt_Edit_Cell::global_width = 0;
void * Flvt_Edit_Cell::ptr = NULL;


Flvt_Edit_Cell::~Flvt_Edit_Cell() {
#ifdef DEBUG
	printf("Flvt_Edit_Cell::~Flvt_Edit_Cell()\n");
#endif
	delete [] v;
}

Flvt_Edit_Cell::Flvt_Edit_Cell() {
#ifdef DEBUG
	printf("Flvt_Edit_Cell::Flvt_Edit_Cell()\n");
#endif
	edited = false;
	callback = NULL;
	if (global_width == 0) {
		w = 0;
		v = NULL;
	}
	else {
		w = global_width;
		v = new char[1 + w];
		v[0] = '\0';
	}
}

Flvt_Edit_Cell::Flvt_Edit_Cell(const char * value, int maxWidth) {
#ifdef DEBUG
	printf("Flvt_Edit_Cell::Flvt_Edit_Cell(const char * value, int maxWidth)\n");
#endif
	edited = false;
	callback = NULL;
	if (maxWidth == 0) {
		if (global_width == 0) {
			if (value == NULL) {
				w = 0;
				v = NULL;
				return;
			}
			else
				w = strlen(value);
		}
		else
			w = global_width;
	}
	else
		w = maxWidth;

	v = new char[1 + w];
	if (value == NULL)
		v[0] = '\0';
	else
		(void) strcpy(v, value);
}

void Flvt_Edit_Cell::set(const char * value) {
#ifdef DEBUG
	printf("Flvt_Edit_Cell::set()\n");
#endif
	if (value == NULL) {
		if (v == NULL)
			return;
		else
			v[0] = '\0';
	}
	else {
		if (v == NULL) {
			if (global_width == 0) {
				w = strlen(value);
				v = new char[1 + w];
				(void) strcpy(v, value);
			}
			else {
				w = global_width;
				v = new char[1 + w];
				(void) strncpy(v, value, global_width);
				v[global_width] = '\0';
			}
		}
		else {
			(void) strncpy(v, value, w);
			v[w] = '\0';
		}
	}
}

void Flvt_Edit_Cell::get(char * value) {
	if (v == NULL)
		value[0] = '\0';
	else
		(void) strcpy(value, v);
}

void Flvt_Edit_Cell::setCommonWidth(int maxWidth) {
	// This permits the global width to be changed as the cells are allocated.
	global_width = maxWidth;
}

void Flvt_Edit_Cell::setCellChangeCallback(void (* fn)(void *, int), int tag) {
	callback = fn;
	data = tag;
}

bool Flvt_Edit_Cell::changed() {
	if (edited) {
		edited = false;
		return true;
	}
	else
		return false;
}
