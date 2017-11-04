
//  Program for testing the Registry class

//  To use it, build the paul_code library using paul_code.dsw,
//  and then build this program using testRegistry.dsw.


#include <string>
#include <iostream>
#include "Registry.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

using namespace std;

#define NFIGS   5
#define NATOMS  4

int main(int argc, char **argv)
{
	Registry r(4, 2);
	int i;

	r.setIntValue("tree.item",0);
	r.setIntValue("model.figureCount", NFIGS + NFIGS);
	r.setIntValue("model.size", 33);
	cout << "Initial level 0 size = " << r.mainSize() << '\n';
	r.deleteKey("model.size");
	r.deleteKey("tree.item");
	cout << "Final level 0 size = " << r.mainSize() << '\n';
    cout << "First dump:\n"; r.dump();

	for (i = 0; i < NFIGS; i++) {
		r.setIntValue("model.figure[%d].primitiveCount", NATOMS, i);
		r.setStringValue("model.figure[%d].Removal", "hello", i);
		for(int j = 0; j < NATOMS; j++) {
			r.setDoubleValue("model.figure[%d].primitive[%d].x",((double)rand()) / RAND_MAX,i,j);
			r.setDoubleValue("model.figure[%d].primitive[%d].y",((double)rand()) / RAND_MAX,i,j);
			r.setDoubleValue("model.figure[%d].primitive[%d].s",((double)rand()) / RAND_MAX,i,j);
		}
	}
	cout << "Initial full size = " << r.size() << '\n';
	for (i = 0; i < NFIGS; i++) {
        string s;
        s = "model.figure[";
        s += (char) ('0' + i); 
        s += "].Removal";
        r.deleteKey(s.data());
    }

	for (i = NFIGS; i < NFIGS + NFIGS; i++) {
		r.setIntValue("model.figure[%d].primitiveCount", NATOMS, i);
		for(int j = 0; j < NATOMS; j++) {
			r.setDoubleValue("model.figure[%d].primitive[%d].x",((double)rand()) / RAND_MAX,i,j);
			r.setDoubleValue("model.figure[%d].primitive[%d].y",((double)rand()) / RAND_MAX,i,j);
			r.setDoubleValue("model.figure[%d].primitive[%d].s",((double)rand()) / RAND_MAX,i,j);
		}
	}
	r.setIntValue("bush",1);
	r.setIntValue("shrub",2);
    r.deleteKey("model");

    int *ia = new int[2];
    ia[0] = 1;
    ia[1] = 2;
    r.setIntArray("array", 2, ia);
    r.deleteKey("array");
    r.deleteKey("@array");

    cout << "Final full size = " << r.size() << '\n';
	if (r.writeToFile("folder.dat"))
        cout << "Wrote folder.dat\n";
	else
        cout << "Could not write output file\n";
    cout << "Major keys are:\n";
	r.printMajorKeys();
//    cout << "Final dump:\n";
//    r.dump();
	return 0;
}

