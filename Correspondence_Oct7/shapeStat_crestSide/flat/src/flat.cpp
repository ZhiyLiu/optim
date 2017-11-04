#include <iostream>
#include <fstream>
//#include "registry.h"
#include "Registry.h"
#include <string.h>

#define LINESIZE 8192        // Big enough for an array

using namespace std;

void usage(char * reason)
{
	char * usageString =
	    "\nUsage: flat [-h] [-s] [-u] [-v] registryFile flatFile\n\n\
Convert a Pablo registry file to unnested, key-value (flat) text, or\n\
vice-versa.\n\
    -h             Print this help message.\n\
    -s             Print a summary of the program's actions.\n\
    -u             Inflate a flat file to a registry file.\n\
    -v             Verbose; list the keys as they are processed.\n\
                       May only be used when no file is written to the\n\
                       standard output.\n\
    registryFile   The registry file.  This will typically be a Pablo\n\
                       model (.m3d) file.\n\
    flatFile       The flat file or `-'.  If a minus sign is used, the\n\
                       flat file will be written to the standard output\n\
                       or read from the standard input.\n\n\
    The default action is to flatten the registry file.  To go the other\n\
    direction, the -u option must be used.\n\n";

    if (reason && *reason)
        cout << "*** " << reason << endl;
    cout << usageString;
}

void reportSummary(char * regFile, int plain, int fancy, bool direction)
{
	cout << "Summary:\n";
	cout << "   " << plain << " regular keys and ";
	cout << fancy << " array keys were ";
	if (direction)
		cout << "flattened";
	else
		cout << "made hierachical";
	cout << ".\n";
}

int main(int argc, char **argv)
{
    Registry reg;
	int len, i;
	int plain, fancy;

    bool unflatten = false;
    bool summary = false;
    bool verbose = false;
    bool error = false;
    char *regFile = NULL;
    char *flatFile = NULL;
	char *order[3] = {"pabloVersion", "model", NULL};
	char so[] = "stdIO";	// Just an indicater, not an actual stream name
    string opt;
	std::ofstream ofs;
	std::ifstream ifs;
	char * space;


	if (argc < 2) {
		usage(NULL);
		return 0;
	}

    for (i = 1; i < argc; i++)
    {
        opt = argv[i];
		if (opt[0] != '-') {
			if (regFile == NULL)
				regFile = argv[i];
			else
				flatFile = argv[i];
			continue;
		}

        if (opt == "-")
        {
            if (regFile == NULL)
				regFile = so;
			else
				flatFile = so;
            continue;
        }
        if (opt == "-h" || opt == "-help")
        {
            usage(NULL);
			return 1;
        }
        else if (opt == "-s")
        {
            summary = true;
            continue;
        }
        else if (opt == "-u")
        {
            unflatten = true;
            continue;
        }
        else if (opt == "-v")
        {
            verbose = true;
            continue;
        }
		else {
			cout << "Error: unknown option " << opt << '\n';
			usage(NULL);
			return 1;
		}
    }
	if (regFile == so && flatFile == so) {
		usage("Error: at least one file name must be provided");
		return 1;
	}
	if (regFile == so) {
		usage("Error: the registry file cannot be written to the standard output");
		return 1;
	}
	if (flatFile == so && ! unflatten) {
		if (verbose || summary) {
			usage("Error: cannot use -s or -v with a conversion to the standard output");
			return 1;
		}
	}

	reg.ordering(const_cast<const char **>(order));
	plain = 0;
	fancy = 0;

    if (! unflatten)    // Registry -> flat file (stdout)
    {
		char ** list;

        // Read registry, print each key-value pair
		try {
	        reg.readFromFile(regFile);
		}
		catch (RException excp) {
			excp.print(std::cout);
			return 1;
		}
        list = reg.collectAllKeys(len, true);
		if (verbose) {
			for (int k = 0; k < len; k++) {
				space = strchr(list[k], ' ');
				*space = '\0';
				cout << list[k] << '\n';
				*space = ' ';
			}
		}

		if (summary) {
			for (int k = 0; k < len; k++) {
				if (strchr(list[k], '{') == NULL)
					plain++;
				else
					fancy++;
			}
			reportSummary(regFile, plain, fancy, true);
		}

		if (flatFile == so) {
			for (int k = 0; k < len; k++)
				cout << list[k] << endl;
		}
		else {
			ofs.open(flatFile);
			for (int k = 0; k < len; k++)
				ofs << list[k] << endl;
			ofs.close();
		}
    }
    else   // Flat file -> registry
    {
        char key[LINESIZE];

		if (flatFile != so) {
			ifs.open(flatFile);
			if (! ifs.is_open()) {
				cout << "Error: could not open file " << flatFile << endl;
				return 1;
			}
		}
        while (true)
		{
			bool fail = false;
			if (flatFile != so) {
				ifs.getline(key, LINESIZE);
				if (ifs.eof())
					break;
				if (ifs.fail())
					fail = true;
			}
			else {
				cin.getline(key, LINESIZE);
				if (cin.eof())
					break;
				if (cin.fail())
					fail = true;
			}
			if (ifs.fail()) {
                cout << "Error: value too long for key " << key << '\n'
					<< "Program must be recompiled with a larger buffer to process this file\n";
				return 1;
			}

            space = strchr(key, ' ');
			if (space == NULL) {
                cout << "Warning: no space char found.  Skipping line:\n\t" << key << '\n';
				continue;
			}
            *space = '\0';   // null-terminate the key

            if (verbose)
				cout << key << '\n';

            // Insert key and value into registry
			char * val = space + 1;
			if (val[0] == '{') {
				int length, type;
				char * array;

				fancy++;
				array = reg.parseArray(val + 1, length, type);
				if (array == NULL) {
					cout << "Error processing array " << key << endl;
					return 1;
				}
				switch (type) {
					case 0:
						reg.setIntArray(key, length, (int *) array);
						break;
					case 1:
						reg.setFloatArray(key, length, (float *) array);
						break;
					case 2:
						reg.setDoubleArray(key, length, (double *) array);
						break;
					default:
						cout << "Error processing array " << key << endl;
						return 1;
				}
			}
			else {
				plain++;
				reg.setStringValue(key, val, NULL);
			}
        }
		reg.writeToFile(regFile);
		if (summary)
			reportSummary(regFile, plain, fancy, false);
    }

    return 0;
}

