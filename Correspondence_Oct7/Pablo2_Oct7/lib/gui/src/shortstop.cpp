// the shortstop application: short utilities.
// DO NOT CHECK LICENSE KEYS in shortstop

//#define DEBUG

#include "pablo_version.h"
#include "P3DControl.h"
#include "ControlParms.h"
#include "Registry.h"
#include "globalBuildParms.h"
#include "Tuning.h"

using namespace std;

extern char * headerString;
extern char * longChangesString;
extern char * shortChangesString;


extern bool getarg(char * shortOption, char * longOption, bool takesArg, string & opt,
    int & argNum, bool & error, bool verbose, int argc, char **argv);


static void printInputsAsTDTable()
{
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
		<< "[BYU FILE] simple marching cubes tiling of binary image. Always uses Quads" << endl;
	cout << "outModelImage = ;" << "\t" << "-omi" << "\t"
		<< "[RAW3 FILE] scan-convert model to image" << endl;

	// undocumented
	//	cout << "distanceMap = ;" << "\t" << "-id" << "\t"
	//		<< "[raw3 FILE] distance map image. Read if exists; written if not" << endl;

}


int shortstop(int argc, char *argv[])
{
	int i;
	bool halt, error;
	bool convertFormat;
	bool optimize;
	bool interactive;
	char * scriptFile;
	bool ignore;
	bool debug;		// Unpublished internal option: super-verbose, for debugging
	bool verbose;	
	bool silent;
	bool versionShown;

	debug = false;
	verbose = false;
	silent = false;
	versionShown = false;

	halt = false;
	error = false;
	convertFormat = false;
	optimize = false;
	scriptFile = NULL;
	ignore = false;
	bool nonInteractive = false;
	char * landmarkModelFile = NULL;
	interactive = false;

	string imageFile;
	string modelFile;

	Registry scriptParms(20, 10);

	// -----------------------------  Start of options parsing -----------------------------

	for (i = 1; i < argc; i++) {
		string opt = argv[i];

		if (getarg("ci", "ignore", false, opt, i, error, verbose, argc, argv)) {
			ignore = true;
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

		if (getarg("ix", "simTransModel", true, opt, i, error, verbose, argc, argv)) {
			scriptParms.setStringValue("simTransModel", argv[i]);
			continue;
		}

		if (getarg("ih", "hist", true, opt, i, error, verbose, argc, argv)) {
			scriptParms.setStringValue("Hist", argv[i]);
			continue;
		}

		//if (getarg("dqf", NULL, true, opt, i, error, verbose, argc, argv)) {
		//	scriptParms.setStringValue("Dqf", argv[i]);
		//	continue;
		//}

		if (getarg("dqfConfig", NULL, true, opt, i, error, verbose, argc, argv)) {
			scriptParms.setStringValue("DqfConfig", argv[i]);
			continue;
		}

		if (getarg("dqfOut", NULL, true, opt, i, error, verbose, argc, argv)) {
			scriptParms.setStringValue("DqfOut", argv[i]);
			continue;
		}

		if (getarg("mpConfig", NULL, true, opt, i, error, verbose, argc, argv)) {
			scriptParms.setStringValue("MpConfig", argv[i]);
			continue;
		}

		if (getarg("mpOut", NULL, true, opt, i, error, verbose, argc, argv)) {
			scriptParms.setStringValue("MpOut", argv[i]);
			continue;
		}

		if (getarg("histConfig", NULL, true, opt, i, error, verbose, argc, argv)) {
			scriptParms.setStringValue("HistConfig", argv[i]);
			continue;
		}

		if (getarg("histOut", NULL, true, opt, i, error, verbose, argc, argv)) {
			scriptParms.setStringValue("HistOut", argv[i]);
			continue;
		}

		if (getarg("histPartitionObject", NULL, true, opt, i, error, verbose, argc, argv)) {
			scriptParms.setStringValue("HistPartitionObject", argv[i]);
			continue;
		}

		if (getarg("histPartitionData", NULL, true, opt, i, error, verbose, argc, argv)) {
			scriptParms.setStringValue("HistPartitionData", argv[i]);
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

			sprintf(key, "dqf%d", f);
			if (getarg(key, NULL, true, opt, i, error, verbose, argc, argv)) {
				sprintf(key, "Dqf%d", f);
				scriptParms.setStringValue(key, argv[i]);
				break;
			}
		}
		if (f < 20)
			continue;


		// The option is unpublished, for use at UNC only
		if (getarg("lD", "debug", false, opt, i, error, verbose, argc, argv)) {
			debug = true;
			continue;
		}

		if (getarg("lv", "verbose", false, opt, i, error, verbose, argc, argv)) {
			verbose = true;
			continue;
		}

		if (getarg("lq", "quiet", false, opt, i, error, verbose, argc, argv)) {
			silent = true;
			continue;
		}

		if (getarg("I", "interactive", false, opt, i, error, verbose, argc, argv)) {
			interactive = true;
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
							tuneVals.setString(w, newStr);
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
							tuneVals.setWeight(w, newWt);
							if (verbose)
								cout << "    " << key << " = " << newWt << '\n';
						}
					}
				}
			}
			continue;
		}

		if (getarg("il", "landmarkModel", true, opt, i, error, verbose, argc, argv)) {
			landmarkModelFile = argv[i];
			scriptParms.setStringValue("LandmarkModel", argv[i]);
			continue;
		}

		if (getarg("id", "distanceMap", true, opt, i, error, verbose, argc, argv)) {		// undocumented
			scriptParms.setStringValue("DistanceMap", argv[i]);
			continue;
		}

		if (getarg("om", "outModel", true, opt, i, error, verbose, argc, argv)) {
			scriptParms.setStringValue("OutModel", argv[i]);
			continue;
		}

		if (getarg("ot", "outTile", true, opt, i, error, verbose, argc, argv)) {
			scriptParms.setStringValue("OutTile", argv[i]);
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

		if (getarg("pd", NULL, false, opt, i, error, verbose, argc, argv)) {
			tuneVals.printDesc();
			halt = true;
			continue;
		}

		if (getarg("ph", NULL, true, opt, i, error, verbose, argc, argv)) {
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
			P3DControl control(0);

			imageFile = argv[i];
			if (! control.loadImage(imageFile.data(), false, true))		// load header only for speed
				return 0;

			// Print like print_image_info() but format like a registry
			Image3D * image = control.getImagePtr();
			cout << "image = " << imageFile << ";" << endl;
			cout << "# voxel counts" << endl;
			cout << "xDim = " << image->getXDim() << ";" << endl;
			cout << "yDim = " << image->getYDim() << ";" << endl;
			cout << "zDim = " << image->getZDim() << ";" << endl;

			cout << "# voxel center to voxel center distance (cm)" << endl;
			cout << "xSpacing = " << image->getXSpacing() << ";" << endl;
			cout << "ySpacing = " << image->getYSpacing() << ";" << endl;
			cout << "zSpacing = " << image->getZSpacing() << ";" << endl;

			Vector3D origin = image->getWorldOrigin();
			cout << "# center of first voxel (cm)" << endl;
			cout << "xOrigin = " << origin.getX() << ";" << endl;
			cout << "yOrigin = " << origin.getY() << ";" << endl;
			cout << "zOrigin = " << origin.getZ() << ";" << endl;

			cout << "# center of last voxel (cm)" << endl;
			Vector3D opposite = origin + Vector3D(
				(image->getXDim() - 1) * image->getXSpacing(),
				(image->getYDim() - 1) * image->getYSpacing(),
				(image->getZDim() - 1) * image->getZSpacing());
			cout << "xOpposite = " << opposite.getX() << ";" << endl;
			cout << "yOpposite = " << opposite.getY() << ";" << endl;
			cout << "zOpposite = " << opposite.getZ() << ";" << endl;

			cout << "# distance from first voxel to last voxel (cm)" << endl;
			Vector3D extent = opposite - origin;
			cout << "xExtent = " << extent.getX() << ";" << endl;
			cout << "yExtent = " << extent.getY() << ";" << endl;
			cout << "zExtent = " << extent.getZ() << ";" << endl;

			cout << "# max distance in any dimension (cm)" << endl;
			double mex = extent.getX();
			if (extent.getY() > mex)
				mex = extent.getY();
			if (extent.getZ() > mex)
				mex = extent.getZ();
			cout << "maxExtent = " << mex << ";" << endl;

			continue;
		}

		if (getarg("pm", "extent", true, opt, i, error, verbose, argc, argv)) {
			double extent;
			double maxExtent;

			maxExtent = 0.0;
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
			P3DControl control(0);

			// Could consume up to the next switch...
			int iCount = 0;		// image count
			while (i < argc) {
				imageFile = argv[i];
				extent = control.calcMaxExtent(imageFile.c_str());
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

		if (getarg("po", NULL, false, opt, i, error, verbose, argc, argv)) {
			optimize = true;
			continue;
		}

		if (getarg("pr", NULL, false, opt, i, error, verbose, argc, argv)) {
			tuneVals.printAsRegistry();
			halt = true;
			continue;
		}

		if (getarg("pt", NULL, false, opt, i, error, verbose, argc, argv)) {
			// Undocumented; for MIDAG internal use only
			// make sure that this header line and printInputsAsTDTable
			// and printAsTDTable all print the same number of columns
			cout << "Registry Name and Default (case independent)"
				 << "\t"<< "Switch Alias"
				 << "\t"<< "Description"
				 << endl;
			printInputsAsTDTable();
			tuneVals.printAsTDTable();
			halt = true;
			continue;
		}

		// Any tuning value can be a switch (case independent)
		if (opt[0] == '-')
		{
			int nparm;

			nparm = (int) tuneVals.indexByLongName(opt.data() + 1);
			if (nparm < 0) {
				cout << "Error: Unknown option: " << opt << "\nTry using -help" << endl;
				error = true;
				continue;
			}
			if (++i >= argc || argv[i][0] == '-')	// Allow paths beginning with '/'
			{
				cout << "Error: Option " << opt << " requires an argument" << endl;
				error = true;
				return false;
			}
			else
				if (verbose)
					cout << opt << " " << argv[i] << endl;

			const char * tag = opt.data() + 1;	// Tuning parameter's name
			if (! tuningLocked(nparm)) {
				// Silently ignore changes to locked values
				double lastWt = tuningWt(nparm);
				tuneVals.setWeight(nparm, atof(argv[i]));

				// Print when non-default values are specified
				if (verbose && lastWt != tuningWt(nparm))
					cout << "  switch changes " << tuningLongName(nparm)
						<< " = " << tuningWt(nparm) << endl;
			}
			if (nparm == ScriptVersion) {
				// Verify value is compatible; see top of file
				double proposedVersion = atof(argv[i]);
				if ((int) proposedVersion != (int) tuningWt(nparm)) {
					cout << "  switch sets " << tuningLongName(nparm)			// AGG: This message doesn't make sense
						<< " = " << tuningWt(nparm)
						<< " -- Warning: expecting "
						<< tuningWt(nparm) << ".xx ... proceed with caution" << endl;
				}
			}
			continue;
		}

		// didn't understand an option
		cout << "Error: Unknown option: " << opt << "\nTry using -help" << endl;
		error = true;
	}

	// -----------------------------  End of options parsing -----------------------------

	// Here's where the normal processing begins, so introduce me
	if (! silent && ! versionShown)
		cout << headerString << endl;


	if (interactive && ! optimize) {
		cout << "Option -I requires using option -po" << endl;
		error = true;
	}


	if (halt || error) {
		if (halt)
			return 0;
		else {
			cout << "Use -h for help" << endl;
			return -1;
		}
	}

	// Set the global controls, including the default for OutputVerbosity
	if (optimize && ! interactive)	// See if this is a scripted optimization run
	{
		bool ret;

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
		globalControl->setDefault(ScaleImages, true);
		globalControl->setDefault(SimTransformSeparate, false);
		globalControl->setDefault(SimTransformMatrix, false);

		// undoLen not significant for scripting, so set it to 1
		P3DControl control(0);

#ifdef BINARY
		ret = control.runBinaryPablo(scriptParms);
#else
		ret = control.runPablo(scriptParms);
#endif
		delete globalControl;
		return (ret ? 0 : 1);
	}


	return 0;
}


