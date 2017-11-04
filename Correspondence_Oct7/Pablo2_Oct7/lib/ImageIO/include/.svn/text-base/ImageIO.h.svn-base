#ifndef IMAGE_IO_H
#define IMAGE_IO_H


#ifndef MAXPATHLEN 
#define MAXPATHLEN 256
#endif

using std::string;



typedef unsigned short u_short;
typedef unsigned char  u_char;

struct ImageStruct;

class ImageIO
{
	public:

		enum ImageType { raw3, analyze, gipl, meta, planIM, dicom, unknown};


		// Guess the format of an image given the path
		ImageType guessImageFormat(string & filename);

		// Load an image guessing its format
		void loadThisImage(string & filename, ImageStruct & imStr,
			bool headerOnly = false, ImageType extension = unknown)
			throw (BasicException);

		// Save an image guessing its format
		void saveImage(string & filename, ImageStruct & imStr);

	protected:

		// Analyze Images
		void loadAnalyze(string & filename, ImageStruct & imStr, bool headerOnly = false);
		void saveAnalyze(string & filename, ImageStruct & imStr);
		void readAnalyzeHeader(const char * hdrname, ImageStruct & imStr);

		// Meta Images
		void loadMeta(string & filename, ImageStruct & imStr, bool headerOnly = false);
		void saveMeta(string & filename, ImageStruct & imStr);

		// Gipl Images
		void loadGipl(string & filename, ImageStruct & imStr, bool headerOnly = false);
		void saveGipl(string & filename, ImageStruct & imStr);

		// PLUNC Images
		void loadPlanim(string & filename, ImageStruct & imStr, bool headerOnly = false);
		void savePlanim(string & filename, ImageStruct & imStr);

};



#endif

