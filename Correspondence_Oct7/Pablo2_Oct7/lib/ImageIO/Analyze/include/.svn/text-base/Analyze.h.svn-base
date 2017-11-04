#ifndef ANALYZE_H
#define ANALYZE_H


using std::string;


class Analyze
{
	public:

		Analyze() {}
		~Analyze() {}

		// Load Analyze image files
		bool load(string filename, ImageStruct & image3D);


		// Write image into Analyze files
		bool save(string filename, ImageStruct & image3D);

	protected :

		void constructAnalyzeNames(const char * filename, char * database,
			char * hdr_path, char * img_path);

		Analyze & operator=(const Analyze &) const;	// Not implemented
		Analyze(const Analyze &);	// Not implemented
};



#endif

