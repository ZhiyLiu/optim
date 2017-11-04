/********************************************************************************/
/*																				*/
/*  	File	:  Samplestats.h												*/
/*																				*/
/*	Description:  statistical analysis class---stores an array of n samples,	*/
/*				  returning sample mean, variance, st.dev. and higher moments.	*/
/*																				*/
/*	Project :  Seurat															*/
/*																				*/
/*	Author  :  A. Thall															*/
/*																				*/
/*	Date	:  16. April 2002													*/
/********************************************************************************/

// Having some stupid VC++ problems with include files
//double mysqrt(double thisnum)
//{

class Samplestats
{
    double *samples;
	long sampsize;

	long samploc;

	double mean, samplevariance, samplestddev;
	bool meancomputed;
	bool varcomputed;
	bool stddevcomputed;

public:
    Samplestats();
	Samplestats(long ssize);
	~Samplestats();

	void init(long ssize);

	// This class doesn't check for overrunning end of sample-array, so don't screw up
	void addsample(double sampval) { samples[samploc++] = sampval; }

	void clearstats() { meancomputed = false; varcomputed = false; stddevcomputed = false; }
	void reinit() { samploc = 0; clearstats(); }

	// Get each of the following statistical measures
	//   the variance computed is the sample variance (Sum[(X_i - X_bar)^2]/(n - 1)),
	//   and the stddev is the sqrt. of the sample variance
	double getmean();
	double getvariance();
	double getstddev();

	// clear any previous and call the three routines above
	void computestats(bool printout = false);

    void printvals(char *message = NULL);
};

/*
 * Constructors
 */
inline Samplestats::Samplestats()
{
	samples = NULL;
	clearstats();
}

inline Samplestats::Samplestats(long ssize)
{
	sampsize = ssize;
	samples = new double[sampsize];
	samploc = 0;
	clearstats();
}
/*
 * Destructor
 */
inline Samplestats::~Samplestats()
{
	if (samples != NULL)
		delete []samples;
}

/*
 * init() -- clear old value
 */
inline void Samplestats::init(long ssize)
{
	sampsize = ssize;
	if (samples != NULL) 
		delete []samples;
	samples = new double[sampsize];
	samploc = 0;
	clearstats();
}

/*
 * computestats() -- compute mean, variance and stddev based on input array
 */
inline void Samplestats::computestats(bool printout)
{
	clearstats();
	getmean();
	getvariance();
	getstddev();

	if (printout)
		printvals("Output statistics are:");
}

/*
 * printvals()
 */
inline void Samplestats::printvals(char *message)
{
    if (message != NULL)
		std::cerr << message << '\n';
	if (samples == NULL) 
		std::cerr << "non-initialized Samplestats class\n";
	else {
		std::cerr << "sample size = " << samploc << " out of space = " << sampsize << ".\n";
		if (meancomputed) {
			std::cerr << "   mean = " << mean << '\n';
			if (varcomputed) {
				std::cerr << "   sampvar = " << samplevariance << '\n';
				if (stddevcomputed)
					std::cerr << "   stddev =" << samplestddev << '\n';
			}
		}
	}
}

