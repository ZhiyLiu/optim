
#include <iostream>
#include <assert.h>
#if defined(_MSC_VER) || defined(__BORLANDC__)
#include <sstream>
#include <time.h>
#else
#include <sstream>
#include <sys/types.h>
#include <unistd.h> 	/* To get the definition of _SC_CLK_TCK */
#include <stdlib.h>
#include <sys/times.h>
#endif

#include "IntervalTimer.h"


static double clock_counts;

using namespace std;


IntervalTimer::IntervalTimer() {	// Constructor
#if defined(_MSC_VER) || defined(__BORLANDC__)
	buffer = (void *) new double;
	clock_counts = CLOCKS_PER_SEC;
#else
	buffer = (void *) new tms;
	clock_counts = sysconf(_SC_CLK_TCK);
#endif
	(void) this->start();
	stop_time = start_time;
}


IntervalTimer::~IntervalTimer() {	// destructor
#if defined(_MSC_VER) || defined(__BORLANDC__)
	delete (double *) buffer;
#else
	delete (tms *) buffer;
#endif
}


double IntervalTimer::start() {
#if defined(_MSC_VER) || defined(__BORLANDC__)
	*((double *)buffer) = clock();
#else
	times((tms *)buffer);
#endif
	assert(clock_counts != 0.0);
#if defined(_MSC_VER) || defined(__BORLANDC__)
	start_time = *((double *) buffer)/clock_counts;
#else
	start_time = double(((tms *) buffer)->tms_utime)/clock_counts;
#endif
	is_running = true;
	return start_time;
}


double IntervalTimer::stop() {
#if defined(_MSC_VER) || defined(__BORLANDC__)
	*((double *)buffer) = clock();
#else
	times((tms *)buffer);
#endif
	assert(clock_counts != 0.0);
#if defined(_MSC_VER) || defined(__BORLANDC__)
	stop_time = *((double *) buffer)/clock_counts;
#else
	stop_time = double(((tms *) buffer)->tms_utime)/clock_counts;
#endif
	is_running = false;
	return stop_time;
}


double IntervalTimer::elapsed() const {
#if defined(_MSC_VER) || defined(__BORLANDC__)
	*((double *)buffer) = clock();
#else
	times((tms *)buffer);
#endif
	assert(clock_counts != 0.0);
	if (is_running)
#if defined(_MSC_VER) || defined(__BORLANDC__)
		return *((double *) buffer)/clock_counts - start_time;
#else
		return double(((tms *) buffer)->tms_utime)/clock_counts
			- start_time;
#endif
	else
		return stop_time - start_time;
}


double IntervalTimer::current_time() const {
#if defined(_MSC_VER) || defined(__BORLANDC__)
	*((double *)buffer) = clock();
#else
	times((tms *)buffer);
#endif
	assert(clock_counts != 0.0);
#if defined(_MSC_VER) || defined(__BORLANDC__)
	return *((double *) buffer)/clock_counts;
#else
	return double(((tms *) buffer)->tms_utime)/clock_counts;
#endif
}


double IntervalTimer::total() {		// get total time
	if (is_running) (void) stop();
	return stop_time - start_time;
}


// Dump processed time to the global debugbuffer Debug as 
// an application message.
void IntervalTimer::process_time(char * prefix, char * suffix) {
	ostringstream message;
	double tmp;


	tmp = elapsed();
	int iHours = int(tmp/3600);
	tmp -= 3600 * iHours;
	int iMinutes = int(tmp/60);
	tmp -= 60 * iMinutes;
	// tmp now contains the seconds

	if (iHours < 10)
		message << "0";
	message << iHours << "h ";

	if (iMinutes < 10)
		message << "0";
	message << iMinutes << "m ";

	if (tmp < 10)
		message << "0";
	message << tmp << "s";
	message << ends;

	const char *str = message.str().c_str();	// Convert ostrstream to string
	if (prefix != NULL)
		std::cout << prefix;		// Send prefix
	std::cout << str;		// Send string to Debugger
	if (suffix != NULL)
		std::cout << suffix;		// Send suffix
//	delete str;

	return;
}


// Dump processed time of a IntervalTimer to an ostream
std::ostream & operator << (std::ostream & ostr, const IntervalTimer & timer) {

	double tmp;


	tmp = timer.elapsed();
	int iHours = int(tmp/3600);
	tmp -= 3600 * iHours;
	int iMinutes = int(tmp/60);
	tmp -= 60 * iMinutes;
	// tmp now contains the seconds

	if (iHours < 10)
		ostr << "0";
	ostr << iHours << "h ";

	if (iMinutes < 10)
		ostr << "0";
	ostr << iMinutes << "m ";

	if (tmp < 10)
		ostr << "0";
	ostr << tmp << "s";

	return ostr;
}


