#ifndef _GPFLIB_TIMER_H_
#define _GPFLIB_TIMER_H_

/*
        UNC/3DCV Image Processing Software Packages


Copyright (c) 1990, 1993
                   3D Computer Vision Research Group,
                   Utrecht University in the Netherlands,
                   Utrecht, The Netherlands.

Copyright (c) 1993, 1996, 1997
                   Computer Science Department,
                   University of North Carolina at Chapel Hill,
                   Chapel Hill, North Carolina, USA.


Permission is given to distribute these files, as long as the copyright
messages are not removed, the file named COPYRIGHT (this file) is
retained unaltered in all distributions and no fees are charged for the
files.  A fee may be charged for reproduction and shipment of the files.

This software is provided strictly on an "as is" basis without warranty
of any kind.  Neither the 3D Computer Vision Research Group nor the
Computer Science Department of the University of North Carolina at
Chapel Hill nor anyone else who has been involved in the creation,
production or delivery of this software shall be liable for any direct,
indirect, consequential or incidental damages arising out of the use or
inability to use this software even if the 3D Computer Vision Research
Group or the Computer Science Department of the University of North
Carolina at Chapel Hill have been advised of the possibility of such
damages.

For a complete statement of the conditions pertaining to the ownership,
distribution and use of this software, refer to the Software License
Agreement provided with the UNC/3DCV Image Processing Software Packages.
If you do not have a copy, you may obtain one from the persons listed
below or from the institutions named above.  By having, retaining or
using a copy of this software, you agree to be subject to the conditions
of the Software License Agreement.

*/


/*  General Timer class.  This class can measures time in seconds in a application.

	Usage:
		The contructor starts the timing automatically.

		double start()
		    Restart the timing and return the current (starting) time in seconds.

		double stop()
		    Stop the timing and return the measured time in seconds.

		double elapsed() const
		    Returns the current time minus the starting time in seconds.  If it
			is running, the timer continues running.

		double total()
		    Returns the total time of the timer in seconds.  If the timer is
			running it is stopped.  Function elapsed() is similar, but does
			not stop the timer.

		double current_time() const
		    Returns the current time of the timer in seconds.  If it is running,
		    the timer continues running.

		void process_time(char * prefix, char * suffix)
		    Generates a message on the standard output, giving the total time of
			the timer in the following format:
			    14h 23m 30.543245s.

		If the timer is running it will continue to run.  The optional strings
		will be inserted in the output before and after the time.  If they are
		not provided, no prefix will be output and the suffix will be a newline.


    IMPLEMENTATION DETAILS:

    bool is_running
    If this variable contains true we are measuring.

    double start_time
    The start time of the measurement.  This variable counts
    seconds in a highly system dependent manner.  The difference
    between start_time and stop_time is system independent.

    double stop_time
    The stop time of the measurement.  This variable counts
    seconds in a highly system dependent manner.  The difference
    between start_time and stop_time is system independent.

    tms * buffer
    Buffer with system time.

*/


class IntervalTimer {

public:
	IntervalTimer();
	~IntervalTimer();

	double start();
	double stop();

	double elapsed() const;
	double total();
	double current_time() const;
	void process_time(char * prefix = NULL, char * suffix = "\n");
	bool running() const;

	friend std::ostream & operator << (std::ostream & ostr,
		const IntervalTimer & timer);


private:

	void * buffer;		// Buffer with system time.
	double start_time;	// Start time of measurement
	double stop_time;	// Stop time of measurement
	bool is_running;	// Running flag.
};


bool inline IntervalTimer::running() const {
	return is_running;
}



#endif	/* _GPFLIB_TIMER_H_ */

