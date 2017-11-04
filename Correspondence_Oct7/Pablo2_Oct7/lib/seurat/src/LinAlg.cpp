
#include <iostream>
#define D_LINALG
#include "Shapedepend.h"

using namespace std;
using namespace ThallCode;


/* NOTICE:

	The error reports in this file exist for very good reasons.
	If you get one of them, then the way to fix your problem is
	to modify your code, not to comment out these error reports.

	In short:  if you do not own this file, do not touch it!
*/

static int errorCount2 = 10;
static int errorCount3 = 10;

void DbVector2::error(bool errorNumber) const
{
	if (--errorCount2 > 0) {
		std::cerr << "Error in DbVector2: ";
		if (errorNumber)
			std::cerr << "self";
		std::cerr << "normalize() of vector with length 0" << endl;
	}
}

void DbVector3::error(bool errorNumber) const
{
	if (--errorCount3 > 0)
		std::cerr << "Fatal error in DbVector3; errorNumber = "
			<< errorNumber << endl;
}

void DbVector3::printvals(char *message)
{
    if (message != NULL)
		std::cerr << message << '\n';
    std::cerr << '[' << x() << ", " << y() << ", " << z() << "]\n";
}


