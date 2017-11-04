/******************************************************************
 * PAUL'S Libraries                                               *
 ******************************************************************
 * Author:						Paul Yushkevich
 *
 * Date:							Apr 1, 1999
 *
 * Description					All my libraries
 *									
 ******************************************************************
 * libs.h
 *	--------
 * Included my all of my libraries and contains namespace declarations
 * and such
 ******************************************************************/
#ifndef PAULY_LIBS
#define PAULY_LIBS


#ifndef NO_NAMESPACE

namespace pauly {

#define NAMESPACE_PAULY_START namespace pauly {
#define NAMESPACE_PAULY_END };
#define USING_NAMESPACE_PAULY using namespace pauly;
#define NAMESPACE_PAULY pauly

}

#else

#define NAMESPACE_PAULY_START ;
#define NAMESPACE_PAULY_END ;
#define USING_NAMESPACE_PAULY ;
#define NAMESPACE_PAULY 

#endif

// Go ahead and let everyone use the namespace.  In case of some conflict, they may use
// NAMESPACE_PAULY:: to prefix the conflicting name
USING_NAMESPACE_PAULY

#endif

