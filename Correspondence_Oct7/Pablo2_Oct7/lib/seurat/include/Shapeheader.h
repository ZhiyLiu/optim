/********************************************************************************/
/*										*/
/*  	File	:  Shapeheader.H						*/
/*										*/
/*	Description:  Type declarations, constants, and standard includes	*/
/*			for Shapemonger classes.				*/
/*										*/
/*	Project :  Shapemonger							*/
/*										*/
/*	Author  :  A. Thall							*/
/*										*/
/*	Date	:  4. November 1996						*/
/*										*/
/*	Modifications:  															*/
/*		30. Oct 00 -- added MIN() and MAX() inlines for ints and doubles		*/
/*		22. Oct 01 -- added INTPOW()											*/
/*																				*/
/********************************************************************************/

#define myself (*this)

// declare constants
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.7182818284590452354
#endif

// declare OFF and ON, for toggles, etc.
#define OFF 0
#define ON 1

// added W to worldaxes----will be used in Quat.H for componenents
typedef enum {X, Y, Z, W} worldaxes;

// All rasterized objects will have distinct names so they can be toggled
//   visible/invisible or deleted entirely from the display list.
#define NULLNAME 0
typedef int objNAME;

typedef enum {ORDINARY, EXTRAORDINARY} ordinarity;

/************************************************************************/
/* INTMOD() -- dumb function to return modular remainder (x MOD modulus)*/
/*      ---use for small loop-counter wraparound functions ONLY.        */
/************************************************************************/
inline int INTMOD(int x, int modulus)
{
    while (x < 0)
        x += modulus;
    while (x >= modulus)
        x -= modulus;
    return x;
}

/********************************************************************************/
/* INTPOW() -- dumb function to raise small integer to small integral power		*/
/*          without checking for overflow or being particularly efficient.		*/
/********************************************************************************/
inline int INTPOW(int x, int power)
{
	int retval = 1;
	while (power-- > 0)
		retval *= x;

	return retval;
}

/********************************************************************************/
/* MIN() and MAX() -- inlines to compare two ints or doubles					*/
/********************************************************************************/
inline int MIN(int a, int b)
{
	return (a < b ? a : b);
}
inline double MIN(double a, double b)
{
	return (a < b ? a : b);
}
inline int MAX(int a, int b)
{
	return (a > b ? a : b);
}
inline double MAX(double a, double b)
{
	return (a > b ? a : b);
}

