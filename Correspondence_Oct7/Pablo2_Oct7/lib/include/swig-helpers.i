/* vim: set syn=cpp:
 * 
 * A swig helper routines files, put in some helper routines,
 * global renames and such here.
 *
 * Put header file specific options in the header itself. Do not
 * pollute this file.
 *
 * - Rohit Saboo
 */

%rename (display) print;
%rename (__add__) operator +;
%rename (change) operator =;
%rename (__sub__) operator -;
%rename (__neg__) operator -();
%rename (__getitem__) operator [];
%rename (__str__) operator char *;
%rename (__call__) operator ();
%rename (__mul__) operator *;
%rename (__div__) operator /;
%rename (__eq__) operator ==;
%rename (__ne__) operator !=;

%feature("autodoc","1");

%define TYPECAST(A, B)
%inline %{
const A* const_ ## A ## _ ## B( const B* b ) {
	return dynamic_cast<const A*>(b);
}
A* A ## _ ## B( B* b ) {
	return dynamic_cast<A*>(b);
}
%}
%enddef

%{

// These two global variables are declared in ControlParms.h
extern ControlParms * globalControl;	// Read the user's preferences file
extern int globalVerbosity;			// Current verbosity level of Pablo

%}

