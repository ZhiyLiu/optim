/** @file gen.h
    This header file contains data type and record declarations of
    general use in the radiation treatment planning (RTP) system.
*/

// HISTORY
// 21-july-2004 zsx: Add doxygen-style comments.

#ifndef _PLUNC_HAVE_GEN_H
#define _PLUNC_HAVE_GEN_H

#ifdef PLAN_WINNT
#include "libmisc.h"
#endif


#ifndef VAR_GREYBAR_ENTRIES
#define GREYBAR_ENTRIES	(100)
#endif

/** This 4-digit number is to be displayed in the form of x.x.xx as
PLUNC version.

    This version number ties the struct's members together so that
    all that needs to be known is this number.

    When any struct in src/include changes, the PLUNC version number
    MUST be changed as well. It is meant to be compared using < or >
    to hardcoded values, as in #if PLUNC_VERSION > 6600
*/
#define PLUNC_VERSION 6700


#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif
#ifndef PI
#define PI (3.141592654)
#endif

/** @name Fixed point math scaling values

Fixed point numbers store real
numbers in an integer using a fixed number of bits (SBITS) for the
fraction. You convert a float into a fixed point representation by
multiplying the float by SCALE and then taking the integer part.
You convert a fixed point number back to a float by multiplying the
fixed point number by INV_SCALE. You can get the fraction bits by
ANDing a fixed point number with MASK. You can get the integer
portion of the fixed point number by shifting it right by SBITS.
*/
//@{
#define SBITS		14
#define SCALE		16384
#define INV_SCALE	6.1035156e-5
#define MASK		0x3fff
//@}

/**  @name Some convenience macros

Macros for the minimum, maximum, integer absolute value,
and floating-point absolute value.
*/
//@{
#define PLAN_MIN(a,b)		(((a) < (b)) ? (a) : (b))
#define PLAN_MAX(a,b)		(((a) > (b)) ? (a) : (b))
#define PLAN_ABS(a)		(((a) < 0) ? -(a) : (a))
#define PLAN_FABS(a)		(((a) < 0.0) ? -(a) : (a))

/** TRUE when the difference between \a a and \a b is
smaller than or equal to \a t. Do not use if the values of
\a a or \a b or \a t change with each invocation of the macro.
*/
#define PLAN_CLOSE_TO(a,b,t)	(PLAN_ABS((a) - (b)) <= (t))

/** Return the number of items in an array. */
#define PLAN_NITEMS(array)      (sizeof(array) / sizeof(*(array)))
//@}

/** Maximum number of intersections between a line and a polygon. */
#define FI_LIMIT		(40)

#define TABLET			(0)  //!< Not currently used.
#define KEYBOARD		(1)  //!< Not currently used.

/** Defines the offset between the value ranges of [-1024, 3072] and
[0, 4096]. */
#define CT_TO_ENCODED (1024)

/** The maximum number of characters for the name string of anything.
*/
#define NAME_LENGTH		(100)

typedef int Boolean;
#ifndef Bool
#define Bool int
#endif

/* Solaris 2.4 would not work with "signed char". It was changed
   in Solaris 2.6. */
#ifdef PLAN_SOLARIS2_4
typedef char SBYTE;
#else
typedef signed char SBYTE;
#endif
typedef unsigned char UBYTE;

#ifndef _PLUNC_HAVE_LIB3D_H

/** PNT3D_struct is a data structure for a 3D-point
    whose coordinates are floats in centimeters (cm). The coordinate
    system is defined by the context in which the PNT3D is used.
*/
typedef struct PNT3D_struct
{
    float x; //!< x-coordinate of the 3D-point in cm.
    float y; //!< y-coordinate of the 3D-point in cm.
    float z; //!< z-coordinate of the 3D-point in cm.
} PNT3D;

/** Data type for a 4x4 matrix that transforms 3D-points
    between two coordinate systems.

    This is a row-vector -- the translation (x,y,z) elements
    are in [3][0], [3][1], and [3][2].

    PLUNC utilizes a number of coordinate systems (CS):
    namely Image CS, Grid CS, Table-Top CS (also as Patient CS),
    Beam CS (left-handed), Unit CS, Filter CS,
    registered Subimage CS, and various Window/view CSs.
*/
typedef float LIB3D_TRANSFORM[4][4];
#endif

/** SPNT_struct is a data structure for a 3D-point whose
    coordinates are of data type SBYTE.
*/
typedef struct SPNT_struct
{
    SBYTE x; //!< x-coordinate of the 3D-point.
    SBYTE y; //!< y-coordinate of the 3D-point.
    SBYTE z; //!< z-coordinate of the 3D-point.
} SPNT;

/** T_STRIP_struct (Triangle-STRIP) records a sequence of vertices of two
contours in a specially-arranged order for the purpose of constructing a
strip of connected triangles to describe the surface formed between
these two contours.

- Since only two edges are needed to define a triangle, the two edges
whose 1st end-point is a vertex of the 1st contour and 2nd end-point
is a vertex of the 2nd contour are recorded in the *pt array. Also
note that one of the two end-points is reused for the 2nd edge.

- Since the next connected triangle shares its first edge with the 2nd
edge of the triangle before it, this shared edge is only recorded once
in the *pt array. */
typedef struct T_STRIP_struct
{
    /** Total number of vertices or end-points of triangle edges in
    this T-strip. */
    int		count;
    int		con_num1; //!< Index of the 1st contour.
    int		con_num2; //!< Index of the 2nd contour.
    /** The indices of the points in *pt. */
    int		*vert;
    /** The vertices from the two contours that serve as the end-points
    of the edges of the triangles in the T-strip. */
    PNT3D	*pt;
    /** Same as *pt for 3D planning, or deformed version of *pt for 4D
    planning. */
    PNT3D	*dpt;
    /** The normal vector of the vertices in *pt. */
    PNT3D	*norm;
} T_STRIP;

/** T_STRIPS_struct is a data structure containing all T_STRIPs
of an anstruct.
*/
typedef struct T_STRIPS_struct
{
    int		count; //!< Number of T_STRIPs in the T_STRIPS.
    T_STRIP	*strip; //!< Pointer to the individual T_STRIPs.
} T_STRIPS;

/** Data type for color in its RGB values and alpha value.
    Element 0-2 for Red, Green, and Blue values respectively.
    Element 3 is the alpha value for specifying the transparency.
*/
typedef float COLOR4[4];

/** PDATE_struct is a data structure for date in 4 integers for
    the day, month, year and day-of-the-week.
*/
typedef struct PDATE_struct
{
    int day;   //!< one based: 1 = first day of month
    int month; //!< one based: 1 = Jan
    int year;  //!< year-1900 if value < 100; else 4-digit year
    int dow;   //!< zero based: 0 = Sunday
} PDATE;

/** PTIME_struct is a data structure for time in 3 integers for
    the hour, minute, and second.
*/
typedef struct PTIME_struct
{
    int hour; /* 24 hour clock */
    int minute;
    int second;
} PTIME;

/** WEIGHT_struct is a data structure to record the weight
components of an external beam or a brachytherapy object.

Weights can be products of several components -
for instance a brachytherapy object can be weighted
by a time and an activity, and, maybe, a decay factor
whereas a beam weight is a beam renormalization factor
times a relative weight factor, etc.
*/
typedef struct WEIGHT_struct
{
    int factors; //!< Number of weight components - brachy=2 beam=4.
    union
    {
	float wt[5]; //!< 5-element array for the generic weight components.
	struct
	{
	    float time; //!< Time duration that an object is present.
	    float activity; //!< The activity of the object.
	    float decay; //!< The decay constant of the object.
	} brachy; //!< Weight component structure for a brachytherapy object.
	struct
	{
	    /** The # of MUs to deliver 1 cGy to the norm point, a parameter
	    for renormalizing the beam. */
	    float renorm;
	    /** The relative weight of this beam defined by the user in the
	    Spreadsheet dialog. */
	    float relative;
	    /** The prescribed dose at norm point per fraction per relative
	    weight - a multiplier to embody target dose, target percent,
	    sum of the beam relative weights in the norm group, and the # of
	    dose fractions.  */
	    float plan_norm;
	    /** The # of dose fractions. */
	    float fractions;
	} tele; //!< Weight component structure for an external beam.
    }f; //!< Union to hold any one of the three weight component types.
} WEIGHT;

/** WEIGHTS_struct is a data structure for multiple WEIGTHs in
    an array.
    A weight file consists of a WEIGHTS structure followed
    by an array of WEIGHT.
*/
typedef struct WEIGHTS_struct
{
    int count; //!< Number of WEIGHTs.
    WEIGHT *weights; //!< Array of WEIGHTs.
} WEIGHTS;

/** TWOD_GRID_struct is a data structure for 2D dose grid on a cm plane.
    The grid is located in the GRID coordinate system. It is used to
    record dose on slice planes.
    A 2D dose grid file is prefaced by a TWOD_GRID and includes:\n
    float dose_matrix[y_count][x_count];\n
    Boolean inside[y_count][x_count];\n
    The units of the dose value matrix are defined by context.
*/
typedef struct TWOD_GRID_struct
{
    float x_start; //!< x-coordinate of the first point.
    float y_start; //!< y-coordinate of the first point.
    float x_inc;   //!< x-distance between neighboring 2D points.
    float y_inc;   //!< y-distance between neighboring 2D points.
    int x_count;   //!< Total number of points in x-direction.
    int y_count;   //!< Total number of points in y-direction.
    float grid_to_pat_T[4][4]; /**< Matrix that transforms a
                                    GRID coordiante point to PATIENT
                                    coordinates. */
    float pat_to_grid_T[4][4]; /**< Inverse of grid_to_pat. */
    float max;                 //!< Maximum dose in this 2D grid.
    float min;                 //!< Minimum dose in this 2D grid.
    float *matrix;             //!< 2D dose values stored in a matrix.
    Boolean *inside;           //!< Not currently used.
} TWOD_GRID;

/** Macro for referencing dose value at location (i, j) of
    a 2D dose grid named GRID_ptr.
*/
#define TWOD_GRID_VALUE(GRID_ptr, i, j) ((GRID_ptr)->matrix[(i) +\
			((j) * (GRID_ptr)->x_count)])

/** GRID_struct is a data structure for 3D dose grid.
    The grid is located in the GRID coordinate system.
    A 3D dose grid file is prefaced by a GRID and includes:\n
    float dose_matrix[z_count][y_count][x_count];\n
    The units of the dose value matrix are defined by context.
*/
typedef struct GRID_struct
{
    PNT3D start; //!< The first point in the GRID.
    PNT3D inc;   //!< The increment between two neighboring 3D points.
    int x_count; //!< Total number of points in x-direction.
    int y_count; //!< Total number of points in y-direction.
    int z_count; //!< Total number of points in z-direction.
    float grid_to_pat_T[4][4]; /**< Matrix that transforms a
                                    GRID coordiante point to PATIENT
                                    coordinates. */
    float pat_to_grid_T[4][4]; /**< Inverse of grid_to_pat. */
    float max; //!< Maximum dose value in this 3D grid.
    float min; //!< Minimum does value in this 3D grid.
    float *matrix; //!< 3D dose values stored in an array.
} GRID;

/** Macro for referencing dose value at integer location (i, j, k) of
    a 3D dose grid named GRID_ptr.
*/
#define GRID_VALUE(GRID_ptr, i, j, k)\
    ((GRID_ptr)->matrix[(i) +\
			(GRID_ptr)->x_count *\
			 ((j) + ((k) * (GRID_ptr)->y_count))])

/** CALC_POINT_struct is a data structure for a point of interest
    in the patient.
*/
typedef struct CALC_POINT_struct
{
    char label[NAME_LENGTH]; //!< Label of the point.
    PNT3D p;                 //!< The point of interest.
    PNT3D def_p; //!< The deformed point (deformation under development).

} CALC_POINT;

/** CALC_POINTS_struct is a data structure for a series of CALC_POINTs
    and their corresponding dose values.

    The point file contains a CALC_POINTS followed by an array of
    CALC_POINT.\n

    The point dose file contains a CALC_POINTS followed by an array of
    CALC_POINT followed by the point doses (float dose[count]).
*/
typedef struct CALC_POINTS_struct
{
    int count;
    CALC_POINT *points;
    float *dose;
} CALC_POINTS, POINT_DOSES;

/** GRIDS_struct is a data structure containing a dose grid and a
file descriptor for the related dose grid file.
    GRIDS is used in both photon and electron dose code.
*/
typedef struct GRIDS_struct
{
    GRID	grid3; //!< Dose grid associated with the file descriptor.
    int		out_fdes; //!< Dose grid file descriptor.
} GRIDS;

/** C_POINTS_struct is a data structure for the file descriptor and
    associated CALC_POINTS.
    C_POINTS are used in both photon and electron dose code.
*/
typedef struct C_POINTS_struct
{
    CALC_POINTS points; //!< Points of interest.
    float *ssd;   //!< Corresponding SSDs of the points.
    float *batho; //!< Corresponding Batho coefficients of the points.
    int	out_fdes; //!< Points file descriptor.
} C_POINTS;

/** CONTOUR_struct is a data structure for a planar polygon delineating
    a cm area/volume of interest. A contour can have a Z thickness and is
    then called a prism (for example a contour on an image slice).

    A 2D contour file is prefaced by a CONTOUR and includes:

    float vertices[2][vertex_count];

    A stack of related CONTOURs form a volume of interest
    represented by data type ANASTRUCT.
*/
typedef struct CONTOUR_struct
{
    int vertex_count; //!< Number of vertices in this contour.
    int slice_number; //!< Index of the associated CT slice.
    PNT3D max;        /**< Maximum extent of prism
                           represented by this contour. */
    PNT3D min;        /**< Minimum extent of prism
                           represented by this contour. */
    float density;	  //!< gm/cc of tissue, filter, etc.
    float z;          //!< Nominal slice position.
    float *x; //!< Pointer to the x-positions of the vertices.
    float *y; //!< Pointer to the y-positions of the vertices.
    PNT3D *def_vert;  //!< For deformation in 4D-RTP, ignore it for now.
} CONTOUR;

/** ANASTRUCT_struct is a data structure for an ANAtomical STRUCTure
    that is referred as an anastruct in PLUNC.

    It corresponds roughly to a DICOM Structure Set.

    There is one ANASTRUCT per anatomical structure, although whether
    a pair of lungs is one structure or two is defined by context. The
    ANASTRUCT contains a set of CONTOURS (zero, one, or more contours
    per image slice). Overlapping anastructs (tumor vs. non-tumor),
    nested contours and negative spaces are all allowed.

    An anastruct file is prefaced by an ANASTRUCT and is the
    concatenation of contour_count CONTOUR files.
*/
typedef struct ANASTRUCT_struct
{
    char label[NAME_LENGTH]; //!< Anastruct name.
    int contour_count;       //!< Number of contours consisted.
    PNT3D max; //!< Maximum extent of all prisms in this anastruct.
    PNT3D min; //!< Minimum extent of all prisms in this anastruct.
    CONTOUR *contours;       //!< Pointer to the contours.
} ANASTRUCT;

/** PLAN_struct is a data structure that contains the orientation
of a view plane in terms of the patient-to-plan transformation matrix,
the isodose curves of interest in the view plane, and the reference
dose represented by 100% isodose curve of a treatment plan.

    A plan file consists of a PLAN_struct, followed by a vector of
    isodoses.

    Notice that PLAN_struct overloads the term "plan" that generally
    refers to a radiotherapy treatment plan, which consists of, in
    addition to the plan as defined here as PLAN_struct, the beams,
    norm points, beam weights, points of interest, dose grid,
    visibility of objects, etc.
*/
typedef struct PLAN_struct
{
    float pat_to_plan_T[4][4]; /**< If 2D, z=0.0 plane is
                               selected; if 3D, this is just a view.
                               */
    int isodose_count;         //!< Number of isodose values.
    float ref_dose;    //!< The isodose is % of this reference dose.
    float *isodoses;   //!< The isodose values of interest.
} PLAN;

/** THREED_MAP_struct is a data structure containing multiple 2D maps
    of data.

    Multiple unrelated data may occupy multiple maps in Z. The data is
    defined by context.
*/
typedef struct THREED_MAP_struct
{
   int		x_dim; //!< The map dimension in x-direction.
   float	x_start; //!< Start position in x-direction.
   float    x_end; //!< End position in x-direction.
   float	x_scale; //!< Step size in x-direction.
   float    inv_x_scale; //!< Inverse of x_scale.
   int		y_dim;  //!< The map dimension in y-direction.
   float	y_start; //!< Start position in y-direction.
   float    y_end; //!< End position in y-direction.
   float	y_scale; //!< Step size in y-direction.
   float    inv_y_scale; //!< Inverse of y_scale.
   int		z_dim; //!< The map dimension in z-direction.
   float	z_start; //!< Start position in z-direction.
   float    z_end; //!< End position in z-direction.
   float	z_scale; //!< Step size in z-direction.
   float    inv_z_scale; //!< Inverse of z_scale.
   int		slicecnt; //!< # of pixels in a slice.
   float	*data; //!< Linear array containing the data of the 3D-map.
   /* Only used for variable zdim map */
   /** The number of data members in each data group. */
   int		depth_inc;
   /** A 2D map that records the number of data groups for each
   pixel. */
   UBYTE	*cnt_map;
   /** A 2D map that records the maximum number of data groups that
   have ever been for each pixel. */
   UBYTE	*max_map;
   /** Linear array of pointers to all of the data groups of a variable
   zdim map. */
   float	**vdata;
} THREED_MAP;

/** NORM_struct is a data structure containing the
    dose normalization points of the treatment plan.
*/
typedef struct NORM_struct
{
    int		num;     //!< Number of norm groups.
    int		current; //!< Index of the current norm group.
    /** An array of the index to the plan's point list for the norm
    point of the norm groups. */
    int		*point;
    int		*points;   //!< Not currently used.
    /** The index of the first beam in the next norm group. */
    int		*beam;
    /** An array of Rx doses, one for each norm group at the isodose
    level given in the corrsponding *percent, while the Norm Pt of
    each group is at 100% level. */
    float	*dose;
    /** The number of fractions to deliver each of the above Rx doses.
    */
    float	*fraction;
    /** An array of isodose level, in %, at which each of the above
    doses is prescribed. */
    float	*percent;
    char	**label;   //!< Label of the norm points.
} NORM;

/** USER_DVH_struct is a data structure for user-defined DVH
    of a chosen anastruct that is used as the objective for
    dose optimization.
*/
typedef struct USER_DVH_struct
{
    int		target; //!< Index of the anastruct for this USER_DVH.
    int		type; /**< 1: if the anastruct is the treatment target;
                       0: if the anastruct is a critical structure
                       of the treatment target. */
    float	weight; //!< Weighting factor for this anastruct.
    int		count;  //!< Number of points in the DVH.
    float	*dose;  //!< Dose values of the points in the DVH.
    float	*vol;   //!< % volume of the points in the DVH.
    float	goal;   //!< % volume meeting the asked dose requirement.
    float	total;  //!< Not currently used.
    int		use_eud;//!< Use EUD based optimization.
    float	A_val;  //!< The exponent in the Generalized EUD equation.
} USER_DVH;

/** DVH_struct is a data structure for differential
    dose-volume histograms (DVH).

    The differential DVH is a plot of the volume receiving a dose
    within a specified dose bin as a function of dose. The other
    form of DVH is cumulative DVH, a plot of the volume that
    receives the indicated dose or higher as a function of dose.
*/
typedef struct DVH_struct
{
    int		num_entries; /**< Number of doses in a DVH plot.
                              Dose starts from 0. */
    float	bin_size;    //!< Size of the dose bins.
    float	*vol;        //!< Pointer to the volumes in the dose-bins.
    float	*weight;     //!< Optimization weighting factor.
} DVH;

/** TEMPLATE_struct is a data structure for beam templates.
*/
typedef struct TEMPLATE_struct
{
    char	name[NAME_LENGTH]; //!< Name of the beam template.
    int		numbeams;          //!< Number of beams in the template.
    PNT3D	*beam_vectors;     //!< Pointer to the beam vectors.
} TEMPLATE;

/** RGB_PIXEL_struct is a data structure for a pixel's color.
*/
typedef struct RGB_PIXEL_struct
{
    unsigned char	red;   //!< RED component of color.
    unsigned char	green; //!< GREEN component of color.
    unsigned char	blue;  //!< BLUE component of color.
} RGB_PIXEL;

/** MOSAIQ_DATA_struct is a data structure for MOSAIQ Tolerance
    Tables for PHOTON and ELECTRON machines.

    The values are stored in UP/site/mosaiq_data file.
*/
typedef struct MOSAIQ_DATA_struct
{
    int		tolerance_photon; //!< Tolerance for PHOTON machine.
    int		tolerance_electron; //!< Tolerance for ELECTRON machine.
} MOSAIQ_DATA;

/* Version constant of the deformation methodology in 4D-RTP. */
#define DEFORMATION_VERSION     4
#define DEFORMATION_SINGLE      1
#define DEFORMATION_MULTI       2

typedef struct DEFORMATION_struct
{
    int		version;
    int		num_x, num_y, num_z, num_t;
    int		type;
    int		roi_flag;
    int		dim_x, dim_y, dim_z;
    int		off_x, off_y, off_z;
    int		*tx, *ty, *tz;
    int		compressed;
    float	scale;
    PNT3D	*v;
    SPNT	*sv;
} DEFORMATION;

#endif
