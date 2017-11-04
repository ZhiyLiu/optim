/** @file libplanio.h
    This header file contains function prototypes for data input
    and output provided in the library \b libplanio.

    XDR is an architecture-neutral file format that is independent of
    CPU, memory byte ordering, OS and vendor. As each new PLUNC data
    structure is defined, it's I/O routines are defined in parallel by
    defining it's XDR primitive.

    We use the word \b record below to indicate an impression of a
    data structure when recorded in a file, that is, the composite
    form of multiple XDR data fields. It is unusual in that it is
    defined proceedurally, not by any data structure. Routines named
    "*_ll_*" (low level) are packages of these low-level XDR calls;
    these should be considered to be private to libplanio, that is,
    not called directly from applications.

    See RFC 1014 "XDR: External Data Representation Standard" from
    Sun Microsystems, Inc, June 1987 for more information.
*/

// HISTORY
// 8-Sep-2004 zsx: Add doxygen-style comments.

#ifndef _PLUNC_HAVE_LIBPLANIO_H
#define _PLUNC_HAVE_LIBPLANIO_H

#ifdef PLAN_WINNT

#include "rpc_winnt/xdr.h"
#include "getopt.h"
struct passwd {
    char *pw_dir;
};
int geteuid();
passwd *getpwuid(int);

#else

#include <pwd.h>

#endif


#ifndef WIN32
// This is defined on WIN32 system.  For other systems
// defining it as 0 will have no effect, but is necessary,
// since it is built into the code.
#ifndef O_BINARY
#define O_BINARY	0
#endif
#endif


/*
#ifdef __cplusplus
extern "C" {
#endif
*/


#ifdef _PLUNC_HAVE_EXTBM_H
/** Get the filter with the specified id.
\param phys_dat_dir string, the directory where the filter file is.
\param filter_id an integer, the id of the filter to get.
\param filter a FILTER pointer, to the storage of the obtained filter
information.
\return 0 if successful, -1 otherwise. */
int get_filter(char *phys_dat_dir, int filter_id, FILTER *filter);
#endif

/** Get the directory containing the physics and dosimetry data
of the treatment machines from the predefined environment variable
PHY_DAT_DIR.
\return a string, the directory defined as environmental variable
PHY_DAT_DIR, or the current directory if the variable is not defined
in the user's environment. */
char *get_phys_dat_dir();

/** Get the site specific PLUNC environment variables from text file. */
int  get_plunc_env();

/** Get the root directory of PLUNC from the predefined environment
variable UP (User Planning was the old name of PLUNC). All directories
that are not specific to a patient may be found within the root directory.
\return a string, the root directory defined in UP, or the current
directory as a default. */
char *get_plunc_root_dir();

#ifdef _PLUNC_HAVE_PLAN_IM_H
/** Read a specified slice of plan_im image from file
\param fdes an integer, the plan_im file descriptor.
\param scan a array of short, which stores the scan pixels just read.
\param h a pointer to plan_im_header, which contains the header information.
\param slice an integer, the slice of plan_im to read.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int read_scan_xy(int fdes, short *scan, plan_im_header *h, int slice);

/** Write a specified slice of plan_im image to file
\param fdes an integer, the plan_im file descriptor.
\param scan a array of short, which stores the scan pixels to be written.
\param h a pointer to plan_im_header, which contains the header information.
\param slice an integer, the slice of plan_im to write.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise.*/
int write_scan_xy(int fdes, short *scan, plan_im_header *h, int slice,
		  bool in_place = true);

/** Get a subset of slices (image scans) from a plan_im file.
\param fd an integer, the given plan_im file descriptor.
\param im_header a plan_im_header pointer, to the image header info.
\param start_slice an integer, the index of the starting slice index of
interest within per_scan[].
\param end_slice an integer, the index of the ending slice of
interest.
\param buffer an array of PIXELTYPE to store the pixels of the
slices of interest. The caller must allocate enough space for
(end_slice - start_slice + 1) * x_dim * y_dim pixels.
\return 0 if successful, -1 otherwise. */
int
get_slice(int fd, plan_im_header *im_header,
	  int start_slice, int end_slice, PIXELTYPE *buffer);

/** Overwrite a subset of slices in a plan_im file with what is stored
in the buffer.
\param fd an integer, the given plan_im file descriptor.
\param im_header a plan_im_header pointer, to the image header info.
\param start_slice an integer, the index of the starting slice of
interest.
\param end_slice an integer, the index of the ending slice of
interest.
\param buffer an array of PIXELTYPEs containing the pixels of the
slices to replace what are currently in the plan_im file.
\return 0 if successful, -1 otherwise. */
int
put_slice(int fd, plan_im_header *im_header,
	   int start_slice, int end_slice, PIXELTYPE *buffer);
#endif

#ifdef _PLUNC_HAVE_EXTBM_H
/** Get the time calculation (time_calc) information of the specified
unit.
\param phys_dat_dir string, the directory containing the \b time_calc file.
\param unit_id an integer, the id of the unit whose time_calc to get.
\param time_calc a TIME_CALC_PARAMETERS pointer, to the storage of
the obtained time_calc information.
\return 0 if successful, -1 otherwise. */
int
get_time_calc(char *phys_dat_dir, int unit_id,
	      TIME_CALC_PARAMETERS *time_calc);

/** Get the accessory information with the specified id.
\param phys_dat_dir string, the directory where the accessory file is.
\param accessory_id an integer, the id of the accessory to get.
\param accessory an ACCESSORY pointer, to the storage of
the obtained accessory information.
\return 0 if successful, -1 otherwise. */
int
get_accessory(
    char        *phys_dat_dir,
    int         accessory_id,
    ACCESSORY	*accessory);

/** Get the tray information with the specified id.
\param phys_dat_dir string, the directory where the tray file is.
\param tray_id an integer, the id of the tray to get.
\param tray a TRAY pointer, to the storage of the obtained tray
information.
\return 0 if successful, -1 otherwise. */
int
get_tray(
    char        *phys_dat_dir,
    int         tray_id,
    TRAY        *tray);

/** Get the unit information with the specified id.
\param phys_dat_dir string, the directory where the unit file is.
\param unit_id an integer, the id of the unit to get.
\param unit a UNIT pointer, to the storage of the obtained unit
information.
\return 0 if successful, -1 otherwise. */
int
get_unit(
    char        *phys_dat_dir,
    int         unit_id,
    UNIT        *unit);
#endif

#ifdef _PLUNC_HAVE_PLAN_XDR_H
/** XDR primitive for ANASTRUCT, excluding contours.
\param xdr an XDR handle.
\param an a pointer to ANASTRUCT.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int
xdr_ANASTRUCT(
    XDR         *xdr,
    ANASTRUCT   *an);

/** XDR primitive for the header of CONTOUR.
\param xdr an XDR handle.
\param contour a pointer to CONTOUR.
\return TRUE if successful, FALSE otherwise. */
int
xdr_CONTOUR(
    XDR *xdr,
    CONTOUR *contour);

/** XDR primitive for the vertices of CONTOUR. Needs to be
called twice, once for c->x and then once for c->y.
\param xdr an XDR handle.
\param xy the array of x or y position of the vertices of CONTOUR.
\param count the number of vertices in the contour.
\return TRUE if successful, FALSE otherwise. */
int
xdr_CONTOUR_X_Y(
    XDR *xdr,
    float *xy,
    int count);

/** XDR primitive for PNT3D.
\param xdr an XDR handle.
\param point a pointer to PNT3D.
\return TRUE if successful, FALSE otherwise. */
int
xdr_POINT(
    XDR *xdr,
    PNT3D *point);

/** XDR primitive for Boolean.
\param xdr an XDR handle.
\param b a pointer to Boolean.
\return TRUE if successful, FALSE otherwise. */
int
xdr_Boolean(
    XDR *xdr,
    Boolean *b);

#ifdef _PLUNC_HAVE_EXTBM_H
int
/** XDR primitive for EXT_BEAM, excluding the beam_outline,
collimated_beam_outline, and custom_filter[2] in EXT_BEAM.
\param xdr an XDR handle.
\param eb a pointer to EXT_BEAM.
\return TRUE if successful, FALSE otherwise. */
xdr_EXT_BEAM(
    XDR         *xdr,
    EXT_BEAM    *eb);

/** XDR primitive for EXT_BEAM in the old format where the filters
are not included.
\param xdr an XDR handle.
\param eb a pointer to EXT_BEAM.
\return TRUE if successful, FALSE otherwise. */
int
xdr_EXT_OLD_BEAM(
    XDR         *xdr,
    EXT_BEAM    *eb);
#endif

#ifdef _PLUNC_HAVE_BRACHY_H
/** XDR primitive for sourc_count and seed_count in IMPLANT_DESC.
\param xdr an XDR handle.
\param i a pointer to IMPLANT_DESC.
\return TRUE if successful, FALSE otherwise. */
int
xdr_IMPLANT_DESC(
    XDR *xdr,
    IMPLANT_DESC *i);

/** XDR primitive for SOURCE.
\param xdr an XDR handle.
\param s a pointer to SOURCE.
\return TRUE if successful, FALSE otherwise. */
int
xdr_SOURCE(
    XDR *xdr,
    SOURCE *s);

/** XDR primitive for SEED.
\param xdr an XDR handle.
\param s a pointer to SEED.
\return TRUE if successful, FALSE otherwise. */
int
xdr_SEED(
    XDR *xdr,
    SEED *s);
#endif

/** XDR primitive for TWOD_GRID, excluding matrix and inside.
\param xdr an XDR handle.
\param g a pointer to TWOD_GRID.
\return TRUE if successful, FALSE otherwise. */
int
xdr_TWOD_GRID(
    XDR *xdr,
    TWOD_GRID *g);

/** XDR primitive for GRID, excluding matrix.
\param xdr an XDR handle.
\param g a pointer to GRID.
\return TRUE if successful, FALSE otherwise. */
int
xdr_GRID(
    XDR *xdr,
    GRID *g);

/** XDR primitive for count in WEIGHTS. This function is not
currently used.
\param xdr an XDR handle.
\param w a pointer to WEIGHTS.
\return TRUE if successful, FALSE otherwise. */
int
xdr_WEIGHTS(
    XDR *xdr,
    WEIGHTS *w);

/** XDR primitive for WEIGHT.
\param xdr an XDR handle.
\param w a pointer to WEIGHT.
\return TRUE if successful, FALSE otherwise. */
int
xdr_WEIGHT(
    XDR *xdr,
    WEIGHT *w);

/** XDR primitive for CALC_POINT.
\param xdr an XDR handle.
\param c a pointer to CALC_POINT.
\return TRUE if successful, FALSE otherwise. */
int
xdr_CALC_POINT(
    XDR *xdr,
    CALC_POINT *c);

/** XDR primitive for PLAN, excluding isodoses.
\param xdrs an XDR handle.
\param p a pointer to PLAN.
\return TRUE if successful, FALSE otherwise. */
int
xdr_PLAN(
    XDR *xdrs,
    PLAN *p);

#ifdef _PLUNC_HAVE_PLAN_IM_H
/** A combined routine for reading or writing plan_im_header.
\param xdr an XDR handle.
\param p a pointer to plan_im_header.
\param write_flag an integer, TRUE for write, or FLASE for read.
\return TRUE if successful, FALSE otherwise. */
int
xdr_PLAN_IM_HEADER(
    XDR                 *xdr,
    plan_im_header      *p,
    int                 write_flag);

/** XDR primitive for writing plan_im_header. Handles extended header.
\param xdr an XDR handle.
\param p a pointer to plan_im_header.
\return TRUE if successful, FALSE otherwise. */
int
xdr_write_PLAN_IM_HEADER(
    XDR                 *xdr,
    plan_im_header      *p);

/** XDR primitive for reading plan_im_header. Handles extended header.
\param xdr an XDR handle.
\param p a pointer to plan_im_header.
\return TRUE if successful, FALSE otherwise. */
int
xdr_read_PLAN_IM_HEADER(
    XDR                 *xdr,
    plan_im_header      *p);

/** XDR primitive for reading or writing plan_im_header.
\param xdr an XDR handle.
\param p a pointer to plan_im_header.
\param write_flag an integer, TRUE for write, or FLASE for read.
\return TRUE if successful, FALSE otherwise. */
int
xdr_read_write_PLAN_IM_HEADER(
    XDR                 *xdr,
    plan_im_header      *p,
    int			write_flag);

/** XDR primitive for per_scan_info.
\param xdr an XDR handle.
\param s a pointer to per_scan_info.
\return TRUE if successful, FALSE otherwise. */
int
xdr_PER_SCAN_INFO(
    XDR *xdr,
    per_scan_info *s);
#endif

/** XDR primitive for PDATE.
\param xdr an XDR handle.
\param date a pointer to PDATE.
\return TRUE if successful, FALSE otherwise. */
int
xdr_DATE(
    XDR *xdr,
    PDATE  *date);

/** XDR primitive for scan data. This function is not currently used.
\param xdr an XDR handle.
\param scan an array of short, which contains the scan's pixel values.
\param resolution an integer, the number of pixels in one dimension.
\return TRUE if successful, FALSE otherwise. */
int
xdr_SCAN( XDR *xdr, short *scan, int resolution);

#ifdef _PLUNC_HAVE_BRACHY_H
int
xdr_ll_brachy_object(XDR_fdes *xdr_fdes, BRACHY_OBJECT *o1, int read_flag);

/** XDR primitive for SEED_SPEC. This function is not currently used.
\param xdr an XDR handle.
\param s a pointer to SEED_SPEC.
\return TRUE if successful, FALSE otherwise. */
int
xdr_SEED_SPEC( XDR *xdr, SEED_SPEC *s);
#endif

#endif

/** Read one TWOD_GRID record from a file at the current file position.
\param fdes an integer, the 2D grid file descriptor.
\param g a pointer to TWOD_GRID, which stores the 2D grid just read.
\param desc_only an integer, if TRUE, only the grid description
(header) is read; otherwise, the grid matrix is also read.
\return 0 if successful, -1 otherwise. */
int
read_2d_grid(int fdes, TWOD_GRID *g, int desc_only);

/** Write one TWOD_GRID record to a file at the current file position.
\param fdes an integer, the 2D grid file descriptor.
\param g a pointer to TWOD_GRID, which stores the 2D grid to be written.
\param desc_only an integer, if TRUE, only the grid description
(header) is written; otherwise, the grid matrix is also written.
\return 0 if successful, -1 otherwise. */
int
write_2d_grid(int fdes, TWOD_GRID *g, int desc_only);

/** Read one ANASTRUCT record from a file at the current file position.
\param fdes an integer, the anastruct file descriptor.
\param anastruct_desc a pointer to ANASTRUCT, which stores the
anastruct just read.
\return 0 if successful, -1 otherwise. */
int
read_anastruct(int fdes, ANASTRUCT *anastruct_desc);

/** Write one ANASTRUCT record to a file at the current file position.
\param fdes an integer, the anastruct file descriptor.
\param anastruct_desc a pointer to ANASTRUCT, which stores the
anastruct to be written.
\return 0 if successful, -1 otherwise. */
int
write_anastruct(int fdes, ANASTRUCT *anastruct_desc);

#ifdef _PLUNC_HAVE_EXTBM_H
/** Read one EXT_BEAM record from a file at the current file position.
\param fdes an integer, the beam file descriptor.
\param eb a pointer to EXT_BEAM, which stores the beam
just read.
\return 0 if successful, -1 otherwise. */
int
read_beam(int fdes, EXT_BEAM *eb);

/** Write one EXT_BEAM record to a file at the current file position.
\param fdes an integer, the beam file descriptor.
\param eb a pointer to EXT_BEAM, which stores the beam
to be written.
\return 0 if successful, -1 otherwise. */
int
write_beam(int fdes, EXT_BEAM *eb);

/** Write one EXT_BEAM record of an old format to a file.
\param fdes an integer, the beam file descriptor.
\param eb a pointer to EXT_BEAM, which stores the beam
to be written.
\return 0 if successful, -1 otherwise. */
int
write_old_beam(int fdes, EXT_BEAM *eb);

/** Checks to see if the beam being read in has a different number
of leaf pairs than what the unit currently specifies and adjusts
the number to match the unit. This was added to ease the transition
from our 58 leaf head to 82 leaf head and only works for that case.
\param eb a pointer to the beam to check. */
void check_mlc_head(EXT_BEAM *eb);
#endif

/** @name Routines for reading and writing basic scalars and arrays of
    scalars.

    Why this layer exists: 1) call these instead of the XDR lib
    directly to deal with bugs in (some vendor's) XDR lib's in a
    uniform fashion. 2) These routine names are explicit as to whether
    the user wants to read or write the stream, whereas XDR contains
    read/write as state.
*/
//@{
int read_byte(int fdes, char *ptr);
int write_byte(int fdes, char *ptr);
int read_bytes(int fdes, char *ptr, int num);
int write_bytes(int fdes, char *ptr, int num);
int read_write_byte(int fdes, char *ptr, char mode);
int read_int(int fdes, int * ptr);
int write_int(int fdes, int *ptr);
int read_ints(int fdes, int *ptr, int num);
int write_ints(int fdes, int *ptr, int num);
int read_write_ints(int fdes, int *ptr, int num, char mode);
int read_short(int fdes, short *ptr);
int write_short(int fdes, short *ptr);
int read_shorts(int fdes, short *ptr, int num);
int write_shorts(int fdes, short *ptr, int num);
int read_write_short(int fdes, short *ptr, char mode);
int read_float(int fdes, float *ptr);
int write_float(int fdes, float *ptr);
int read_floats(int fdes, float *ptr, int num);
int write_floats(int fdes, float *ptr, int num);
int read_write_floats(int fdes, float *ptr, int num, char mode);
int read_byte_array(int fdes, char *ptr, int num);
int write_byte_array(int fdes, char *ptr, int num);
//@}

/** Read one CONTOUR record from a file at the current file position.
\param fdes an integer, the contour file descriptor.
\param c a pointer to CONTOUR, which stores the contour
just read.
\return 0 if successful, -1 otherwise. */
int read_contour(int fdes, CONTOUR *c);

/** Write one CONTOUR record to a file at the current file position.
\param fdes an integer, the contour file descriptor.
\param c a pointer to CONTOUR, which stores the contour to be written.
\return 0 if successful, -1 otherwise. */
int write_contour(int fdes, CONTOUR *c);

#ifdef _PLUNC_HAVE_EXTBM_H
/** Read one FILTER record from a file at the current file position.
\param fdes an integer, the filter file descriptor.
\param filter a pointer to FILTER, which stores the filter just read.
\return 0 if successful, -1 otherwise. */
int read_filter(int fdes, FILTER *filter);

/** Write one FILTER record to a file at the current file position.
\param fdes an integer, the filter file descriptor.
\param filter a pointer to FILTER, which stores the filter to be
written.
\return 0 if successful, -1 otherwise. */
int write_filter(int fdes, FILTER *filter);
#endif

/** Read one GRID record from a file at the current file position.
\param fdes an integer, the grid file descriptor.
\param g a pointer to GRID, which stores the grid just read.
\param desc_only an integer, if TRUE, only the grid description is
read; otherwise, the grid matrix is also read.
\return 0 if successful, -1 otherwise. */
int read_grid(int fdes, GRID *g, int desc_only);

/** Write one GRID record to a file at the current file position.
\param fdes an integer, the grid file descriptor.
\param g a pointer to GRID, which stores the grid to be written.
\param desc_only an integer, if TRUE, only the grid description is
written; otherwise, the grid matrix is also written.
\return 0 if successful, -1 otherwise. */
int write_grid(int fdes, GRID *g, int desc_only);

#ifdef _PLUNC_HAVE_PLAN_IM_H
/** Read plan_im header from a file at the current file position.
\param fdes an integer, the plan_im file descriptor.
\param header a pointer to plan_im_header, which stores the image
header.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int read_image_header(int fdes, plan_im_header *header);

/** Write plan_im header to a file at the current file position.
\param fdes an integer, the plan_im file descriptor.
\param header a pointer to plan_im_header, which stores the image
header to be written.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int write_image_header(int fdes, plan_im_header *header);

/** Sets the static int image_compression variable to determine
what compression level (0: no compression, 1: lossless compression)
to use when saving a plan_im file.
\param val an integer, the plan_im compression level.
\return 0 if successful, -1 otherwise. */
int set_image_compression(int val);

#endif

#ifdef _PLUNC_HAVE_BRACHY_H
/** Read one IMPLANT_DESC record from a file at the current file
position.
\param fdes an integer, the implant file descriptor.
\param implant_desc a pointer to IMPLANT_DESC, which stores the
implant just read.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int read_implant(int fdes, IMPLANT_DESC *implant_desc);

/** Not currently used. */
int read_object(int fdes, BRACHY_OBJECT *object_desc);

/** Read one BRACHY_OBJECTS record from a file at the current file
position.
\param fdes an integer, the brachy objects file descriptor.
\param objects_desc a pointer to BRACHY_OBJECTS, which stores the
brachy objects just read.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int read_objects(int fdes, BRACHY_OBJECTS *objects_desc);
#endif

#ifdef _PLUNC_HAVE_EXTBM_H
/** Read one ACCESSORY record from a file at the current file
position.
\param fdes an integer, the accessory file descriptor.
\param accessory a pointer to ACCESSORY, which stores the accessory
just read.
\return 0 if successful, -1 otherwise. */
int read_accessory(int fdes, ACCESSORY *accessory);

/** Write one ACCESSORY record to a file at the current file
position.
\param fdes an integer, the accessory file descriptor.
\param accessory a pointer to ACCESSORY, which stores the accessory
to be written.
\return 0 if successful, -1 otherwise. */
int write_accessory(int fdes, ACCESSORY *accessory);

/** Read all ACCESSORY records from the accessory file in phys_dat_dir.
\param count a pointer to integer, which store the total number of
accessories just read.
\param a a pointer to an array of ACCESSORY, which stores all of the
accessories just read.
\return 0 if successful, -1 otherwise. */
int read_accessories(int *count, ACCESSORY **a);

/** Read all FILTER records from the filter file in phys_dat_dir.
\param count a pointer to integer, which store the total number of
filters just read.
\param f a pointer to an array of FILTER, which stores all of the
filters jsut read.
\return filter count if successful, -1 otherwise. */
int read_filters(int *count, FILTER **f);

/** Read all TRAY records from the tray file in phys_dat_dir.
\param count a pointer to integer, which store the total number of
trays just read.
\param t a pointer to an array of TRAY, which stores all of the
trays just read.
\return 0 if successful, -1 otherwise. */
int read_trays(int *count, TRAY **t);

/** Read all UNIT records from the unit file in phys_dat_dir.
\param count a pointer to integer, which store the total number of
units just read.
\param u a pointer to an array of UNIT, which stores all of the
units just read.
\return 0 if successful, -1 otherwise. */
int read_units(int *count, UNIT **u);
#endif

/** Read a specified number of slices of a THREED_MAP record from a
file.
\param name a string, the name of the file storing the map to be read.
\param map a pointer to THREED_MAP, the storage for the map just read.
\param slice an integer, the index of the first slice to be read.
\param num_slices an integer, the total number of slices to read.
\return an integer, 0 if file can't be open, XDR_NOERROR if successful,
XDR_ERROR otherwise. */
int read_map(char *name, THREED_MAP *map, int slice, int num_slices);

/** Write a specified number of slices of a THREED_MAP record to a
file.
\param name a string, the name of the file to store the map.
\param map a pointer to THREED_MAP, the storage for the map to be
written.
\param slice an integer, the index of the first slice to be written.
\param num_slices an integer, the total number of slices to written.
\return an integer, 0 if file can't be open, XDR_NOERROR if successful,
XDR_ERROR otherwise. */
int write_map(char *name, THREED_MAP *map, int slice, int num_slices);

/** Read in a transformation matrix.
\param name a pointer to string, the name of the matrix file.
\param mat a LIB3D_TRANSFORM, the transformation matrix to set.
\return 0 if successful, non-zero otherwise. */
int read_matrix(char *name, LIB3D_TRANSFORM mat);

/** Write out a transformation matrix.
\param name a pointer to string, the name of the matrix file.
\param mat a LIB3D_TRANSFORM, the transformation matrix.
\return 0 if successful, non-zero otherwise. */
int write_matrix(char *name, LIB3D_TRANSFORM mat);

/** Read one PLAN record from a file at the current file position.
\param fdes an integer, the plan file descriptor.
\param p a pointer to PLAN, which stores the plan just read.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int read_plan(int fdes, PLAN *p);

/** Write one PLAN record to a file at the current file position.
\param fdes an integer, the plan file descriptor.
\param p a pointer to PLAN, which stores the plan to be written.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int write_plan(int fdes, PLAN *p);

/** Read all PLAN records from a file at current file position.
\param fdes an integer, the plan file descriptor.
\param plan_desc a pointer to an array of PLAN, which stores all of
the plans just read.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int read_plans(int fdes, PLAN *plan_desc);

/** Write a number of PLAN records to a file at current file position.
\param fdes an integer, the plan file descriptor.
\param plan_desc a pointer to an array of PLAN, which stores all of
the plans to be written.
\param num in integer, the specified number of PLAN records to write.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int write_plans(int fdes, PLAN *plan_desc, int num);

/** Read one CALC_POINTS record from a file at current file position.
\param fdes an integer, the points file descriptor.
\param p a pointer to CALC_POINTS, which stores the points just read.
\param desc_only an integer, if TRUE, only the point count and positions
are read; otherwise, the point doses are also read.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int read_points(int fdes, CALC_POINTS *p, int desc_only);

/** Write one CALC_POINTS record to a file at current file position.
\param fdes an integer, the points file descriptor.
\param p a pointer to CALC_POINTS, which stores the points to be written.
\param desc_only an integer, if TRUE, only the point count and positions
are written; otherwise, the point doses are also written.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int write_points(int fdes, CALC_POINTS *p, int desc_only);

#ifdef _PLUNC_HAVE_EXTBM_H
/** Read one SAR_TABLE record from a file at current file position.
\param fdes an integer, the SAR table file descriptor.
\param sar a pointer to SAR_TABLE, which stores the SAR table just read.
\return 0 if successful, XDR_ERROR otherwise. */
int read_sar(int fdes, SAR_TABLE *sar);

/** Write one SAR_TABLE record to a file at current file position.
\param fdes an integer, the SAR table file descriptor.
\param sar a pointer to SAR_TABLE, which stores the SAR table to be
written.
\return 0 if successful, XDR_ERROR otherwise. */
int write_sar(int fdes, SAR_TABLE *sar);
#endif


/** SWAp Byte - exchange every two adjacent bytes up to a specified
total number of bytes in the source array, and save the result in the
destination array. The source and destination may overlap but must be
the same size. "my_" in the name is intended to reduce namespace
collisions.
\param src an array of bytes, the source array.
\param dest an array of bytes, the destination array.
\param num an integer, the total number of bytes to process. */
void my_swab(char *src, char *dest, int num);

#ifdef _PLUNC_HAVE_EXTBM_H
/** Read one TIME_CALC_PARAMETERS record from a file at the current
file position.
\param fdes an integer, the time_calc file descriptor.
\param time_calc a pointer to TIME_CALC_PARAMETERS, which stores the
time_calc just read.
\return 0 if successful, XDR_ERROR otherwise. */
int read_time_calc(int fdes, TIME_CALC_PARAMETERS *time_calc);

/** Write one TIME_CALC_PARAMETERS record to a file at the current
file position.
\param fdes an integer, the time_calc file descriptor.
\param time_calc a pointer to TIME_CALC_PARAMETERS, which stores the
time_calc to be written.
\return 0 if successful, XDR_ERROR otherwise. */
int write_time_calc(int fdes, TIME_CALC_PARAMETERS *time_calc);

/** Read one TRAY record from a file at the current file position.
\param fdes an integer, the tray file descriptor.
\param tray a pointer to TRAY, which stores the tray just read.
\return 0 if successful, XDR_ERROR otherwise. */
int read_tray(int fdes, TRAY *tray);

/** Write one TRAY record to a file at the current file position.
\param fdes an integer, the tray file descriptor.
\param tray a pointer to TRAY, which stores the tray to be written.
\return 0 if successful, XDR_ERROR otherwise. */
int write_tray(int fdes, TRAY *tray);

/** Read one UNIT record from a file at the current file position.
\param fdes an integer, the unit file descriptor.
\param unit a pointer to UNIT, which stores the unit just read.
\return 0 if successful, XDR_ERROR otherwise. */
int read_unit(int fdes, UNIT *unit);

/** Write one UNIT record to a file at the current file position.
\param fdes an integer, the unit file descriptor.
\param unit a pointer to UNIT, which stores the unit to be written.
\return 0 if successful, XDR_ERROR otherwise. */
int write_unit(int fdes, UNIT *unit);
#endif

/** Read or write one NORM record at the current file position.
\param fdes an integer, the norm info file descriptor.
\param norm a pointer to NORM, which stores the norm info to be
written.
\param mode a char, XDR_READ_MODE for read, and XDR_WRITE_MODE for
write.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int read_write_norm(int fdes, NORM *norm, char mode);

/** Read one NORM record from a file at the current file position.
\param fdes an integer, the norm info file descriptor.
\param norm a pointer to NORM, which stores the norm info just read.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int read_norm(int fdes, NORM *norm);

/** Write one NORM record to a file at the current file position.
\param fdes an integer, the norm info file descriptor.
\param norm a pointer to NORM, which stores the norm info to be
written.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int write_norm(int fdes, NORM *norm);

/** Read one WEIGHTS record from a file at the current file position.
\param fdes an integer, the weights file descriptor.
\param w a pointer to WEIGHTS, which stores the weights just read.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int read_weights(int fdes, WEIGHTS *w);

/** Write one WEIGHTS record to a file at the current file position.
\param fdes an integer, the weights file descriptor.
\param w a pointer to WEIGHTS, which stores the weights to be written.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int write_weights(int fdes, WEIGHTS *w);

/** A convenience routine that calls either calloc() or
realloc(). Calls realloc() when ptr != NULL and size != 0. If size =
0, it defaults to size = sizeof(int); If ptr == NULL, ptr =
calloc(size,1). */
void *Realloc(void *ptr, int size);

/** Read buff_size bytes from a socket or file.
\param fdes an integer, the file descriptor.
\param buff a pointer to chars, the array of bytes read into.
\param buff_size an integer, the number of bytes to read.
\return 0 if successful, -1 otherwise */
int sock_read(int fdes, char *buff, unsigned int buff_size);

/** Write buff_size bytes to a socket or file.
\param fdes an integer, the file descriptor.
\param buff a pointer to chars, the array of bytes to write.
\param buff_size an integer, the number of bytes to write.
\return 0 if successful, -1 otherwise */
int sock_write(int fdes, char *buff, unsigned int buff_size);

/** Not currently used. */
int runlen_encode(char *filename, RGB_PIXEL *buf, int xdim, int ydim);

#ifdef _PLUNC_HAVE_BRACHY_H
/** Write one IMPLANT_DESC record to a file at the current file
position.
\param fdes an integer, the implant file descriptor.
\param implant_desc a pointer to IMPLANT_DESC, which stores the
implant to be written.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int write_implant(int fdes, IMPLANT_DESC *implant_desc);

/** Not currently used. */
int write_object(int fdes, BRACHY_OBJECT *object_desc);

/** Write one BRACHY_OBJECTS record to a file at the current file
position.
\param fdes an integer, the brachy objects file descriptor.
\param objects_desc a pointer to BRACHY_OBJECTS, which stores the
brachy objects to be written.
\return XDR_NO_ERROR if successful, XDR_ERROR otherwise. */
int write_objects(int fdes, BRACHY_OBJECTS *objects_desc);
#endif

/** Check to see if a socket is ready for reading.
\param fd an integer, the socket descriptor.
\param msec an integer, the number of milliseconds to wait.
\return 0 if no data present, non-zero otherwise. */
int fd_poll(int fd, int msec);

#ifdef _PLUNC_HAVE_PLAN_XDR_H
/** Low level routine for reading bytes to an XDR stream.
\param xdrs a pointer to an XDR stream.
\param dest a pointer to a character array to read into.
\param count an integer, the number of bytes to read.
\return 1 if successful, 0 otherwise. */
bool_t xdrfd_get_bytes(XDR *xdrs, caddr_t dest, int count);

/** Low level routine for writing bytes to an XDR stream.
\param xdrs a pointer to an XDR stream.
\param src a pointer to the character array to write from.
\param count an integer, the number of bytes to write.
\return 1 if successful, 0 otherwise. */
bool_t xdrfd_put_bytes(XDR *xdrs, caddr_t src, int count);

/** Determine whether the end of file is reached.
\param fdes an integer, the file descriptor.
\return TRUE if end of file is reached, FALSE if not. */
bool_t xdr_eof(int fdes);

/** Initialize a stdio xdr stream.
\param xdrs an XDR, the xdr stream handle for use on the stream file.
\param fd an integer, the stream file descriptor.
\param op an xdr_op, set as the xdr operation flag. */
void xdrfd_create(XDR *xdrs, int fd, enum xdr_op op);

/* Determine the current position in an XDR stream.
\param xdrs a pointer to an XDR stream.
\return an integer, the current byte position. */
int xdr_pos(XDR *xdrs);

/* Force the XDR stream to be written to disk (or socket).
\param xdrs a pointer to an XDR stream.
\return 1 on success, 0 otherwise. */
bool_t xdrfd_flush_buff(XDR *xdrs);

/** Open a stdio xdr stream.
\param fdes an integer, the file descriptor.
\param mode a char, the file mode: XDR_READ_MODE or XDR_WRITE_MODE.
\return an pointer to an XDR_fdes, or NULL if an error occurs. */
XDR_fdes * xdr_ll_open_stream(int fdes, char mode);

/** Close a stdio xdr stream, writing the buffer first if in write mode.
\param xdr_fdes a pointer to an XDR_fdes, the xdr file descriptor,
containing the file descriptor and the xdr stream. */
void xdr_ll_close_stream(XDR_fdes *xdr_fdes);

/** Close a stdio xdr stream.
\param fdes an integer, the file descriptor. */
void xdr_close(int fdes);

/** Print an error message.
\param str a string, contains the error message to be printed. */
void XDR_ReportError(char *str);

/** On an output stream, causes any buffered but unwritten data to
be written.
\param xdr_fdes a pointer to an xdr file descriptor, containing the
file descriptor and the xdr stream. */
void xdr_ll_flush(XDR_fdes *xdr_fdes);

#endif

/** Read Lantis export site specific data from $UP\\site\\lantis_data
file. */
int get_mosaiq_data(MOSAIQ_DATA *mosaiq);

/** Set the global is_sock[fd] to reflect that this fd is for a socket.
\param fd an integer, the socket descriptor. */
void set_socket(int fd);

/** Clear the global is_sock[fd] to reflect that this fd is not for a socket.
\param fd an integer, the file or socket descriptor. */
void clear_socket(int fd);

/** Checks to see if a descriptor is for a socket or file.
\param fd an integer, the descriptor in question.
\returns 1 if fd is for a socket, 0 otherwise. */
int is_socket(int fd);

/** Proper method to close a socket to guarantee that is_socket() works.
\param fd an integer, the socket descriptor.
\returns 0 if successfully closed, -1 otherwise. */
int close_socket(int fd);

/*
#ifdef __cplusplus
}
#endif
*/

#endif

