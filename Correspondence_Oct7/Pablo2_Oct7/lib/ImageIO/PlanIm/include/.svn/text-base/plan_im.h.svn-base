// HISTORY
// 2-Aug-2004 zsx: Add doxygen-style comments.

/** @file plan_im.h
    This header file describes the image storage format used in the
    PLUNC system.

  The conventions for image geometry are as follows: \n
  - Scans should always be viewed from the foot of the patient. In the
     simple case where the patient is supine and head first, the
     patient's left is on the right side of the image and the
     pixel_to_patient[0][0] is positive. For supine, feet first, the
     patient's left is still on the right side of the image but
     the pixel_to_patient[0][0] is negative. \n
  - PLanUNC's convention is different than the old PLUNC convention
     where the scans were always (?) viewed from the foot of the
     scanner couch. \n
  - In PIXEL COORDINATES, \f$p_x=0\f$ is at the left and \f$p_y=0\f$
     at the top. Pixels are 0-indexed. \n
  - Scans are conceived as being collections of SAMPLES.  This
     implies that the overall size of the imaged area is (resolution -
     1) x pixel_size. \n
  - In PATIENT COORDINATES - a RIGHT-HANDED system - when looking
     towards the isocenter of the gantry from the foot of the couch
     when it is at its norminal position: \n
     positive x is to the right of the observer horizontally; \n
     positive y is straight up; \n
     positive z is toward the foot of the couch; \n
     the x origin is couch center; \n
     the y origin is tabletop; \n
     the z origin is an arbitrary clinically-locatable table position; \n
     the "table_height" field is the number of pixels above the bottom
     of the image that the top of the table is seen. \n
  - In PIXEL COORDINATES: \n
     x = 0 is px0 = (resolution - 1) / 2.0 in FLOATING POINT - for a scan
     with an even number of pixels; \n
     px0 will have a fractional part of 0.5, for instance the
     x = 0 patient coordinate is found in a 256x256 image at
     px0 = 127.5 (that's midway between pixels 127 and 128); \n
     y = 0 is py0 = (resolution - 1 - table_height).
     For example if the table height for a 256x256 image is 25,
     the y = 0 patient coordinate occurs at py0 = 230.0; \n
     The z = 0 patient coordinate can have any fractional pixel
     value since it is not constrained to lie in a sample plane. \n
  - The "pixel_to_patient_TM" transform is a 4x4 transform that takes
     pixel coordinates from slice 0 and turns them into patient
     coordinates.  To obtain a transformation matrix for any other slice,
     replace pixel_to_patient_TM[3][2] with the z_position of that slice.
     This assumes we are restricted to table and gantry angles of 0.  If
     we want to be able to handle non-zero values for these fields later
     on, it will be necessary to have a unique pixel_to_patient_TM for
     each slice (i.e., put it in per_scan record). \n
  - The "first slice", or "slice 0" is defined as that record to which
     plan_im_header.per_scan[0] points. \n

  Note that this may not be the same as the scanner's inherent
  conventions.  The transformation is made by some scanner-specific magic
  which lives in the routines which create planning image files from
  tape, floppy, or network.  The original sense of the scan lines is
  thrown away.  The original table angle is demagicated and stored
  here as the table angle for the 'standard scanner' - the counterclockwise
  rotation of the table as viewed from above.  Similarly the gantry
  tilt is rationalized to be the counterclockwise rotation of the
  gantry as viewed from the right side of the machine (from +x in
  patient coordinates).  No magic remains in this specification.

\note

  - We are currently restricting table and gantry tilt to be zero in
  PLUNC software.  The corresponding fields (with implied
  ability to have different values for each scan) are for future
  expansion and are provided here in the interest of saving our having
  to think this spec through again later.  That is probably wishful
  thinking.

  - It would probably pay to create a scanner description file (much like
  the treatment unit description file) that rationalizes the magic and
  further allows users to pick the scanner type from a menu.  This
  will certainly be required if it turns out that a single scanner
  uses different view-from conventions for different types of scans
  (head-cephalic vs. body-caudal for instance).

*/

#ifndef _PLUNC_HAVE_PLAN_IM_H
#define _PLUNC_HAVE_PLAN_IM_H

/** Old maximum number of scans allowed in a plan_im image. */
#define OLD_ABSURDLY_LARGE_MAX_SCANS (500)

/** Maximum number of scans allowed in a plan_im image. */
#define ABSURDLY_LARGE_MAX_SCANS (4096)

/** @name PLAN_IM extended header opcodes

   Extended header opcodes. These are four byte integer opcodes.
   The most significant byte is 01 for a PLUNC recognizable opcode
   (other application specific categories can use something else).
   The next two bytes give 65536 possible opcodes. The last byte
   indicates what type of thing follows: \n
	0x00	no data \n
	0x01	xdr int \n
	0x02	xdr float \n
	0x03	xdr string \n
	0x04	xdr bytes \n
	0x05	xdr array \n
	0xff	unknown \n
	Others will be defines as needed.\n
   The four bytes in the file is an xdr int indicating how many
   bytes there are for the data field, it is the application's
   writer routine to compute this appropriately (for example, a
   xdr string of length 9 would have 4 bytes for the length, and
   then 12 bytes (due to paddind to a four byte multiple) for
   the string with the last 3 bytes being 0's, so the total length
   of the xdr string data would be 16 bytes.
*/
//@{
#define PLAN_IM_FIRST_OPCODE	0x01000000
#define PLAN_IM_UNIT_NUMBER	0x01000103
#define PLAN_IM_Y_DIM		0x01000201
#define PLAN_IM_Y_SIZE		0x01000302
#define PLAN_IM_TIME_COUNT	0x01000401
#define PLAN_IM_CUT		0x01000501
#define PLAN_IM_COMPRESSED	0x01000601
#define PLAN_IM_END		0x01ffff00
//@}

typedef short PIXELTYPE;

/** @name Enumerators and descriptions

   The "enumerator description" array (*_desc[]) assigns a unique
   textual name to each enum that should be used across all
   applications. An application queries the list of all enum's of a
   specified type simply by looking through all the elements of the
   array.

   These are used in user interfaces which need to present choices in
   textual form. This eliminates the problem of finding & updating all
   the applications that duplicate this info as hard-coded constants
   every time the enum's change.

   Make sure the enums and the enum_desc arrays are always
   synchronized.
*/
//@{
enum position_list
{
    bogus_position, prone, supine, on_left_side, on_right_side
};

#ifdef PLAN_IM_EXTERN
int position_desc_count = 5;
/** The textual names of the patient's table positions. */
char *position_desc[] =
{
  "in a bogus position",
  "prone",
  "supine",
  "on left side",
  "on right side"
};
#else
extern int position_desc_count;
extern char *position_desc[];
#endif

enum entry_list
{
    bogus_entry, feet, head
};

#ifdef PLAN_IM_EXTERN
int entry_desc_count = 3;
/** The textual names of the entries (what goes first into the
    gantry). */
char *entry_desc[] =
{
  "bogus",
  "feet",
  "head"
};
#else
extern int entry_desc_count;
extern char *entry_desc[];
#endif

enum cut_list
{
    bogus, axial, sagital, coronal 
};

#ifdef PLAN_IM_EXTERN
int cut_desc_count = 4;
char *cut_desc[] =
{
  "bogus",
  "axial",
  "sagital",
  "coronal"
};
#else
extern int cut_desc_count;
extern char *cut_desc[];
#endif

enum scanners
{
    bogus_scanner, somatom, delta_scanner, ge9800,
    picker, simulix, film_scanner, philips_brilliance, emotion
};

#ifdef PLAN_IM_EXTERN
int scanners_desc_count = 9;
/** The textual names of the scanners. */
char *scanners_desc[] =
{
  "bogus scanner",
  "somatom",
  "delta scanner",
  "ge9800",
  "picker",
  "simulix",
  "film scanner",
  "philips brilliance",
  "emotion"
};
#else
extern int scanners_desc_count;
extern char *scanners_desc[];
#endif

/** Pixels in a plan_im image is one of the following types
    related to the imaging modalities. */
enum pixel_types
{
    bogus_pixel, ct_number, scout_view, mri_number, pet_number
};

#ifdef PLAN_IM_EXTERN
int pixel_desc_count = 5;
/** The textual names of the pixel types. */
char *pixel_desc[] =
{
  "bogus pixel",
  "CT number",
  "scout view",
  "MRI number",
  "PET number"
};
#else
extern int pixel_desc_count;
extern char *pixel_desc[];
#endif
//@}

/** _per_scan_info is a data structure describing the essential
    parameters assiciated with each scan of a plan_im image.
*/
typedef struct _per_scan_info
{
    /** z position of the center of this slice, related to table
    position for standard scanner. */
    float z_position;
    /** Offset to the data of this slice, including header offset. As
	returned from lseek() */
    unsigned int   offset_ptrs;
    /** Number of bytes data for slice */
    int nbytes;		
    /** Gantry angle as interpreted for standard scanner. */
    float gantry_angle;
    /** Table angle as interpreted for standard scanner. */
    float table_angle;
    int	  scan_number;	//!< The scanner's concept. -1 for interpolated scans
} per_scan_info;

/** _plan_im_header is a data structure describing the header
    information of a plan_im image.

    A plan_im file consists of a plan_im_header, then optional
    extended header codes, then the pixel intensities. \b Extensible
    means that the field may have extra content past it's normal
    length (eg, more characters or extra dimensions) in the header
    extension. \b Extended means that the value exists only in the
    header extension.

    The scans are stored in the file in arbitrary order - probably
    the order in which they are read. The per_scan array is kept in
    ascending z order.
*/
typedef struct _plan_im_header
{
    char  unit_number[20];	//!< Hospital id #. Extensible.
    char  patient_name[100]; //!< Patient name.
    PDATE date;	//!< When the image was first created by the scanner.
    char  comment[100]; //!< Comment.
    enum scanners machine_id;	//!< Currently of only potential use.
    enum position_list patient_position; //!< Patient table position.
    enum entry_list whats_first;//!< What goes first into the gantry.
    enum cut_list cut;			//!< Cut type: axial, sagital, coronal
    enum pixel_types pixel_type;//!< Pixel type: ct, mri, pet number
    int   time_count;  //!< For 4D planning, number of time sequences. Extended
    int   slice_count; //!< Total number of slices in the image.
    int   resolution;  /**< The number of pixels in x and y
			                dimension. Extensible, but
			                unused if x_dim != y_dim. */
    float pixel_size;  /**< In absolute value of cm. The pixels are
                            assumed to be squares unless x_size != y_size
			    in which case pixel_size should be ignored. */
    int   x_dim;  //!< Defaults to x_dim = resolution.
    int   y_dim;  //!< Defaults to y_dim = resolution.
    int   compression;	/**< 0:no compression, 1-9:zlib compression 
    				         default: 1 */
    float x_size; //!< Defaults to x_size = pixel_size.
    float y_size; //!< Defaults to y_size = pixel_size.
    int   table_height;	/**< Height of table in image in pixels from
			     the bottom of the image;
			     Defines the Patient CS y-origin. */
    int   min;	  //!< Min pixel value in study.
    int   max;	  //!< Max pixel value in study.
    /** Transformation matrix from pixel to patient coordinates */
    float pixel_to_patient_TM[4][4];
    /** Array of scan-specific info for each scan of the plan_im.
	Only the first slice_count values are valid. */
    per_scan_info per_scan[ABSURDLY_LARGE_MAX_SCANS];
} plan_im_header;

#endif
