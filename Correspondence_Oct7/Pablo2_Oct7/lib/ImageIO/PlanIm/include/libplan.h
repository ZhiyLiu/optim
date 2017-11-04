/** @file libplan.h
    This header file contains function prototypes for planning
    a radiation treatment provided in the library \b libplan.
*/

// HISTORY
// 13-Sept-2004 zsx: Add doxygen-style comments.

#ifndef _PLUNC_HAVE_LIBPLAN_H
#define _PLUNC_HAVE_LIBPLAN_H

#include <stdio.h>

#define CCW_SENSE	1
#define CW_SENSE	-1

/*
#ifdef __cplusplus
extern "C" {
#endif
*/

#ifdef PLAN_WINNT
#define rint(x) (((x) >= 0.0) ? ((int)((x) + 0.5)) : ((int)((x) - 0.5)))
#endif

/** Find approximate points of intersection between a 2D line segment
and a polygon specified by a set of vertex coordinates of a contour.
If the line and a polygon edge are collinear, the ray is moved a tiny
distance towards the origin to simplify the calcuation.
\param Contour a pointer to CONTOUR, the contour of interest.
\param Point1 an array of 2 floats, contains the x and y coordinates
of the start point of the line segment.
\param Point2 an array of 2 floats, contains the x and y coordinates
of the end point of the line segment.
\param Intersects an array of float pairs (order is x1,y1,x2,y2,...),
stores the points of intersection of the line segment and the polygon
up to the first FI_LIMIT points.
\param Distances an array of floats, stores the distances
from Point1 to the intersections up to the first FI_LIMIT points.
\return an integer, the total number of intersections.
\note If the return value > FI_LIMIT, not all the intersections could
be return and the calling routine needs to cope. */
int
FindNtrsc(
    CONTOUR *Contour,
    float Point1[2],
    float Point2[2],
    float Intersects[2][FI_LIMIT],
    float Distances[FI_LIMIT]);

/** Find approximate points of intersection between a bounded ray (2D
line segment) and an anastruct.
\param obj a pointer to ANASTRUCT, the anastruct of interest.
\param p1 a PNT3D, the start point of the line segment.
\param p2 a PNT3D, the end point of the line segment.
\param alpha an array of float, stores the fractional distances
to each intersection from p1.\
alpha(intersect) = d(intersect - p1)/d(p2 - p1).
\param debug an integer, not currently used.
\return an integer, the number of intersections found.
\note Based on Bob Siddon's Prism intersection paper.
\note If the number of intersections of the ray with the anastruct
in 2D is greater than the preset limit FI_LIMIT, the function
complains to stderr and returns.*/
int
ana_ntrsc(
    ANASTRUCT   *obj,
    PNT3D       p1,
    PNT3D       p2,
    float       *alpha,
    int         debug);

/** Compute the thickness of the compensator (not including the base
plate) for the line going from the beam source through the point pt
defined in beam coordinates.
\param pt a PNT3D, a beam coordinate point at the compensator definition distance.
\param ana a pointer to ANASTRUCT, the anastruct representing the compensator.
\return a float, the thickness of the compensator. */
float
compensator_ntrsc(PNT3D pt, ANASTRUCT *ana);

/** Find the cm extent of a given anastruct, and save the result
as max.x, min.x, max.y, min.y, max.z, and min.z of the anastruct.
\param ana a pointer to ANASTRUCT, the anastruct of interest.
\return an integer, always 0. */
int anastruct_min_max(ANASTRUCT *ana);

/** Not currently used. */
int anastruct_normal(ANASTRUCT *ana, PNT3D pt, PNT3D *norm);

/** Computes the volume of an anastruct in cc's.
\param ana a pointer to ANASTRUCT, the anastruct whose volume is
computed.
\return a float, the volume of the anastruct in cc's. */
float anastruct_volume(ANASTRUCT *ana);

#ifdef _PLUNC_HAVE_EXTBM_H

/** Create a bolus for a given anastruct in a given beam.
\param beam a pointer to EXT_BEAM, the beam for which the bolus
is to be created.
\param ana a pointer to ANASTRUCT, the anastruct for which the
bolus is to be created (should always be skin).
\param bolus a pointer to ANASTRUCT, the storage for the created
bolus.
\param thickness a float, the thickness of the bolus in cm.
\param resolution a float, not currently used.
\param border a float, the margin around the collimated beam
outline to extend the bolus beyond. */
void
beam_bolus(
    EXT_BEAM    *beam,
    ANASTRUCT   *ana,
    ANASTRUCT   *bolus,
    float       thickness,
    float       resolution,
    float       border);

#ifdef _PLUNC_HAVE_PLAN_IM_H
/* Not currently used. */
int export_laser_iso(plan_im_header *h, EXT_BEAM *beam, int beam_count);
#endif
#endif

/** Determine the element in an array (sorted in ascending order)
whose value is just below the value of interest.
\param value a float, the value of interest.
\param vector an array of float, the known array in ascending order.
\param vec_size an integer, the length of the array.
\return an integer, the index of the array element whose value is
just below the value of interest. */
int bin_search(float value, float *vector, int vec_size);

#ifdef _PLUNC_HAVE_EXTBM_H

/** Compare two beams for all of the data members of EXT_BEAM except
serial_number, name, and version. Used in xvsim to determine if a
dose calculation needs updatding.
\param bm1 a pointer to EXT_BEAM, the first beam to be compared.
\param bm2 a pointer to EXT_BEAM, the second beam to be compared.
\return 0 if the two beams have the same value for every parameter
compared; 1 if otherwise. */
int cmp_beams(EXT_BEAM *bm1, EXT_BEAM *bm2);

/** Determine the collimated beam outline.
\param beam a pointer to EXT_BEAM, the beam whose collimated outline
is to be determined.
\return an integer, always 0. */
int collimate(EXT_BEAM *beam);

/** Determine the collimator jaw positions that will enclose the beam
outline.
\param beam a pointer to EXT_BEAM, the beam upon which the collimator
jaw positions are based.
\param fix_max an integer, if set to 1, the collimator jaws intersect
the beam outline at the middle of the jaws' widths; if set to 0, the
collimator jaws are in contact with the beam outline without blocking
any part of it.
\return an integer, 0 if successfully collimated, or 1 if the MLC can't
define the beam shape. */
int collimator_enclosure(EXT_BEAM *beam, int fix_max);

/** Check whether MLC collimation matches the field shape.
\param beam a pointer to EXT_BEAM, the beam whose collimation is being
checked.
\return an integer,
	- 0 if MLC collimation matches the field shape
	- 1 if the beam outline has more than 2 in/out crossings for a leaf pair
	- 2 if a left leaf doesn't fit the beam outline
	- 3 if a right leaf doesn't fit the beam outline
*/
int check_ML_collimator(EXT_BEAM *beam);

/** Check whether there are any open X leaf pairs that are fully
covered by a Y jaw (illegal for current Siemens MLC).
\param beam a pointer to EXT_BEAM, the beam to be checked.
\retun an integer, 1 for yes, and 0 for no. */
int check_ML_y_coverage(EXT_BEAM *beam);
#endif

/** Find the extent of a given contour, and save
the result as max.x, min.x, max.y, and min.y of the contour.
\param c a pointer to CONTOUR, the contour of interest.
\return an integer, always 0. */
int contour_min_max(CONTOUR *c);

#ifdef _PLUNC_HAVE_PLAN_IM_H
/** Convert the bounding box of a contour from the PATIENT coordinates
to the PIXEL coordinates.
\param h a pointer to plan_im_header, the plan_im header.
\param con a pointer to CONTOUR, the contour.
\param min_i a pointer to integer, the minimum x-index of the pixels
forming the contour's bounding box in the PIXEL coordinate system.
\param max_i a pointer to integer, the maximum x-index of the pixels
forming the contour's bounding box in the PIXEL coordinate system.
\param min_j a pointer to integer, the minimum y-index of the pixels
forming the contour's bounding box in the PIXEL coordinate system.
\param max_j a pointer to integer, the maximum y-index of the pixels
forming the contour's bounding box in the PIXEL coordinate system. */
void
con_image_bounding_box(
    plan_im_header *h,
    CONTOUR *con,
    int *min_i,
    int *max_i,
    int *min_j,
    int *max_j);

/** Set the most significant bit of the image voxel value if it, or any of
its neighbors, exceeds a threshold. The neighbors of voxel 0 in a
scan are point 1 to 7 in the same scan and the next higher scan
(point 4 is directly above point 0):
. . . .            . . . .
. 0 1 .            . 4 5 .
. 2 3 .            . 6 7 .
. . . .            . . . .
Also note that the image voxel values and the threshold are offset
by certain amount to make all of them non-negative.
\param header a pointer to plan_im_header, the plan_im header.
\param image an array of PIXELTYPE, the image data set.
\param threshold an integer, the threshold above which the voxel is
flagged. */
void
flag_above_threshold(
    plan_im_header      *header,
    PIXELTYPE           *image,
    int                 threshold);

/** Flag the image voxels that are within an anastruct.
\param header a pointer to plan_im_header, the header of the plan_im.
\param image an array of PIXELTYPE, the image data.
\param ana a pointer to ANASTRUCT, the anastruct of interest. */
void
flag_inside_ana(
    plan_im_header      *header,
    PIXELTYPE           *image,
    ANASTRUCT		*ana);

/** Set the specified bit of the image voxel value when it, or any of
its neighbors, exceeds a threshold. The neighbors of voxel 0 in a
scan are point 1 to 7 in the same scan and the next higher scan
(point 4 is directly above point 0):
. . . .            . . . .
. 0 1 .            . 4 5 .
. 2 3 .            . 6 7 .
. . . .            . . . .
\param header a pointer to plan_im_header, the plan_im header.
\param image an array of PIXELTYPE, the image data set.
\param threshold an integer, the threshold above which the pixel is
flagged.
\param bit an integer, the specified bit of the voxel value to be
set, starting from the most significant (bit=0). */
void
set_threshold_bit(
    plan_im_header      *header,
    PIXELTYPE           *image,
    int                 threshold,
    int                 bit);

/** Generate an indirect indexing array to allow a non-uniformly spaced
image to be accessed quickly as though it had uniform spacing. The resulting
array will have scale*header->slice_count entries, each entry is a fixed point
number with scale bits of fraction. The integer portion will be the index
of the desired slice, the fraction portion tells how far between this slice
and the next slice the point of interest will be. For example, if you want
to get the z position that is 1/5th of the way between the 1st and last slices
you would use the entry ptr[0.2*scale*header->slice_count], the integer
portion of that fixed point number would give the base slice and the fraction
portion would tell how far between that base slice and the next slice the
point of interest lies.
\param header a pointer to a plan_im_header, the plan_im header.
\param ptr for the indirect index array that will be created.
\param scale a float, the scaling factor for the fixed point numbers.
\return an integer, the number of elements in the indexing array created. */
int
fake_uniform_spacing(
    plan_im_header      *header,
    int                 **ptr,
    float               scale);
#endif

/** Generate a slice of 2D dose grid for a viewing plane from the 3D
dose grid.
\param grid a pointer to GRID, the storage of the 3D dose grid.
\param pat_to_view_T a 4x4 array of float, the transformation matrix
from PATIENT system to the view plane.
\param x_inc a float, the increment in x-dimension for the 2D grid.
\param y_inc a float, the increment in y-dimension for the 2D grid.
\param grid2 a pointer to GRID, the storge for the 2D grid to be
generated.
\param alloc_flag an integer, if set to FALSE, for example, in case
of locking grid, the matrix of grid2 is not populated while the other
parameters of GRID are. */
void
new_grid_slicer(
    GRID	*grid,
    float	pat_to_view_T[4][4],
    float	x_inc,
    float	y_inc,
    GRID	*grid2,
    int		alloc_flag);

/** If any dimension of a 3D dose grid has only one element, replace
it with two elements of the original value and positioned 0.05 mm
before and after the original position.
\param old_grid a pointer to GRID, the original dose grid.
\param new_grid a pointer to GRID, the new, and possibly sized-up,
dose grid.
\return an integer, always 0. */
int flesh_out_grid(GRID *old_grid, GRID *new_grid);

/** Get the current local date information in the PDATE format from
the system clock.
\param date a pointer to PDATE, output storage for the current local
date. */
void get_date(PDATE *date);

/** Get the current local time information in the PTIME format from
the system clock.
\param time a pointer to PTIME, output storage for the current local
time. */
void get_time(PTIME *time);

/** Gets a list of subdirectories (and files if dirs_only is FALSE)
in the directory "name" that end in the suffix "suffix".
\param name a pointer to the directory to look in.
\param suffix a pointer to the suffix string to match.
\param dirs_only an integer, flag that determines if files are included.
\return a list of directories (and files). */
char ** get_dirs(char *name, char *suffix, int dirs_only);

/** Find all of the files with a specified extension in a given
subdirectory.
\param name a string, the name of the subdirectory.
\param suffix a string, the specified file extension (including .),
or NULL for any extension.
\return an array of string, the file names that meet the requirement.
*/
char ** get_dirs(char *name, char *suffix);

/** Not currently used. */
char ** get_dirs2(char *name, char *suffix);

/** Free the memory used by the file name strings.
\param dir_names an array of string, the file names whose memory is to be freed.
*/
void free_dirs(char **dir_names);

/** Returns TRUE if "dir" exists and is a directory.
\param dir a string, the name of the directory to check.
\return an integer, TRUE if dir exists and is a directory. */
int is_dir(char *dir);

/** Returns TRUE if "file" exists and is a file.
\param file a string, the name of the file to check.
\return an integer, TRUE if file exists and is a file. */
int is_file(char *file);

/** Generates a CRC (Cyclical Redundancy Check) value for the input
string, used for Mosaiq .rtp files.
\param str a string, the input string to the CRC algorithm.
\return an integer, CRC value for that string. */
int rtpcrc(char *str);

/** Determine whether a file is present in the current directory.
\param file a string, the name of the file being checked for existence.
\return 1 if TRUE, or 0 if FALSE. */
int is_file(char *file);

#ifdef _PLUNC_HAVE_EXTBM_H

/** Determine the sense and offset for the collimator, gantry, and
table of a treatment unit.
\param unit a pointer to UNIT, the unit of interest.
\param c_sense a pointer to float, storage for the sense for the
collimator.
\param c_off a pointer to float, storage for the offset for the
collimator.
\param g_sense a pointer to float, storage for the sense for the
gantry.
\param g_off a pointer to float, storage for the offset for the
gantry.
\param t_sense a pointer to float, storage for the sense for the
table.
\param t_off a pointer to float, storage for the offset for the
table.
*/
void
get_machine_senses(
    UNIT        *unit,
    float       *c_sense,
    float       *c_off,
    float       *g_sense,
    float       *g_off,
    float       *t_sense,
    float       *t_off);
#endif

/** Get PLUNC version in the form of x.x.xx (major.minor.revision).
\returns a string, the PLUNC version. */
char *get_plunc_version();

/** Generate a serial number from the system time, the number of
seconds elapsed since midnight (00:00:00), January 1, 1970,
coordinated universal time, according to the system clock.
\param num a pointer to unsigned long, storage location for the
serial number. */
void get_serial_num(unsigned long *num);

/** Compute the max and min dose values in a given grid.
\param g a pointer to GRID, containing the grid values.\ The max and
min values are assigned to the corresponding data members of g. */
void grid_min_max(GRID *g);

/** Determine whether a 3D point is inside or outside an anastruct.
\param ana a pointer to ANASTRUCT, storage for the anastruct of
interest.
\param point a PNT3D, the point whose location is to be determined. */
int in_ana(ANASTRUCT *ana, PNT3D point);

/** Determine the dose scaling factors for each dose bin that will
bring the actual DVH to meet the user's desired DVH.
\param hist a pointer to DVH, the actual DVH of the current plan.
\param target a pointer to USER_DVH, the user's desired DVH.
\param index_list a pointer to DVH, whose *vol data member is used
to record the reciprocal of the user desired dose that has the same
accumulative volume for that dose bin.
\param ref_dose a float, the reference dose of the plan.
\param sensitive an integer, set to TRUE for a sensitive structure
and the values in *index_list are negated. */
void
set_index_hist2(
    DVH		*hist,
    USER_DVH	*target,
    DVH		*index_list,
    float	ref_dose,
    int		sensitive);

/** Determine whether a point is inside or outside a contour.
\param con a pointer to CONTOUR, storage for the anastruct contour
of interest.
\param point a PNT3D, the point whose location is to be determined. */
int inout(CONTOUR *con, PNT3D point);

/** Linear interpolation with mode control, that is, if mode is 1 and
the point to be interpolated is outside of the given line segment,
the point is clipped at the closest endpoint of the given segment.
\param mode an integer, 0 for normal linear interpolation without
point clipping, and 1 for clipping the interpolated point to the
closest endpoint of the given line segment.
\param x1 a float, the x-coordinate of the start point of the given
line segment.
\param x2 a float, the x-coordinate of the endpoint of the given
line segment.
\param x a float, the x-coordinate of the point to be interpolated.
\param y1 a float, the y-coordinate of the start point of the given
line segment.
\param y2 a float, the y-coordinate of the endpoint of the given
line segment.
\param fx a pointer to float, the signed ratio between the distance
from the point to be interpolated to the start point and the distance
from the endpoint to the start point.
\return a float, the y-coordinate of the interpolated point.
*/
float
interp(
    int mode,
    float x1,
    float x2,
    float x,
    float y1,
    float y2,
    float *fx);

/** Not currently used. */
void
least_squares_1D(
    float       *x,
    float       *y,
    int         dim,
    float       *output,
    int         degree);

#ifdef _PLUNC_HAVE_PLAN_IM_H

/** Determine the bounding box (xo,yo) to (xf,yf) in the given viewing plane
specified by pat_to_plan_T that will be large enough to hold a 2D resampling
of the image.
\param im a pointer to plan_im_header, the plan_im header.
\param pat_to_plan_T a LIB3D_TRANSFORM, the PATIENT to PLAN transform.
\param inv_image_T a LIB3D_TRANSFORM, the transform matrix from image
n to image 0 (used for secondary images registered to the base image).
\param xres an integer, not currently used.
\param yres an integer, not currently used.
\param xo a pointer to float, the minimum x value of the image in
the view plane.
\param yo a pointer to float, the minimum y value of the image in
the view plane.
\param xf a pointer to float, the maximum x value of the image in
the view plane.
\param yf a pointer to float, the maximum y value of the image in
the view plane. */
void
plan_im_sizer(
    plan_im_header	*im,
    LIB3D_TRANSFORM	pat_to_plan_T,
    LIB3D_TRANSFORM	inv_image_T,
    int			xres,
    int			yres,
    float		*xo,
    float		*yo,
    float		*xf,
    float		*yf);

/** Generate a 2D image for the given view plane by slicing the plan_im at
that plane given the (xo,yo) to (xf,yf) bounding box.
\param im a pointer to plan_im_header, plan_im header.
\param pat_to_plan_T a LIB3D_TRANSFORM, the PATIENT to PLAN transform.
\param inv_image_T a LIB3D_TRANSFORM, the transform matrix from image
n to image 0 (used for secondary images registered to the base image).
\param in an array of PIXELTYPE, the input plan_im.
\param out an array of PIXELTYPE, the generated 2D image slice for the
given view plane.
\param skin a pointer to ANASTRUCT, all voxels outside of this
anastruct, if not NULL, are set to WHITE.
\param xres an integer, number of pixels in x direction of the generated
2D image in the view plane.
\param yres an integer, number of pixels in y direction of the generated
2D image in the view plane.
\param xo a float, the minimum x value of the image in the view plane.
\param yo a float, the minimum y value of the image in the view plane.
\param xf a float, the maximum x value of the image in the view plane.
\param yf a float, the maximum y value of the image in the view plane.
\param z_interp an integer, if set to TRUE, the voxel values are
tri-linearly interpolated; otherwise, the interpolation is only
done bi-linearly within the nearest image slice. */
void
plan_im_slicer(
    plan_im_header	*im,
    LIB3D_TRANSFORM	pat_to_plan_T,
    LIB3D_TRANSFORM	inv_image_T,
    PIXELTYPE		*in,
    PIXELTYPE		*out,
    ANASTRUCT		*skin,
    int			xres,
    int			yres,
    float		xo,
    float		yo,
    float		xf,
    float		yf,
    int			z_interp);

#endif

/** Not currently used. */
void
print_anastruct(ANASTRUCT *a,
                int debug,
                int surface,
                int volume);

/** Print a date in a file in this format: Thursday, September 16, 2004.
\param stream a pointer to FILE, the file to write the date string.
\param date a pointer to PDATE, the date to be written to the file. */
void
print_date(FILE *stream, PDATE *date);

/** Generate triangle-strips (tstrips) from the vertices of the
neighbouring contours of an anastruct for rendering the anastruct as a
shaded-surface.
\param ana a pointer to ANASTRUCT, the anastruct to generate tstrips for.
\param tstrips a pointer to T_STRIPS, the storage for the generated
tstrips. */
int generate_tiles(ANASTRUCT *ana, T_STRIPS *tstrips);

/** Find the pathlength of a ray intersecting an anastruct. The anastruct
is treated as a stack of prisms. A prism is a three dimensional object that
can be created by translating a plane polygon (contour) along its normal by
an amount equal to the thickness of the contour.
\param obj a pointer to ANASTRUCT, the anastruct of interest.
\param p1 a PNT3D, the start point of the ray.
\param p2 a PNT3D, the end point of the ray.
\param inside a pointer to integer, TRUE if p1 is inside the anastruct.
\param nearest a pointer to float, the distance from p1 to the first
intersection with the anastruct, or 1e10 if there is no intersection.
\param distance a pointer to float, the rayline thickness, ignoring
the density of anastruct. This is the total distance along the ray that
is inside the anastruct.
\param pathlength a pointer to float, the pathlength is the rayline
thickness weighted by anastruct density. */
void
prism(ANASTRUCT   *obj,
      PNT3D       p1,
      PNT3D       p2,
      int         *inside,
      float       *nearest,
      float       *distance,
      float       *pathlength);

/** Find all intersections of a 3D line described by two points and
a contour (treated as a prism).
\param con a pointer to CONTOUR, the contour of interest.
\param p1 a PNT3D, the first point defining the line.
\param p2 a PNT3D, the second point defining the line.
\param intersects an array of points that stores the coordinates of
the intersects, up to FI_LIMIT points.
\param int_alpha an array of float, stores the signed distances from
p1 to the intersections as fraction of distance from p1 to p2, up to
FI_LIMIT intersects.
\param inside a pointer to integer, stores whether p1 is inside or
outside the contour prism.
\return an integer, the total number of intersects.
\note If the return value is greater than FI_LIMIT, not all the
intersections could be return and the calling routine needs to
cope. */
int
prism_ntrsc(
    CONTOUR     *con,
    PNT3D       p1,
    PNT3D       p2,
    float       intersects[2][FI_LIMIT],
    float       int_alpha[FI_LIMIT],
    int         *inside);

#ifdef _PLUNC_HAVE_EXTBM_H
/** Bring beam description into agreement with the \b pat_to_beam
    matrix, that is, set the beam's gantry_angle, coll_angle and
    table_angle to correspond to the angles implied by the matrix.
\param bm a pointer to EXT_BEAM, the beam of interest.
\param force an integer, if set to FALSE, no need to reconcile when
the gantry is within 0.5 degrees from vertical orientations. */
int reconcile_beam(EXT_BEAM *bm, int force);
#endif

/** Smooth a contour by deleting vertices: continue deleting vertices
    as long as the distances from the removed vertices to the new cord
    (remaining after the candidate vertex would have been deleted) are
    smaller than the criterion distance (cm).
\param criterion a float, the maximum allowed distance from the
removed vertices to the new line segment of the contour.
\param con a pointer to CONTOUR, the address of the contour.
\param alloc_flag an integer, flag whether the resultant contour
occupy the original or the new allocated memory. */
void reduce_con(float criterion, CONTOUR *con, int alloc_flag);

/** Determine the minimum distance from a point to the line segments
of a contour.
\param pt a PNT3D, the point of interest, assumed to be on the same
slice as the contour is.
\param con a pointer to CONTOUR, the contour of interest.
\return a float, the minimum distance from pt to con. */
float point_to_contour_dist(PNT3D pt, CONTOUR *con);

/** Determine the distance from a point to a line segment on the same
slice. The distance is defined as the perpendicular distance between
the line segment and the point if the perpendicular intersects the
line segment between the endpoints. Otherwise the distance is the
distance between the point and the nearest endpoint of the line
segment.
\param pt a PNT3D, the point of interest.
\param l1 a PNT3D, the start endpoint of the line segment.
\param l2 a PNT3D, the end endpoint of the line segment.
\return a float, the distance from pt to the line segment formed by
l1 and l2. */
float point_to_line_segment_dist(PNT3D pt, PNT3D l1, PNT3D l2);

/** Find the dose at a point by triple-linearly interpolating
(trilerp'ing) the dose grid.
\param grid a pointer to GRID, the dose grid to interpolate.
\param pt a pointer to PNT3D, the point to find dose for.
\return a float, the dose of the point, or 0.0 if the point is
outside the grid. */
float sample_grid(GRID *grid, PNT3D *pt);

/** Replace characters *~/?,&()-\ in a string with _ to make sure it
can be used as a legal file name. The chars to replace are OS-specific
-- this is the union of any illegal char in any OS.
\param name a string, storage for the file name. */
void safe_filename(char *name);

#ifdef _PLUNC_HAVE_PLAN_IM_H
/** Not currently used. */
int scale_image(
    int         inres,
    PIXELTYPE   *inptr,
    int         outres,
    PIXELTYPE   *outptr,
    int         min,
    float       scale,
    float       offx,
    float       offy);

/** Not currently used. */
int
scale_image_xy(
    int         inx,
    int		iny,
    PIXELTYPE   *inptr,
    int         outx,
    int         outy,
    PIXELTYPE   *outptr,
    int         min,
    float       scale,
    float       offx,
    float       offy);
#endif

/** Print \b str (if non-NULL) followed by the time elapsed in seconds
(to microsecond precision) since this function was last called.
\param str a string, the provided message to be output; don't output
if it is NULL. */
void spittime(char *str);

/** Same as \b spittime() but returns time.
\param str a string, the provided message to be output; ignore if NULL.
\return  a double, the time elapsed in seconds (to microsecond
precision) since this function was last called. */
double dspittime(char *str);

/** Return the time in seconds with microsecond precision since Epoch,
and if the provided message is not NULL, output it to the stderr
followed by the time.
\param str a string, the provided message to be output; don't output
if it is NULL.
\return a double, the time in seconds with microsecond precision
since Epoch. */
double dmytime(char *str);

/** Not currently used. */
double dcputime(char *str);

/** This function intercepts and rewrite system calls on commands \b
    cp (copy files or dirs), \b rm (delete files) and \b rd (delete
    dir) so that they will run on the deployed OS.
\param command a string, the instruction with it's parameters, eg "rm tempFile" */
int system1(char *command);

/** Same as \b system1() but presents a modal dialog that instructs
the user to "please wait" until the instruction has completed.
\param command a string, the command-line instruction. */
int system2(char *command);

/** Add a set of weighted same-size dose grids, sample-wise.
\param mcount an integer, the number of dose grids to add.
\param wts an array of float, the weights of the dose grids.
\param nmat an array of GRID, the dose grids to sum.
\param result a pointer to GRID, the storage allocated by the caller
for the resultant dose grid.
\return an integer, always 0. */
int sum_mat(int mcount, float *wts, GRID *nmat, GRID *result);

/** Add the dose of corresponding weighted points in multiple
    same-size point lists.
\param mcount an integer, the number of point lists to add.
\param wts an array of float, the weights of the point lists.
\param nvec array of CALC_POINTS, the point lists to sum.
\param result a pointer to CALC_POINTS, storage for the resultant
point list.
\return an integer, always 0. */
int sum_vec(int mcount, float *wts, CALC_POINTS *nvec, CALC_POINTS *result);

/** Sleep the thread for some time (microseconds). Sleep is highly
    OS-specific so be sure to see all the variations on this routine.
\param usec an unsigned integer, the number of microseconds to sleep.
\return an integer, always 0. */
int my_usleep(unsigned int usec);

#ifdef PLAN_WINNT
/** Sleep some number of seconds.
\param usec an unsigned integer, the number of econds to sleep.
\return an integer, always 0. */
int sleep(int sec);
#endif

/** Interpolate a new value from a sorted array: find, in a sorted
array, two adjacent elements that bound a given value x, and then
interpolate these two bounding numbers to find the linear
corresponding value of x in the dependent array.
\param mode an integer, if 0 allow extrapolation if outside the array bounds,
if 1 don't extrapolate beyond the array (clip to first or last entry).
\param vec_size an integer, the size of the data arrays, which
affects the internal searching method.
\param xvec an array of float, the array of x values in ascending
order.
\param x a float, the given value whose dependent value is to be
interpolated.
\param yvec an array of float, the array that contains the dependent
values of xvec.
\param index a pointer to integer, the index of the number in xvec
that is just smaller than x.
\param fx a pointer to float, the distance from xvec[*index] to x as
a fraction of the distance from xvec[*index] to xvec[*index + 1].
\return a float, the interpolated dependent value of x. */
float
v_interp(
    int         mode,
    int         vec_size,
    float       xvec[],
    float       x,
    float       yvec[],
    int         *index,
    float       *fx);

#ifdef _PLUNC_HAVE_EXTBM_H
/** Check and fix leaf positions of each MLC segment of the beam;
conditions are the leaf pairs outside maximum jaw openings are closed,
the first open pair has greater than 0.3cm gap, closed leaf pairs
moved to a specified location (1.9cm) -- add other conditions as
needed.
\param beam a pointer to EXT_BEAM, the beam to be validate.
\param unit a pointer to UNIT, the unit to whom the beam belongs. */
void validate_beam_segments(EXT_BEAM *beam, UNIT *unit);
#endif

/** A quick version of the standard error function erf() defined as
2/sqrt(pi) * integral from 0 to x of exp(-t*t)dt.
This is used when computing Gaussian blur functions.
\param x a float, the value you want to compute the erf() of.
\return a float, the error function value. */
float my_erf(float x);

/** Reset the Gaussian-distributed random number generator with
a given seed.
\param seed an integer, the seed for random number generation. */
void reset_rand48(int seed);

/** Gaussian-distributed random number generator.
\param sigma a float, the required sigma of the Gaussian-distributed
random number to be generated.
\return a float, the generated Gaussian-distributed random number. */
float grand48(float sigma);

#ifdef PLAN_WINNT
/** Generate uniformly distributed pseudo-random numbers between 0
and RAND_MAX.
\param seed an integer, the seed for the generator. */
void srand48(int seed);
#endif

/** Detect and fix a contour's sense (CW or CCW).
\param con a pointer to CONTOUR, the contour in question.
\param sense_flag a integer, the required contour sense, CW_SENSE or
CCW_SENSE. */
void fix_contour_sense(CONTOUR *con, int sense_flag);

/* deformation routines, ignore these for now */
#ifdef _PLUNC_HAVE_PLAN_IM_H
void init_deformation(DEFORMATION *deform);
int read_deformation(int fd, DEFORMATION *deform);
int write_deformation(int fd, DEFORMATION *deform);
int
deform_image(plan_im_header *im, short *data,
	     DEFORMATION *deform, int it, float dt);
int
deform_grid(GRID *in, GRID *out,
	    float T[4][4], float inv_T[4][4],
	    DEFORMATION *deform, int it, float dt);
int
deform_grid_delta(GRID *in, GRID *out,
		  float T[4][4], float inv_T[4][4],
		  DEFORMATION *deform, int i0, int i1, float d0, float d1);
int
deform_ana(ANASTRUCT *ana,
	   float T[4][4], float inv_T[4][4],
	   DEFORMATION *deform, int it, float dt);
int
deform_calc_points(CALC_POINTS *p,
		   float T[4][4], float inv_T[4][4],
		   DEFORMATION *deform, int it, float dt);
int deform_pnt(PNT3D p1, PNT3D *p2, DEFORMATION *deform, int it, float dt);
int
deform_pnt_delta(PNT3D p1, PNT3D *p2, DEFORMATION *deform,
		 int i0, int i1, float d0, float d1);

short sample_image(plan_im_header *im, short *d, float x, float y, float z);
int
deform_tstrip(T_STRIPS *tstrip,
	      float T[4][4], float inv_T[4][4],
	      DEFORMATION *deform, int it, float dt);
#endif

/*
#ifdef __cplusplus
}
#endif
*/

#endif

