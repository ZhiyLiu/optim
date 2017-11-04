
#include <stdio.h>
#include "plan_file_io.h"
#include "plan_xdr.h"
#include "gen.h"
#include "extbm.h"
#include "libplanio.h"
#include "libplan.h"

/*
   ***NOTE***:
   The version variable has been moved from plan_xdr_defs.cxx's
   xdr_EXT_BEAM() to here in order to modify version 6 of the
   EXT_BEAM structure to be able to be forward compatible with
   the upcoming version 7.
   This was only done because of a mistake that put a copy of
   version 7 beams into the clinic that the current version 6
   had to be able to read.
*/

/*
 ---------------------------------------------------------------------------
 
   NAME
 	read_beam.c
 
   SYNOPSIS
        read_beam(fdes, ext_beam_desc)
	int fdes;
	EXT_BEAM *ext_beam_desc;

        write_beam(fdes, ext_beam_desc)
	int fdes;
	EXT_BEAM *ext_beam_desc;

 
   DESCRIPTION
        read or write an EXT_BEAM structure to file
 
   RETURN VALUE
        0  -> ok
	-1 -> error
 
   SEE ALSO
        plan_xdr_defs.c, read_contour.c
  
 ---------------------------------------------------------------------------
*/

int
read_beam(int fdes, EXT_BEAM *eb)
{
    XDR_fdes	*xdr_fdes;
    CONTOUR	dummy_con;

    xdr_fdes = xdr_ll_open_stream(fdes, XDR_READ_MODE);
    if (xdr_fdes == NULL) {
	XDR_ReportError("Cannot open XDR_fdes");
	return(XDR_ERROR);
    }

    if (! xdr_EXT_BEAM(xdr_fdes->xdrs, eb)) {
	XDR_ReportError("Cannot read external beam");
	return(XDR_ERROR);
    }

    xdr_ll_close_stream(xdr_fdes); 

    if (read_contour(fdes, &(eb->beam_outline)))
	return(XDR_ERROR);
    if (read_contour(fdes, &(eb->collimated_beam_outline)))
	return(XDR_ERROR);

    /* To allow growth of the arrays */
    eb->beam_outline.x =
	(float *)Realloc(eb->beam_outline.x, 500*sizeof(float));
    eb->beam_outline.y =
	(float *)Realloc(eb->beam_outline.y, 500*sizeof(float));
    eb->collimated_beam_outline.x =
	(float *)Realloc(eb->collimated_beam_outline.x, 500*sizeof(float));
    eb->collimated_beam_outline.y =
	(float *)Realloc(eb->collimated_beam_outline.y, 500*sizeof(float));

/* ***NOTE*** This version 7 contour needs to be read and thown away */
    if (eb->version == -7) read_contour(fdes, &dummy_con);

    if (eb->filter_type[0] == COMPENSATOR_FILTER) {
	if (read_anastruct(fdes, &(eb->custom_filter[0])))
	    return(XDR_ERROR);
    }
    else {
	eb->custom_filter[0].contours = NULL;
	eb->custom_filter[0].contour_count = 0;
    }
    if (eb->filter_type[1] == BOLUS_FILTER) {
	if (read_anastruct(fdes, &(eb->custom_filter[1])))
	    return(XDR_ERROR);
    }
    else {
    	eb->custom_filter[1].contours = NULL;
    	eb->custom_filter[1].contour_count = 0;
    }

    /* See if the MLC has to be converted */
    check_mlc_head(eb);

    return(0);
}

int
write_beam(int fdes, EXT_BEAM *eb)
{
    XDR_fdes *xdr_fdes;
    int status = XDR_NO_ERROR;

    xdr_fdes = xdr_ll_open_stream(fdes, XDR_WRITE_MODE);
    if (xdr_fdes == NULL) {
	XDR_ReportError("Cannot open XDR_fdes");
	return(XDR_ERROR);
    }

    if (! xdr_EXT_BEAM(xdr_fdes->xdrs, eb)) {
	XDR_ReportError("Cannot write external beam");
	status = XDR_ERROR;
	return(status);
    }

    xdr_ll_close_stream(xdr_fdes);

    if (write_contour(fdes, &(eb->beam_outline)))
	return(XDR_ERROR);

    if (write_contour(fdes, &(eb->collimated_beam_outline)))
	return(XDR_ERROR);

    if (eb->filter_type[0] == COMPENSATOR_FILTER) {
	if (write_anastruct(fdes, &(eb->custom_filter[0])))
	    return(XDR_ERROR);
    }
    if (eb->filter_type[1] == BOLUS_FILTER) {
	if (write_anastruct(fdes, &(eb->custom_filter[1])))
	    return(XDR_ERROR);
    }

    return(0);
}

int
write_old_beam(int fdes, EXT_BEAM *eb)
{
    XDR_fdes *xdr_fdes;
    int status = XDR_NO_ERROR;

    xdr_fdes = xdr_ll_open_stream(fdes, XDR_WRITE_MODE);
    if (xdr_fdes == NULL) {
	XDR_ReportError("Cannot open XDR_fdes");
	return(XDR_ERROR);
    }
    if (eb->filter_type[0] == COMPENSATOR_FILTER)
	eb->filter_type[0] = NO_FILTER;
    if (eb->filter_type[0] == BOLUS_FILTER)
	eb->filter_type[0] = NO_FILTER;

    if (! xdr_EXT_OLD_BEAM(xdr_fdes->xdrs, eb)) {
	XDR_ReportError("Cannot write external beam");
	status = XDR_ERROR;
	return(status);
    }

    xdr_ll_close_stream(xdr_fdes);

    if (write_contour(fdes, &(eb->beam_outline)))
	return(XDR_ERROR);

    if (write_contour(fdes, &(eb->collimated_beam_outline)))
	return(XDR_ERROR);

    return(0);
}

void
check_mlc_head(EXT_BEAM *eb)
{
    int		i, j, k, missing;
    float	x, Y1, Y2;
    UNIT	unit;

    /* Code to check number of leafs against the unit and adjust as needed */
    /* Note: the jaw[][[] arrays were malloc'd properly in xdr_EXT_BEAM() */
    if (get_unit(get_phys_dat_dir(), eb->unit_id, &unit) == 0) {
	/* This is to handle our Primus's going from 29 to 41 leaf pairs */
	if (eb->x_jaw_count < unit.x_jaw_count) {
	    missing = (unit.x_jaw_count - eb->x_jaw_count)/2;
	    for (j = 0; j < eb->segment_count; j++) {
		/* Y jaw */
		eb->jaw[j][unit.x_jaw_count] = eb->jaw[j][eb->x_jaw_count];
		/* replicate the last X jaw */
		k = unit.x_jaw_count - missing;
		for (i = 0; i < missing; i++) {
		    eb->jaw[j][i+k] = eb->jaw[j][eb->x_jaw_count-1];
		}
		/* Move jaws x_jaw_count-1  thru 1 */
		k = missing;
		for (i = eb->x_jaw_count-1; i > 0; i--) {
		    eb->jaw[j][i+k] = eb->jaw[j][i];
		}
		/* replicate jaw 0 */
		for (i = 0; i < missing; i++) {
		    eb->jaw[j][i+1] = eb->jaw[j][0];
		}
		/* Fix all min_w and max_w */
		for (i = 0; i < unit.x_jaw_count; i++) {
		    eb->jaw[j][i].min_w = unit.jaw[i].min_w;
		    eb->jaw[j][i].max_w = unit.jaw[i].max_w;
		}
	    }
	    eb->x_jaw_count = unit.x_jaw_count;
	    /* Close any open leaf pairs beyond the Y jaws */
	    for (j = 0; j < eb->segment_count; j++) {
		if (eb->jaw[j][0].independent < 2) break;
		Y1 = eb->jaw[j][eb->x_jaw_count].pos_1 + 0.05f;
		Y2 = eb->jaw[j][eb->x_jaw_count].pos_2 - 0.05f;
		for (i = 0; i < eb->x_jaw_count; i++) {
		    if (eb->jaw[j][i].max_w <= Y1 ||
			eb->jaw[j][i].min_w >= Y2) {
			x = 0.5f*(eb->jaw[j][i].pos_1 + eb->jaw[j][i].pos_2);
			eb->jaw[j][i].pos_1 = x;
			eb->jaw[j][i].pos_2 = x;
		    }
		}
	    }
	}
	/* This is to handle our Primus's going from 41 to 29 leaf pairs */
	/* Will lose information if missing jaws are open */
	if (eb->x_jaw_count > unit.x_jaw_count) {
	    missing = (eb->x_jaw_count - unit.x_jaw_count)/2;
	    for (j = 0; j < eb->segment_count; j++) {
		/* Populate jaws thru 1 thru x_jaw_count - 1 */
		for (i = 1; i < unit.x_jaw_count-1; i++) {
		    eb->jaw[j][i] = eb->jaw[j][i+missing];
		}
		/* move the last X jaw */
		eb->jaw[j][unit.x_jaw_count-1] = eb->jaw[j][eb->x_jaw_count-1];
		/* Y jaw */
		eb->jaw[j][unit.x_jaw_count] = eb->jaw[j][eb->x_jaw_count];
		/* Fix all min_w and max_w */
		for (i = 0; i < unit.x_jaw_count; i++) {
		    eb->jaw[j][i].min_w = unit.jaw[i].min_w;
		    eb->jaw[j][i].max_w = unit.jaw[i].max_w;
		}
	    }
	    eb->x_jaw_count = unit.x_jaw_count;
	}
    }
}
