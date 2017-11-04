
#include <stdio.h>
#include "plan_file_io.h"
#include "plan_xdr.h"
#include "gen.h"
#include "extbm.h"
#include "libplanio.h"

int read_write_sar(int fdes, SAR_TABLE *sar, char mode);


int
read_sar(int fdes, SAR_TABLE *sar)
{
    return(read_write_sar(fdes, sar, XDR_READ_MODE));
}

int
write_sar(int fdes, SAR_TABLE *sar)
{
    return(read_write_sar(fdes, sar, XDR_WRITE_MODE));
}

int
read_write_sar(int fdes, SAR_TABLE *sar, char mode)
{
    int		loop, jloop;
    XDR_fdes	*xdr_fdes;
    XDR		*xdr;

    xdr_fdes = xdr_ll_open_stream(fdes, mode);
    if (xdr_fdes == NULL) {
	XDR_ReportError("Cannot open XDR_fdes");
	return(XDR_ERROR);
    }
    xdr = xdr_fdes->xdrs;

    if (!(xdr_int(xdr, &sar->version))) {
	XDR_ReportError("Cannot read sar");
	return(XDR_ERROR);
    }
    if (sar->version < 0) {
	if (!(xdr_int(xdr, &sar->unit_id))) {
	    XDR_ReportError("Cannot read unit");
	    return(XDR_ERROR);
	}
	sar->version = -sar->version;
    }
    else {
	sar->unit_id = sar->version;
	sar->version = 0;
    }

    if (!(xdr_int(xdr, &sar->type) &&
	  xdr_float(xdr, &sar->tran) &&
	  xdr_int(xdr, &sar->depth_count) &&
	  xdr_int(xdr, &sar->radius_count) &&
	  xdr_int(xdr, &sar->flatness_depth_count) &&
	  xdr_int(xdr, &sar->flatness_radius_count))) {
	XDR_ReportError("Cannot read sar");
	return(XDR_ERROR);
    }
    for (loop = 0; loop < sar->depth_count; loop++) {
	if (!(xdr_float(xdr, &sar->depth_vec[loop]))) {
	    XDR_ReportError("Cannot read/write sar");
	    return(XDR_ERROR);
	}
    }
    for (loop = 0; loop < sar->radius_count; loop++) {
	if (!(xdr_float(xdr, &sar->radius_vec[loop]))) {
	    XDR_ReportError("Cannot read/write sar");
	    return(XDR_ERROR);
	}
    }
    for (loop = 0; loop < sar->flatness_depth_count; loop++) {
	if (!(xdr_float(xdr, &sar->flatness_depth_vec[loop]))) {
	    XDR_ReportError("Cannot read/write sar");
	    return(XDR_ERROR);
	}
    }
    for (loop = 0; loop < sar->flatness_radius_count; loop++) {
	if (!(xdr_float(xdr, &sar->flatness_radius_vec[loop]))) {
	    XDR_ReportError("Cannot read/write sar");
	    return(XDR_ERROR);
	}
    }

    if (sar->version >= 5) {
	if (!xdr_int(xdr, &sar->TAR0_depth_count)) {
	    XDR_ReportError("Cannot read/write sar");
	    return(XDR_ERROR);
	}
	for (loop = 0; loop < sar->TAR0_depth_count; loop++) {
	    if (!(xdr_float(xdr, &sar->TAR0_depth_vec[loop]))) {
		XDR_ReportError("Cannot read/write sar");
		return(XDR_ERROR);
	    }
	}
    }
    else {
	sar->TAR0_depth_count = sar->depth_count;
	for (loop = 0; loop < sar->depth_count; loop++) {
	    sar->TAR0_depth_vec[loop] = sar->depth_vec[loop];
	}
    }
    for (jloop = 0; jloop < sar->depth_count; jloop++) {
	for (loop = 0; loop < sar->radius_count; loop++) {
	    if (!(xdr_float(xdr, &sar->SAR[jloop][loop]))) {
		XDR_ReportError("Cannot read/write sar");
		return(XDR_ERROR);
	    }
	}
    }
    for (loop = 0; loop < sar->TAR0_depth_count; loop++) {
	if (!(xdr_float(xdr, &sar->TAR0[loop]))) {
	    XDR_ReportError("Cannot read/write sar");
	    return(XDR_ERROR);
	}
    }
    for (jloop = 0; jloop < sar->flatness_radius_count; jloop++) {
	for (loop = 0; loop < sar->flatness_depth_count; loop++) {
	    if (!(xdr_float(xdr, &sar->flatness[jloop][loop]))) {
		XDR_ReportError("Cannot read/write sar");
		return(XDR_ERROR);
	    }
	}
    }
    for (loop = 0; loop < sar->radius_count; loop++) {
	if (!(xdr_float(xdr, &sar->Sp[loop]))) {
	    XDR_ReportError("Cannot read/write sar");
	    return(XDR_ERROR);
	}
    }
    if (!(xdr_float(xdr, &sar->Sp0))) {
	XDR_ReportError("Cannot read/write sar");
	return(XDR_ERROR);
    }
    for (loop = 0; loop < sar->radius_count; loop++) {
	if (!(xdr_float(xdr, &sar->Sp_prime[loop]))) {
	    XDR_ReportError("Cannot read/write sar");
	    return(XDR_ERROR);
	}
    }

    xdr_ll_close_stream(xdr_fdes);
    return(0);
}

