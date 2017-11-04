
#include <stdio.h>
#include <stdlib.h>
#include "plan_file_io.h"
#include "plan_xdr.h"
#include "gen.h"
#include "extbm.h"
#include "libplanio.h"

int read_write_accessory(int fdes, ACCESSORY *accessory, char mode);

/*
 ---------------------------------------------------------------------------
 
   NAME
 	read_accessory - read a accessory description from file
 
   SYNOPSIS
 
 
   DESCRIPTION
 	The hole_list space will be malloced;
 
   RETURN VALUE
 	0 if OK, -1 otherwise.
 
 ---------------------------------------------------------------------------
*/

int
read_accessory(int fdes, ACCESSORY *accessory)
{
    return(read_write_accessory(fdes, accessory, XDR_READ_MODE));
}

int
write_accessory(int fdes, ACCESSORY *accessory)
{
    return(read_write_accessory(fdes, accessory, XDR_WRITE_MODE));
}

int
read_write_accessory(int fdes, ACCESSORY *accessory, char mode)
{
    char	*cp;
    XDR_fdes    *xdr_fdes;
    XDR         *xdr;

    xdr_fdes = xdr_ll_open_stream(fdes, mode);
    if (xdr_fdes == NULL) {
	XDR_ReportError("Cannot open XDR_fdes");
	return(XDR_ERROR);
    }
    xdr = xdr_fdes->xdrs;
    cp = &(accessory->name[0]);
    if (!(xdr_int(xdr, &accessory->version) &&
	  xdr_int(xdr, &accessory->unit_id) &&
	  xdr_int(xdr, &accessory->accessory_id) &&
	  xdr_int(xdr, &accessory->type) &&
	  xdr_string(xdr, &cp, NAME_LENGTH) &&
	  xdr_float(xdr, &accessory->accessory_factor))) {
	XDR_ReportError("Cannot read/write accessory");
	return(XDR_ERROR);
    }
    xdr_ll_close_stream(xdr_fdes);
    return(0);
}

