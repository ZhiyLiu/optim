
#include <stdio.h>
#include "plan_file_io.h"
#include "plan_xdr.h"
#include "gen.h"
#include "libbrachy.h"
#include "libplanio.h"

int XDR_WriteImplant(int fdes, IMPLANT_DESC *id);


int
write_implant(int fdes, IMPLANT_DESC *implant_desc)
{
    return(XDR_WriteImplant(fdes, implant_desc));
}


int
XDR_WriteImplant(int fdes, IMPLANT_DESC *id)
{
    char *cp;
    XDR_fdes *xdr_fdes;
    int status = XDR_NO_ERROR;
    unsigned int lp;

    xdr_fdes = xdr_ll_open_stream(fdes, XDR_WRITE_MODE);
    if (xdr_fdes == NULL) {
	XDR_ReportError("Cannot open XDR_fdes");
	return(XDR_ERROR);
    }

    if (! xdr_IMPLANT_DESC(xdr_fdes->xdrs, id))
    {
	XDR_ReportError("Cannot write implant desc header");
	status = XDR_ERROR;
	goto cleanup;
    }

    if (id->source_count)
    {
	lp = id->source_count;
	cp = (char *)id->source;
	if (! xdr_array(xdr_fdes->xdrs, &cp, &lp, id->source_count,
			sizeof(SOURCE), (xdrproc_t)xdr_SOURCE))
	{
	    XDR_ReportError("Cannot write implant sources");
	    status = XDR_ERROR;
	    goto cleanup;
	}
    }

    if (id->seed_count)
    {
	lp = id->seed_count;
	cp = (char *)id->seed;
	if (! xdr_array(xdr_fdes->xdrs, &cp, &lp, id->seed_count,
			sizeof(SEED), (xdrproc_t)xdr_SEED))
	{
	    XDR_ReportError("Cannot write implant seeds");
	    status = XDR_ERROR;
	    goto cleanup;
	}
    }

	
 cleanup:
    xdr_ll_close_stream(xdr_fdes);
    return(status);
}

