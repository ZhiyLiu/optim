
#include <stdio.h>
#include <stdlib.h>
#include "plan_file_io.h"
#include "plan_xdr.h"
#include "gen.h"
#include "libbrachy.h"
#include "libplanio.h"

int XDR_ReadImplant(int fdes, IMPLANT_DESC *id);

int
read_implant(int fdes, IMPLANT_DESC *implant_desc)
{
    return(XDR_ReadImplant(fdes, implant_desc));
}


int
XDR_ReadImplant(int fdes, IMPLANT_DESC *id)
{
    char *cp;
    XDR_fdes *xdr_fdes;
    int status = XDR_NO_ERROR;
    unsigned int lp;

    xdr_fdes = xdr_ll_open_stream(fdes, XDR_READ_MODE);
    if (xdr_fdes == NULL) {
	XDR_ReportError("Cannot open XDR_fdes");
	return(XDR_ERROR);
    }

    if (! xdr_IMPLANT_DESC(xdr_fdes->xdrs, id)) {
	XDR_ReportError("Cannot read implant desc header");
	status = XDR_ERROR;
	goto cleanup;
    }

    if (id->source_count) {
	id->source = (SOURCE *) malloc(id->source_count * sizeof(SOURCE));
	if (id->source == NULL) {
	    XDR_ReportError("Cannot malloc implant sources");
	    status = XDR_ERROR;
	    goto cleanup;
	}

	cp = (char *)id->source;
	if (! xdr_array(xdr_fdes->xdrs, &cp, &lp,
			id->source_count, sizeof(SOURCE),
			(xdrproc_t)xdr_SOURCE)) {
	    XDR_ReportError("Cannot read implant sources");
	    status = XDR_ERROR;
	    goto cleanup;
	}
    }
    else
	id->source = (SOURCE *) NULL;

    if (id->seed_count) {
	id->seed = (SEED *) malloc(id->seed_count * sizeof(SEED));
	if (id->seed == NULL) {
	    XDR_ReportError("Cannot malloc implant seeds");
	    status = XDR_ERROR;
	    goto cleanup;
	}

	cp = (char *)id->seed;
	if (! xdr_array(xdr_fdes->xdrs, &cp, &lp,
			id->seed_count, sizeof(SEED),
			(xdrproc_t)xdr_SEED)) {
	    XDR_ReportError("Cannot read implant seeds");
	    status = XDR_ERROR;
	    goto cleanup;
	}
    }
    else
	id->seed = (SEED *) NULL;

	
 cleanup:
    xdr_ll_close_stream(xdr_fdes);
    return(status);
}

