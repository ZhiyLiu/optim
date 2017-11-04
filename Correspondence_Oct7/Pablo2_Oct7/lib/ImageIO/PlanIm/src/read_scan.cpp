
#include <stdio.h>
#include <stdlib.h>
#include "plan_file_io.h"
#include "plan_strings.h"
#include "plan_xdr.h"
#include "plan_sys.h"
#include "gen.h"
#include "plan_im.h"
#include "libplanio.h"
#include "zlib.h"

#define ABS(x) (((x) < 0) ? -(x) : x)
static int compress(char *buf, int nbytes);
static int uncompress(char *buf, int cbytes, int nbytes);

/*
 ---------------------------------------------------------------------------
 
   NAME
 	read_scan.c
 
   SYNOPSIS
        read_scan(fdes, scan_desc, offset, resolution)
	int fdes;
	short *scan_desc;
	int offset;
	int resolution;
 
 
   DESCRIPTION
        read a scan of the specified resolution from the file
	starting at the given offset
 
   RETURN VALUE
        0  -> ok
	-1 -> error
 
   SEE ALSO
        plan_xdr_defs.c
 
 
 ---------------------------------------------------------------------------
*/

int
read_scan_xy(int fdes, short *scan_desc, plan_im_header *h, int slice)
{
    int nbytes;
    int ret;
    int status = XDR_NO_ERROR;
    int x_dim = h->x_dim;
    int y_dim = h->y_dim;
    int offset = h->per_scan[slice].offset_ptrs;
    int cbytes = h->per_scan[slice].nbytes;
    XDR_fdes *xdr_fdes;

    xdr_fdes = xdr_ll_open_stream(fdes, XDR_READ_MODE);
    if (xdr_fdes == NULL) {
	XDR_ReportError("sock_read_scan_1 : Cannot open scan xdr_fdes");
	return(XDR_ERROR);
    }

    if (!is_socket(fdes)) {
	ret = xdr_setpos(xdr_fdes->xdrs, offset);
	if (! ret) {
	    XDR_ReportError("sock_read_scan_1 : cannot setpos");
	    status = XDR_ERROR;
	    xdr_ll_close_stream(xdr_fdes);
	    return(status);
	}
    }

    nbytes = x_dim*y_dim*sizeof(short);
    if (cbytes == 0) cbytes = nbytes;
    if (! xdrfd_get_bytes(xdr_fdes->xdrs, (caddr_t)scan_desc, cbytes)) {
	XDR_ReportError("sock_read_scan_1 : cannot read scan");
	status = XDR_ERROR;
    }

    if (cbytes != nbytes) uncompress((char *)scan_desc, cbytes, nbytes);

#ifndef PLAN_BIG_ENDIAN
    /* swap bytes - we are running on a little-endian machine */
    if (status != XDR_ERROR)
    	my_swab((char *)scan_desc, (char *)scan_desc, nbytes);
#endif
    xdr_ll_close_stream(xdr_fdes);
    return(status);
}

int
write_scan_xy(int fdes, short *scan_desc, plan_im_header *h, int slice, bool in_place)
{
    int x_dim = h->x_dim;
    int y_dim = h->y_dim;
    int offset;
    int status = XDR_NO_ERROR;
    int nbytes;
    char *tmp;
    XDR_fdes *xdr_fdes;

    h->per_scan[slice].offset_ptrs = lseek(fdes, 0, SEEK_CUR);
    offset = h->per_scan[slice].offset_ptrs;
 
    xdr_fdes = xdr_ll_open_stream(fdes, XDR_WRITE_MODE);
    if (xdr_fdes == NULL) {
	XDR_ReportError("sock_write_scan_1 : Cannot open scan xdr_fdes");
	return(-1);
    }

    if (x_dim == 0) {
	h->x_dim = h->y_dim = h->resolution;
	x_dim = y_dim = h->resolution;
    }
    nbytes = x_dim*y_dim*sizeof(short);

#ifndef PLAN_BIG_ENDIAN
    /* we are running on a little-endian machine */
    /* malloc space for swap buffer */
    tmp = (char *)malloc(nbytes);
    if (tmp == NULL) {
	perror("Cannot malloc for scan write");
	return(-1);
    }
    /* swap bytes */
    my_swab((char *)scan_desc, tmp, nbytes);
#else
    /* we are running on a big-endian machine */
    if (h->compression && ! in_place) {
	/* copy to avoid corrupting the provided buffer */
	tmp = (char *)malloc(nbytes);
	if (tmp == NULL) {
	    perror("Cannot malloc for scan write");
	    return(-1);
	}
	(void) memcpy(tmp, scan_desc, nbytes);
    }
    else
	tmp = (char *)scan_desc;
#endif

    if (h->compression) nbytes = compress(tmp, nbytes);

    if (status == XDR_NO_ERROR) {
	if (! xdrfd_put_bytes(xdr_fdes->xdrs, tmp, nbytes)) {
	    XDR_ReportError("sock_write_scan_1 : cannot read scan");
	    status = XDR_ERROR;
	}
    }

    if (tmp != (char *)scan_desc)
	free(tmp);

    xdr_ll_close_stream(xdr_fdes);
    h->per_scan[slice].nbytes = nbytes;
    return(0);
}

void
my_swab(char *src, char *dest, int num)
{   int		i;
    char	temp;

    for (i = 0; i < num; i += 2) {
	temp = src[i];
	dest[i] = src[i+1];
	dest[i+1] = temp;
    }
}

static int
compress(char *buf, int nbytes)
{   int			ret;
    unsigned char	*out;
    z_stream		strm;

    out = (unsigned char *)malloc(nbytes);
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = deflateInit(&strm, 1);
    if (ret != Z_OK) {
	return(-1);
    }
    strm.avail_in = nbytes;
    strm.next_in = (unsigned char *)buf;

    strm.avail_out = nbytes;
    strm.next_out = out;
    ret = deflate(&strm, Z_FINISH);

    nbytes = nbytes - strm.avail_out;
    memcpy(buf, out, nbytes);

    deflateEnd(&strm);
    free(out);
    return(nbytes);
}

static int
uncompress(char *buf, int cbytes, int nbytes)
{   int			ret;
    unsigned char	*out;
    z_stream		strm;

    out = (unsigned char *)calloc(nbytes, 1);
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit(&strm);
    if (ret != Z_OK) return(-1);

    strm.avail_in = cbytes;
    strm.next_in = (unsigned char *)buf;

    strm.avail_out = nbytes;
    strm.next_out = out;
    ret = inflate(&strm, Z_NO_FLUSH);

    memcpy(buf, out, nbytes);

    inflateEnd(&strm);
    free(out);
    return(nbytes);
}

