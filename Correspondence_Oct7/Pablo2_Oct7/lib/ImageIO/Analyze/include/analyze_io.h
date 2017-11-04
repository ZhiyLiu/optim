#ifndef ANALYZE_IO_H
#define ANALYZE_IO_H


// Low-level function definitions

struct dsr * read_analyze_header(char * hdr_filename, short int datatype,
	bool & swapBytes);

bool read_analyze_image(char * img_filename, short int datatype,
	ImageStruct * im, bool swapBytes);

void swap_hdr(struct dsr *pntr);

bool write_analyze(char * hdr_name, char * img_name, char * database,
	ImageStruct * im);


#endif

