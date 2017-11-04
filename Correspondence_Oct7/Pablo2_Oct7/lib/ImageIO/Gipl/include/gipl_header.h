/****************************************************************
 *								*
 *	GUY'S IMAGE PROCESSING LIBRARY INCLUDE FILE		*
 *								*
 ****************************************************************/

#ifndef GIPL_HEADER
#define GIPL_HEADER

/*
  The original image data files provided by the clinical partners and other
  collaborators were in either of the following formats: GE Genesis, Philips
  ACR-NEMA, Analyze. Because of their complexity and large differences
  between these formats, it was decided to use a simpler format (GIPL)
  consisting of a fixed size header (256 bytes) that precedes the raw
  image data. This format was developed as part of the Guy's Image Processing
  Libraries by the Image Processing Group at the United Medical and Dental
  Schools of Guy's and St Thomas' Hospitals, London.

  UKC has made the following modifications / additions:

  + The float field user_def1 is now used to store the value of
  the inter-slice gap. As the dimensions given in the vector
  vector pixdim include this value, the true slice thickness
  must be obtained as follows:

  sl_thickness = pixdim[3] - user_def1

  Therefore, if there is no gap between slices, then sl_thickness =
  pixdim[3].
*/

/*		IMAGE TYPE DEFINITIONS				*/

#define	IT_BINARY		1
#define	IT_CHAR			7
#define	IT_U_CHAR		8
#define	IT_SHORT		15
#define	IT_U_SHORT		16
#define	IT_U_INT		31
#define	IT_INT			32
#define	IT_FLOAT		64
#define	IT_DOUBLE		65
#define	IT_C_SHORT		144
#define	IT_C_INT		160
#define	IT_C_FLOAT		192
#define	IT_C_DOUBLE		193
#define IT_SURFACE              200
#define IT_POLYGON              201


/*	        ORIENTATION DEFINITIONS (flag1)		       */

#define UNDEFINED 0
#define UNDEFINED_PROJECTION 1
#define AP_PROJECTION 2
#define LATERAL_PROJECTION 3
#define OBLIQUE_PROJECTION 4
#define UNDEFINED_TOMO 8
#define AXIAL 9
#define CORONAL 10
#define SAGITTAL 11
#define OBLIQUE_TOMO 12


/*		FORMAT DEFINITIONS				*/

#define FORMAT_GIPL		0
#define FORMAT_GIPL_STRING	"Gipl"
#define FORMAT_MAYO		1
#define FORMAT_MAYO_STRING	"Mayo"
#define FORMAT_NM_IGE		2
#define FORMAT_NM_IGE_STRING	"Starcam"

#define GIPL_MAGIC_NUMBER	0xefffe9b0

/*		GIPL HEADER STRUCTURE DEFINITION		*/

/* Thanks to the great designer of this header gipl, file headers 
   have to be parsed entry by entry, as they are not word 
   (4 Bytes on 32 bit machines) aligned !!!

   msturm, 05/11/96
*/

typedef struct
{                                  /* Offset Size         Description        */

  unsigned short   dim[4];   /*    0    8  X,Y,Z,T image dimensions    */
  unsigned short   image_type;/*   8    2  Type of image data (below)  */
  float   pixdim[4];         /*   10   16  X,Y,Z,T pixel dimensions mm */
  char    line1[80];         /*   26   80  Patient / Text field        */
  float   matrix[20];        /*  106   80                              */ 
  char    flag1;             /*  186    1  Orientation flag (below)    */
  char    flag2;             /*  187    1                              */
  double  min;               /*  188    8  Minimum voxel value         */
  double  max;               /*  196    8  Maximum voxel value         */
  double  origin[4];         /*  204   32  X,Y,Z,T offset              */
  float   pixval_offset;     /*  236    4                              */
  float   pixval_cal;        /*  240    4                              */
  float   user_def1;         /*  244    4  Inter-slice Gap             */
  float   user_def2;         /*  248    4  User defined field          */
  unsigned int magic_number; /*  252    4                              */
  /* TOT = 256                              */
} GIPL_HD;

// Constructor implemented in giplio.cc
GIPL_HD *ipNewGiplHeader(void);

typedef struct {
  FILE	*hd,*im;
  GIPL_HD	header;
  long	image_pos;
  int	rwmode;
  short	format;
  char	format_string[32];
  char	mayo_orient;
} IFILE;

typedef struct {
  double	min_in,max_in;
  double	min_out,max_out;
} SCALE_LIMITS;

typedef struct {
  double	min,max;
  double	pixval_cal,pixval_offset;
} SCALE_FACTORS;

/*
#define		FALSE		0
#define		TRUE		1

extern	IFILE	*image_open();
extern	IFILE	*image_create();
extern	FILE	*io_error_message();

#define get_io_error_status()				io_error_status(0,0)
#define set_io_error_stream(stream)			io_error_message(0,stream,1)
#define get_io_error_stream()				io_error_message(0,0,0)
#define get_io_error_message(string)			io_error_message(string,0,-1)
*/

#endif

