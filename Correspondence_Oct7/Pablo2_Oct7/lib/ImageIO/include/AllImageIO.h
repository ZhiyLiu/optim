#ifndef ALL_IMAGE_IO_H
#define ALL_IMAGE_IO_H




/*  Class AllImageIO.

	This class provides multi-format Image I/O capability to Pablo.
	Although class RAWImageFile is in the m3d library, it is actually
	used via this class.

	There are only four usable functions in this class, besides the
	constructor and destructor.

	1.  Image3D * read(const char * filename, bool stacked, bool
			headerOnly)

		This function is used to read an image from a specified file.
		The filename should be a complete path, including the file's
		extension.  It is assumed that the extension correctly indicates
		the format of the file.  If the image is to be treated as a
		stacked (i.e. bit-plane) image after reading, then the second
		argument should be true.  This only changes the way the
		intensity maximum is set and is used to make it easier to
		adjust the image intensities for viewing stacked images.  If
		headerOnly is true, then the user should assume that no
		intensities will be read from the file, although for some file
		formats this may not be true.  All other parameters will be
		set in the image returned.  However such images should only be
		used for obtaining header values for printing; they should
		then be discarded.

	2.	bool write(const char * filename, const Image3D & image,
			int minIntensity, int maxIntensity, image_t ifmt)

		This function is used to write a regular image to disk.  It is
		not intended for the production of stacked images, which Pablo
		is not equipped to write.  The file format used will be
		determined by the extension provided in the filename.  The
		second argument is the image itself.  The next two arguments
		may indicate the intensity range to be recorded in the image;
		if absent they will be determined by the write() function.
		The final argument may be used to specify that a particular
		file format be used, in which case any existing extension will
		be ignored.  This is primarily of use when no file name
		extension is provided.  If the ifmt argument is unused and no
		file format can be determined from the file name, then a Raw3
		file will be produced.

	3.  void setImageScaling(bool scaleInput)

		This function is used to store the user's ScaleImages preference
		for use by the read() function.  A value of true indicates that
		the intensities read from the image will be mapped to the full
		range supported by Pablo.  This is the normal behavior.

	4.  void setImageShifting(bool shiftInput, int shiftAmount)

		This function may be used to specify an alternative intensity to
		be subtracted from the voxels during input remapping.  Normally,
		the minimum intensity is subtracted.  When used in cooperation
		with setImageScaling(), an image may be mapped to any allowable
		range.  If shiftInput is false, the shiftAmount is not used and
		the minimum intensity is used.
        
    Note: Images with particular modalities (presently CT and SHIFTED_CT
    images) are assigned the predefined intensity range of the modality
    before any remapping occurs, regardless of what intensity range is
    found in the image header.

        The following functions are only provided for special applications
        and should never be called from inside Pablo.

	5.  void setImageBitLength(int nbits)

        Instead of using functions setImageScaling() and setImageShifting(),
        function setImageBitLength() may be used to specify the number of
        bits to which the image should be scaled and shifted.  If the range
        of voxels is already within nbits, then no remapping will be done
        when the image is read.  If a remapping occurs, the image will be
        stored to cover the full range of [0, (nbits << 1) - 1].
        
	6.  void setImageMapActual(bool yesNo)

        Normally, CT data is considered to be in Houndsfield units, and the
        image is mapped using the full range of possible CT intensities.  In
        the case of SHIFTED_CT, for which 1024 has been added to the CT range,
        the full range is also retained.  In both cases, the file intensities
        may not fully occupy the range.  Calling this function with a value of
        true, eliminates this feature, so that whatever intensity range is
        present in the file will be mapped to the GreyValue range stored by
        class Image3D.
*/


class Image3D;
struct ImageStruct;


class AllImageIO
{
	public:

		// This must correspond to enum ImageType in ImageIO.h
		enum image_t {
			Raw3, Analyze, Gipl, Meta, PlanIM, Dicom, NoImageFormat
		};

		AllImageIO() {}
		~AllImageIO() {}

		Image3D * read(const char * filename, bool stacked = false,
			bool headerOnly = false);

		bool write(const char * filename, const Image3D & image,
			int minIntensity = 0, int maxIntensity = -1,
			image_t ifmt = NoImageFormat);

        // The type of the last image loaded or saved; or NoImageFormat,
		// if the read or write failed.
        image_t lastImageType() const { return last_format; }

		static void setImageScaling(bool scaleInput);
		static void setImageShifting(bool shiftInput, int shiftAmount = 0);
		static void setImageBitLength(int nbits);
		static void setImageMapActual(bool yesNo);

	private:

		Image3D * read_raw3(const char * filename, bool headerOnly = false);
		bool write_raw3(const char * filename, const Image3D * image,
			unsigned short min, unsigned short max);
		Image3D * convertInputImage(ImageStruct & image_struct, bool stacked,
			bool verbose = false);
		void remapImageForOutput(ImageStruct & image_struct, unsigned short min,
			unsigned short max, unsigned short * outVoxels, bool verbose = false);
		void extractOutputImageInfo(ImageStruct & image_struct, 
			unsigned short * min, unsigned short * max, const Image3D & image3d);

		AllImageIO(const AllImageIO &);					// Undefined
		AllImageIO & operator=(const AllImageIO &);		// Undefined

		static bool scale_input;
		static bool shift_input;
		static int shift;
		static int bits;
		static bool map_actual;
                static image_t last_format;
};


#endif

