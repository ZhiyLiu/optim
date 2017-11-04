
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#if ! defined(_MSC_VER) && ! defined(__BORLANDC__)
#include <unistd.h>
#include <ctype.h>
#endif
#include "RAWImageFile.h"
#include "DQFImage.h"
#include "zlib.h"
#include "ControlParms.h"
#include "utility.h"


//#define DEBUG

#ifndef S_ISREG
#define S_ISREG(m)  (((S_IFREG & m) != 0) ? true : false)
#endif


const int UNKNOWN_GREY_VALUE = 99999;

const int MAX_UNCOMPRESS_TRIES = 3;
const float COMPRESS_FACTOR = 1.1f;
const float UNCOMPRESS_INITIAL_FACTOR = 2.0f;
const float UNCOMPRESS_GROWTH_FACTOR = 1.25f;
const int CT_MIN_INTENS = -1024;
const int CT_MAX_INTENS = 3071;


enum order_t { ORDERERR = 0, NORMAL = 1, BYTEREV,
		PAIRREV, REVERSE, NATIVE, RETAIN };

order_t Order()
{
	union {
	    unsigned int	i;
	    unsigned char	c[4];
	}	word;

	order_t	Machine_Byte_Order;

	for (register int j = 0; j < 4; j++) word.c[j] = j + 1;
		if (word.i == 000100401404) Machine_Byte_Order = NORMAL;
		else if (word.i == 00400601001) Machine_Byte_Order = REVERSE;
		else if (word.i == 00200202003) Machine_Byte_Order = BYTEREV;
		else Machine_Byte_Order = PAIRREV;
	return Machine_Byte_Order;
}

int cpuByteOrder()	// Global function
{
	switch (Order()) {
		case NORMAL:	return 0;

		case REVERSE:	return 1;

		default:	return -1;
	}
}

order_t preferredOrder()
{
	switch (globalControl->readInt(ByteOrder)) {
		case 0:		return REVERSE;

		case 1:		return NATIVE;

		case 2:		return NORMAL;

		case 3:		return RETAIN;

		default:	return ORDERERR;
	}
}

#ifdef DEBUG
const char * printOrder(order_t bo)
{
	switch (bo) {
		case REVERSE:	return "REVERSE";

		case NATIVE:	return "NATIVE";

		case NORMAL:	return "NORMAL";

		case RETAIN:	return "RETAIN";

		default:		return "ORDERERR";
	}
}
#endif

/*  Determine the byte order for writing a file.   Argument current_byte_order
    must be either REVERSE or NORMAL.  The last argument will be set to indicate
    whether or not a byte swap is needed.  The return value will be the new byte
    order, either NORMAL or REVERSE, or ORDERERR if a problem occurs.
*/
order_t newOrder(order_t current_byte_order, order_t desired_byte_order,
				 bool & swapBytes)
{
	order_t bo;

	if (current_byte_order != REVERSE && current_byte_order != NORMAL)
		return ORDERERR;

	swapBytes = false;

	switch (desired_byte_order) {
		case NATIVE:
			bo = Order();	// Only NORMAL or REVERSE are allowed
			desired_byte_order = bo;
			if (bo == NORMAL) {
				if (current_byte_order == REVERSE)
					swapBytes = true;
				break;
			}
			else if (bo == REVERSE) {
				if (current_byte_order == NORMAL)
					swapBytes = true;
				break;
			}
			else {
				std::cout << "Machine has an unsupported byte order" << std::endl;
				return ORDERERR;
			}

		case RETAIN:
			desired_byte_order = current_byte_order;
			break;

		case REVERSE:
			desired_byte_order = REVERSE;
			if (current_byte_order != REVERSE)
				swapBytes = true;
			break;

		case NORMAL:
			desired_byte_order = NORMAL;
			if (current_byte_order != NORMAL)
				swapBytes = true;
			break;

		default:
			return ORDERERR;
	}
	return desired_byte_order;
}

void swapLong(unsigned long & word)
{
	unsigned char * w = (unsigned char *) &word;
	unsigned char c;

	c = w[0];
	w[0] = w[3];
	w[3] = c;
	c = w[1];
	w[1] = w[2];
	w[2] = c;
}

void swapShortArray(int size, short * input, short * output)
{
    register int i;
	register GreyValue val0;

	for (i = 0; i < size; i++) {
		val0 = ((GreyValue *) input)[i];
		val0 = ((val0 >> 8) & 0x00FF) | (val0 << 8);
		output[i] = (short) val0;
	}
}


using namespace std;


bool RAWImageFile::scale_input = true;
bool RAWImageFile::shift_input = true;  // In Pablo, this never changes 
int RAWImageFile::shift = 0;
int RAWImageFile::bits = 0;
bool RAWImageFile::map_actual = false;

const char unknown[] = "Image modality is unknown\n";

void RAWImageFile::setImageScaling(bool scaleInput)
{
	scale_input = scaleInput;
    bits = 0;
}

void RAWImageFile::setImageShifting(bool shiftInput, int shiftAmount)
{
	shift_input = shiftInput;
	shift = shiftAmount;
    bits = 0;
}

void RAWImageFile::setImageBitLength(int nbits)
{
    bits = nbits;
    shift_input = false;
    shift = 0;
    scale_input = false;
}

void RAWImageFile::setImageMapActual(bool yesNo)
{
    map_actual = yesNo;
}

void writeExtendedHeader(FILE * fp, ExtendedHeader * eh)
{
	fprintf(fp, "1, ");		// Version number
	if (eh->haveModality) {
		if (eh->modality == CT)
			fprintf(fp, "CT,");
		else if (eh->modality == SHIFTED_CT)
			fprintf(fp, "SHIFTED_CT,");
		else if (eh->modality == MRI)
			fprintf(fp, "MRI,");
		else if (eh->modality == DQF)
			fprintf(fp, "DQF,");
	}
	else
		fprintf(fp, ",");
	fprintf(fp, "\n");

	// For some modalities, a second line is present
	if (eh->modality == DQF) {
		// Get next field, depending on the modality observed
		fprintf(fp, "%d %d %d %d %lf %d %d %d %d\n", eh->sWidth, eh->xOrigin, 
			eh->yOrigin, eh->zOrigin, eh->scale, eh->shift, eh->separateXDim,
			eh->separateYDim, eh->separateZDim);
	}
}

int readExtendedHeader(FILE * fp, ExtendedHeader * eh)
{
    char str[64];	// Entire line, including newline
	char * token;
	char * p;
	int i;

	i = 0;
	do {
		fscanf(fp, "%c", str + i);
	} while (str[i++] != '\n');
	str[--i] = '\0';	// Replace newline

	// Get first field, the record type.
	// This is presently enforced, although unused.
	// It exists to facilitate future reviions of the format.
	token = strtok(str, ",");
	if (token == NULL)
		return false;	// Invalid format
	if (token[0] != '\0') {
		p = token;
		while (*p == ' ')	// Skip leading spaces
			p++;
		if (*p != '1')		// Record type of '1'
			return false;	// Invalid format
	}
	else
		return false;	// Invalid format

	token = strtok(NULL, ",");
	if (token == NULL)
		return true;

	if (token[0] == '\0') {
		eh->modality = UNKNOWN_MODALITY;
		eh->haveModality = false;
		if (globalVerbosity >= 1)
			cout << unknown;
	}
	else {
		// Process the image's modality
		char * m_str = token;
		while (*m_str == ' ')	// Skip leading spaces
			m_str++;
		p = m_str;
		while (*p) {	// Convert string to upper case
			if (*p == ' ') {
				// Modality is a single word - discard final space(s)
				*p = '\0';
				break;
			}
			*p = toupper(*p);
			p++;
		}
		if (0 == strcmp(m_str, "CT")) {
			eh->modality = CT;
			eh->haveModality = true;
			if (globalVerbosity >= 1)
				cout << "Image modality is CT\n";
		}
		else if (0 == strcmp(m_str, "SHIFTED_CT")) {
			eh->modality = SHIFTED_CT;
			eh->haveModality = true;
			if (globalVerbosity >= 1)
				cout << "Image modality is SHIFTED_CT\n";
		}
		else if (0 == strcmp(m_str, "MRI")) {
			eh->modality = MRI;
			eh->haveModality = true;
			if (globalVerbosity >= 1)
				cout << "Image modality is MRI\n";
		}
		else if (0 == strcmp(m_str, "DQF")) {
			eh->modality = DQF;
			eh->haveModality = true;
			if (globalVerbosity >= 1)
				cout << "Image contains DQF data\n";
		}
		else {
			eh->modality = UNKNOWN_MODALITY;
			eh->haveModality = false;
			if (globalVerbosity >= 1)
				cout << unknown;
		}
	}

	// Get next field
/*	token = strtok(NULL, ",");
	if (token == NULL)
		return true;
	if (token[0] == '\0') {
	}
	else {
		cout << "token is " << token << endl;
	}
*/

	// For some modalities, a second line is present
	if (eh->modality == DQF) {
		// Get next field, depending on the modality observed
		if (9 != fscanf(fp, "%d %d %d %d %lf %d %d %d %d", &eh->sWidth,
				&eh->xOrigin, &eh->yOrigin, &eh->zOrigin, &eh->scale, &eh->shift,
				&eh->separateXDim, &eh->separateYDim, &eh->separateZDim))
			return false;	// Invalid format
	}
	return true;
}

Image3D * RAWImageFile::read(const char * filename, bool headerOnly)
{
    Image3D * image;
	ExtendedHeader * eh;

	bool compressed;
	bool swapBytes;
	order_t	file_byte_order;

	bool failed;
	bool convertCompression, convertFormat;
	bool hasComments;
    bool xFlip, yFlip, zFlip;

    char lsbStr[4];
    int xDim, yDim, zDim;
    int size;
    bool findMin, findMax;
    int minIntensity, maxIntensity;

    GreyValue * voxelArray;
	short * tempArray;	// Not GreyValue: intensities in files can be negative
	char c;

    double xSpacing, ySpacing, zSpacing;
    double xOrig, yOrig, zOrig;
    Vector3D worldOriginPos;

    FILE * fp;
    struct stat buf;


    if (stat(filename, &buf) != 0 || S_ISREG(buf.st_mode) == false) {
#ifdef BINARY
		cout << "WARNING: ";
#endif
		cout << filename << " is not a valid image file" << endl;
		return NULL;
    }
    fp = fopen(filename, "rb");
    if (fp == NULL) {
		cout << "Could not open file " << filename << endl;
        return NULL;
	}
    hasComments = skipComments(fp);

    // Determine file byte format LSB or MSB first
    fscanf(fp, "%3s%c", lsbStr, &c);
	if (c == '\n')
		eh = new ExtendedHeader;
	else
		eh = NULL;
    lsbStr[3] = '\0';
    if (strcmp(&(lsbStr[0]), "lsb") == 0) {
		file_byte_order = REVERSE;
#ifdef DEBUG
		cout << "File byte order is REVERSE (lsb)\n";
#endif
	}
	else if (strcmp(&(lsbStr[0]), "msb") == 0) {
		file_byte_order = NORMAL;
#ifdef DEBUG
		cout << "File byte order is NORMAL (msb)\n";
#endif
	}
    else {  // Not a raw3 file
#ifdef BINARY
		cout << "WARNING: ";
#endif
		cout << filename << " is not a valid image file" << endl;
        return NULL;
    }
	(void) newOrder(file_byte_order, NATIVE, swapBytes);

	if (eh != NULL) {
		// Read extended header information: a series of space-separated
		// strings terminated by a space-separated semicolon.
		if (! readExtendedHeader(fp, eh)) {
#ifdef BINARY
			cout << "WARNING: " << filename << " is not a valid image file"
				<< ": failed to parse extended header" << endl;
#else
			cout << filename << " is not a valid image file" << endl;
#endif
			return NULL;
		}
	}

    fscanf(fp, "%d %d %d ", &xDim, &yDim, &zDim);
    // If any of our dimensions are negative, we make them positive,
	// and flip the image (below) in the indicated direction.  The
	// result of a flip is that, for the specified axis, the model
	// space origin and bound will be exchanged with respect to the
	// image (and thus world) coordinate systems of that axis.  However,
	// since the model space is considered to always be right handed,
	// this effectively changes the handedness of the world and image
	// coordinate systems.
    if (xDim < 0) {
        xFlip = true;
        xDim = -xDim;
		if (globalVerbosity >= 1)
			cout << "Input image will be flipped along the X axis" << endl;
    }
	else
		xFlip = false;

    if (yDim < 0) {
        yFlip = true;
        yDim = -yDim;
		if (globalVerbosity >= 1)
			cout << "Input image will be flipped along the Y axis" << endl;
    }
	else
		yFlip = false;

    if (zDim < 0) {
        zFlip = true;
        zDim = -zDim;
		if (globalVerbosity >= 1)
			cout << "Input image will be flipped along the Z axis" << endl;
    }
	else
		zFlip = false;

    size = xDim * yDim * zDim;
#ifdef DEBUG
	cout << "Image dimensions: " << xDim << " x " << yDim << " x " << zDim << '\n';
#endif
    // Read in the world coordinate transformation parameters
    fscanf(fp, "%lg %lg %lg %lg %lg %lg ",
        &xSpacing, &ySpacing, &zSpacing, &xOrig, &yOrig, &zOrig);
#ifdef DEBUG
	cout << "Image spacings: " << xSpacing << ", " << ySpacing << ", " << zSpacing << '\n';
	cout << "Image origins: " << xOrig << ", " << yOrig << ", " << zOrig << '\n';
#endif

    // Determine the minimum and maximum values from the file, or find
    // them later if they are unknown (signified by UNKOWN_PIXVAL).
	// Note that no check is made to assure that values provided are
	// actually correct.
    fscanf(fp, "%d %d", &minIntensity, &maxIntensity);
	if (globalVerbosity > 1)
		cout << "Image file's intensity range = [" << minIntensity << ", "
			<< maxIntensity << "]\n";
    if(minIntensity == UNKNOWN_GREY_VALUE)
        findMin = true;
    else
        findMin = false;
    if(maxIntensity == UNKNOWN_GREY_VALUE)
        findMax = true;
    else
        findMax = false;

	compressed = false;
	if (eh == NULL) {
		// Old format (unextended header)
		fscanf(fp, "%c", &c);
		if (c == '\t')	// Indicator that data is compressed
			compressed = true;
	}
	else {
		// New format with extended header
		c = ' ';
		while (c == ' ')
			fscanf(fp, "%c", &c);
		if (c == 'Z')	// Indicator that data is compressed
			compressed = true;
		else if (c != 'u')
			cout << "Error: invalid compression indicator\n";
		char c2;
		fscanf(fp, "%c%c", &c, &c2);
		if (c != '\n' || c2 != '\n')
			cout << "Error: header incorrectly terminated\n";
	}

    worldOriginPos.set(xOrig, yOrig, zOrig);

	if (headerOnly)
		voxelArray = NULL;
	else {

		// Read in the binary data
		if (compressed) {
			int len, s;
			int count, ret;

			// Uncompress the voxels using gzip
			fscanf(fp, "%d ", &len);
			char * compressedVoxels = new char [len];
			s = safeFread(compressedVoxels, sizeof(char) * len, fp);
			if (s != (sizeof(char) * len)) {
				cout << filename << " is truncated or could not be read" << endl;
				(void) fclose(fp);
				return NULL;
			}
			s = (int) (UNCOMPRESS_INITIAL_FACTOR*size);
			tempArray = NULL;
			count = 0;
			while (count < MAX_UNCOMPRESS_TRIES) {
				// With a value of 2 for UNCOMPRESS_INITIAL_FACTOR, testing showed
				// only one pass through this loop was required, so the loop probably is
				// not necessary.  A value of less than 2 was invariably too small.  The
				// zlib documentation, however, does not give any explanation of this
				// phenomenon.
				delete [] tempArray;
				tempArray = new short[s];
#ifdef BINARY
				if (globalVerbosity > 0)
#endif
					cout << "Uncompressing image .." << flush;
				ret = uncompress((Bytef *) tempArray, (uLongf *) &s, (const Bytef *) compressedVoxels,
					(uLongf) len);
#ifdef BINARY
				if (globalVerbosity > 0)
#endif
					cout << ".. done" << endl;
				if (ret != Z_OK) {
					if (ret == Z_DATA_ERROR)
						break;
					if (s > 0)	// Z_MEM_ERROR or Z_BUF_ERROR
						s = (int) (UNCOMPRESS_GROWTH_FACTOR*s);
				}
				else {
					unsigned long file_crc;
					unsigned long crc = 0;
					crc = crc32(crc, (const Bytef *)tempArray, s);
					ret = safeFread(&file_crc, sizeof(unsigned long) * 1, fp);
#ifdef DEBUG
					cout << "CRC  on " << s << " bytes= " << crc << endl;
					cout << "CRC in file = " << file_crc << endl;
#endif
					if (ret != (sizeof(unsigned long)*1))
						cout << "Error reading CRC from input file: the image may be corrupted" << endl;
					else {
						if (swapBytes)
							swapLong(file_crc);
#ifdef DEBUG
						cout << "CRC in file after possible byte swap = " << file_crc << endl;
#endif
						if (crc != file_crc)
							cout << "Mismatched CRC encountered: the image may be corrupted" << endl;
					}

					ret = Z_OK;
					break;
				}
				count++;
			}
			if (ret != Z_OK) {
#ifdef BINARY
				cout << "WARNING: " << filename << ": ";
#endif
				if (ret == Z_MEM_ERROR)
					cout << "Uncompress did not have enough memory" << endl;
				else if (ret == Z_BUF_ERROR)
					cout << "Uncompress did not have room in its output buffer" << endl;
				else {	// Z_DATA_ERROR
					s = -1;
					cout << "Uncompress received corrupted data" << endl;
				}
			}

			delete [] compressedVoxels;
			if (s != size*sizeof(GreyValue)) {
#ifdef BINARY
				cout << "WARNING: " << filename << ": ";
#endif
				cout << "Image did not uncompress properly: received " << s << " bytes; expecting "
					<< size << " bytes" << endl;
			}
		}
		else {
			// Read uncompressed voxels
			tempArray = new short[size];
			safeFread(&(tempArray[0]), sizeof(short) * size, fp);
		}
		(void) fclose(fp);

		// Determine whether or not the file must be (un)compressed and rewritten
		convertCompression = false;
		if (globalControl->readBool(ConvertImages) == true) {
			if (compressed) {
				if (globalControl->readBool(CompressImages) == false)
					convertCompression = true;
			}
			else {
				if (globalControl->readBool(CompressImages) == true)
					convertCompression = true;
			}
		}

		// Determine whether or not the file's format must be converted
		convertFormat = false;
		if (globalControl->readBool(ConvertImageFormat) == true) {
			if (eh == NULL) {
				if (globalControl->readInt(ImageFormat) == 1)
					convertFormat = true;
			}
			else {
				if (globalControl->readInt(ImageFormat) == 0)
					convertFormat = true;
			}
		}

		// Allocate space for the voxel array
		voxelArray = new GreyValue[size];

		// Convert old file to be compressed/uncompressed depending on user's preference
		failed = false;
		if (convertCompression || convertFormat) {
			bool swap;

#ifdef DEBUG
			cout << "Converting image file\n";
#endif
			char * old_filename;
			old_filename = new char[5 + strlen(filename)];
			strcpy(old_filename, filename);
			strcat(old_filename, ".tmp");
			if (rename(filename, old_filename) != 0) {
				cout << "Error: could not write file " << old_filename << '\n';
				failed = true;
			}
			else {
				order_t ord = newOrder(file_byte_order, preferredOrder(), swap);
#ifdef DEBUG
				cout << "New byte order = " << printOrder(ord) << '\n';
#endif
				if (ord == ORDERERR)
					failed = true;
				else {
					if (swap && swapBytes) {
						swapShortArray(size, tempArray, tempArray);
						swap = false;
						swapBytes = false;
					}
					if (swap)
						// Temporarily use voxelArray for storing byte-swap voxels
						swapShortArray(size, tempArray, (short *) voxelArray);

					if (convertFormat) {
						if (eh == NULL) {
							eh = new ExtendedHeader;
							// Make an educated guess as to the modality
							if (minIntensity < 0 && minIntensity >= CT_MIN_INTENS &&
								maxIntensity > 0 && maxIntensity <= CT_MAX_INTENS)
							{
								eh->haveModality = true;
								eh->modality = CT;
								cout << "Warning: Assuming the modality of file\n    " << filename << '\n'
									<< "    is CT; if incorrect, manual editing of the file will be needed.\n";
							}
							else if (minIntensity == 0 &&
								maxIntensity > 0 && maxIntensity <= CT_MAX_INTENS - CT_MIN_INTENS)
							{
								eh->haveModality = true;
								eh->modality = SHIFTED_CT;
								cout << "Warning: Assuming the modality of file\n    " << filename << '\n'
									<< "    is SHIFTED_CT; if incorrect, manual editing of the file will be needed.\n";
							}
							else {
								eh->haveModality = false;
								eh->modality = UNKNOWN_MODALITY;
								if (globalVerbosity >= 1)
									cout << unknown;
							}
						}
						else {
							delete eh;
							eh = NULL;
						}
					}

					/*	The logic below may look wrong, but is not.  It is based on the
						following logic table.

						convertImages	compressImages	compressed	convertCompression	compress
						-------------	--------------	----------	------------------	--------
							  1				   1			0				1			    1
							  1				   0			0				0			    0
							  0				   1			0				0			    0
							  0				   0			0				0			    0

							  1				   1			1				0			    1
							  1				   0			1				1			    0
							  0				   1			1				0			    1
							  0				   0			1				0			    1
					*/
					bool compress = convertCompression;
					if (convertCompression) {
						if (compressed)
							compress = ! compress;
					}

					if (write(filename, xDim, yDim, zDim, xSpacing, ySpacing, zSpacing,
						xOrig, yOrig, zOrig, minIntensity, maxIntensity,
						 (swap ? (short *) voxelArray : tempArray),
						ord, compress, eh, (hasComments ? old_filename : NULL)))
					{
						if (unlink(old_filename) != 0)
							cout << "Could not delete file " << old_filename << endl;
					}
					else
						failed = true;
				}
			}
			delete [] old_filename;
			if (! failed) {
				cout << "Converted file " << filename << "\n    to ";
				if (convertCompression) {
					cout << (compressed ? "un" : "") << "compressed";
					if (convertFormat)
						cout << " and ";
				}
				if (convertFormat)
					cout << (eh == NULL ? "old" : "new") << " format";
				cout  << endl;
			}
		}
		if (failed) {
			cout << "Error: could not convert file " << filename << "\n    to ";
			if (convertCompression) {
				cout << (compressed ? "un" : "") << "compressed";
				if (convertFormat)
					cout << " and ";
			}
			if (convertFormat)
				cout << (eh == NULL ? "old" : "new") << " format";
			cout << endl;
		}
	}	// ! headerOnly

#ifdef DEBUG
	if (swapBytes)
		cout << "Bytes will be swapped\n";
	else
		cout << "Bytes will not be swapped\n";
#endif

    int xIndex, yIndex, zIndex;
    register int i;

    // Find the min or the max, as needed.  Swap bytes, as needed,
	// if either the min or max is calculated.
    if (eh && eh->haveModality && ! map_actual) {
		if (eh->modality == CT) {
			minIntensity = CT_MIN_INTENS;
			maxIntensity = CT_MAX_INTENS;
		}
		else if (eh->modality == SHIFTED_CT) {
			minIntensity = 0;
			maxIntensity = CT_MAX_INTENS - CT_MIN_INTENS;
		}
		// Presently MRI values are not fixed to a specific range
		//else if (eh->modality == MRI) {
		//}
	}
	else {
		if (! headerOnly && (findMin || findMax)) {
#ifdef DEBUG
			if (swapBytes)
				cout << "Swapping bytes during range computation\n";
#endif
			calc_rangeWithSwap(findMin, minIntensity, findMax, maxIntensity, swapBytes,
				tempArray, size);
			swapBytes = false;
		}
	}

	if (globalVerbosity >= 1) {
		if (minIntensity == UNKNOWN_GREY_VALUE)
			cout << "Voxel intensity range is unknown\n";
		else if (! eh || eh->modality != DQF)
			cout << "Voxel intensity range = [" << minIntensity << ", "
				<< maxIntensity << "]\n";
	}

    double scale;
	int diff;

	// At this point, the data's intensity range will be [minIntensity, maxIntensity],
	// unless the image was CT (when a fixed range is used), or the image's header
	// contained incorrect values for the range.
	if (! headerOnly) {
        bool scaling;
        bool shifting;

        if (bits > 0) {
            int max = (1 << bits) - 1;
            // Compare current and future ranges
            if (maxIntensity - minIntensity > max) {
                scaling = true;
                diff = minIntensity;
                shifting = (diff != 0);
                scale = max/(double) (maxIntensity - minIntensity);
            }
            else {
                // Only need possibly to shift
                scaling = false;
                scale = 1.0;
                if (maxIntensity > max) {
                    diff = max - maxIntensity;
                    shifting = true;
                }
                else if (minIntensity < 0) {
                    diff = minIntensity;
                    shifting = true;
                    if (diff == 0)
                        shifting = false;
                }
				else {
					shifting = false;
					diff = 0;
				}
            }
            if (globalVerbosity >= 1)
                cout << "Voxels will be mapped to fall within [0, " << max << "], " << bits << " bits\n";
        }
        else {
            if (scale_input)
			    scale = (double) MAX_GREY_VALUE / (maxIntensity - minIntensity);
		    else
			    scale = 1.0;

		    if (! scale_input && shift_input && shift != 0)
			    diff = shift;
		    else
			    diff = minIntensity;

		    if (globalVerbosity >= 1) {
				if (eh && (eh->modality == DQF || eh->modality == SHIFTED_CT))
				    cout << "Voxels will be read without remapping\n";
			    else if (scale_input)
				    cout << "Voxels will be mapped to [0, " << MAX_GREY_VALUE << "]\n";
			    else if (shift_input)
				    cout << "Voxels will be mapped to [0, " << maxIntensity - diff << "]\n";
			    else
				    cout << "Voxels will be read without remapping\n";
		    }

		    scaling = scale_input;
		    shifting = shift_input;
		    if (scale == 1.0)
			    scaling = false;
		    if (diff == 0)
			    shifting = false;
        }

#ifdef DEBUG
		cout << "scale: " << scale << ",  scale_input: " << scale_input
			<< ",  shift: " << shift << ",  shift_input: " << shift_input << endl;
#endif

		// Fit values into the GreyValue range; flip axis orientations, as necessary
#ifdef DEBUG
		if (swapBytes)
			cout << "Swapping bytes during intensity remapping\n";
#endif

		bool adjust = scaling || shifting;
		i = 0;
		for (int slicenum = 0; slicenum < zDim; slicenum++)
		{
			for (int rownum = 0; rownum < yDim; rownum++)
			{
				for (int colnum = 0; colnum < xDim; colnum++)
				{
					register GreyValue finalpixval;
					register int xlatepixval;
					register GreyValue val0;

					val0 = (GreyValue) tempArray[i];
					if (swapBytes)
						val0 = ((val0 >> 8) & 0x00FF) | (val0 << 8);

					if (adjust) {
						xlatepixval = ((short) val0) - diff;	// short cannot be replaced by int here
						finalpixval = (GreyValue) (scale*xlatepixval + 0.5);
					}
					else
						finalpixval = val0;

					if (xFlip)
						xIndex = xDim - colnum - 1;
					else
						xIndex = colnum;

					if (yFlip)
						yIndex = yDim - rownum - 1;
					else
						yIndex = rownum;

					if (zFlip)
						zIndex = zDim - slicenum - 1;
					else
						zIndex = slicenum;

					voxelArray[xIndex + xDim * (yIndex + yDim * zIndex)] = finalpixval;
					i++;
				}
			}
		}

		delete [] tempArray;

		if (adjust) {
		    minIntensity = (int) ((minIntensity - diff)*scale + 0.5);
		    maxIntensity = (int) ((maxIntensity - diff)*scale + 0.5);
		}
		// Now the data's current range is [minIntensity, maxIntensity]

	}	// ! headerOnly

	// Construct the image object
	if (eh && eh->modality == DQF) {
		image = new DQFImage(xDim/eh->sWidth, yDim, zDim, eh->sWidth,
			eh->separateXDim, eh->separateYDim, eh->separateZDim);
		((DQFImage *) image)->setVoxels(voxelArray, xDim/eh->sWidth, yDim, zDim,
			eh->sWidth);
		((DQFImage *) image)->setMapping(eh->scale, eh->shift);
	}
	else {
		image = new Image3D;
		image->setVoxels(voxelArray, xDim, yDim, zDim);
	}

	if (maxIntensity > MAX_GREY_VALUE || minIntensity < MIN_GREY_VALUE) {
	    cout << "Warning: this image may have loaded incorrectly.  It has an apparent intensity\n";
	    cout << "    range of [" << minIntensity << ", " << maxIntensity
		<< "], after the voxels were remapped.  You may want to\n";
	    cout << "    change the program's image input parameters to shift or rescale the image.\n";
		if (maxIntensity > MAX_GREY_VALUE) {
			cout << "Adjusting intensity maximum from " << maxIntensity << " to "
				<< MAX_GREY_VALUE << endl;
			maxIntensity = MAX_GREY_VALUE;
		}
		if (minIntensity < MIN_GREY_VALUE) {
			cout << "Adjusting intensity minimum from " << minIntensity << " to "
				<< MIN_GREY_VALUE << endl;
			minIntensity = MIN_GREY_VALUE;
		}
	}
    image->setRange((GreyValue) minIntensity, (GreyValue) maxIntensity);    // Range of voxels in image
#ifdef DEBUG
	cout << "setRange: min = " << minIntensity << ",  maxIntensity = "
		<< maxIntensity << endl;
#endif
    image->setSpacingAndOrigin(xSpacing, ySpacing, zSpacing, &worldOriginPos);
    if (eh)
		image->setModality(eh->modality);

    if (! headerOnly) {
		image->setIntensityMapping(scale, diff);
#ifdef DEBUG
		cout << "setIntensityMapping: scale = " << scale << ",  diff = "
			<< diff << endl;
#endif
	}
	if (eh && eh->modality == DQF)
		((DQFImage *) image)->roi(eh->xOrigin, eh->yOrigin, eh->zOrigin,
			eh->xOrigin + xDim/eh->sWidth - 1, eh->yOrigin + yDim - 1,
			eh->zOrigin + zDim -1);

	if (eh) {
		delete eh;
	}

    return image;
}

// Write an image to the specified file.
bool RAWImageFile::write(const char * filename, const Image3D & image,
	GreyValue min, GreyValue max)       // AGG: These args should be ints
{
    int xDim, yDim, zDim;
    double xSpacing, ySpacing, zSpacing;
    Vector3D originPos;
    GreyValue * voxels;
    short * tempArray;
	ExtendedHeader eh;
	ExtendedHeader * peh;
	bool compress;


    xDim = image.getXDim();
    yDim = image.getYDim();
    zDim = image.getZDim();

    xSpacing = image.getXSpacing();
    ySpacing = image.getYSpacing();
    zSpacing = image.getZSpacing();

    originPos = image.getWorldOrigin();

    voxels = image.getVoxels();

	bool swapBytes;
	order_t ord = newOrder(Order(), preferredOrder(), swapBytes);
#ifdef DEBUG
	cout << "New byte order = " << printOrder(ord) << '\n';
#endif
	if (ord == ORDERERR)
		return false;	// Should never happen

	if (swapBytes) {
		int s;

#ifdef DEBUG
		cout << "Swapping bytes\n";
#endif
		s = xDim*yDim*zDim;
		tempArray = new short[s];
		swapShortArray(s, (short *) voxels, tempArray);
	}
	else
		tempArray = (short *) voxels;

	eh.modality = image.modality();
	if (eh.modality == UNKNOWN_MODALITY)
		eh.haveModality = false;
	else
		eh.haveModality = true;
	compress = globalControl->readBool(CompressImages);
	if (eh.modality == DQF) {
		Vector3D origin = ((DQFImage &) image).getROIOrigin(); 
		eh.xOrigin = (int) origin.getX();
		eh.yOrigin = (int) origin.getY();
		eh.zOrigin = (int) origin.getZ();
		eh.sWidth = ((DQFImage &) image).getWidth();
		((DQFImage &) image).mapping(eh.scale, eh.shift);
		eh.separateXDim = ((DQFImage &) image).getXDim();
		eh.separateYDim = ((DQFImage &) image).getYDim();
		eh.separateZDim = ((DQFImage &) image).getZDim();
		peh = &eh;    // New format 1
	}
	else if (globalControl->readInt(ImageFormat) == 0)
		peh = NULL;   // Old format 0
	else
		peh = &eh;    // New format 1

	// Note: there was originally no way to set the intensity range recorded
	// upon output; the min and max arguments were always MIN_GREY_VALUE
	// and MAX_GREY_VALUE.  See pablo_format.txt for details.  Either min
	// or max may equal 99999.
    return write(filename, xDim, yDim, zDim, xSpacing, ySpacing, zSpacing,
		originPos.getX(), originPos.getY(), originPos.getZ(),
		min, max, tempArray, (int) ord, compress, peh);
}

/*  Low-level write function.  This does not swap the bytes of the image.
    The voxels should be pre-swapped to match the byte order specified for
    the output.  If this is NATIVE, then the byte order should correspond
    to the order of the machine.  The old_filename should be provided if
    comments are to be copied from it.  If eh is NULL, format 0 is used,
    otherwise format 1 is written.
*/ 
bool RAWImageFile::write(const char * filename, int xDim, int yDim, int zDim,
						 double xSpacing, double ySpacing, double zSpacing,
						 double xOrig, double yOrig, double zOrig, 
						 int minVal, int maxVal, const short * voxels,
						 int byte_order, bool useGzip, ExtendedHeader * eh,
						 const char * old_filename)
{
    FILE * fp;
	order_t file_byte_order;
	bool swapBytes;

    order_t bo = Order();
	if (bo != NORMAL && bo != REVERSE) {
		cout << "Machine has an unsupported byte order: image file cannot be written" << endl;
		return false;
	}

    fp = fopen(filename, "wb");
    if (fp == NULL)
        return false;

	if (old_filename)
		copyComments(old_filename, fp);

	// Decide how to write the CRC; this is not used to swap the bytes of the image
	file_byte_order = (order_t) byte_order;
	swapBytes = false;
	if (file_byte_order == REVERSE) {
		if (bo != REVERSE)
			swapBytes = true;
		fprintf(fp, "lsb");
	}
	else if (file_byte_order == NORMAL) {
		if (bo != NORMAL)
			swapBytes = true;
		fprintf(fp, "msb");
	}
	else {
		cout << "Unsupported byte ordering specified: image file cannot be written" << endl;
		(void) fclose(fp);
		(void) unlink(filename);
		return false;
	}

	if (eh != NULL) {
		fprintf(fp, "\n");
		writeExtendedHeader(fp, eh);
	}
	else
		fprintf(fp, " ");

    fprintf(fp, "%d %d %d ", xDim, yDim, zDim);
    fprintf(fp, "%g %g %g ", xSpacing, ySpacing, zSpacing);
    fprintf(fp, "%g %g %g ", xOrig, yOrig, zOrig);
    fprintf(fp, "%d %d", minVal, maxVal);

	if (useGzip) {
		unsigned long crc;

		// Compress the voxels using gzip
#ifdef BINARY
		if (globalVerbosity > 0)
#endif
			cout << "Compressing image .." << flush;
		int size = xDim * yDim * zDim * sizeof(GreyValue);
		int len = ((int) (COMPRESS_FACTOR*size)) + 16;
		char * compressedVoxels = new char [len];

		crc = 0;
		crc = crc32(crc, (const Bytef *)voxels, size);
#ifdef DEBUG
		cout << "CRC  on " << size << " bytes= " << crc << endl;
#endif

		int ret = compress((Bytef *) compressedVoxels, (uLongf *) &len, (const Bytef *) voxels,
			(uLongf) size);
#ifdef BINARY
		if (globalVerbosity > 0)
#endif
			cout << ".. done" << endl;
		if (ret != Z_OK) {
			if (ret == Z_MEM_ERROR)
				cout << "Compress did not have enough memory" << endl;
			else
				cout << "Compress did not have room in its output buffer" << endl;
			useGzip = false;
			cout << "Writing image data uncompressed" << endl;
		}
		else {
			if (eh == NULL)	// Old format
				fprintf(fp, "\t%d ", len);	// TAB indicates that the data is compressed
			else	// With an extended header, regular header ends with newline
				fprintf(fp, " Z\n\n%d ", len);	// Z indicates that the data is compressed
#ifdef DEBUG
			cout << "Compression ratio achieved: " << (float)len/(float)size << ".\n";
#endif
			ret = fwrite(compressedVoxels, sizeof(char), len, fp);
			if (swapBytes)
				swapLong(crc);
			ret += fwrite(&crc, sizeof(unsigned long), 1, fp);
			if (ret != len + 1) {
				cout << "Error: writing of output image failed" << endl;
				(void) fclose(fp);
				(void) unlink(filename);
				return false;
			}
		}
		delete [] compressedVoxels;
	}

	if (! useGzip) {
		int size = xDim * yDim * zDim;
		fprintf(fp, " ");
		if (eh != NULL)	// New format
			fprintf(fp, "u\n\n");		// Indicator that the data is not compressed
		int ret = fwrite(voxels, sizeof(short), size, fp);
		if (ret != size) {
			cout << "Error: writing of " << filename << " failed" << endl;
			(void) fclose(fp);
			(void) unlink(filename);
			return false;
		}
	}

    (void) fclose(fp);
	return true;
}

bool RAWImageFile::skipComments(FILE * fp)
{
    char c;
	bool ret;

    c = getc(fp);
	ret = false;
    if (c == '#')
		ret = true;

    while (c == '#')
    {
        while (c != '\n' && c != EOF)
            c = getc(fp);

        if (c != EOF)
            c = getc(fp);
    }

    ungetc(c, fp);
	return ret;
}

void RAWImageFile::copyComments(const char * old_filename, FILE * new_fp)
{
    char c;
	FILE * old_fp;

    old_fp = fopen(old_filename, "rb");
	if (old_fp  == NULL)
		return;

    c = getc(old_fp);
    while (c == '#')
    {
		(void) putc(c, new_fp);

        while (c != '\n' && c != EOF) {
            c = getc(old_fp);
			(void) putc(c, new_fp);
		}

        if (c != EOF)
            c = getc(old_fp);
    }
	(void) fclose(old_fp);
}



