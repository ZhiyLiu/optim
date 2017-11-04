
#include "utility.h"

#include <sys/stat.h>
#include <fcntl.h>
#if ! defined(_MSC_VER) && ! defined(__BORLANDC__)
#include <unistd.h>
#endif

#ifndef S_ISREG
#define S_ISREG(m)  (((S_IFREG & m) != 0) ? true : false)
#endif


using namespace std;


double inter_pt_distance(const Point & p1, const Point & p2) {
	return (p1 - p2).norm();
}


/*
	Function to project a point to a line to get a new point.  The
	point must be provided in p.  The line is specified by the other
	two points.  The parameter t corresponds to the intersection.
*/
Point project_point_to_line(const Point & p, const Point & l0, const Point & l1, double & t) {
	Vector3D u = l1 - l0;
	Vector3D v = p - l0;
	double dp = u*v;
	t = dp/(u*u);
	Point p_proj = t*l1 + (1.0 - t)*l0;
	return p_proj;
}

/*
	Function to compute the distance of a point from a line defined by
	two points on the line.  The point must be provided in p and the
	line by the other two points.
*/
double pt_to_line_distance(const Point & p, const Point & l0, const Point & l1) {
	Point p_proj = project_point_to_line(p, l0, l1);
	Vector3D v = p - p_proj;
	double d = v.norm();
	return d;
}

/*
	Function to project a point to a plane to get a new point.  The
	point must be provided in p.  The plane, ax + by + cz = d, is
	specified by its normal, n = ai + bj + ck, and the offset d,
	the distance of the plane from the origin.
*/
Point project_point_to_plane(const Point & p, const Vector3D & n, double d) {
	double len = n.norm();
	len *= len;
	double t = (d - n*p)/len;
	Point p_proj = t*n + p;
	return p_proj;
}

/*
	Function to project a point to a plane to get a new point.  The
	point must be provided in p.  The plane is specified by the other
	three verticies.
*/
Point project_point_to_plane(const Point & p, const Point & v0,
							 const Point & v1, const Point & v2)
{
	Vector3D u = v1 - v0;
	Vector3D v = v2 - v0;
	Vector3D n = v.cross(u);
	n.normalize();	// AGG: needed?
	double d = n*v0;
	return project_point_to_plane(p, n, d);
}

/*
	Function to find the point of intersection of a line and a plane.  The
	line must be specified by the two points l0 and l1.  The plane is specified
	by the other three verticies.
*/
Point intersect_line_with_plane(const Point & l0, const Point & l1,
								const Point & v0, const Point & v1, const Point & v2)
{
	Vector3D u = v1 - v0;
	Vector3D v = v2 - v0;
	Vector3D n = v.cross(u);
	n.normalize();	// AGG: needed?
	double d = n*v0;

	double nl0 = n*l0;
	double nl1 = n*l1;
	double t = (d - nl0)/(nl1 - nl0);
	Point q = t*l1 + (1.0 - t)*l0;
	return q;
}

/*	Function to see if point p is on or "above" a designated line, specified
	by points l0 and l1.  Thus, a point inside a triangle, will be "above"
	or "below" all of its sides, as long as the points are similarly ordered.
*/
bool above_line(const Point & p, const Point & l0, const Point & l1) {
	Vector3D u, v, w;

	u = l1 - l0;
	v = p - l0;
	w = u.cross(v);
	if (w.getZ() < 0.0) return false;
	else return true;
}

/*	Function to see if point p is inside the triangle specified by v0,
	v1 and v2.
*/
bool inside_triangle(const Point & p, const Point & v0, const Point & v1,
					 const Point & v2)
{
	int i = (int) above_line(p, v0, v1);
	i += (int) above_line(p, v1, v2);
	i += (int) above_line(p, v2, v0);
	if (i == 0 || i == 3)
		return true;
	else
		return false;
}


/*
	Diagnostic output
*/

void print_image_info(Image3D & im, const char * name) {
	if (name != NULL)
		cout << name << ":\n";
	else
		cout << "Image:\n";
	cout << "\tDimensions: " << im.getXDim() << " " 
		<< im.getYDim() << " " << im.getZDim() << " pixels\n";
	cout << "\tPixel size (cm): " << im.getXSpacing() << " "
		<< im.getYSpacing() << " " << im.getZSpacing() << '\n';
	cout << "\tImage size (cm): " << im.getXExtent() << " "
		<< im.getYExtent() << " " << im.getZExtent() << '\n';
	Point o = im.getWorldOrigin();
	cout << "\tVoxel origin from file (cm): " << o.getX() << " "
		<< o.getY() << " " << o.getZ() << "\n\n";
	o = im.getModelOrigin();
	cout << "\tVoxel coordinates of model origin (cm): " << o.getX() << " "
		<< o.getY() << " " << o.getZ() << "\n\n";
}


/*
	Diagnostic output
	The average inter-atom distances at the atom centers are reported
*/

void list_interatom_dist(vector<M3DPrimitive *> & set) {
	Point p, pn;

	cout << "Inter-atom distances:\n";
	double avg_interatom_dist = 0;
	p = set[0]->getX();
	int natom;
	for (natom = 1; natom < set.size(); natom++) {
		pn = set[natom]->getX();
		double dist = inter_pt_distance(pn, p);
		cout << '\t' << dist << '\n';
		avg_interatom_dist += dist;
		p = pn;
	}
	avg_interatom_dist /= (natom - 1);

	cout << "\taverage = " << avg_interatom_dist << '\n';
}


// Compute the intensity range of a floating point array
void calc_intensity_range(unsigned char * pixels, int size, unsigned char & min,
    unsigned char & max)
{
	register int k;

	min = pixels[0];
	max = pixels[0];

	int start = size & 0x01;

	for(k = start; k < size; k += 2)
	{
		register unsigned char val0 = pixels[k];
		register unsigned char val1 = pixels[k + 1];

		if(val0 < val1) {
			if(val0 < min)
				min = val0;
			if(val1 > max)
				max = val1;
		}
		else {
			if(val1 < min)
				min = val1;
			if(val0 > max)
				max = val0;
		}
	}
}


// Compute the intensity range of a floating point array
void calc_intensity_range(float * pixels, int size, float & min, float & max)
{
	register int k;

	min = pixels[0];
	max = pixels[0];

	int start = size & 0x01;

	for(k = start; k < size; k += 2)
	{
		register float val0 = pixels[k];
		register float val1 = pixels[k + 1];

		if(val0 < val1) {
			if(val0 < min)
				min = val0;
			if(val1 > max)
				max = val1;
		}
		else {
			if(val1 < min)
				min = val1;
			if(val0 > max)
				max = val0;
		}
	}
}


// Compute the intensity range of a GreyValue array
void calc_intensity_range(GreyValue * pixels, int size, GreyValue & min, GreyValue & max)
{
	register int k;

	min = pixels[0];
	max = pixels[0];
#ifdef AE2_BUILD
	min &= IMAGE_BITS;
	max &= IMAGE_BITS;
#endif

	int start = size & 0x01;

	for(k = start; k < size; k += 2)
	{
		register GreyValue val0 = pixels[k];
		register GreyValue val1 = pixels[k + 1];
#ifdef AE2_BUILD
	    val0 &= IMAGE_BITS;
	    val1 &= IMAGE_BITS;
#endif

		if(val0 < val1) {
			if(val0 < min)
				min = val0;
			if(val1 > max)
				max = val1;
		}
		else {
			if(val1 < min)
				min = val1;
			if(val0 > max)
				max = val0;
		}
	}
}


// Compute the intensity range of an Image3D object
void calc_intensity_range(Image3D & image, GreyValue & min, GreyValue & max) {
	int xdim, ydim, zdim;

	xdim = image.getXDim();
	ydim = image.getYDim();
	zdim = image.getZDim();

	GreyValue * pixels = image.getVoxels();
	int size = xdim*ydim*zdim;

	min = pixels[0];
	max = pixels[0];

	calc_intensity_range(pixels, size, min, max);
}


/*  Find the min, the max, or both of a voxel array.  Also, if necessary,
    swap the byte order of the data before the calculation.  If both the
	min and max are needed, a faster algorithm than simply checking all
	values is used.  If neither findMin nor findMax are true, no action
	is taken.  The size argument is the number of voxels.  Note that the
	min and max are integers to allow for negative values.
*/
void calc_rangeWithSwap(bool findMin, int & min, bool findMax, int & max,
		   bool swapBytes, short * data, int size)
{
	GreyValue val0, val1;
	int i;
	int start;

    if (findMin && findMax)
    {
		start = size & 0x01;
		if (start) {
            val0 = (GreyValue) data[0];
            if (swapBytes) {
                val0 = ((val0 >> 8) & 0x00FF) | (val0 << 8);
	            data[0] = val0;
			}

            min = val0;
            max = val0;
        }
		else {
			min = MAX_GREY_VALUE;
			max = MIN_GREY_VALUE;
		}

        for (i = start; i < size; i += 2)
        {
            val0 = (GreyValue) data[i];
            if (swapBytes) {
                val0 = ((val0 >> 8) & 0x00FF) | (val0 << 8);
	            data[i] = val0;
			}

            val1 = (GreyValue) data[i + 1];
            if (swapBytes) {
                val1 = ((val1 >> 8) & 0x00FF) | (val1 << 8);
	            data[i + 1] = val1;
			}

            if (val0 < val1) {
				if (val0 < min)
					min = val0;
				if (val1 > max)
					max = val1;
			}
            else {
				if (val1 < min)
					min = val1;
				if (val0 > max)
					max = val0;
			}
        }
    }
    else {
		if (findMin)
		{
			min = MAX_GREY_VALUE;

			for (i = 0; i < size; i++)
			{
				val0 = (GreyValue) data[i];
				if (swapBytes) {
					val0 = ((val0 >> 8) & 0x00FF) | (val0 << 8);
					data[i] = val0;
				}

				if (val0 < min)
					min = val0;
			}
		}
		if (findMax)
		{
			max = MIN_GREY_VALUE;

			for (i = 0; i < size; i++)
			{
				val0 = (GreyValue) data[i];
				if (swapBytes) {
					val0 = ((val0 >> 8) & 0x00FF) | (val0 << 8);
					data[i] = val0;
				}

				if (val0 > max)
					max = val0;
			}
		}
	}
}


void report_intensity_range(Image3D & image, char * label)
{
	GreyValue min;
	GreyValue max;

	calc_intensity_range(image, min, max);
	cout << "Intensity range";
	if (label)
		cout << " of " << label;
	cout << " = (" << min << ", " << max <<")" << endl;
}

void report_intensity_range(GreyValue * pixels, int size, char * label)
{
	GreyValue min;
	GreyValue max;

	calc_intensity_range(pixels, size, min, max);
	cout << "Intensity range";
	if (label)
		cout << " of " << label;
	cout << " = (" << min << ", " << max <<")" << endl;
}

void verify_squareness(const Point & origin, const Point & corner,
					   const Point & end, char * label)
{
	Vector3D u, v;
	const double fuzz = 1.0e-15;

	u = corner - origin;
	v = end - corner;
	double d = u*v;
	if (d < 0.0) d = -d;
	if (d > fuzz)
		cout << label << " IS NOT ORTHOGONAL (dot product = " << d << ")\n";
}


#ifdef WIN32
#include <windows.h>
#include <winuser.h>
#endif

void beep()
{
#ifdef WIN32
      MessageBeep(MB_ICONERROR);
#else
	  cout << '\007';
#endif
}


bool fileExists(const char * filename)
{
    struct stat buf;

	bool fileMissing = (stat(filename, &buf) != 0 || S_ISREG(buf.st_mode) == false);
	return ! fileMissing;
}


/*  In some cases, fread() fails; this does not.  The first argument
    is a pointer to the memory buffer into which the data will be read.
	The second argument is the number of bytes (note: not items) to be
	read, starting at the current position in the file.  Finally, fp is
	the file from which data is to be read.  This must be open before
	calling this function.
*/
int safeFread(void * buffer, unsigned long size, FILE * fp)
{
	// fread() fails when all these conditions are met:
	// (1) the file is on AFS (a problem specific to Radonc)
	// (2) the file contains many zero's
	// (3) the entire file is read in one pass.
	// --> we can't control #1 & #2, so read in 8K byte chunks

	unsigned long count = 0;		// Bytes read so far
	char * charBuf = (char *) buffer;
	while (count != size && ! feof(fp) && ! ferror(fp))
	{
		size_t n = 8*1024;
		if (size - count < n)
			n = size - count;
		size_t justRead = fread(charBuf + count, 1, n, fp);

		if (justRead == 0)		// Premature undetected EOF?
			break;
		count += justRead;
	}
	return count;
}

#ifdef _WIN32

void reverseBackSlashes(char * str)
{
	while (*str != '\0') {
		if (*str == '\\')
			*str = '/';
		str++;
	}
}

#endif


