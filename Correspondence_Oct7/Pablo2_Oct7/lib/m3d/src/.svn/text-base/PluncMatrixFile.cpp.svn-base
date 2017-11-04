// pluncMatrixFile -- read/write PLUNC matrix file.
/*
 file format is optional comment lines followed by 4 lines of
 4 floating point numbers per line, and unix end-of-line:

[comments]\n
sxx sxy sxz 0\n
syx syy syz 0\n
szx szy szz 0\n
tx  ty  tz  1\n

 (double *mat) is in order [sxx sxy sxz 0 syx ... tx ty tz 1].
 right-most column is [0 0 0 1] for homegeneous mat.
*/

#include <iostream>

#include <fstream>

using namespace std;

// read matrix file 
// Return # of elements found in file -- should be 16.
// Optional comment is printed and discarded.
// Caller alloc/free's the 16 elements of array 'mat.'
int readPluncMatrixFile(const char* matrixFilename, double *mat)
{
	bool ret;
	int nElements = 0;	// # of array elements read
	
	// try to read as [comment]\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n
	FILE *fp = fopen(matrixFilename, "r");
	if (!fp) {
		std::cout << "WARNING: readPluncMatrix: cannot open file " << matrixFilename << endl;
		return 0;
	}
	
	int loop, num;
	char dummy[1000];		// discarded comment line(s)
	ret = true;
	for (loop = 0; loop < 4; loop++)
	{
		num = fscanf(fp, "%lf%lf%lf%lf", 
			mat + loop*4 + 0,
			mat + loop*4 + 1,
			mat + loop*4 + 2,
			mat + loop*4 + 3);
		if (num == 4)
			nElements += 4;
		else 
		{
			// comment or end-of-file?
			if (num == EOF)
			{
				cout << "readPluncMatrix: file " << matrixFilename << ": EOF found reading element "
					<< nElements << endl;
				break;
			}
			fscanf(fp, "%[^\n]", dummy);	// read up to newline; unix line-term only
			cout << "readPluncMatrix: discarding comment string: " << dummy << endl;
			loop--;
			continue;			
		}
	}
	fclose(fp);
	return nElements;
}

int writePluncMatrixFile(const char* matrixFilename, const double *mat)
{
	FILE *fp = fopen(matrixFilename,"wb");
	if (!fp)
	{
		std::cout << "WARNING: writePluncMatrix: could not open matrix file " << matrixFilename << endl;
		return false;
	}

	fprintf(fp, "# matrix\n");
	for (int i = 0; i < 16; i++)
	{
		fprintf(fp, " %f", mat[i]);
		if (0 == ((i+1) % 4)) fprintf(fp, "\n");
	}
	fclose(fp);
	return true;
}
