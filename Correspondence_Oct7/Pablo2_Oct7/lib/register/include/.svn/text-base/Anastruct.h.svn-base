#ifndef ANASTRUCT_LIST_H
#define ANASTRUCT_LIST_H


#ifdef AE2_BUILD

struct Contour2D {
    int pointCount;
    int slice_number;
    float min[3];
    float max[3];
    // These points are in world coordinates
    float z;
    float * x;
    float * y;
    // These points are in model coordinates
    double modelZ;
    double * modelX;
    double * modelY;
};


class Anastruct3D
{
    public:
        int contour_count;
        Contour2D * contours;
        char * name;
        int color;
        float glColor[3];
        bool visible;

        Anastruct3D();
        ~Anastruct3D();

	const char * getName() { return name; }
	int contourCount() { return contour_count; }
	int pointCount(int contourId) { return contours[contourId].pointCount; }
	int getColor() { return color; }
	int isVisible() { return visible; }
};


class AnastructList
{
	std::vector<Anastruct3D *> ana;

    public:

        AnastructList() { }
        ~AnastructList();
	void addAnastruct(Anastruct3D * a);
        void deleteAnastruct(int index);
        bool convertAnastructCoords(Image3D * im);

	int anastructCount() const { return ana.size(); }
	Anastruct3D * getAnastruct(int index) const;
	const char * getName(int index) const;
	int contourCount(int index) const;
	int pointCount(int index, int contourId) const;
	int getColor(int index) const;
};

#endif



#endif  /* ANASTRUCT_LIST_H */

