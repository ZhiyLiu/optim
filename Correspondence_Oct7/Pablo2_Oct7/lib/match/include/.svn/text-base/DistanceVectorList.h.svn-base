#ifndef DISTANCE_VECTOR_LIST_H
#define DISTANCE_VECTOR_LIST_H


/*
 * This file defines the class DistanceVectorList, which is used to
 * show the distance vectors.
 *
 */

struct DistanceVector {
	Vector3D p;
	Vector3D n;
	Vector3D grad;
	double dist;
	int method;		// 0 for cache, 1 for grad, 2 for normal
};

inline std::ostream & operator << (std::ostream & out, const DistanceVector & p )
{
	out << p.p.getX() << ' ' << p.p.getY() << ' ' << p.p.getZ() << ' '
		<< p.n.getX() << ' ' << p.n.getY() << ' ' << p.n.getZ() << ' '
		<< p.grad.getX() << ' ' << p.grad.getY() << ' ' << p.grad.getZ() << ' '
		<< p.dist << ' ' << p.method << std::endl;
	return out;
}

inline std::istream & operator >> (std::istream & in, DistanceVector & p )
{
	double x,y,z;
	in >> x >> y >> z;
	p.p.set(x,y,z);
	in >> x >> y >> z;
	p.n.set(x,y,z);
	in >> x >> y >> z;
	p.grad.set(x,y,z);
	in >> p.dist;
	in >> p.method;
	return in;
}


class DistanceVectorList
{
	public:
		typedef std::vector<DistanceVector> DistVectorPList;
		DistVectorPList plist;
};


#endif	/* DISTANCE_VECTOR_LIST_H */

