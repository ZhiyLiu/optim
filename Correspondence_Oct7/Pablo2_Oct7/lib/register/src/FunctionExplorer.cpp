
#include "globalBuildParms.h"
#ifdef OPTIMIZATION_VISUALIZER


#include <iostream>
#include "LogManager.h"
#include "FunctionExplorer.h"


using std::cout;
using std::endl;



const Vector FunctionExplorer::createOrthogonalVector(const Vector & v)
{
	int n = v.size();
	Vector o(n);
	for (int i = 0; i + 1 < 2*(n/2); i += 2) {
		o(i) = v(i+1);
		o(i+1) = -v(i);
	}
	if (n % 2 == 1) {
		o(n - 1) = 0;
	}
	return o;
}

void FunctionExplorer::explorePlane(Function * f,  const Vector & v, const Vector & dir1, const Vector & dir2)
{	
	cout << "Before exploring plane" << endl;
	cout << "Center point: "; 
	v.print(); cout << endl;
	//v.print();
	cout << "Dir 1: ";
	dir1.print(); 
	cout << endl;
	cout << "Dir 2: ";
	dir2.print(); 
	cout << endl;
//	o.print();
	Vector* point;
	globalLogManager.clearCategory(2);
	globalLogManager.setCategory(2);
	if (!f) { return; }
	for (int x = -10; x <= 10; x++) {
		for (int y = -10; y <= 10; y++) {
			point = new Vector(v + (0.1 * x * dir1) + (0.1 * y * dir2));
//			globalLogManager.beginEvent(point);
			f->evaluate(*point);
		}
	}
	cout << "done exploring plane" << endl;
}


#endif	/* OPTIMIZATION_VISUALIZER */

