

#include <iostream>
#include <vector>
#include "InterfiguralConstraints.h"



InterfiguralConstraints::InterfiguralConstraints(const InterfiguralConstraints & ifc) {
	figures = ifc.figures;
	dists = ifc.dists;
}

InterfiguralConstraints & InterfiguralConstraints::operator= (const
					InterfiguralConstraints & ifc) {
	figures = ifc.figures;
	dists = ifc.dists;
	return *this;
}

int InterfiguralConstraints::addFigure(int figureID, float dist) {
	figures.push_back(figureID);
	dists.push_back(dist);
	return figures.size() - 1;
}

bool InterfiguralConstraints::deleteFigure(int figureID) {
	for (int i = 0; i < size(); i++) {
		if (figures[i] == figureID) {
			figures.erase(figures.begin() + i);
			dists.erase(dists.begin() + i);
			return true;
		}
	}
	return false;
}

bool InterfiguralConstraints::updateFigure(int figureID, float dist) {
	for (int i = 0; i < size(); i++) {
		if (figures[i] == figureID) {
			dists[i] = dist;
			return true;
		}
	}
	return false;
}

bool InterfiguralConstraints::changeFigureId(int oldID, int newID) {
	for (int i = 0; i < size(); i++) {
		if (figures[i] == oldID) {
			figures[i] = newID;
			return true;
		}
	}
	return false;
}

void InterfiguralConstraints::remapFigureIds(int * map, int len) {
	int i, oldID;

	for (i = 0; i < size(); i++) {
		for (oldID = 0; oldID < len; oldID++) {
			if (figures[i] == oldID) {
				figures[i] = map[oldID];
				break;
			}
		}
	}
}

void InterfiguralConstraints::clear() {
	figures.clear();
	dists.clear();
}

void InterfiguralConstraints::print(const char * msg) {

	using namespace std;

	if (msg == NULL)
		cout << "InterfiguralConstraints:\n";
	else
		cout << "InterfiguralConstraints " << msg << ':';
	cout << "	size = " << size() << '\n';
	for (int i = 0; i < size(); i++)
		// The "this" below is needed on Solaris
		cout << '	' << figure(i) << ".  " << this->distance(i) << '\n';
	cout << "	figures = " << &figures << "    dists = " << &dists << endl;
}



