
#include <queue>
#include <set>

#define D_DIATOM
#define D_XFERLIST
#define D_POINTLIST_SERVER2
#include "Shapedepend.h"
#include "SelectedPartialFigures.h"

//#define DEBUG 1

#if defined(VERIFY) || defined(DEBUG)
#include <iostream>
using namespace std;
#endif


using namespace ThallCode;

SelectedPartialFigures::SelectedPartialFigures()
{
	figures = NULL;
	point_dists = NULL;
	quad_dists = NULL;
	plists = NULL;
	maxFigures = 0;
	numFigures = 0;
	disp_fig_num = -1;
	recomputed = false;
#ifdef VERIFY
	point_sizes = NULL;
	tile_sizes = NULL;
#endif
}

SelectedPartialFigures::~SelectedPartialFigures() {
	if (figures != NULL) {
		for (int i = 0; i < numFigures; i++) {
			delete [] point_dists[i];
			delete [] quad_dists[i];
		}
		delete [] figures;
	}
	if (plists != NULL) {
		for (int i = 0; i < numFigures; i++)
			delete ((Pointlist_server2 **) plists)[i];
		delete [] (Pointlist_server2 **) plists;
	}
#ifdef VERIFY
	delete [] point_sizes;
	delete [] tile_sizes;
#endif
}


void SelectedPartialFigures::initialize(int nfigs) {
	int i;

#ifdef DEBUG
	cout << "SelectedPartialFigures::initialize(" << nfigs << ") called\n";
#endif
	if (figures != NULL) {
		for (i = 0; i < numFigures; i++) {
			delete [] point_dists[i];
			delete [] quad_dists[i];
		}
		delete [] figures;
		figures = NULL;
		point_dists = NULL;
		quad_dists = NULL;
	}
	if (plists != NULL) {
		for (i = 0; i < numFigures; i++)
			delete ((Pointlist_server2 **) plists)[i];
		delete [] (Pointlist_server2 **) plists;
		plists = NULL;
	}
#ifdef VERIFY
	delete [] point_sizes;
	point_sizes = NULL;
	delete [] tile_sizes;
	tile_sizes = NULL;
#endif
	recomputed = false;

	maxFigures = nfigs;
	numFigures = 0;
	if (nfigs <= 0) return;

	figures = new int[nfigs];
#ifdef VERIFY
	point_sizes = new int[nfigs];
	tile_sizes = new int[nfigs];
#endif
	point_dists = new double *[nfigs];
	quad_dists = new double *[nfigs];
	for (i = 0; i < nfigs; i++) {
		point_dists[i] = NULL;
		quad_dists[i] = NULL;
	}
	plists = (void **) new Pointlist_server2 *[nfigs];
	for (i = 0; i < nfigs; i++)
		plists[i] = (void *) new Pointlist_server2;
	recomputed = true;
}

void SelectedPartialFigures::figure(int fig_num, int figureID) {
	if (fig_num >= maxFigures)
		return;
	if (fig_num + 1 > numFigures)
		numFigures = fig_num + 1;
	figures[fig_num] = figureID;
}

void SelectedPartialFigures::figureSize(int fig_num, int num_boundary_pts, int num_tiles) {
	if (fig_num >= maxFigures)
		return;
	if (fig_num + 1 > numFigures)
		numFigures = fig_num + 1;

	point_dists[fig_num] = new double[num_boundary_pts];
	quad_dists[fig_num] = new double[4*num_tiles];
#ifdef DEBUG
	cout << "SelectedPartialFigures: for file figure " << fig_num
		<< ", set number of boundary points to \n    "
		<< num_boundary_pts << " and number of tile verticies to "
		<< 4*num_tiles << endl;
#endif
#ifdef VERIFY
	point_sizes[fig_num] = num_boundary_pts;
	tile_sizes[fig_num] = 4*num_tiles;
#endif
}

void SelectedPartialFigures::point_distance(int fig_num, int boundary_pt, double dist) {
	if (fig_num >= maxFigures)
		return;
	if (fig_num + 1 > numFigures)
		numFigures = fig_num + 1;

#ifdef VERIFY
	if (boundary_pt >= point_sizes[fig_num])
		cout << "Invalid storage of point distance" << endl;
	else
#endif
		point_dists[fig_num][boundary_pt] = dist;
}

void SelectedPartialFigures::quad_distance(int fig_num, int vertex, double dist) {
	if (fig_num >= maxFigures)
		return;
	if (fig_num + 1 > numFigures)
		numFigures = fig_num + 1;

#ifdef VERIFY
	if (vertex >= tile_sizes[fig_num])
		cout << "Invalid storage of tile distance" << endl;
	else
#endif
		quad_dists[fig_num][vertex] = dist;
}

int SelectedPartialFigures::figureID(int fig_num) {
	if (fig_num < 0 || fig_num >= numFigures)
		return -1;
	else
		return figures[fig_num];
}

double SelectedPartialFigures::point_distance(int fig_num, int boundary_pt) {
	if (fig_num >= numFigures)
		return 0.0;
	else {
#ifdef VERIFY
		if (boundary_pt >= point_sizes[fig_num]) {
			cout << "Invalid request for point distance" << endl;
			return 0.0;
		}
		else
#endif
			return point_dists[fig_num][boundary_pt];
	}
}

double SelectedPartialFigures::quad_distance(int fig_num, int vertex) {
	if (fig_num >= numFigures)
		return 0.0;
	else {
#ifdef VERIFY
		if (vertex >= tile_sizes[fig_num]) {
			cout << "Invalid request for tile distance" << endl;
			return 0.0;
		}
		else
#endif
			return quad_dists[fig_num][vertex];
	}
}

double * SelectedPartialFigures::point_distances(int fig_num) {
	if (fig_num >= numFigures)
		return NULL;
	else
		return point_dists[fig_num];
}

double * SelectedPartialFigures::quad_distances(int fig_num) {
	if (fig_num >= numFigures)
		return NULL;
	else
		return quad_dists[fig_num];
}

void SelectedPartialFigures::color(const float * c) {
	for (int i = 0; i < 3; i++)
		render_color[i] = c[i];
}

using namespace std;

void SelectedPartialFigures::print() {
	cout << "SelectedPartialFigures:\n";
	cout << "	size = " << size() << '\n';
	for (int i = 0; i < size(); i++)
		cout << "	" << i << ".  " << figureID(i) << '\n';
	cout << "	display figure = " << displayFigure() << endl;
}

