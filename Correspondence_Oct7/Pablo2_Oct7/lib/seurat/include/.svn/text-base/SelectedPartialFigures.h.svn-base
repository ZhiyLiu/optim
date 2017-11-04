#ifndef SELECTEDPARTIALFIGURES_H
#define SELECTEDPARTIALFIGURES_H


/*  A class to store a list of figures for interpenetration avoidance.
    With the figures are stored the distances of the four verticies of
	each tile from the boundary of the marked figure.  Because the
	length of the distances lists are not stored, the using program
	must be consistent. 

    Function displayFigure() is used to store the figure currently
	being rendered.  An argument of -1 to displayFigure() indicates
	that all surfaces should be displayed.
*/

//#define VERIFY 1


class SelectedPartialFigures
{
public:
	SelectedPartialFigures();
	~SelectedPartialFigures();

	void initialize(int nfigs);

	// Functions for storing parameters
	void figure(int fig_num, int figureID);
	void figureSize(int fig_num, int num_boundary_pts, int num_tiles);
	void point_distance(int fig_num, int boundary_pt, double dist);
	void quad_distance(int fig_num, int vertex, double dist);
	void color(const float * c);
	void pointListsChangeReset() { recomputed = false; }
	void displayFigure(int figureID) { disp_fig_num = figureID; }
	void markedFigure(int figureID) { marked_fig = figureID; }

	// Information functions
	int size() { return maxFigures; }
	int number() { return numFigures; }
	int figureID(int fig_num);
	double point_distance(int fig_num, int boundary_pt);
	double quad_distance(int fig_num, int vertex);
	double * point_distances(int fig_num);
	double * quad_distances(int fig_num);
	const float * color() { return render_color; }
	void ** pointLists() { return plists; }
	bool pointListsChanged() { return recomputed; }
	int displayFigure() { return disp_fig_num; }
	int markedFigure() { return marked_fig; }

	void print();

private:
	int * figures;
	float render_color[3];
	int disp_fig_num;
	int marked_fig;
	int numFigures;
	int maxFigures;
	double ** point_dists;
	double ** quad_dists;
	void ** plists;
	bool recomputed;
#ifdef VERIFY
	int * point_sizes;
	int * tile_sizes;
#endif

};

#endif

