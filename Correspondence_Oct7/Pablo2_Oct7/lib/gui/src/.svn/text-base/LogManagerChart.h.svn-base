#ifndef _LOG_MANAGER_CHART_H_
#define _LOG_MANAGER_CHART_H_

#include <FL/Fl_Widget.H>
#include "LogManager.h"



class P3DControl;
class Plot3DWindowWrapper;


class LogManagerChart : public Fl_Widget
{

public:

	enum SelectionType {SelectCenter, SelectX1, SelectX2, SelectDiffX1,
		SelectDiffX2, SelectProjection, SelectUnknown };

	LogManagerChart(int x, int y, int w, int h);
	void setControl(P3DControl * c) { control = c; }

	void setThreedView(Plot3DWindowWrapper * p3w) { threed = p3w; }
	void setCategory(int cat) { category = cat; }
	void setStepSize(int step) { stepSize = step; }
	void setDataRange(int first, int last) { 
			firstIndex = first;
			lastIndex = last;
	}

	virtual void draw();
	virtual int handle(int event);

	void setAttributeIsVisible(int idx, int val) { 
		if ((idx >= 0) && (idx < MAX_ATTRIBUTES)) {
			attributeIsVisible[idx] = val;
		}
	}

	void calculateBounds();
	void setSelectionType(SelectionType st) { selectionType = st; }

	double minY, maxY;

private:

	int category; 
	// Map from x in chart coordinates to log manager data index
	int dataIndices[MAX_EVENTS];
	int stepSize;	// Draw every nth data point
	int stp;		// AGG: Why is this variable needed?
	int firstIndex, lastIndex;	// Range of data to draw
	int selectionX;		// To select at item
	int n;	// The number of items drawn

	int dataToScreen(int dataX);
	int screenToData(int screenX);
	int attributeIsVisible[MAX_ATTRIBUTES];

	SelectionType selectionType;
	P3DControl * control;	// To access the model slideshow
	Plot3DWindowWrapper * threed;

	const Vector getDerivative(Vector * point);
};




#endif	/* _LOG_MANAGER_CHART_H_ */

