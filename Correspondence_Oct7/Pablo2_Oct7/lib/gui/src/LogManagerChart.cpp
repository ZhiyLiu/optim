#include "globalBuildParms.h"
#ifdef OPTIMIZATION_VISUALIZER

#include <FL/Fl.H>
#include <FL/Enumerations.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Widget.H>
#include <iostream>
#include "LogManagerChart.h"
#include "P3DUserInterfaceCallback.h"
#include "Plot3DWindowWrapper.h"
#include "FunctionExplorer.h"


using namespace std;




inline int LogManagerChart::dataToScreen(int dataX)
{
	int width = (w() - 16) / n;
	return (8 + dataX * width);
}

inline int LogManagerChart::screenToData(int screenX)
{
	int width = (w() - 16) / n;
	double d ((screenX - 8) * 1.0);
	d = ((d / width)  + 0.5);
	return d;
}

LogManagerChart::LogManagerChart(int x, int y, int w, int h) :
	Fl_Widget(x, y, w, h, NULL)
{
		stepSize = 17;
		category = 0;
		firstIndex = 0;
		lastIndex = -1;
		n = 0;
		selectionX = 0;
		selectionType = SelectCenter;	// Matching Origin button in OptVisualizerUI.fl
		threed = NULL;
}

void LogManagerChart::draw()
{
	fl_color(FL_BACKGROUND_COLOR);
	fl_rectf(x(), y(), w(), h());
	draw_box();
	draw_label();
	if ((globalLogManager.numEventsPerCategory[category] > 0) &&
		(n > 0))
	{		
		// Next calculate the range to draw one
		int attr;
		double yScale;
		int i;

		// Assume calculateBounds is called before we are invalidated
		yScale = (1.0 * (h() - 16)) / (maxY - minY);

		// Now draw a line plot
		int width = (w() - 16) / n;		
		for (attr = 0; attr < globalLogManager.numAttributes; attr++) {
			if (attributeIsVisible[attr]) {
				fl_color(attr);	

				int prevX = 8;
				int prevY = h() - 8 - (globalLogManager.eventsByCategory[category][firstIndex][attr]
					- minY) * yScale;
				int currX, currY;
				for (i = 1; i < n; i++) { 
					currX = 8 + i * width;
					currY = h() - 8 - (globalLogManager.eventsByCategory[category][dataIndices[i]][attr]
						- minY) * yScale;
					fl_line(prevX, prevY, currX, currY);
					fl_line(prevX, prevY - 1, currX, currY - 1);
					fl_line(prevX - 3, prevY - 3, prevX + 3, prevY + 3);
					fl_line(prevX - 3, prevY + 3, prevX + 3, prevY - 3);
					prevX = currX;
					prevY = currY;
				}
				fl_line(prevX - 3, prevY - 3, prevX + 3, prevY + 3);
				fl_line(prevX - 3, prevY + 3, prevX + 3, prevY - 3);
			}
		}	
	}
	else
		n = 0;

	fl_color(0);
	fl_line(selectionX, 2, selectionX, h() - 2);
}

int LogManagerChart::handle(int event)
{
	if (event == FL_PUSH)
	{
		if (Fl::event_button() == FL_LEFT_MOUSE) {
			selectionX  = Fl::event_x();
			//		cout << "Left button click: selectionX = " << selectionX<< endl;	
			if (n > 0) {
				int idx = ((selectionX  * n) - 8) / (w()-16); 
				idx = screenToData(selectionX);
				//			cout << "click is related to data item # " << idx << endl;
				if(idx < n ) {
					cout << "data index is " << dataIndices[idx] << endl;
					for (int j = 0; j < globalLogManager.numAttributes; j++) {
						cout << ">>> " << globalLogManager.attributeNames[j] << ": " <<
							globalLogManager.eventsByCategory[category][dataIndices[idx]][j] << endl;
					}
					int width = (w() - 16) / n;		
					selectionX = 8 + idx * width;
					Vector * v = globalLogManager.allParameters[(category * MAX_EVENTS)
						+ dataIndices[idx]];

					switch (selectionType) {
						case SelectCenter: threed->setCenterVector(v); break;
						case SelectX1: threed->setX1Vector(getDerivative(v)); break;
						case SelectX2: threed->setX2Vector(getDerivative(v)); break;
						case SelectDiffX1: threed->setDiffX1Vector(v); break;
						case SelectDiffX2: threed->setDiffX2Vector(v); break;
						case SelectProjection: threed->setPointVector(v); break;
						case SelectUnknown: return 0;
					}
					damage(-1);
					return 1; // Now wait for the ui to request a drawing

					if (v) {
						Vector * other = (idx + 1 < n ?
							globalLogManager.allParameters[(category * MAX_EVENTS) + dataIndices[idx +1]] :
							globalLogManager.allParameters[(category * MAX_EVENTS) + dataIndices[idx - 1]]);

	  					Vector dirA((*other) - (*v));
						if (Fl::event_alt()) {
							cout << "*** alt key is pressed, other point will be origin" << endl;
							dirA = - (*v);
						}
						dirA *= 0.125;

						FunctionExplorer::explorePlane(globalLogManager.getProblem(), *v, dirA);
						if (threed)
							threed->damage(-1);
					}
				}
			}
			damage(1); 
			return 1;
		}
		else if (Fl::event_button() == FL_RIGHT_MOUSE)
		{
			int attr;
			int event;
			cout << endl << endl << "#";
			for (attr = 0; attr < globalLogManager.numAttributes; attr++) {
				if (attributeIsVisible[attr]) {
					cout << "[" << globalLogManager.attributeNames[attr] << "]\t";
				}
			}
			cout << endl;
			for (event = 0; event < n; event++) {
				for (attr = 0; attr < globalLogManager.numAttributes; attr++) {
					if (attributeIsVisible[attr]) {
						cout << globalLogManager.eventsByCategory[category][dataIndices[event]][attr]
							<< "\t";
					}
				}
				cout << endl;
			}
			cout << endl;
			return 1;
		}
		else
			return 0;
	}
	else
		return 0;
}

void LogManagerChart::calculateBounds()
{
	int i, attr;
	minY = 1; // globalLogManager.eventsByCategory[category][firstIndex][0]; // choose attribute later
	maxY = -1; // minY;

	// First calculate the range of events to draw
	// Skip events if needed to make the ones we draw be at least one pixel wide
	int end = ((lastIndex == -1) || (globalLogManager.numEventsPerCategory[category]
		- 0 < lastIndex)) ? globalLogManager.numEventsPerCategory[category] - 1 : lastIndex;

	stp = stepSize;
	while ((n = (end - firstIndex)/stp) > (w() - 16))
		stp *= 2;

	for (i = 0; i < n; i++) {
		for (attr = 0; attr < globalLogManager.numAttributes; attr++) {
			if (attributeIsVisible[attr]) {		

				double temp;

				temp = globalLogManager.eventsByCategory[category][firstIndex + (i * stp)][attr];
				if (temp < minY) {
					minY = temp;
				}
				else if (temp > maxY) {
					maxY = temp;
				}
			}
		}
		dataIndices[i] = firstIndex + i*stp;
	}
}

const Vector LogManagerChart::getDerivative(Vector * point)
{
	DifferentiableFunction * f;
	if ((point && globalLogManager.getOptimizer()) &&
		(f = globalLogManager.getOptimizer()->getDifferentiableProblem()))
	{
		Vector v(point->size());
		f->computeOneJet(*point, v);
		return v;
	}
	else
		return Vector(0);
}



#endif	/* OPTIMIZATION_VISUALIZER */


