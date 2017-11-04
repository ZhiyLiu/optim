#include "M3DMainFigureProblem.h"
#include "M3DFigureElongater.h"
#include "Tuning.h"

//#define DISPLAY_PENALTIES	/* Uncomment to print details of penalty calculations */


const double M3DMainFigureProblem::MAIN_FIGURE_MAX_RETURN = 1.7e308;
const double M3DMainFigureProblem::PENALTY_SCALE_FACTOR = 1.0;
const double M3DMainFigureProblem::CONSTRAINT_PENALTY_SCALE_FACTOR = 0.02;


M3DMainFigureProblem::M3DMainFigureProblem()
{
    match = NULL;
    refObject = NULL;
	addToResidueObject	= NULL;
    tree_index = 0;
    lastPenalty = 0.0;
    mainFigProblemSimTransform.setToIdentity();
}

void M3DMainFigureProblem::initialize(Match * m, M3DObject * referenceObject,
	int treeIndex)
{
    match = m;
    refObject = referenceObject;
    tree_index = treeIndex;
#ifndef BINARY
	mainFigProblemSimTransform = match->getInitialTransform();
#endif

#ifdef DISPLAY_PENALTIES
    std::cout << "Main figure penalties:  \n";
#endif
}

#ifdef OPTIMIZATION_VISUALIZER

void M3DMainFigureProblem::copyResultsToVisualizer(const Vector * x, double totalPenalty) {
	if (x && match) {
		globalLogManager.beginEvent(x);
		const Match::matchResult * results = match->getFigureStageResults();
		if (results) {
			double wtMatch;
			for (int i = 0; i < MAX_NUM_FIG_MATCH_RESULTS; i++) {
				if (results[i].useIt)
					wtMatch = tuningWt(results[i].tuningParm)*results[i].value;
				else
					wtMatch = 0.0;
				globalLogManager.setValue(1 + i, wtMatch);
			}
			globalLogManager.setValue(0, totalPenalty);
		}
	}
}

#endif	/* OPTIMIZATION_VISUALIZER */

