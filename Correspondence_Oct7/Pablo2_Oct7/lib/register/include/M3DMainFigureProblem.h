#ifndef M3D_MAIN_FIGURE_PROBLEM_H
#define M3D_MAIN_FIGURE_PROBLEM_H

#include "Match.h"
#include "SimilarityTransform3D.h"
#include "optima.h"
#include "LogManager.h"

#ifdef BINARY
extern SimilarityTransform3D GLOBALMOMTransform;
#endif

class M3DMainFigureProblem : public Function
{
	public:
		M3DMainFigureProblem();
		~M3DMainFigureProblem() {}

		virtual void initialize(Match * m, M3DObject * referenceObject,
			int treeIndex);
		bool setAddToResidueObject( const M3DObject* object ) {
			addToResidueObject	= object;
			return true;
		}
		virtual void setParameterCount(int count) { nparms = count; }

		Match * getMatch() { return match; }

		virtual double evaluate(const Vector &x) = 0;

		double getLastPenalty() { return lastPenalty; }
		virtual int parameterCount() const { return nparms; }

		// CL: predict turns on/off prediction of other figures;
		// it is only used in Adaptive Pablo
		virtual M3DObject * createTargetObject(const Vector & x, int & figureId,
			bool predict = false) = 0;

		static const double PENALTY_SCALE_FACTOR;
		static const double CONSTRAINT_PENALTY_SCALE_FACTOR;
		static const double MAIN_FIGURE_MAX_RETURN;

	protected:

#ifdef OPTIMIZATION_VISUALIZER
		virtual void copyResultsToVisualizer(const Vector * x, double totalPenalty);
#endif

		Match * match;
		M3DObject * refObject;
		const M3DObject * addToResidueObject;

		SimilarityTransform3D mainFigProblemSimTransform;

		double lastPenalty;
		int tree_index;
		int nparms;
};

#endif

