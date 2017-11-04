#ifndef M3D_REGISTRATION_PROBLEM_H
#define M3D_REGISTRATION_PROBLEM_H

#include "Match.h"
#include "SimilarityTransform3D.h"
#include "optima.h"


#define PARAMETER_SIZE 7
#define MAX_RETURN 1.7e308



//#define EUCLIDEAN_ONLY	/* restrict to euclidean transform only */

class M3DRegistrationProblem : public Function
{
	public:

		M3DRegistrationProblem();
		M3DRegistrationProblem(Match * _match, double _penaltyWeight = 0.0);
		virtual void initialize(Match * _match, double _penaltyWeight = 0.0,
			double _constraintPenaltyWeight = 0.0);

		~M3DRegistrationProblem() {}

		Match * getMatch() { return match; }

		void setPenaltyWeight(double w) { penaltyWeight = w; }
		void setConstraintsPenaltyWeight(double w) { constraintsPenaltyWeight = w; }

		// The vector parameter is evaluated as a 7-tuple:
		// (translation in x,y,z; rotation around x,y,z; scale)
		virtual double evaluate(const Vector &x);

		virtual M3DObject * createTargetObject(const Vector & x) {
			return NULL; }

		double computePenalty(const Vector &x);

		void computeTransformation(SimilarityTransform3D & matrix, const Vector &x);

		double getPenalty() { return lastPenalty; }

		static const double ENSEMBLE_OBJECTIVE_SCALE_FACTOR;
		static const double ENSEMBLE_CONSTRAINT_PENALTY_SCALE_FACTOR;
		static const double ENSEMBLE_MAX_RETURN;
		static const double ENSEMBLE_PGA_SCALE_FACTOR;

	protected:

		Match * match;

		double penaltyWeight;
		double constraintsPenaltyWeight;
		double lastPenalty;

		static const double ENSEMBLE_TRANSLATION_FACTOR;
		static const double ENSEMBLE_ROTATION_FACTOR;
		static const double ENSEMBLE_SCALE_FACTOR;
};


#endif	/* M3D_REGISTRATION_PROBLEM_H */

