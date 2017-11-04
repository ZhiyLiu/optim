#ifndef _LOG_MANAGER_H_
#define _LOG_MANAGER_H_

#include "globalBuildParms.h"	/* Get definition of OPTIMIZATION_VISUALIZER */


#ifdef OPTIMIZATION_VISUALIZER


#include <utility>
#include "optima.h"
//#pragma warning(disable : 4786)


class OptimizerBase;

using std::pair;


#define MAX_ATTRIBUTES  (11)
#define MAX_EVENTS (10000)
 //for now categories could be "gradient calculation" and "search"
// and 8 parameters of of similarity transform
#define MAX_CATEGORIES (5)


class LogManager
{

public:

	typedef double record[MAX_ATTRIBUTES];
	typedef record recordList[MAX_EVENTS];

	int numAttributes;
	int numCategories;
	int numTotalEvents;

	LogManager();
	virtual ~LogManager() {}
	void setOptimizer(OptimizerBase * opt) { optimizer = opt; }
	OptimizerBase * getOptimizer() const { return optimizer; }
	Function * getProblem() const { return problem; }

	// Call these addX... methods once at setup time
	// to get indices to be used in inner loops
	int addCategory(char * name);
	int addAttribute(char * name);

	// Conjugate gradient method will call begin event before calls to 
	// f->evaluate()
	void setCategory(int newCat);
	void beginEvent(const Vector * parameter = NULL);
	int getCurrentCategory() { return currentCategory; }

	// f->evaluate() will call setValue
	void setValue(int idx, double value);

	void clearCategory(int cat);

	char * categoryNames[MAX_CATEGORIES];
	const char * attributeNames[MAX_ATTRIBUTES];
	recordList eventsByCategory[MAX_CATEGORIES];

	int numEventsPerCategory[MAX_CATEGORIES];

	// 1D search is the last event in a category
	int lastEventIndices[MAX_EVENTS][MAX_CATEGORIES]; 
	int numLastEvents[MAX_CATEGORIES];

	Vector * allParameters[MAX_CATEGORIES * MAX_EVENTS];

	// this sets up a bunch of categories and attributes that
	// I think are interesting.
	//
	// in the future, I think this should go away and individual
	// optimizers should call
	// addCategory
	// addAttribute
	void initialize();

private:

	OptimizerBase * optimizer;
	Function * problem;

	// It might be useful to see all events in chronological order
	std::pair<int, int> allEventIndex[MAX_CATEGORIES * MAX_EVENTS];

	int currentCategory;
};

extern LogManager globalLogManager;



#endif	/* OPTIMIZATION_VISUALIZER */

#endif	/* _LOG_MANAGER_H_ */


