
#include "globalBuildParms.h"
#ifdef OPTIMIZATION_VISUALIZER


#include "matrix.h"
#include "LogManager.h"
#include <cstring>

LogManager globalLogManager;


	// call these addX... methods once at setup time
	// to get indices to be used in inner loops
int LogManager::addCategory(char * name)
{
	if (numCategories < MAX_CATEGORIES) {
		categoryNames[numCategories] = name;
		numEventsPerCategory[numCategories] = 0;
		numCategories++;
		return (numCategories -1);
	}
	return -1;
}


int LogManager::addAttribute(char * name)
{
	if (numAttributes < MAX_ATTRIBUTES) {
		attributeNames[numAttributes] = name;
		numAttributes++;
		return (numAttributes - 1);
	}
	return -1;
}

// Conjugate gradient method will call begin event before calls to 
// f->evaluate()
void LogManager::setCategory(int newCat)
{
	// This doesn't seem to matter
	if (currentCategory != newCat)  {
		if (numEventsPerCategory[currentCategory]> 0) {
			lastEventIndices[currentCategory][numLastEvents[currentCategory]]
				= numEventsPerCategory[currentCategory] - 1;
			numLastEvents[currentCategory]++;
		}
		currentCategory = newCat;
	}

}

void LogManager::beginEvent(const Vector * parameter)
{
	Vector * x;

	if ((currentCategory >= 0) && (currentCategory < numCategories)) {
		if (numTotalEvents < MAX_EVENTS) {
			if (parameter == NULL)
				x = NULL;
			else
				x = new Vector(*parameter);
			allParameters[(currentCategory * MAX_EVENTS) + numEventsPerCategory[currentCategory]] = x;
			allEventIndex[numTotalEvents] =
				pair<int, int>(currentCategory, numEventsPerCategory[currentCategory]);
			numTotalEvents++;
			numEventsPerCategory[currentCategory]++;
		} // What to do on array out of bounds?
	}
}

// f->evaluate() will call setValue
void LogManager::setValue(int idx, double value)
{
	if ((currentCategory >= 0) && (currentCategory < numCategories))
		eventsByCategory[currentCategory][numEventsPerCategory[currentCategory] - 1][idx]
			= value;
}

void LogManager::clearCategory(int cat)
{
	if ((cat >= 0) && (cat < numCategories)) {
		for (int i = 0; i < numEventsPerCategory[cat]; i++) {
			Vector * v = allParameters[(cat * MAX_EVENTS) + numEventsPerCategory[cat]];
			allParameters[(cat * MAX_EVENTS) + numEventsPerCategory[cat]] = 0;
			delete v;
		}
		numEventsPerCategory[cat] = 0;
		while ((numTotalEvents > 0) && (allEventIndex[numTotalEvents - 1].first == cat)) {
			numTotalEvents--;
		}
		// What about the allEventIndex... 
	}
}

LogManager::LogManager()
{
	int i;

	for (i = 0; i < MAX_CATEGORIES; i++) {
		categoryNames[i] = 0;
		numEventsPerCategory[i] = 0;
		numLastEvents[i] = 0;
	}
	for (i = 0; i < MAX_ATTRIBUTES; i++)
		attributeNames[i] =  0;

	// Do not initialize these guys since their couters are allready at zero
	// recordList[MAX_CATEGORIES] eventsByCategory;
	// int[MAX_EVENTS][MAX_CATEGORIES] lastEventIndices; 

	int numAttributes = 0;
	int numCategories = 0;
	int numTotalEvents = 0;

	problem = NULL;
	optimizer = NULL;

	memset(allParameters, 0, sizeof(Vector *) * MAX_EVENTS * MAX_CATEGORIES);
}

void LogManager::initialize()
{
	addCategory("Iterations");
	addCategory("Search");
	addCategory("Exploration");
	addCategory("Projection");
	addCategory("Staggered");

	// Attributes are hard-coded for figure stage binary...  better would
	// be to initialize them based on tuning values
	addAttribute("Objective Function");
	addAttribute("Image Match");
	addAttribute("Landmarks");
	addAttribute("Model");
	addAttribute("Interpenetration");
	addAttribute("Mahalanobis Distance");
	addAttribute("Param 4");
	addAttribute("Param 5");
	addAttribute("Param 6");
	addAttribute("Param 7");
	addAttribute("Param 8");
}




#endif	/* OPTIMIZATION_VISUALIZER */

