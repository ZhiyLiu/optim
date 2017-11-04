#include "P3DUndoList.h"

P3DUndoList::P3DUndoList(int _maxLevel)
{
	// P3DUndoList::addUndo() requires at least one level of undo
	if (_maxLevel > 1)
		maxLevel = _maxLevel;
	else
		maxLevel = 1;
}

P3DUndoList::~P3DUndoList()
{
	clearRedoList();
	clearUndoList();
}

void P3DUndoList::addUndo(M3DObject * objPtr, int markedPrimitiveId)
{
	model_t model;

    if(objPtr == NULL)
        return;

	// Clear the redo list
	clearRedoList();

	// Remove the end of the list if it is full
	if(undoObjList.size() >= maxLevel)
	{
		model = undoObjList.back();
		if(model.object != NULL)
			delete model.object;
		undoObjList.pop_back();
	}

	// Make a copy of the object
	model.object = objPtr->clone();
	model.markedPrimitive = markedPrimitiveId;

	// Add our new object to the list
	undoObjList.push_front(model);
}

M3DObject * P3DUndoList::undoObject(M3DObject * objPtr, int & markedPrimitiveId)
{
	model_t model;

	// Don't change the current object, if undo list is empty
	// or we are passed a NULL ptr
	if(undoObjList.empty() || objPtr == NULL)
		return objPtr;

	// Move current object to the redo list
	model.markedPrimitive = markedPrimitiveId;
	model.object = objPtr;
	redoObjList.push_front(model);

	// Replace object with saved object
	model = undoObjList.front();

	// Remove saved object from the undo list
	undoObjList.pop_front();

	markedPrimitiveId = model.markedPrimitive;
	return model.object;
}

M3DObject * P3DUndoList::deleteUndo(int & markedPrimitiveId)
{
	// Don't perform delete if undo list is empty
	if(undoObjList.empty())
		return NULL;

	// Get saved object for return
	model_t model = undoObjList.front();

	// Remove saved object from the undo list
	undoObjList.pop_front();

	markedPrimitiveId = model.markedPrimitive;
	return model.object;
}

M3DObject * P3DUndoList::redoObject(M3DObject * objPtr, int & markedPrimitiveId)
{
	model_t model;

	// Don't change the current object, if undo list is empty
	// or we are passed a NULL ptr
	if(redoObjList.empty() || objPtr == NULL)
		return objPtr;

	// Move current object to the undo list
	model.markedPrimitive = markedPrimitiveId;
	model.object = objPtr;
	undoObjList.push_front(model);

	// Get (and remove) front of redo list
	model = redoObjList.front();
	redoObjList.pop_front();

	markedPrimitiveId = model.markedPrimitive;
	return model.object;
}

void P3DUndoList::clearRedoList()
{
	int size = redoObjList.size();
	int i;

	model_t model;

	for(i = 0; i < size; i++)
	{
		model = redoObjList.back();
		if(model.object != NULL)
			delete model.object;

		redoObjList.pop_back();
	}
}

void P3DUndoList::clearUndoList()
{
	int size = undoObjList.size();
	int i;

	model_t model;

	for(i = 0; i < size; i++)
	{
		model = undoObjList.back();
		if(model.object != NULL)
			delete model.object;

		undoObjList.pop_back();
	}
}

