#ifndef P3D_UNDO_LIST_H
#define P3D_UNDO_LIST_H

#include <list>
#include "M3DObject.h"

#define DEFAULT_MAX_UNDO_LEVEL 15

class P3DUndoList
{
	public:

		P3DUndoList(int _maxLevel = DEFAULT_MAX_UNDO_LEVEL);
		~P3DUndoList();

		void addUndo(M3DObject * objPtr, int markedPrimitiveId = -1);
		M3DObject * undoObject(M3DObject * objPtr, int & markedPrimitiveId);
		M3DObject * redoObject(M3DObject * objPtr, int & markedPrimitiveId);
		M3DObject * deleteUndo(int & markedPrimitiveId);

		void clearRedoList();
		void clearUndoList();
		int size() { return undoObjList.size(); }
		int maxSize() { return maxLevel; }

	private:

		struct model_t {
			M3DObject * object;
			int	markedPrimitive;
		};

		std::list<model_t> undoObjList;
		std::list<model_t> redoObjList;
		int maxLevel;
};

#endif

