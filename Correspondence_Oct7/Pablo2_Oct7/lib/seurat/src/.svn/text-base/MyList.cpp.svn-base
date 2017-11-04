#ifndef MYLIST_CPP
#define MYLIST_CPP

#include "MyList.h"

template<class Type> inline void MyList<Type>::Nullify() { 

		// dibyendu : checking if list already has some memory occupied

		if( list != NULL )
			delete [] list ;

		list = NULL; length = 0; max_length = 0; 
}

template<class Type> inline void MyList<Type>::Destroy() {

        if (list != NULL)
            delete [] list;

        list = NULL; length = 0; max_length = 0;
    }


/********************************************************************************/
/* operator[] -- return reference to indexed element of list			*/
/********************************************************************************/
/*template<class Type> inline Type& MyList<Type>::operator[](int index)
{
#ifdef DEBUG
    if(index < length && index >= 0)
        return list[index];
    else {
        cerr << "MyList Warning: indexing wrong into the MyList[]: index = "
             << index << " length = " << length << endl;
        return list[0];
    }
#else
    return list[index];
#endif
}*/

/********************************************************************************/
/* operator() -- return pointer to indexed element of list			*/
/********************************************************************************/
/*template<class Type> inline Type* MyList<Type>::operator()(int index)
{
#ifdef DEBUG
    if (index < length && index >= 0) 
        return list+index;
    else {
        cerr << "MyList Warning: indexing wrong into the MyList(): index = "
             << index << " length = " << length << endl;
        return list;
    }
#else
    return list+index;
#endif
}*/

/********************************************************************************/
/* create() -- allocate new list of desired size, if one doesn't already exist 	*/
/********************************************************************************/
/*template<class Type> inline void MyList<Type>::Create(int n)
{
#ifdef DEBUG
    if (max_length != 0) {
        cerr << "Warning: creating MyList when one already exists" << endl;
        //delete list;
	self.Destroy();
    }
#endif
    if (n != 0) 
        list = new Type[n]; length = n; max_length = n;
}*/

#endif
