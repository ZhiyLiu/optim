#ifndef MYLIST_H
#define MYLIST_H

/********************************************************************************/
/*																				*/
/*  	File	:  MyList.H														*/
/*																				*/
/*	Description:  class declaration for MyList template class					*/
/*				This has similar functionality to std::Vector<>, but is a bit	*/
/*				simpler.														*/
/*																				*/
/*	Project :  Rakshasa															*/
/*	Author  :  Stephen Ehmann (from Subdivision code project)					*/
/*				rewritten for Rakshasa by Thall									*/
/*	Date	:  5. Sept 1999														*/
/*																				*/
/*	Modifications:  															*/
/********************************************************************************/
#include <cstddef>

template<class Type>
class MyList
{
    int length;
    int max_length;
    Type* list;
public:
    MyList() { length = 0; max_length = 0; list = NULL; }
    ~MyList() {
        //cerr << "Starting to delete list" << endl;
		Destroy();
 //       delete [] list;
        //cerr << "Done deleting list" << endl;
    }

    int Length() { return length; }
    int Max_Length() { return max_length; }
    Type& operator[](int index);   // return reference
    Type* operator()(int index);   // return pointer
    void Create(int n);

    void Set_Length(int len) { length = len; }
    void Nullify();// { list = NULL; length = 0; max_length = 0; }
    void Destroy();/* {
        //cerr << "Inside destroy: list = " << list << endl;
        if (list != NULL)
            delete []list;
        list = NULL; length = 0; max_length = 0;
    }*/
};

/********************************************************************************/
/* operator[] -- return reference to indexed element of list			*/
/********************************************************************************/
template<class Type> inline Type& MyList<Type>::operator[](int index)
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
}

/********************************************************************************/
/* operator() -- return pointer to indexed element of list			*/
/********************************************************************************/
template<class Type> inline Type* MyList<Type>::operator()(int index)
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
}

/********************************************************************************/
/* create() -- allocate new list of desired size, if one doesn't already exist 	*/
/********************************************************************************/
template<class Type> inline void MyList<Type>::Create(int n)
{
#ifdef DEBUG
    if (max_length != 0) {
        cerr << "Warning: creating MyList when one already exists" << endl;
        delete list;
	self.Destroy();
    }
#endif
	
    if (n != 0) 
    {
        list = new Type[n]; 
        length = n; 
        max_length = n;
	}
}

#include "MyList.cpp"

#endif
