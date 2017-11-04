/********************************************************************************/
/*																				*/
/*  	File	:  PseudoSet.h													*/
/*																				*/
/*	Description:  Implements a fast set class (no union or difference operators)*/
/*			to substitute for slow Std. Template Set class.						*/
/*																				*/
/*			All function 
/*				hasmember()
/*				insert()
/*				del()
/*				clear()
/*			are O(1) time.														*/
/*																				*/
/*			There is an initial O(numkeys) initialization, and the routine		*/
/*			is O(numkeys) space, as it's just a smarter version of the old		*/
/*			boolean-vector-based set.  Uses integer keys in range				*/
/*			[0,...,numkeys - 1].												*/
/*																				*/
/*			Also an O(numelements) operator=() copy function.					*/
/*			Also an O(numelements) operator==() test for equality				*/
/*																				*/
/*	Project :  Seurat															*/
/*																				*/
/*	Dependencies:  iostreams for output.										*/
/*																				*/
/*	Author  :  A. Thall															*/
/*																				*/
/*	Date	:  8. May 2002														*/
/*	Modifications:																*/
/********************************************************************************/
const int NULLSETNAME = -1;

class elnode {
public:
   int setname;
   int el_location;
};

class PseudoSet
{
   elnode *thislist;
   int *thisset;
   int numkeys;
   int elcount;
   int currentsetname;
public:
   PseudoSet();
   PseudoSet(int nkeys);
   ~PseudoSet();

   void init(int nkeys);
   bool hasmember(int setel) const;
   bool insert(int setel);
   void del(int setel);
   void clear();

   int size() const { return elcount; }
   int maxsize() const { return numkeys; }
   bool empty() const { return (elcount == 0); }
   bool full() const { return (elcount == numkeys); }

   int operator[](int index) const { return thisset[index]; }
   PseudoSet& operator=(const PseudoSet& cpset);
   bool operator==(const PseudoSet& cpset) const;

   void printvals(char *message = NULL);
};

/*
 * constructors() -- initialize arrays and counters, set setnames to NULLSETNAME
 *		--if nkeys = zero, set pointers to NULL and do all
 *		allocations with a separate init command.
 */
inline PseudoSet::PseudoSet()
{
	thislist = NULL;
	thisset = NULL;
	elcount = 0;
	currentsetname = 0;
	numkeys = 0;
}

inline PseudoSet::PseudoSet(int nkeys)
{
	numkeys = nkeys;

	thislist = new elnode[numkeys];
	thisset = new int[numkeys];
	elcount = 0;
	currentsetname = 0;

	for (int i = 0; i < numkeys; i++)
		thislist[i].setname = NULLSETNAME;
}

/*
 * destructor()
 */
inline PseudoSet::~PseudoSet()
{
	if (thislist != NULL)
		delete []thislist;
	if (thisset != NULL)
		delete []thisset;
}

/*
 * init() -- initialize arrays to nkeys.  if numkeys == nkeys, do nothing
 */
inline void PseudoSet::init(int nkeys)
{
	if (numkeys != nkeys) {
		numkeys = nkeys;

		if (thislist != NULL)
			delete []thislist;
		if (thisset != NULL)
			delete []thisset;
		thislist = new elnode[numkeys];
		thisset = new int[numkeys];
		elcount = 0;
		currentsetname = 0;

		for (int i = 0; i < numkeys; i++)
			thislist[i].setname = NULLSETNAME;
	}
	else
		clear();
}

/*
 * operator=() -- copy operator ( O(elcount) )
 */
inline PseudoSet& ThallCode::PseudoSet::operator=(const PseudoSet& cpset)
{
	if (numkeys != cpset.numkeys)
		init(cpset.numkeys);
	else {
		if (currentsetname == INT_MAX) {
			currentsetname = 0;
			for (int i = 0; i < numkeys; i++)
				thislist[i].setname = NULLSETNAME;
		}
		else
			currentsetname += 1;
	}

	for (int i = 0; i < cpset.size(); i++) {
		thisset[i] = cpset[i];
		thislist[thisset[i]].el_location = i;
		thislist[thisset[i]].setname = currentsetname;
	}
	elcount = cpset.size();

	return (*this);
}

/*
 * operator==() -- comparision operator ( O(elcount) )
 */
inline bool PseudoSet::operator==(const PseudoSet& cpset) const
{
	if (numkeys != cpset.numkeys || elcount != cpset.elcount)
		return false;
	else {
		for (int i = 0; i < cpset.elcount; i++) {

			if (!hasmember(cpset[i]) || !cpset.hasmember((*this)[i]))
				return false;
		}
	}

	return true;
}

/*
 * find() -- return true if setel is an element of the set, else false
 */
inline bool PseudoSet::hasmember(int setel) const
{
	return (thislist[setel].setname == currentsetname);
}

/*
 * insert() -- add the setel to the set, returning false if it is already
 *            a set element.
 */
inline bool PseudoSet::insert(int setel)
{
	if (thislist[setel].setname == currentsetname)
		return false;
	else {
		thislist[setel].setname = currentsetname;
		thislist[setel].el_location = elcount;
		thisset[elcount++] = setel;
		return true;
	}
}

/*
 * del() -- delete the element from the set.  If not a member, does nothing
 */
inline void PseudoSet::del(int setel)
{
	if (thislist[setel].setname == currentsetname) {

		int el_loc = thislist[setel].el_location;

		// put the last element in the set (at elcount-1) in position el_loc,
		//   change its pointer in thislist to el_loc,
		//   and set thislist[setel].setname to the NULLSETNAME
		thislist[setel].setname = NULLSETNAME;

		if (elcount > 1) {
			int lastelement = thisset[elcount - 1];
			thisset[el_loc] = lastelement;
			thislist[lastelement].el_location = el_loc;
		}
		elcount--;
	}
}

/*
 * clear() -- delete all elements from the set
 *      ---reinitializes if currentsetname == REALLYBIGINT, to avoid name overlap
 */
inline void PseudoSet::clear()
{
	elcount = 0;
	if (currentsetname != INT_MAX)
		currentsetname++;
	else {
		currentsetname = 0;
		for (int i = 0; i < numkeys; i++)
			thislist[i].setname = NULLSETNAME;
	}
}

/*
 * printvals()
 */
inline void PseudoSet::printvals(char *message)
{
	if (message != NULL)
		std::cerr << message << '\n';

	std::cerr << '[';
	for (int i = 0; i < elcount; i++) {
		std::cerr << thisset[i];
		if (i != elcount - 1)
			std::cerr << ", "; 
	}
	std::cerr << "]\n";
}


