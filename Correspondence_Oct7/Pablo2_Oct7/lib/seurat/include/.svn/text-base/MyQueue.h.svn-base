#ifndef _MYQUEUE_H_
#define _MYQUEUE_H_


#define MY_DATATYPE int

class MyQueue
{
	//MY_REAL x, y, z;
	int size, head, tail;
	int *contents;

public:
	//////////////////////////////////////////////////////////////
	// Constructors
	MyQueue()
	{
		head=tail=0;
		size=0;
		contents=NULL;
	}
	MyQueue(int s)
	{
		contents=NULL;
		Initialize(s);
	}
	~MyQueue()
	{
		if(contents!=NULL)
			delete []contents;
	}

	void Initialize(int s)
	{
		size=s;
		head=tail=0;
		if(contents!=NULL)
			delete []contents;
		contents=new int[size];
	}
	void Clear()
	{
		head=tail=0;
	}

	//////////////////////////////////////////////////////////////
	// En/De-que
	MY_DATATYPE Deque()
	{
		MY_DATATYPE t;
		if(!Empty())
		{
			t=contents[head];
			head=(head+1)%size;
			return t;
		}
		else
			return -1;
	}
	bool Enque(int e)
	{
		if(Full())
			return false;
		else
		{
			contents[tail]=e;
			tail=(tail+1)%size;
			return true;
		}
	}
	bool Empty()
	{
		return (head==tail);
	}
	bool Full()
	{
		return (tail==head-1 || (tail==size-1 && head==0));
	}

	//////////////////////////////////////////////////////////////
	// Debugging tools
/*
	void print()
	{
		int i;
		int end;	// Uninitialized variable!

		printf("Quaue\tMaxsize:%d\tsize:%d\n", size, end-head);
		if(!Empty())
		{
			if(head<tail)
				end=tail;
			else
				end=size+tail;
			for(i=head; i<end; i++)
				printf("%d ", contents[i%size]);
		}
		printf("\n");
	}
*/
};

#endif

