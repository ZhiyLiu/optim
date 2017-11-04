include ../makeinclude

VPATH = src

TARGET = ../m3d/include/pablo_version.h

INCLUDES = -I.

CPPFILES = version.cpp

OBJECTS = $(CPPFILES:.cpp=.o)


libs:	all

libs_b:	all

regular:	all

binary:	all


all:	bin/version
	rm -f $(TARGET)
	bin/version $(TARGET) ../gui/src/.svn/entries \
	    ../m3d/include/.svn/entries ../m3d/src/.svn/entries \
	    ../match/include/.svn/entries ../match/src/.svn/entries \
	    ../paul_code/include/.svn/entries ../paul_code/src/.svn/entries \
	    ../planes/include/.svn/entries ../planes/src/.svn/entries \
	    ../register/include/.svn/entries ../register/src/.svn/entries \
	    ../seurat/include/.svn/entries ../seurat/src/.svn/entries

bin/version:	$(OBJECTS)
	$(CXX) $(LDFLAGS) -o bin/version $(OBJECTS)

clean:
	rm -f *.o
	rm -f bin/version
	rm -f $(TARGET)
	cp src/pablo_version.h $(TARGET)

