include ../../makeinclude

VPATH = src

INCLUDES = -Iinclude -I../include -I../utilities/include

CPPFILES = read_analyze.cpp write_analyze.cpp Analyze.cpp

OBJECTS = $(CPPFILES:.cpp=.o)

LIB_OBJECTS = read_analyze.o write_analyze.o Analyze.o

../lib/libAnalyze.a: $(OBJECTS)
	$(LIBCOMMAND) ../lib/libAnalyze.a $(OBJECTS)

libs:	$(LIB_OBJECTS)
	$(LIBCOMMAND) ../lib/libAnalyze.a $(LIB_OBJECTS)

clean:
	rm -f *.o
	rm -f ../lib/libAnalyze.a

.SUFFIXES: .cpp .h .o

.cpp.o :
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $<
