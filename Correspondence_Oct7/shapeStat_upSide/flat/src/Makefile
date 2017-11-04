PABLO = /NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2
PABLO_BUILD = /NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/binPablo2

INCLUDES = -I. -I$(PABLO)/lib/paul_code/include -I$(PABLO)/lib/m3d/include

LIBS = -lpaul_code -lm3d

GLDLIBS = 

LIBDIRS = -L$(PABLO_BUILD)

CPPFILES = flat.cpp

OBJECTS = $(CPPFILES:.cpp=.o)

DEPEND = $(PABLO_BUILD)/libpaul_code.a $(PABLO_BUILD)/libm3d.a

.SUFFIXES: .cpp .h .o


all:	../bin/flat 

.cpp.o :
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $<

../bin/flat: $(DEPEND) $(OBJECTS)
	$(CXX) -o ../bin/flat $(LINKFLAGS) $(OBJECTS) $(LIBDIRS) $(LIBS) $(GLDLIBS)

clean:
	rm -f *.o
	rm -f ../bin/flat

