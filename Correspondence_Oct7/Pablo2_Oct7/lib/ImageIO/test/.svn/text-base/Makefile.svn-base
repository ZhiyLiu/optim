include ../../makeinclude

VPATH = src

TARGET = testIO

INCLUDES = -I../include -I../../m3d/include

CPPFILES = testIO.cpp

OBJECTS = $(CPPFILES:.cpp=.o)

LIBS = -L../lib -lImageIO -lAnalyze -lGipl -lMeta -lutility \
	-L../../m3d/lib -lm3d -L../../paul_code/lib -lregistry \
	-L../../zlib/lib -lzlib

all:	$(TARGET)

$(TARGET):	$(OBJECTS)
	$(CXX) -o $(TARGET) $(OBJECTS) $(LIBS)

clean:
	rm -f *.o
	rm -f $(TARGET)

.SUFFIXES: .cpp .h .o

.cpp.o :
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $<

