include ../../makeinclude

MEX = /afs/isis.unc.edu/pkg/matlab-2008a/.install/i386_linux26/bin/mex

INCLUDES = -I../include -I../../include -I../Templates -I../utilities/include \
	   -I../Analyze/include -I../Dicom/include -I../Meta/include \
	   -I../Gipl/include -I../PlanIm/include -I../planio/include \
	   -I../PlanIM/include -I../../m3d/include -I../../paul_code/include \

LIBS = -lplanes -lregister -lmatch -lregister -lseurat -loptima \
	-lmatrix -lImageIO -lm3d -lz -lregistry -lGipl -lMeta -lAnalyze \
	-lPlanIm -lioUtilities -lflvw $(FLTK) $(CLAPACK) -lF77 $(SYSLIBS)

LIBS_B = -lplanes_b -lregister_b -lmatch_b -lregister_b -lseurat_b -lm3d_b \
		-loptima_b -lregistry_b -lmatrix_b -lzlib -lImageIO -lGipl \
		-lMeta -lAnalyze -lPlanIm -lioUtilities -lflvw $(FLTK) \
		$(CLAPACK) $(SYSLIBS)

GLDLIBS = -L$(GL) $(GL_LIBS) -L$(X11) -lX11 -lXext -lm

LIBDIRS = -L../../register/lib -L../../seurat/lib -L../../planes/lib \
		-L../../match/lib -L../../m3d/lib -L../../paul_code/lib \
		-L../../zlib/lib -L../../ImageIO/lib \
		-L$(CLAPACKDIR) -L$(CLAPACKDIR)/F2CLIBS \
		-L$(FLVWDIR)/lib -L$(FLTKDIR)/$(FLTKLIBDIR)

# CPPFILES = ImageIO.cpp AllImageIO.cpp AllImageIO2.cpp

# OBJECTS = $(CPPFILES:.cpp=.o)

all: readmedim.mexglx writemedim.mexglx

READ_MEX_CPPFILE = readmedim.cpp

readmedim.mexglx: $(READ_MEX_CPPFILE) 
	$(MEX) $(INCLUDES) $(LIBDIRS) $(LIBS) $(READ_MEX_CPPFILE)

WRITE_MEX_CPPFILE = writemedim.cpp

writemedim.mexglx: $(WRITE_MEX_CPPFILE) 
	$(MEX) $(INCLUDES) $(LIBDIRS) $(LIBS) $(WRITE_MEX_CPPFILE)

clean:
	rm -f *.o
	rm -f ../lib/libImageIO.a

.SUFFIXES: .cpp .h .o

# .cpp.o :
#	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $<
