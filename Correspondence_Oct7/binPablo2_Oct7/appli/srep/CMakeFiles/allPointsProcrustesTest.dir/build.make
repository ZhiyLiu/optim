# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /NIRAL/tools/bin_linux64/cmake

# The command to remove a file.
RM = /NIRAL/tools/bin_linux64/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /NIRAL/tools/bin_linux64/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /work/ltu/Correspondence_Oct7/Pablo2_Oct7

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /work/ltu/Correspondence_Oct7/binPablo2_Oct7

# Include any dependencies generated for this target.
include appli/srep/CMakeFiles/allPointsProcrustesTest.dir/depend.make

# Include the progress variables for this target.
include appli/srep/CMakeFiles/allPointsProcrustesTest.dir/progress.make

# Include the compile flags for this target's objects.
include appli/srep/CMakeFiles/allPointsProcrustesTest.dir/flags.make

appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o: appli/srep/CMakeFiles/allPointsProcrustesTest.dir/flags.make
appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o: /work/ltu/Correspondence_Oct7/Pablo2_Oct7/appli/srep/allPointsProcrustesTest.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /work/ltu/Correspondence_Oct7/binPablo2_Oct7/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o"
	cd /work/ltu/Correspondence_Oct7/binPablo2_Oct7/appli/srep && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o -c /work/ltu/Correspondence_Oct7/Pablo2_Oct7/appli/srep/allPointsProcrustesTest.cpp

appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.i"
	cd /work/ltu/Correspondence_Oct7/binPablo2_Oct7/appli/srep && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /work/ltu/Correspondence_Oct7/Pablo2_Oct7/appli/srep/allPointsProcrustesTest.cpp > CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.i

appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.s"
	cd /work/ltu/Correspondence_Oct7/binPablo2_Oct7/appli/srep && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /work/ltu/Correspondence_Oct7/Pablo2_Oct7/appli/srep/allPointsProcrustesTest.cpp -o CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.s

appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o.requires:
.PHONY : appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o.requires

appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o.provides: appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o.requires
	$(MAKE) -f appli/srep/CMakeFiles/allPointsProcrustesTest.dir/build.make appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o.provides.build
.PHONY : appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o.provides

appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o.provides.build: appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o

# Object files for target allPointsProcrustesTest
allPointsProcrustesTest_OBJECTS = \
"CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o"

# External object files for target allPointsProcrustesTest
allPointsProcrustesTest_EXTERNAL_OBJECTS =

allPointsProcrustesTest: appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkFiltering.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkCommon.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkRendering.a
allPointsProcrustesTest: libSRep.a
allPointsProcrustesTest: libSRepInterpolation.a
allPointsProcrustesTest: libSRepVisualization.a
allPointsProcrustesTest: libCorrespondence.a
allPointsProcrustesTest: libSRepIO.a
allPointsProcrustesTest: libregister.a
allPointsProcrustesTest: libVisualization.a
allPointsProcrustesTest: libCorrespondence.a
allPointsProcrustesTest: libVisualization.a
allPointsProcrustesTest: libSRepIO.a
allPointsProcrustesTest: libcalquantum.a
allPointsProcrustesTest: libSRepVisualization.a
allPointsProcrustesTest: libregister.a
allPointsProcrustesTest: libImageIO.a
allPointsProcrustesTest: libpaul_code.a
allPointsProcrustesTest: /usr/lib64/libblas.so
allPointsProcrustesTest: /usr/lib64/liblapack.so
allPointsProcrustesTest: /usr/lib64/libblas.so
allPointsProcrustesTest: /usr/lib64/liblapack.so
allPointsProcrustesTest: libmatch.a
allPointsProcrustesTest: libSRep.a
allPointsProcrustesTest: libSRepInterpolation.a
allPointsProcrustesTest: libseurat.a
allPointsProcrustesTest: /usr/lib64/libGLU.so
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkGenericFiltering.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkGeovis.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkproj4.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkCharts.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkViews.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkInfovis.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkWidgets.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkVolumeRendering.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkHybrid.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkRendering.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkGraphics.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkverdict.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkImaging.a
allPointsProcrustesTest: /usr/lib64/libXt.so
allPointsProcrustesTest: /usr/lib64/libSM.so
allPointsProcrustesTest: /usr/lib64/libICE.so
allPointsProcrustesTest: /usr/lib64/libX11.so
allPointsProcrustesTest: /usr/lib64/libXext.so
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkIO.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkFiltering.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkCommon.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkDICOMParser.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkNetCDF_cxx.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libLSDyna.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtksys.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkmetaio.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtksqlite.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkpng.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtktiff.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkjpeg.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkexpat.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkftgl.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkfreetype.a
allPointsProcrustesTest: /usr/lib64/libGL.so
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkexoIIc.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkNetCDF.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkhdf5_hl.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkhdf5.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtklibxml2.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkzlib.a
allPointsProcrustesTest: /devel/linux/Pablo/install/lib/vtk-5.10/libvtkalglib.a
allPointsProcrustesTest: libm3d.a
allPointsProcrustesTest: libzlib.a
allPointsProcrustesTest: /usr/lib64/libuuid.so
allPointsProcrustesTest: appli/srep/CMakeFiles/allPointsProcrustesTest.dir/build.make
allPointsProcrustesTest: appli/srep/CMakeFiles/allPointsProcrustesTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../allPointsProcrustesTest"
	cd /work/ltu/Correspondence_Oct7/binPablo2_Oct7/appli/srep && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/allPointsProcrustesTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
appli/srep/CMakeFiles/allPointsProcrustesTest.dir/build: allPointsProcrustesTest
.PHONY : appli/srep/CMakeFiles/allPointsProcrustesTest.dir/build

appli/srep/CMakeFiles/allPointsProcrustesTest.dir/requires: appli/srep/CMakeFiles/allPointsProcrustesTest.dir/allPointsProcrustesTest.cpp.o.requires
.PHONY : appli/srep/CMakeFiles/allPointsProcrustesTest.dir/requires

appli/srep/CMakeFiles/allPointsProcrustesTest.dir/clean:
	cd /work/ltu/Correspondence_Oct7/binPablo2_Oct7/appli/srep && $(CMAKE_COMMAND) -P CMakeFiles/allPointsProcrustesTest.dir/cmake_clean.cmake
.PHONY : appli/srep/CMakeFiles/allPointsProcrustesTest.dir/clean

appli/srep/CMakeFiles/allPointsProcrustesTest.dir/depend:
	cd /work/ltu/Correspondence_Oct7/binPablo2_Oct7 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work/ltu/Correspondence_Oct7/Pablo2_Oct7 /work/ltu/Correspondence_Oct7/Pablo2_Oct7/appli/srep /work/ltu/Correspondence_Oct7/binPablo2_Oct7 /work/ltu/Correspondence_Oct7/binPablo2_Oct7/appli/srep /work/ltu/Correspondence_Oct7/binPablo2_Oct7/appli/srep/CMakeFiles/allPointsProcrustesTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : appli/srep/CMakeFiles/allPointsProcrustesTest.dir/depend

