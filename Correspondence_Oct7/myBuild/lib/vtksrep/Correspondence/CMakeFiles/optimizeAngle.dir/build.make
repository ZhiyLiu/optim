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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /playpen/workspace/newuoa/Correspondence_Oct7/myBuild

# Include any dependencies generated for this target.
include lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/depend.make

# Include the progress variables for this target.
include lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/progress.make

# Include the compile flags for this target's objects.
include lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/flags.make

lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o: lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/flags.make
lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/OptimizeAngle.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/Correspondence && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/OptimizeAngle.cpp

lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/Correspondence && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/OptimizeAngle.cpp > CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.i

lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/Correspondence && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/OptimizeAngle.cpp -o CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.s

lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o.requires:
.PHONY : lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o.requires

lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o.provides: lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o.requires
	$(MAKE) -f lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/build.make lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o.provides.build
.PHONY : lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o.provides

lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o.provides.build: lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o

# Object files for target optimizeAngle
optimizeAngle_OBJECTS = \
"CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o"

# External object files for target optimizeAngle
optimizeAngle_EXTERNAL_OBJECTS =

optimizeAngle: lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o
optimizeAngle: lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/build.make
optimizeAngle: /playpen/software/vtk_build/bin/libvtkCommon.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkFiltering.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkImaging.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkGraphics.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkGenericFiltering.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkIO.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkRendering.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkVolumeRendering.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkHybrid.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkWidgets.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkInfovis.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkGeovis.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkViews.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkCharts.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkFiltering.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkCommon.a
optimizeAngle: libSRep.a
optimizeAngle: libSRepIO.a
optimizeAngle: libSRepInterpolation.a
optimizeAngle: libSRepVisualization.a
optimizeAngle: libcalquantum.a
optimizeAngle: libVisualization.a
optimizeAngle: libm3d.a
optimizeAngle: libregister.a
optimizeAngle: libCorrespondence.a
optimizeAngle: libVisualization.a
optimizeAngle: libCorrespondence.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkGenericFiltering.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkGeovis.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkproj4.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkCharts.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkViews.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkInfovis.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkWidgets.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkVolumeRendering.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkHybrid.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkRendering.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkImaging.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkIO.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkDICOMParser.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkNetCDF_cxx.a
optimizeAngle: /playpen/software/vtk_build/bin/libLSDyna.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkmetaio.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtksqlite.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkpng.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtktiff.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkjpeg.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkexpat.a
optimizeAngle: /usr/lib/x86_64-linux-gnu/libXt.so
optimizeAngle: /playpen/software/vtk_build/bin/libvtkexoIIc.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkNetCDF.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkhdf5_hl.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkhdf5.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtklibxml2.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkzlib.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkalglib.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkftgl.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkfreetype.a
optimizeAngle: libSRepIO.a
optimizeAngle: libSRepVisualization.a
optimizeAngle: libcalquantum.a
optimizeAngle: libregister.a
optimizeAngle: libImageIO.a
optimizeAngle: libpaul_code.a
optimizeAngle: /usr/lib/libblas.so
optimizeAngle: /usr/lib/liblapack.so
optimizeAngle: /usr/lib/libblas.so
optimizeAngle: /usr/lib/liblapack.so
optimizeAngle: libmatch.a
optimizeAngle: libSRep.a
optimizeAngle: libSRepInterpolation.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkGraphics.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkFiltering.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkCommon.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtksys.a
optimizeAngle: /playpen/software/vtk_build/bin/libvtkverdict.a
optimizeAngle: libm3d.a
optimizeAngle: libzlib.a
optimizeAngle: libseurat.a
optimizeAngle: /usr/lib/x86_64-linux-gnu/libGLU.so
optimizeAngle: /usr/lib/x86_64-linux-gnu/libGL.so
optimizeAngle: /usr/lib/x86_64-linux-gnu/libSM.so
optimizeAngle: /usr/lib/x86_64-linux-gnu/libICE.so
optimizeAngle: /usr/lib/x86_64-linux-gnu/libX11.so
optimizeAngle: /usr/lib/x86_64-linux-gnu/libXext.so
optimizeAngle: /usr/lib/x86_64-linux-gnu/libuuid.so
optimizeAngle: lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../../optimizeAngle"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/Correspondence && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/optimizeAngle.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/build: optimizeAngle
.PHONY : lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/build

lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/requires: lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/OptimizeAngle.cpp.o.requires
.PHONY : lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/requires

lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/clean:
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/Correspondence && $(CMAKE_COMMAND) -P CMakeFiles/optimizeAngle.dir/cmake_clean.cmake
.PHONY : lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/clean

lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/depend:
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7 /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence /playpen/workspace/newuoa/Correspondence_Oct7/myBuild /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/Correspondence /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/vtksrep/Correspondence/CMakeFiles/optimizeAngle.dir/depend
