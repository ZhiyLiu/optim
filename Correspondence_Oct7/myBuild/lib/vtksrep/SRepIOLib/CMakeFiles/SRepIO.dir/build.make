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
CMAKE_COMMAND = /playpen/software/cmake-2.8.12-Linux-i386/bin/cmake

# The command to remove a file.
RM = /playpen/software/cmake-2.8.12-Linux-i386/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /playpen/software/cmake-2.8.12-Linux-i386/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /playpen/workspace/newuoa/Correspondence_Oct7/myBuild

# Include any dependencies generated for this target.
include lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/depend.make

# Include the progress variables for this target.
include lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/progress.make

# Include the compile flags for this target's objects.
include lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/flags.make

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o: lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/flags.make
lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepIOLib/vtkreadsrep.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/SRepIOLib && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepIOLib/vtkreadsrep.cpp

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/SRepIOLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepIOLib/vtkreadsrep.cpp > CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.i

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/SRepIOLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepIOLib/vtkreadsrep.cpp -o CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.s

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o.requires:
.PHONY : lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o.requires

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o.provides: lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o.requires
	$(MAKE) -f lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/build.make lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o.provides.build
.PHONY : lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o.provides

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o.provides.build: lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o: lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/flags.make
lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepIOLib/vtkwritesrep.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/SRepIOLib && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepIOLib/vtkwritesrep.cpp

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/SRepIOLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepIOLib/vtkwritesrep.cpp > CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.i

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/SRepIOLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepIOLib/vtkwritesrep.cpp -o CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.s

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o.requires:
.PHONY : lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o.requires

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o.provides: lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o.requires
	$(MAKE) -f lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/build.make lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o.provides.build
.PHONY : lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o.provides

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o.provides.build: lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o

# Object files for target SRepIO
SRepIO_OBJECTS = \
"CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o" \
"CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o"

# External object files for target SRepIO
SRepIO_EXTERNAL_OBJECTS =

libSRepIO.a: lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o
libSRepIO.a: lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o
libSRepIO.a: lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/build.make
libSRepIO.a: lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library ../../../libSRepIO.a"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/SRepIOLib && $(CMAKE_COMMAND) -P CMakeFiles/SRepIO.dir/cmake_clean_target.cmake
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/SRepIOLib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SRepIO.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/build: libSRepIO.a
.PHONY : lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/build

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/requires: lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkreadsrep.cpp.o.requires
lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/requires: lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/vtkwritesrep.cpp.o.requires
.PHONY : lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/requires

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/clean:
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/SRepIOLib && $(CMAKE_COMMAND) -P CMakeFiles/SRepIO.dir/cmake_clean.cmake
.PHONY : lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/clean

lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/depend:
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7 /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepIOLib /playpen/workspace/newuoa/Correspondence_Oct7/myBuild /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/SRepIOLib /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/vtksrep/SRepIOLib/CMakeFiles/SRepIO.dir/depend

