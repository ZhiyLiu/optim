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
include lib/zlib/CMakeFiles/zlib.dir/depend.make

# Include the progress variables for this target.
include lib/zlib/CMakeFiles/zlib.dir/progress.make

# Include the compile flags for this target's objects.
include lib/zlib/CMakeFiles/zlib.dir/flags.make

lib/zlib/CMakeFiles/zlib.dir/src/trees.c.o: lib/zlib/CMakeFiles/zlib.dir/flags.make
lib/zlib/CMakeFiles/zlib.dir/src/trees.c.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/trees.c
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/zlib/CMakeFiles/zlib.dir/src/trees.c.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/src/trees.c.o   -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/trees.c

lib/zlib/CMakeFiles/zlib.dir/src/trees.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/src/trees.c.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/trees.c > CMakeFiles/zlib.dir/src/trees.c.i

lib/zlib/CMakeFiles/zlib.dir/src/trees.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/src/trees.c.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/trees.c -o CMakeFiles/zlib.dir/src/trees.c.s

lib/zlib/CMakeFiles/zlib.dir/src/trees.c.o.requires:
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/trees.c.o.requires

lib/zlib/CMakeFiles/zlib.dir/src/trees.c.o.provides: lib/zlib/CMakeFiles/zlib.dir/src/trees.c.o.requires
	$(MAKE) -f lib/zlib/CMakeFiles/zlib.dir/build.make lib/zlib/CMakeFiles/zlib.dir/src/trees.c.o.provides.build
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/trees.c.o.provides

lib/zlib/CMakeFiles/zlib.dir/src/trees.c.o.provides.build: lib/zlib/CMakeFiles/zlib.dir/src/trees.c.o

lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.o: lib/zlib/CMakeFiles/zlib.dir/flags.make
lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/gzio.c
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/src/gzio.c.o   -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/gzio.c

lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/src/gzio.c.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/gzio.c > CMakeFiles/zlib.dir/src/gzio.c.i

lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/src/gzio.c.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/gzio.c -o CMakeFiles/zlib.dir/src/gzio.c.s

lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.o.requires:
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.o.requires

lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.o.provides: lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.o.requires
	$(MAKE) -f lib/zlib/CMakeFiles/zlib.dir/build.make lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.o.provides.build
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.o.provides

lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.o.provides.build: lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.o

lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.o: lib/zlib/CMakeFiles/zlib.dir/flags.make
lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/uncompr.c
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/src/uncompr.c.o   -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/uncompr.c

lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/src/uncompr.c.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/uncompr.c > CMakeFiles/zlib.dir/src/uncompr.c.i

lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/src/uncompr.c.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/uncompr.c -o CMakeFiles/zlib.dir/src/uncompr.c.s

lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.o.requires:
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.o.requires

lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.o.provides: lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.o.requires
	$(MAKE) -f lib/zlib/CMakeFiles/zlib.dir/build.make lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.o.provides.build
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.o.provides

lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.o.provides.build: lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.o

lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.o: lib/zlib/CMakeFiles/zlib.dir/flags.make
lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/zutil.c
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/src/zutil.c.o   -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/zutil.c

lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/src/zutil.c.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/zutil.c > CMakeFiles/zlib.dir/src/zutil.c.i

lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/src/zutil.c.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/zutil.c -o CMakeFiles/zlib.dir/src/zutil.c.s

lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.o.requires:
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.o.requires

lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.o.provides: lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.o.requires
	$(MAKE) -f lib/zlib/CMakeFiles/zlib.dir/build.make lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.o.provides.build
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.o.provides

lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.o.provides.build: lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.o

lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.o: lib/zlib/CMakeFiles/zlib.dir/flags.make
lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/adler32.c
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/src/adler32.c.o   -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/adler32.c

lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/src/adler32.c.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/adler32.c > CMakeFiles/zlib.dir/src/adler32.c.i

lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/src/adler32.c.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/adler32.c -o CMakeFiles/zlib.dir/src/adler32.c.s

lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.o.requires:
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.o.requires

lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.o.provides: lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.o.requires
	$(MAKE) -f lib/zlib/CMakeFiles/zlib.dir/build.make lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.o.provides.build
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.o.provides

lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.o.provides.build: lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.o

lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.o: lib/zlib/CMakeFiles/zlib.dir/flags.make
lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/inffast.c
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/src/inffast.c.o   -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/inffast.c

lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/src/inffast.c.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/inffast.c > CMakeFiles/zlib.dir/src/inffast.c.i

lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/src/inffast.c.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/inffast.c -o CMakeFiles/zlib.dir/src/inffast.c.s

lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.o.requires:
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.o.requires

lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.o.provides: lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.o.requires
	$(MAKE) -f lib/zlib/CMakeFiles/zlib.dir/build.make lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.o.provides.build
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.o.provides

lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.o.provides.build: lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.o

lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.o: lib/zlib/CMakeFiles/zlib.dir/flags.make
lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/deflate.c
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/src/deflate.c.o   -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/deflate.c

lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/src/deflate.c.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/deflate.c > CMakeFiles/zlib.dir/src/deflate.c.i

lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/src/deflate.c.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/deflate.c -o CMakeFiles/zlib.dir/src/deflate.c.s

lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.o.requires:
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.o.requires

lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.o.provides: lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.o.requires
	$(MAKE) -f lib/zlib/CMakeFiles/zlib.dir/build.make lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.o.provides.build
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.o.provides

lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.o.provides.build: lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.o

lib/zlib/CMakeFiles/zlib.dir/src/infback.c.o: lib/zlib/CMakeFiles/zlib.dir/flags.make
lib/zlib/CMakeFiles/zlib.dir/src/infback.c.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/infback.c
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/zlib/CMakeFiles/zlib.dir/src/infback.c.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/src/infback.c.o   -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/infback.c

lib/zlib/CMakeFiles/zlib.dir/src/infback.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/src/infback.c.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/infback.c > CMakeFiles/zlib.dir/src/infback.c.i

lib/zlib/CMakeFiles/zlib.dir/src/infback.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/src/infback.c.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/infback.c -o CMakeFiles/zlib.dir/src/infback.c.s

lib/zlib/CMakeFiles/zlib.dir/src/infback.c.o.requires:
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/infback.c.o.requires

lib/zlib/CMakeFiles/zlib.dir/src/infback.c.o.provides: lib/zlib/CMakeFiles/zlib.dir/src/infback.c.o.requires
	$(MAKE) -f lib/zlib/CMakeFiles/zlib.dir/build.make lib/zlib/CMakeFiles/zlib.dir/src/infback.c.o.provides.build
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/infback.c.o.provides

lib/zlib/CMakeFiles/zlib.dir/src/infback.c.o.provides.build: lib/zlib/CMakeFiles/zlib.dir/src/infback.c.o

lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.o: lib/zlib/CMakeFiles/zlib.dir/flags.make
lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/inftrees.c
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_9)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/src/inftrees.c.o   -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/inftrees.c

lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/src/inftrees.c.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/inftrees.c > CMakeFiles/zlib.dir/src/inftrees.c.i

lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/src/inftrees.c.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/inftrees.c -o CMakeFiles/zlib.dir/src/inftrees.c.s

lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.o.requires:
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.o.requires

lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.o.provides: lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.o.requires
	$(MAKE) -f lib/zlib/CMakeFiles/zlib.dir/build.make lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.o.provides.build
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.o.provides

lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.o.provides.build: lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.o

lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.o: lib/zlib/CMakeFiles/zlib.dir/flags.make
lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/crc32.c
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_10)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/src/crc32.c.o   -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/crc32.c

lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/src/crc32.c.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/crc32.c > CMakeFiles/zlib.dir/src/crc32.c.i

lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/src/crc32.c.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/crc32.c -o CMakeFiles/zlib.dir/src/crc32.c.s

lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.o.requires:
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.o.requires

lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.o.provides: lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.o.requires
	$(MAKE) -f lib/zlib/CMakeFiles/zlib.dir/build.make lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.o.provides.build
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.o.provides

lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.o.provides.build: lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.o

lib/zlib/CMakeFiles/zlib.dir/src/compress.c.o: lib/zlib/CMakeFiles/zlib.dir/flags.make
lib/zlib/CMakeFiles/zlib.dir/src/compress.c.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/compress.c
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_11)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/zlib/CMakeFiles/zlib.dir/src/compress.c.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/src/compress.c.o   -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/compress.c

lib/zlib/CMakeFiles/zlib.dir/src/compress.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/src/compress.c.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/compress.c > CMakeFiles/zlib.dir/src/compress.c.i

lib/zlib/CMakeFiles/zlib.dir/src/compress.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/src/compress.c.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/compress.c -o CMakeFiles/zlib.dir/src/compress.c.s

lib/zlib/CMakeFiles/zlib.dir/src/compress.c.o.requires:
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/compress.c.o.requires

lib/zlib/CMakeFiles/zlib.dir/src/compress.c.o.provides: lib/zlib/CMakeFiles/zlib.dir/src/compress.c.o.requires
	$(MAKE) -f lib/zlib/CMakeFiles/zlib.dir/build.make lib/zlib/CMakeFiles/zlib.dir/src/compress.c.o.provides.build
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/compress.c.o.provides

lib/zlib/CMakeFiles/zlib.dir/src/compress.c.o.provides.build: lib/zlib/CMakeFiles/zlib.dir/src/compress.c.o

lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.o: lib/zlib/CMakeFiles/zlib.dir/flags.make
lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.o: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/inflate.c
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/CMakeFiles $(CMAKE_PROGRESS_12)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.o"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/src/inflate.c.o   -c /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/inflate.c

lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/src/inflate.c.i"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/inflate.c > CMakeFiles/zlib.dir/src/inflate.c.i

lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/src/inflate.c.s"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/src/inflate.c -o CMakeFiles/zlib.dir/src/inflate.c.s

lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.o.requires:
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.o.requires

lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.o.provides: lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.o.requires
	$(MAKE) -f lib/zlib/CMakeFiles/zlib.dir/build.make lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.o.provides.build
.PHONY : lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.o.provides

lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.o.provides.build: lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.o

# Object files for target zlib
zlib_OBJECTS = \
"CMakeFiles/zlib.dir/src/trees.c.o" \
"CMakeFiles/zlib.dir/src/gzio.c.o" \
"CMakeFiles/zlib.dir/src/uncompr.c.o" \
"CMakeFiles/zlib.dir/src/zutil.c.o" \
"CMakeFiles/zlib.dir/src/adler32.c.o" \
"CMakeFiles/zlib.dir/src/inffast.c.o" \
"CMakeFiles/zlib.dir/src/deflate.c.o" \
"CMakeFiles/zlib.dir/src/infback.c.o" \
"CMakeFiles/zlib.dir/src/inftrees.c.o" \
"CMakeFiles/zlib.dir/src/crc32.c.o" \
"CMakeFiles/zlib.dir/src/compress.c.o" \
"CMakeFiles/zlib.dir/src/inflate.c.o"

# External object files for target zlib
zlib_EXTERNAL_OBJECTS =

libzlib.a: lib/zlib/CMakeFiles/zlib.dir/src/trees.c.o
libzlib.a: lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.o
libzlib.a: lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.o
libzlib.a: lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.o
libzlib.a: lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.o
libzlib.a: lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.o
libzlib.a: lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.o
libzlib.a: lib/zlib/CMakeFiles/zlib.dir/src/infback.c.o
libzlib.a: lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.o
libzlib.a: lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.o
libzlib.a: lib/zlib/CMakeFiles/zlib.dir/src/compress.c.o
libzlib.a: lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.o
libzlib.a: lib/zlib/CMakeFiles/zlib.dir/build.make
libzlib.a: lib/zlib/CMakeFiles/zlib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C static library ../../libzlib.a"
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && $(CMAKE_COMMAND) -P CMakeFiles/zlib.dir/cmake_clean_target.cmake
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/zlib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/zlib/CMakeFiles/zlib.dir/build: libzlib.a
.PHONY : lib/zlib/CMakeFiles/zlib.dir/build

lib/zlib/CMakeFiles/zlib.dir/requires: lib/zlib/CMakeFiles/zlib.dir/src/trees.c.o.requires
lib/zlib/CMakeFiles/zlib.dir/requires: lib/zlib/CMakeFiles/zlib.dir/src/gzio.c.o.requires
lib/zlib/CMakeFiles/zlib.dir/requires: lib/zlib/CMakeFiles/zlib.dir/src/uncompr.c.o.requires
lib/zlib/CMakeFiles/zlib.dir/requires: lib/zlib/CMakeFiles/zlib.dir/src/zutil.c.o.requires
lib/zlib/CMakeFiles/zlib.dir/requires: lib/zlib/CMakeFiles/zlib.dir/src/adler32.c.o.requires
lib/zlib/CMakeFiles/zlib.dir/requires: lib/zlib/CMakeFiles/zlib.dir/src/inffast.c.o.requires
lib/zlib/CMakeFiles/zlib.dir/requires: lib/zlib/CMakeFiles/zlib.dir/src/deflate.c.o.requires
lib/zlib/CMakeFiles/zlib.dir/requires: lib/zlib/CMakeFiles/zlib.dir/src/infback.c.o.requires
lib/zlib/CMakeFiles/zlib.dir/requires: lib/zlib/CMakeFiles/zlib.dir/src/inftrees.c.o.requires
lib/zlib/CMakeFiles/zlib.dir/requires: lib/zlib/CMakeFiles/zlib.dir/src/crc32.c.o.requires
lib/zlib/CMakeFiles/zlib.dir/requires: lib/zlib/CMakeFiles/zlib.dir/src/compress.c.o.requires
lib/zlib/CMakeFiles/zlib.dir/requires: lib/zlib/CMakeFiles/zlib.dir/src/inflate.c.o.requires
.PHONY : lib/zlib/CMakeFiles/zlib.dir/requires

lib/zlib/CMakeFiles/zlib.dir/clean:
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib && $(CMAKE_COMMAND) -P CMakeFiles/zlib.dir/cmake_clean.cmake
.PHONY : lib/zlib/CMakeFiles/zlib.dir/clean

lib/zlib/CMakeFiles/zlib.dir/depend:
	cd /playpen/workspace/newuoa/Correspondence_Oct7/myBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7 /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib /playpen/workspace/newuoa/Correspondence_Oct7/myBuild /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib /playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib/CMakeFiles/zlib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/zlib/CMakeFiles/zlib.dir/depend

