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
CMAKE_SOURCE_DIR = /home/ksb11/GlobalOptimization/GlobalOpt

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ksb11/GlobalOptimization/GlobalOpt/build

# Include any dependencies generated for this target.
include src/LBFGSB/CMakeFiles/LBFGSB.dir/depend.make

# Include the progress variables for this target.
include src/LBFGSB/CMakeFiles/LBFGSB.dir/progress.make

# Include the compile flags for this target's objects.
include src/LBFGSB/CMakeFiles/LBFGSB.dir/flags.make

src/LBFGSB/CMakeFiles/LBFGSB.dir/blas.f.o: src/LBFGSB/CMakeFiles/LBFGSB.dir/flags.make
src/LBFGSB/CMakeFiles/LBFGSB.dir/blas.f.o: ../src/LBFGSB/blas.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ksb11/GlobalOptimization/GlobalOpt/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object src/LBFGSB/CMakeFiles/LBFGSB.dir/blas.f.o"
	cd /home/ksb11/GlobalOptimization/GlobalOpt/build/src/LBFGSB && /global/apps/gcc/4.8.4/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/ksb11/GlobalOptimization/GlobalOpt/src/LBFGSB/blas.f -o CMakeFiles/LBFGSB.dir/blas.f.o

src/LBFGSB/CMakeFiles/LBFGSB.dir/blas.f.o.requires:
.PHONY : src/LBFGSB/CMakeFiles/LBFGSB.dir/blas.f.o.requires

src/LBFGSB/CMakeFiles/LBFGSB.dir/blas.f.o.provides: src/LBFGSB/CMakeFiles/LBFGSB.dir/blas.f.o.requires
	$(MAKE) -f src/LBFGSB/CMakeFiles/LBFGSB.dir/build.make src/LBFGSB/CMakeFiles/LBFGSB.dir/blas.f.o.provides.build
.PHONY : src/LBFGSB/CMakeFiles/LBFGSB.dir/blas.f.o.provides

src/LBFGSB/CMakeFiles/LBFGSB.dir/blas.f.o.provides.build: src/LBFGSB/CMakeFiles/LBFGSB.dir/blas.f.o

src/LBFGSB/CMakeFiles/LBFGSB.dir/lbfgsb.f.o: src/LBFGSB/CMakeFiles/LBFGSB.dir/flags.make
src/LBFGSB/CMakeFiles/LBFGSB.dir/lbfgsb.f.o: ../src/LBFGSB/lbfgsb.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ksb11/GlobalOptimization/GlobalOpt/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object src/LBFGSB/CMakeFiles/LBFGSB.dir/lbfgsb.f.o"
	cd /home/ksb11/GlobalOptimization/GlobalOpt/build/src/LBFGSB && /global/apps/gcc/4.8.4/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/ksb11/GlobalOptimization/GlobalOpt/src/LBFGSB/lbfgsb.f -o CMakeFiles/LBFGSB.dir/lbfgsb.f.o

src/LBFGSB/CMakeFiles/LBFGSB.dir/lbfgsb.f.o.requires:
.PHONY : src/LBFGSB/CMakeFiles/LBFGSB.dir/lbfgsb.f.o.requires

src/LBFGSB/CMakeFiles/LBFGSB.dir/lbfgsb.f.o.provides: src/LBFGSB/CMakeFiles/LBFGSB.dir/lbfgsb.f.o.requires
	$(MAKE) -f src/LBFGSB/CMakeFiles/LBFGSB.dir/build.make src/LBFGSB/CMakeFiles/LBFGSB.dir/lbfgsb.f.o.provides.build
.PHONY : src/LBFGSB/CMakeFiles/LBFGSB.dir/lbfgsb.f.o.provides

src/LBFGSB/CMakeFiles/LBFGSB.dir/lbfgsb.f.o.provides.build: src/LBFGSB/CMakeFiles/LBFGSB.dir/lbfgsb.f.o

src/LBFGSB/CMakeFiles/LBFGSB.dir/linpack.f.o: src/LBFGSB/CMakeFiles/LBFGSB.dir/flags.make
src/LBFGSB/CMakeFiles/LBFGSB.dir/linpack.f.o: ../src/LBFGSB/linpack.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ksb11/GlobalOptimization/GlobalOpt/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object src/LBFGSB/CMakeFiles/LBFGSB.dir/linpack.f.o"
	cd /home/ksb11/GlobalOptimization/GlobalOpt/build/src/LBFGSB && /global/apps/gcc/4.8.4/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/ksb11/GlobalOptimization/GlobalOpt/src/LBFGSB/linpack.f -o CMakeFiles/LBFGSB.dir/linpack.f.o

src/LBFGSB/CMakeFiles/LBFGSB.dir/linpack.f.o.requires:
.PHONY : src/LBFGSB/CMakeFiles/LBFGSB.dir/linpack.f.o.requires

src/LBFGSB/CMakeFiles/LBFGSB.dir/linpack.f.o.provides: src/LBFGSB/CMakeFiles/LBFGSB.dir/linpack.f.o.requires
	$(MAKE) -f src/LBFGSB/CMakeFiles/LBFGSB.dir/build.make src/LBFGSB/CMakeFiles/LBFGSB.dir/linpack.f.o.provides.build
.PHONY : src/LBFGSB/CMakeFiles/LBFGSB.dir/linpack.f.o.provides

src/LBFGSB/CMakeFiles/LBFGSB.dir/linpack.f.o.provides.build: src/LBFGSB/CMakeFiles/LBFGSB.dir/linpack.f.o

src/LBFGSB/CMakeFiles/LBFGSB.dir/timer.f.o: src/LBFGSB/CMakeFiles/LBFGSB.dir/flags.make
src/LBFGSB/CMakeFiles/LBFGSB.dir/timer.f.o: ../src/LBFGSB/timer.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ksb11/GlobalOptimization/GlobalOpt/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object src/LBFGSB/CMakeFiles/LBFGSB.dir/timer.f.o"
	cd /home/ksb11/GlobalOptimization/GlobalOpt/build/src/LBFGSB && /global/apps/gcc/4.8.4/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/ksb11/GlobalOptimization/GlobalOpt/src/LBFGSB/timer.f -o CMakeFiles/LBFGSB.dir/timer.f.o

src/LBFGSB/CMakeFiles/LBFGSB.dir/timer.f.o.requires:
.PHONY : src/LBFGSB/CMakeFiles/LBFGSB.dir/timer.f.o.requires

src/LBFGSB/CMakeFiles/LBFGSB.dir/timer.f.o.provides: src/LBFGSB/CMakeFiles/LBFGSB.dir/timer.f.o.requires
	$(MAKE) -f src/LBFGSB/CMakeFiles/LBFGSB.dir/build.make src/LBFGSB/CMakeFiles/LBFGSB.dir/timer.f.o.provides.build
.PHONY : src/LBFGSB/CMakeFiles/LBFGSB.dir/timer.f.o.provides

src/LBFGSB/CMakeFiles/LBFGSB.dir/timer.f.o.provides.build: src/LBFGSB/CMakeFiles/LBFGSB.dir/timer.f.o

# Object files for target LBFGSB
LBFGSB_OBJECTS = \
"CMakeFiles/LBFGSB.dir/blas.f.o" \
"CMakeFiles/LBFGSB.dir/lbfgsb.f.o" \
"CMakeFiles/LBFGSB.dir/linpack.f.o" \
"CMakeFiles/LBFGSB.dir/timer.f.o"

# External object files for target LBFGSB
LBFGSB_EXTERNAL_OBJECTS =

src/LBFGSB/libLBFGSB.a: src/LBFGSB/CMakeFiles/LBFGSB.dir/blas.f.o
src/LBFGSB/libLBFGSB.a: src/LBFGSB/CMakeFiles/LBFGSB.dir/lbfgsb.f.o
src/LBFGSB/libLBFGSB.a: src/LBFGSB/CMakeFiles/LBFGSB.dir/linpack.f.o
src/LBFGSB/libLBFGSB.a: src/LBFGSB/CMakeFiles/LBFGSB.dir/timer.f.o
src/LBFGSB/libLBFGSB.a: src/LBFGSB/CMakeFiles/LBFGSB.dir/build.make
src/LBFGSB/libLBFGSB.a: src/LBFGSB/CMakeFiles/LBFGSB.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking Fortran static library libLBFGSB.a"
	cd /home/ksb11/GlobalOptimization/GlobalOpt/build/src/LBFGSB && $(CMAKE_COMMAND) -P CMakeFiles/LBFGSB.dir/cmake_clean_target.cmake
	cd /home/ksb11/GlobalOptimization/GlobalOpt/build/src/LBFGSB && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LBFGSB.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/LBFGSB/CMakeFiles/LBFGSB.dir/build: src/LBFGSB/libLBFGSB.a
.PHONY : src/LBFGSB/CMakeFiles/LBFGSB.dir/build

src/LBFGSB/CMakeFiles/LBFGSB.dir/requires: src/LBFGSB/CMakeFiles/LBFGSB.dir/blas.f.o.requires
src/LBFGSB/CMakeFiles/LBFGSB.dir/requires: src/LBFGSB/CMakeFiles/LBFGSB.dir/lbfgsb.f.o.requires
src/LBFGSB/CMakeFiles/LBFGSB.dir/requires: src/LBFGSB/CMakeFiles/LBFGSB.dir/linpack.f.o.requires
src/LBFGSB/CMakeFiles/LBFGSB.dir/requires: src/LBFGSB/CMakeFiles/LBFGSB.dir/timer.f.o.requires
.PHONY : src/LBFGSB/CMakeFiles/LBFGSB.dir/requires

src/LBFGSB/CMakeFiles/LBFGSB.dir/clean:
	cd /home/ksb11/GlobalOptimization/GlobalOpt/build/src/LBFGSB && $(CMAKE_COMMAND) -P CMakeFiles/LBFGSB.dir/cmake_clean.cmake
.PHONY : src/LBFGSB/CMakeFiles/LBFGSB.dir/clean

src/LBFGSB/CMakeFiles/LBFGSB.dir/depend:
	cd /home/ksb11/GlobalOptimization/GlobalOpt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ksb11/GlobalOptimization/GlobalOpt /home/ksb11/GlobalOptimization/GlobalOpt/src/LBFGSB /home/ksb11/GlobalOptimization/GlobalOpt/build /home/ksb11/GlobalOptimization/GlobalOpt/build/src/LBFGSB /home/ksb11/GlobalOptimization/GlobalOpt/build/src/LBFGSB/CMakeFiles/LBFGSB.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/LBFGSB/CMakeFiles/LBFGSB.dir/depend
