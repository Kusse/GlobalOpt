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
CMAKE_SOURCE_DIR = /home/ksb11/GAMIN/GAWALE/Benchmarks/ACSTEM/psedo_images

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ksb11/GAMIN/GAWALE/Benchmarks/ACSTEM/psedo_images/build

# Include any dependencies generated for this target.
include CMakeFiles/Test.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Test.dir/flags.make

CMakeFiles/Test.dir/main.f90.o: CMakeFiles/Test.dir/flags.make
CMakeFiles/Test.dir/main.f90.o: ../main.f90
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ksb11/GAMIN/GAWALE/Benchmarks/ACSTEM/psedo_images/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/Test.dir/main.f90.o"
	/usr/bin/f95  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/ksb11/GAMIN/GAWALE/Benchmarks/ACSTEM/psedo_images/main.f90 -o CMakeFiles/Test.dir/main.f90.o

CMakeFiles/Test.dir/main.f90.o.requires:
.PHONY : CMakeFiles/Test.dir/main.f90.o.requires

CMakeFiles/Test.dir/main.f90.o.provides: CMakeFiles/Test.dir/main.f90.o.requires
	$(MAKE) -f CMakeFiles/Test.dir/build.make CMakeFiles/Test.dir/main.f90.o.provides.build
.PHONY : CMakeFiles/Test.dir/main.f90.o.provides

CMakeFiles/Test.dir/main.f90.o.provides.build: CMakeFiles/Test.dir/main.f90.o

# Object files for target Test
Test_OBJECTS = \
"CMakeFiles/Test.dir/main.f90.o"

# External object files for target Test
Test_EXTERNAL_OBJECTS =

Test: CMakeFiles/Test.dir/main.f90.o
Test: CMakeFiles/Test.dir/build.make
Test: CMakeFiles/Test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking Fortran executable Test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Test.dir/build: Test
.PHONY : CMakeFiles/Test.dir/build

CMakeFiles/Test.dir/requires: CMakeFiles/Test.dir/main.f90.o.requires
.PHONY : CMakeFiles/Test.dir/requires

CMakeFiles/Test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Test.dir/clean

CMakeFiles/Test.dir/depend:
	cd /home/ksb11/GAMIN/GAWALE/Benchmarks/ACSTEM/psedo_images/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ksb11/GAMIN/GAWALE/Benchmarks/ACSTEM/psedo_images /home/ksb11/GAMIN/GAWALE/Benchmarks/ACSTEM/psedo_images /home/ksb11/GAMIN/GAWALE/Benchmarks/ACSTEM/psedo_images/build /home/ksb11/GAMIN/GAWALE/Benchmarks/ACSTEM/psedo_images/build /home/ksb11/GAMIN/GAWALE/Benchmarks/ACSTEM/psedo_images/build/CMakeFiles/Test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Test.dir/depend

