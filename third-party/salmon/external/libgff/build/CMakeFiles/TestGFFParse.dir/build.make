# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /home/linuxbrew/.linuxbrew/Cellar/cmake/3.6.2/bin/cmake

# The command to remove a file.
RM = /home/linuxbrew/.linuxbrew/Cellar/cmake/3.6.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/scratch7/avi/test/shoal-paper/third-party/salmon/external/libgff

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/scratch7/avi/test/shoal-paper/third-party/salmon/external/libgff/build

# Include any dependencies generated for this target.
include CMakeFiles/TestGFFParse.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/TestGFFParse.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TestGFFParse.dir/flags.make

CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o: CMakeFiles/TestGFFParse.dir/flags.make
CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o: ../src/TestGFFParse.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/scratch7/avi/test/shoal-paper/third-party/salmon/external/libgff/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o"
	/home/linuxbrew/.linuxbrew/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o -c /mnt/scratch7/avi/test/shoal-paper/third-party/salmon/external/libgff/src/TestGFFParse.cpp

CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.i"
	/home/linuxbrew/.linuxbrew/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/scratch7/avi/test/shoal-paper/third-party/salmon/external/libgff/src/TestGFFParse.cpp > CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.i

CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.s"
	/home/linuxbrew/.linuxbrew/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/scratch7/avi/test/shoal-paper/third-party/salmon/external/libgff/src/TestGFFParse.cpp -o CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.s

CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o.requires:

.PHONY : CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o.requires

CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o.provides: CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o.requires
	$(MAKE) -f CMakeFiles/TestGFFParse.dir/build.make CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o.provides.build
.PHONY : CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o.provides

CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o.provides.build: CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o


# Object files for target TestGFFParse
TestGFFParse_OBJECTS = \
"CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o"

# External object files for target TestGFFParse
TestGFFParse_EXTERNAL_OBJECTS =

TestGFFParse: CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o
TestGFFParse: CMakeFiles/TestGFFParse.dir/build.make
TestGFFParse: libgff.a
TestGFFParse: CMakeFiles/TestGFFParse.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/scratch7/avi/test/shoal-paper/third-party/salmon/external/libgff/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable TestGFFParse"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestGFFParse.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TestGFFParse.dir/build: TestGFFParse

.PHONY : CMakeFiles/TestGFFParse.dir/build

CMakeFiles/TestGFFParse.dir/requires: CMakeFiles/TestGFFParse.dir/src/TestGFFParse.cpp.o.requires

.PHONY : CMakeFiles/TestGFFParse.dir/requires

CMakeFiles/TestGFFParse.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/TestGFFParse.dir/cmake_clean.cmake
.PHONY : CMakeFiles/TestGFFParse.dir/clean

CMakeFiles/TestGFFParse.dir/depend:
	cd /mnt/scratch7/avi/test/shoal-paper/third-party/salmon/external/libgff/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/scratch7/avi/test/shoal-paper/third-party/salmon/external/libgff /mnt/scratch7/avi/test/shoal-paper/third-party/salmon/external/libgff /mnt/scratch7/avi/test/shoal-paper/third-party/salmon/external/libgff/build /mnt/scratch7/avi/test/shoal-paper/third-party/salmon/external/libgff/build /mnt/scratch7/avi/test/shoal-paper/third-party/salmon/external/libgff/build/CMakeFiles/TestGFFParse.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/TestGFFParse.dir/depend

