# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb2

# Include any dependencies generated for this target.
include exe/CMakeFiles/h2_fci.dir/depend.make

# Include the progress variables for this target.
include exe/CMakeFiles/h2_fci.dir/progress.make

# Include the compile flags for this target's objects.
include exe/CMakeFiles/h2_fci.dir/flags.make

exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.o: exe/CMakeFiles/h2_fci.dir/flags.make
exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.o: ../exe/h2_fci.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.o"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb2/exe && /opt/local/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/h2_fci.dir/h2_fci.cpp.o -c /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/exe/h2_fci.cpp

exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/h2_fci.dir/h2_fci.cpp.i"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb2/exe && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/exe/h2_fci.cpp > CMakeFiles/h2_fci.dir/h2_fci.cpp.i

exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/h2_fci.dir/h2_fci.cpp.s"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb2/exe && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/exe/h2_fci.cpp -o CMakeFiles/h2_fci.dir/h2_fci.cpp.s

exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.o.requires:

.PHONY : exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.o.requires

exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.o.provides: exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.o.requires
	$(MAKE) -f exe/CMakeFiles/h2_fci.dir/build.make exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.o.provides.build
.PHONY : exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.o.provides

exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.o.provides.build: exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.o


# Object files for target h2_fci
h2_fci_OBJECTS = \
"CMakeFiles/h2_fci.dir/h2_fci.cpp.o"

# External object files for target h2_fci
h2_fci_EXTERNAL_OBJECTS =

exe/h2_fci: exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.o
exe/h2_fci: exe/CMakeFiles/h2_fci.dir/build.make
exe/h2_fci: src/libci.a
exe/h2_fci: /usr/local/hf/lib/libhf.a
exe/h2_fci: /usr/local/bmqc/lib/libbmqc.a
exe/h2_fci: /usr/local/cpputil/lib/libcpputil.a
exe/h2_fci: /usr/local/libwint/lib/libwint.a
exe/h2_fci: /usr/local/libint/2.3.1/lib/libint2.a
exe/h2_fci: /usr/local/numopt/lib/libnumopt.a
exe/h2_fci: exe/CMakeFiles/h2_fci.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable h2_fci"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb2/exe && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/h2_fci.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
exe/CMakeFiles/h2_fci.dir/build: exe/h2_fci

.PHONY : exe/CMakeFiles/h2_fci.dir/build

exe/CMakeFiles/h2_fci.dir/requires: exe/CMakeFiles/h2_fci.dir/h2_fci.cpp.o.requires

.PHONY : exe/CMakeFiles/h2_fci.dir/requires

exe/CMakeFiles/h2_fci.dir/clean:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb2/exe && $(CMAKE_COMMAND) -P CMakeFiles/h2_fci.dir/cmake_clean.cmake
.PHONY : exe/CMakeFiles/h2_fci.dir/clean

exe/CMakeFiles/h2_fci.dir/depend:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/exe /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb2 /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb2/exe /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb2/exe/CMakeFiles/h2_fci.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : exe/CMakeFiles/h2_fci.dir/depend

