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
CMAKE_BINARY_DIR = /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb

# Include any dependencies generated for this target.
include src/CMakeFiles/ci.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/ci.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/ci.dir/flags.make

src/CMakeFiles/ci.dir/BaseCI.cpp.o: src/CMakeFiles/ci.dir/flags.make
src/CMakeFiles/ci.dir/BaseCI.cpp.o: ../src/BaseCI.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/ci.dir/BaseCI.cpp.o"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/src && /opt/local/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ci.dir/BaseCI.cpp.o -c /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/src/BaseCI.cpp

src/CMakeFiles/ci.dir/BaseCI.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ci.dir/BaseCI.cpp.i"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/src/BaseCI.cpp > CMakeFiles/ci.dir/BaseCI.cpp.i

src/CMakeFiles/ci.dir/BaseCI.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ci.dir/BaseCI.cpp.s"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/src/BaseCI.cpp -o CMakeFiles/ci.dir/BaseCI.cpp.s

src/CMakeFiles/ci.dir/BaseCI.cpp.o.requires:

.PHONY : src/CMakeFiles/ci.dir/BaseCI.cpp.o.requires

src/CMakeFiles/ci.dir/BaseCI.cpp.o.provides: src/CMakeFiles/ci.dir/BaseCI.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/ci.dir/build.make src/CMakeFiles/ci.dir/BaseCI.cpp.o.provides.build
.PHONY : src/CMakeFiles/ci.dir/BaseCI.cpp.o.provides

src/CMakeFiles/ci.dir/BaseCI.cpp.o.provides.build: src/CMakeFiles/ci.dir/BaseCI.cpp.o


src/CMakeFiles/ci.dir/DOCI.cpp.o: src/CMakeFiles/ci.dir/flags.make
src/CMakeFiles/ci.dir/DOCI.cpp.o: ../src/DOCI.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/ci.dir/DOCI.cpp.o"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/src && /opt/local/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ci.dir/DOCI.cpp.o -c /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/src/DOCI.cpp

src/CMakeFiles/ci.dir/DOCI.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ci.dir/DOCI.cpp.i"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/src/DOCI.cpp > CMakeFiles/ci.dir/DOCI.cpp.i

src/CMakeFiles/ci.dir/DOCI.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ci.dir/DOCI.cpp.s"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/src/DOCI.cpp -o CMakeFiles/ci.dir/DOCI.cpp.s

src/CMakeFiles/ci.dir/DOCI.cpp.o.requires:

.PHONY : src/CMakeFiles/ci.dir/DOCI.cpp.o.requires

src/CMakeFiles/ci.dir/DOCI.cpp.o.provides: src/CMakeFiles/ci.dir/DOCI.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/ci.dir/build.make src/CMakeFiles/ci.dir/DOCI.cpp.o.provides.build
.PHONY : src/CMakeFiles/ci.dir/DOCI.cpp.o.provides

src/CMakeFiles/ci.dir/DOCI.cpp.o.provides.build: src/CMakeFiles/ci.dir/DOCI.cpp.o


# Object files for target ci
ci_OBJECTS = \
"CMakeFiles/ci.dir/BaseCI.cpp.o" \
"CMakeFiles/ci.dir/DOCI.cpp.o"

# External object files for target ci
ci_EXTERNAL_OBJECTS =

src/libci.a: src/CMakeFiles/ci.dir/BaseCI.cpp.o
src/libci.a: src/CMakeFiles/ci.dir/DOCI.cpp.o
src/libci.a: src/CMakeFiles/ci.dir/build.make
src/libci.a: src/CMakeFiles/ci.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libci.a"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/src && $(CMAKE_COMMAND) -P CMakeFiles/ci.dir/cmake_clean_target.cmake
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ci.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/ci.dir/build: src/libci.a

.PHONY : src/CMakeFiles/ci.dir/build

src/CMakeFiles/ci.dir/requires: src/CMakeFiles/ci.dir/BaseCI.cpp.o.requires
src/CMakeFiles/ci.dir/requires: src/CMakeFiles/ci.dir/DOCI.cpp.o.requires

.PHONY : src/CMakeFiles/ci.dir/requires

src/CMakeFiles/ci.dir/clean:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/src && $(CMAKE_COMMAND) -P CMakeFiles/ci.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/ci.dir/clean

src/CMakeFiles/ci.dir/depend:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/src /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/src /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tb/src/CMakeFiles/ci.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/ci.dir/depend
