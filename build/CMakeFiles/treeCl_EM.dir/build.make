# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.2.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.2.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/kgori/Documents/repositories/treeCl_EM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/kgori/Documents/repositories/treeCl_EM/build

# Include any dependencies generated for this target.
include CMakeFiles/treeCl_EM.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/treeCl_EM.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/treeCl_EM.dir/flags.make

CMakeFiles/treeCl_EM.dir/PLL.cpp.o: CMakeFiles/treeCl_EM.dir/flags.make
CMakeFiles/treeCl_EM.dir/PLL.cpp.o: ../PLL.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/kgori/Documents/repositories/treeCl_EM/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/treeCl_EM.dir/PLL.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/treeCl_EM.dir/PLL.cpp.o -c /Users/kgori/Documents/repositories/treeCl_EM/PLL.cpp

CMakeFiles/treeCl_EM.dir/PLL.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/treeCl_EM.dir/PLL.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/kgori/Documents/repositories/treeCl_EM/PLL.cpp > CMakeFiles/treeCl_EM.dir/PLL.cpp.i

CMakeFiles/treeCl_EM.dir/PLL.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/treeCl_EM.dir/PLL.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/kgori/Documents/repositories/treeCl_EM/PLL.cpp -o CMakeFiles/treeCl_EM.dir/PLL.cpp.s

CMakeFiles/treeCl_EM.dir/PLL.cpp.o.requires:
.PHONY : CMakeFiles/treeCl_EM.dir/PLL.cpp.o.requires

CMakeFiles/treeCl_EM.dir/PLL.cpp.o.provides: CMakeFiles/treeCl_EM.dir/PLL.cpp.o.requires
	$(MAKE) -f CMakeFiles/treeCl_EM.dir/build.make CMakeFiles/treeCl_EM.dir/PLL.cpp.o.provides.build
.PHONY : CMakeFiles/treeCl_EM.dir/PLL.cpp.o.provides

CMakeFiles/treeCl_EM.dir/PLL.cpp.o.provides.build: CMakeFiles/treeCl_EM.dir/PLL.cpp.o

CMakeFiles/treeCl_EM.dir/main.cpp.o: CMakeFiles/treeCl_EM.dir/flags.make
CMakeFiles/treeCl_EM.dir/main.cpp.o: ../main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/kgori/Documents/repositories/treeCl_EM/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/treeCl_EM.dir/main.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/treeCl_EM.dir/main.cpp.o -c /Users/kgori/Documents/repositories/treeCl_EM/main.cpp

CMakeFiles/treeCl_EM.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/treeCl_EM.dir/main.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/kgori/Documents/repositories/treeCl_EM/main.cpp > CMakeFiles/treeCl_EM.dir/main.cpp.i

CMakeFiles/treeCl_EM.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/treeCl_EM.dir/main.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/kgori/Documents/repositories/treeCl_EM/main.cpp -o CMakeFiles/treeCl_EM.dir/main.cpp.s

CMakeFiles/treeCl_EM.dir/main.cpp.o.requires:
.PHONY : CMakeFiles/treeCl_EM.dir/main.cpp.o.requires

CMakeFiles/treeCl_EM.dir/main.cpp.o.provides: CMakeFiles/treeCl_EM.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/treeCl_EM.dir/build.make CMakeFiles/treeCl_EM.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/treeCl_EM.dir/main.cpp.o.provides

CMakeFiles/treeCl_EM.dir/main.cpp.o.provides.build: CMakeFiles/treeCl_EM.dir/main.cpp.o

# Object files for target treeCl_EM
treeCl_EM_OBJECTS = \
"CMakeFiles/treeCl_EM.dir/PLL.cpp.o" \
"CMakeFiles/treeCl_EM.dir/main.cpp.o"

# External object files for target treeCl_EM
treeCl_EM_EXTERNAL_OBJECTS =

treeCl_EM: CMakeFiles/treeCl_EM.dir/PLL.cpp.o
treeCl_EM: CMakeFiles/treeCl_EM.dir/main.cpp.o
treeCl_EM: CMakeFiles/treeCl_EM.dir/build.make
treeCl_EM: CMakeFiles/treeCl_EM.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable treeCl_EM"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/treeCl_EM.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/treeCl_EM.dir/build: treeCl_EM
.PHONY : CMakeFiles/treeCl_EM.dir/build

CMakeFiles/treeCl_EM.dir/requires: CMakeFiles/treeCl_EM.dir/PLL.cpp.o.requires
CMakeFiles/treeCl_EM.dir/requires: CMakeFiles/treeCl_EM.dir/main.cpp.o.requires
.PHONY : CMakeFiles/treeCl_EM.dir/requires

CMakeFiles/treeCl_EM.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/treeCl_EM.dir/cmake_clean.cmake
.PHONY : CMakeFiles/treeCl_EM.dir/clean

CMakeFiles/treeCl_EM.dir/depend:
	cd /Users/kgori/Documents/repositories/treeCl_EM/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/kgori/Documents/repositories/treeCl_EM /Users/kgori/Documents/repositories/treeCl_EM /Users/kgori/Documents/repositories/treeCl_EM/build /Users/kgori/Documents/repositories/treeCl_EM/build /Users/kgori/Documents/repositories/treeCl_EM/build/CMakeFiles/treeCl_EM.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/treeCl_EM.dir/depend

