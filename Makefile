# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/chris/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/212.5457.51/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/chris/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/212.5457.51/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/chris/academics/PhD/codebase/phononSqueezing

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/chris/academics/PhD/codebase/phononSqueezing

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/home/chris/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/212.5457.51/bin/cmake/linux/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/home/chris/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/212.5457.51/bin/cmake/linux/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target test
test:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running tests..."
	/home/chris/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/212.5457.51/bin/cmake/linux/bin/ctest --force-new-ctest-process $(ARGS)
.PHONY : test

# Special rule for the target test
test/fast: test
.PHONY : test/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/chris/academics/PhD/codebase/phononSqueezing/CMakeFiles /home/chris/academics/PhD/codebase/phononSqueezing//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/chris/academics/PhD/codebase/phononSqueezing/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named PhononSqueezing.out

# Build rule for target.
PhononSqueezing.out: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 PhononSqueezing.out
.PHONY : PhononSqueezing.out

# fast build rule for target.
PhononSqueezing.out/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/build
.PHONY : PhononSqueezing.out/fast

src/checkSomeStuff.o: src/checkSomeStuff.cpp.o
.PHONY : src/checkSomeStuff.o

# target to build an object file
src/checkSomeStuff.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/checkSomeStuff.cpp.o
.PHONY : src/checkSomeStuff.cpp.o

src/checkSomeStuff.i: src/checkSomeStuff.cpp.i
.PHONY : src/checkSomeStuff.i

# target to preprocess a source file
src/checkSomeStuff.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/checkSomeStuff.cpp.i
.PHONY : src/checkSomeStuff.cpp.i

src/checkSomeStuff.s: src/checkSomeStuff.cpp.s
.PHONY : src/checkSomeStuff.s

# target to generate assembly for a file
src/checkSomeStuff.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/checkSomeStuff.cpp.s
.PHONY : src/checkSomeStuff.cpp.s

src/mainSqueeze.o: src/mainSqueeze.cpp.o
.PHONY : src/mainSqueeze.o

# target to build an object file
src/mainSqueeze.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/mainSqueeze.cpp.o
.PHONY : src/mainSqueeze.cpp.o

src/mainSqueeze.i: src/mainSqueeze.cpp.i
.PHONY : src/mainSqueeze.i

# target to preprocess a source file
src/mainSqueeze.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/mainSqueeze.cpp.i
.PHONY : src/mainSqueeze.cpp.i

src/mainSqueeze.s: src/mainSqueeze.cpp.s
.PHONY : src/mainSqueeze.s

# target to generate assembly for a file
src/mainSqueeze.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/mainSqueeze.cpp.s
.PHONY : src/mainSqueeze.cpp.s

src/matrixOperations.o: src/matrixOperations.cpp.o
.PHONY : src/matrixOperations.o

# target to build an object file
src/matrixOperations.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/matrixOperations.cpp.o
.PHONY : src/matrixOperations.cpp.o

src/matrixOperations.i: src/matrixOperations.cpp.i
.PHONY : src/matrixOperations.i

# target to preprocess a source file
src/matrixOperations.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/matrixOperations.cpp.i
.PHONY : src/matrixOperations.cpp.i

src/matrixOperations.s: src/matrixOperations.cpp.s
.PHONY : src/matrixOperations.s

# target to generate assembly for a file
src/matrixOperations.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/matrixOperations.cpp.s
.PHONY : src/matrixOperations.cpp.s

src/onePhonon/calcGSOnePh.o: src/onePhonon/calcGSOnePh.cpp.o
.PHONY : src/onePhonon/calcGSOnePh.o

# target to build an object file
src/onePhonon/calcGSOnePh.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/calcGSOnePh.cpp.o
.PHONY : src/onePhonon/calcGSOnePh.cpp.o

src/onePhonon/calcGSOnePh.i: src/onePhonon/calcGSOnePh.cpp.i
.PHONY : src/onePhonon/calcGSOnePh.i

# target to preprocess a source file
src/onePhonon/calcGSOnePh.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/calcGSOnePh.cpp.i
.PHONY : src/onePhonon/calcGSOnePh.cpp.i

src/onePhonon/calcGSOnePh.s: src/onePhonon/calcGSOnePh.cpp.s
.PHONY : src/onePhonon/calcGSOnePh.s

# target to generate assembly for a file
src/onePhonon/calcGSOnePh.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/calcGSOnePh.cpp.s
.PHONY : src/onePhonon/calcGSOnePh.cpp.s

src/onePhonon/evalDoccInGSOnePh.o: src/onePhonon/evalDoccInGSOnePh.cpp.o
.PHONY : src/onePhonon/evalDoccInGSOnePh.o

# target to build an object file
src/onePhonon/evalDoccInGSOnePh.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/evalDoccInGSOnePh.cpp.o
.PHONY : src/onePhonon/evalDoccInGSOnePh.cpp.o

src/onePhonon/evalDoccInGSOnePh.i: src/onePhonon/evalDoccInGSOnePh.cpp.i
.PHONY : src/onePhonon/evalDoccInGSOnePh.i

# target to preprocess a source file
src/onePhonon/evalDoccInGSOnePh.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/evalDoccInGSOnePh.cpp.i
.PHONY : src/onePhonon/evalDoccInGSOnePh.cpp.i

src/onePhonon/evalDoccInGSOnePh.s: src/onePhonon/evalDoccInGSOnePh.cpp.s
.PHONY : src/onePhonon/evalDoccInGSOnePh.s

# target to generate assembly for a file
src/onePhonon/evalDoccInGSOnePh.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/evalDoccInGSOnePh.cpp.s
.PHONY : src/onePhonon/evalDoccInGSOnePh.cpp.s

src/onePhonon/evalExpectationOnePh.o: src/onePhonon/evalExpectationOnePh.cpp.o
.PHONY : src/onePhonon/evalExpectationOnePh.o

# target to build an object file
src/onePhonon/evalExpectationOnePh.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/evalExpectationOnePh.cpp.o
.PHONY : src/onePhonon/evalExpectationOnePh.cpp.o

src/onePhonon/evalExpectationOnePh.i: src/onePhonon/evalExpectationOnePh.cpp.i
.PHONY : src/onePhonon/evalExpectationOnePh.i

# target to preprocess a source file
src/onePhonon/evalExpectationOnePh.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/evalExpectationOnePh.cpp.i
.PHONY : src/onePhonon/evalExpectationOnePh.cpp.i

src/onePhonon/evalExpectationOnePh.s: src/onePhonon/evalExpectationOnePh.cpp.s
.PHONY : src/onePhonon/evalExpectationOnePh.s

# target to generate assembly for a file
src/onePhonon/evalExpectationOnePh.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/evalExpectationOnePh.cpp.s
.PHONY : src/onePhonon/evalExpectationOnePh.cpp.s

src/onePhonon/setUpGlobalHamiltonianOnePh.o: src/onePhonon/setUpGlobalHamiltonianOnePh.cpp.o
.PHONY : src/onePhonon/setUpGlobalHamiltonianOnePh.o

# target to build an object file
src/onePhonon/setUpGlobalHamiltonianOnePh.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/setUpGlobalHamiltonianOnePh.cpp.o
.PHONY : src/onePhonon/setUpGlobalHamiltonianOnePh.cpp.o

src/onePhonon/setUpGlobalHamiltonianOnePh.i: src/onePhonon/setUpGlobalHamiltonianOnePh.cpp.i
.PHONY : src/onePhonon/setUpGlobalHamiltonianOnePh.i

# target to preprocess a source file
src/onePhonon/setUpGlobalHamiltonianOnePh.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/setUpGlobalHamiltonianOnePh.cpp.i
.PHONY : src/onePhonon/setUpGlobalHamiltonianOnePh.cpp.i

src/onePhonon/setUpGlobalHamiltonianOnePh.s: src/onePhonon/setUpGlobalHamiltonianOnePh.cpp.s
.PHONY : src/onePhonon/setUpGlobalHamiltonianOnePh.s

# target to generate assembly for a file
src/onePhonon/setUpGlobalHamiltonianOnePh.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/setUpGlobalHamiltonianOnePh.cpp.s
.PHONY : src/onePhonon/setUpGlobalHamiltonianOnePh.cpp.s

src/onePhonon/setupBasicOpertorsOnePh.o: src/onePhonon/setupBasicOpertorsOnePh.cpp.o
.PHONY : src/onePhonon/setupBasicOpertorsOnePh.o

# target to build an object file
src/onePhonon/setupBasicOpertorsOnePh.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/setupBasicOpertorsOnePh.cpp.o
.PHONY : src/onePhonon/setupBasicOpertorsOnePh.cpp.o

src/onePhonon/setupBasicOpertorsOnePh.i: src/onePhonon/setupBasicOpertorsOnePh.cpp.i
.PHONY : src/onePhonon/setupBasicOpertorsOnePh.i

# target to preprocess a source file
src/onePhonon/setupBasicOpertorsOnePh.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/setupBasicOpertorsOnePh.cpp.i
.PHONY : src/onePhonon/setupBasicOpertorsOnePh.cpp.i

src/onePhonon/setupBasicOpertorsOnePh.s: src/onePhonon/setupBasicOpertorsOnePh.cpp.s
.PHONY : src/onePhonon/setupBasicOpertorsOnePh.s

# target to generate assembly for a file
src/onePhonon/setupBasicOpertorsOnePh.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/setupBasicOpertorsOnePh.cpp.s
.PHONY : src/onePhonon/setupBasicOpertorsOnePh.cpp.s

src/onePhonon/timeEvolution.o: src/onePhonon/timeEvolution.cpp.o
.PHONY : src/onePhonon/timeEvolution.o

# target to build an object file
src/onePhonon/timeEvolution.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/timeEvolution.cpp.o
.PHONY : src/onePhonon/timeEvolution.cpp.o

src/onePhonon/timeEvolution.i: src/onePhonon/timeEvolution.cpp.i
.PHONY : src/onePhonon/timeEvolution.i

# target to preprocess a source file
src/onePhonon/timeEvolution.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/timeEvolution.cpp.i
.PHONY : src/onePhonon/timeEvolution.cpp.i

src/onePhonon/timeEvolution.s: src/onePhonon/timeEvolution.cpp.s
.PHONY : src/onePhonon/timeEvolution.s

# target to generate assembly for a file
src/onePhonon/timeEvolution.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/timeEvolution.cpp.s
.PHONY : src/onePhonon/timeEvolution.cpp.s

src/onePhonon/timeStep.o: src/onePhonon/timeStep.cpp.o
.PHONY : src/onePhonon/timeStep.o

# target to build an object file
src/onePhonon/timeStep.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/timeStep.cpp.o
.PHONY : src/onePhonon/timeStep.cpp.o

src/onePhonon/timeStep.i: src/onePhonon/timeStep.cpp.i
.PHONY : src/onePhonon/timeStep.i

# target to preprocess a source file
src/onePhonon/timeStep.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/timeStep.cpp.i
.PHONY : src/onePhonon/timeStep.cpp.i

src/onePhonon/timeStep.s: src/onePhonon/timeStep.cpp.s
.PHONY : src/onePhonon/timeStep.s

# target to generate assembly for a file
src/onePhonon/timeStep.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/onePhonon/timeStep.cpp.s
.PHONY : src/onePhonon/timeStep.cpp.s

src/setupElectronicOperatorsSmall.o: src/setupElectronicOperatorsSmall.cpp.o
.PHONY : src/setupElectronicOperatorsSmall.o

# target to build an object file
src/setupElectronicOperatorsSmall.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/setupElectronicOperatorsSmall.cpp.o
.PHONY : src/setupElectronicOperatorsSmall.cpp.o

src/setupElectronicOperatorsSmall.i: src/setupElectronicOperatorsSmall.cpp.i
.PHONY : src/setupElectronicOperatorsSmall.i

# target to preprocess a source file
src/setupElectronicOperatorsSmall.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/setupElectronicOperatorsSmall.cpp.i
.PHONY : src/setupElectronicOperatorsSmall.cpp.i

src/setupElectronicOperatorsSmall.s: src/setupElectronicOperatorsSmall.cpp.s
.PHONY : src/setupElectronicOperatorsSmall.s

# target to generate assembly for a file
src/setupElectronicOperatorsSmall.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/setupElectronicOperatorsSmall.cpp.s
.PHONY : src/setupElectronicOperatorsSmall.cpp.s

src/twoPhonons/evalDoccInGS.o: src/twoPhonons/evalDoccInGS.cpp.o
.PHONY : src/twoPhonons/evalDoccInGS.o

# target to build an object file
src/twoPhonons/evalDoccInGS.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/twoPhonons/evalDoccInGS.cpp.o
.PHONY : src/twoPhonons/evalDoccInGS.cpp.o

src/twoPhonons/evalDoccInGS.i: src/twoPhonons/evalDoccInGS.cpp.i
.PHONY : src/twoPhonons/evalDoccInGS.i

# target to preprocess a source file
src/twoPhonons/evalDoccInGS.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/twoPhonons/evalDoccInGS.cpp.i
.PHONY : src/twoPhonons/evalDoccInGS.cpp.i

src/twoPhonons/evalDoccInGS.s: src/twoPhonons/evalDoccInGS.cpp.s
.PHONY : src/twoPhonons/evalDoccInGS.s

# target to generate assembly for a file
src/twoPhonons/evalDoccInGS.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/twoPhonons/evalDoccInGS.cpp.s
.PHONY : src/twoPhonons/evalDoccInGS.cpp.s

src/twoPhonons/evalExpectation.o: src/twoPhonons/evalExpectation.cpp.o
.PHONY : src/twoPhonons/evalExpectation.o

# target to build an object file
src/twoPhonons/evalExpectation.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/twoPhonons/evalExpectation.cpp.o
.PHONY : src/twoPhonons/evalExpectation.cpp.o

src/twoPhonons/evalExpectation.i: src/twoPhonons/evalExpectation.cpp.i
.PHONY : src/twoPhonons/evalExpectation.i

# target to preprocess a source file
src/twoPhonons/evalExpectation.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/twoPhonons/evalExpectation.cpp.i
.PHONY : src/twoPhonons/evalExpectation.cpp.i

src/twoPhonons/evalExpectation.s: src/twoPhonons/evalExpectation.cpp.s
.PHONY : src/twoPhonons/evalExpectation.s

# target to generate assembly for a file
src/twoPhonons/evalExpectation.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/twoPhonons/evalExpectation.cpp.s
.PHONY : src/twoPhonons/evalExpectation.cpp.s

src/twoPhonons/setUpGlobalHamiltonian.o: src/twoPhonons/setUpGlobalHamiltonian.cpp.o
.PHONY : src/twoPhonons/setUpGlobalHamiltonian.o

# target to build an object file
src/twoPhonons/setUpGlobalHamiltonian.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/twoPhonons/setUpGlobalHamiltonian.cpp.o
.PHONY : src/twoPhonons/setUpGlobalHamiltonian.cpp.o

src/twoPhonons/setUpGlobalHamiltonian.i: src/twoPhonons/setUpGlobalHamiltonian.cpp.i
.PHONY : src/twoPhonons/setUpGlobalHamiltonian.i

# target to preprocess a source file
src/twoPhonons/setUpGlobalHamiltonian.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/twoPhonons/setUpGlobalHamiltonian.cpp.i
.PHONY : src/twoPhonons/setUpGlobalHamiltonian.cpp.i

src/twoPhonons/setUpGlobalHamiltonian.s: src/twoPhonons/setUpGlobalHamiltonian.cpp.s
.PHONY : src/twoPhonons/setUpGlobalHamiltonian.s

# target to generate assembly for a file
src/twoPhonons/setUpGlobalHamiltonian.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/twoPhonons/setUpGlobalHamiltonian.cpp.s
.PHONY : src/twoPhonons/setUpGlobalHamiltonian.cpp.s

src/twoPhonons/setupBasicOpertors.o: src/twoPhonons/setupBasicOpertors.cpp.o
.PHONY : src/twoPhonons/setupBasicOpertors.o

# target to build an object file
src/twoPhonons/setupBasicOpertors.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/twoPhonons/setupBasicOpertors.cpp.o
.PHONY : src/twoPhonons/setupBasicOpertors.cpp.o

src/twoPhonons/setupBasicOpertors.i: src/twoPhonons/setupBasicOpertors.cpp.i
.PHONY : src/twoPhonons/setupBasicOpertors.i

# target to preprocess a source file
src/twoPhonons/setupBasicOpertors.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/twoPhonons/setupBasicOpertors.cpp.i
.PHONY : src/twoPhonons/setupBasicOpertors.cpp.i

src/twoPhonons/setupBasicOpertors.s: src/twoPhonons/setupBasicOpertors.cpp.s
.PHONY : src/twoPhonons/setupBasicOpertors.s

# target to generate assembly for a file
src/twoPhonons/setupBasicOpertors.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/PhononSqueezing.out.dir/build.make CMakeFiles/PhononSqueezing.out.dir/src/twoPhonons/setupBasicOpertors.cpp.s
.PHONY : src/twoPhonons/setupBasicOpertors.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... test"
	@echo "... PhononSqueezing.out"
	@echo "... src/checkSomeStuff.o"
	@echo "... src/checkSomeStuff.i"
	@echo "... src/checkSomeStuff.s"
	@echo "... src/mainSqueeze.o"
	@echo "... src/mainSqueeze.i"
	@echo "... src/mainSqueeze.s"
	@echo "... src/matrixOperations.o"
	@echo "... src/matrixOperations.i"
	@echo "... src/matrixOperations.s"
	@echo "... src/onePhonon/calcGSOnePh.o"
	@echo "... src/onePhonon/calcGSOnePh.i"
	@echo "... src/onePhonon/calcGSOnePh.s"
	@echo "... src/onePhonon/evalDoccInGSOnePh.o"
	@echo "... src/onePhonon/evalDoccInGSOnePh.i"
	@echo "... src/onePhonon/evalDoccInGSOnePh.s"
	@echo "... src/onePhonon/evalExpectationOnePh.o"
	@echo "... src/onePhonon/evalExpectationOnePh.i"
	@echo "... src/onePhonon/evalExpectationOnePh.s"
	@echo "... src/onePhonon/setUpGlobalHamiltonianOnePh.o"
	@echo "... src/onePhonon/setUpGlobalHamiltonianOnePh.i"
	@echo "... src/onePhonon/setUpGlobalHamiltonianOnePh.s"
	@echo "... src/onePhonon/setupBasicOpertorsOnePh.o"
	@echo "... src/onePhonon/setupBasicOpertorsOnePh.i"
	@echo "... src/onePhonon/setupBasicOpertorsOnePh.s"
	@echo "... src/onePhonon/timeEvolution.o"
	@echo "... src/onePhonon/timeEvolution.i"
	@echo "... src/onePhonon/timeEvolution.s"
	@echo "... src/onePhonon/timeStep.o"
	@echo "... src/onePhonon/timeStep.i"
	@echo "... src/onePhonon/timeStep.s"
	@echo "... src/setupElectronicOperatorsSmall.o"
	@echo "... src/setupElectronicOperatorsSmall.i"
	@echo "... src/setupElectronicOperatorsSmall.s"
	@echo "... src/twoPhonons/evalDoccInGS.o"
	@echo "... src/twoPhonons/evalDoccInGS.i"
	@echo "... src/twoPhonons/evalDoccInGS.s"
	@echo "... src/twoPhonons/evalExpectation.o"
	@echo "... src/twoPhonons/evalExpectation.i"
	@echo "... src/twoPhonons/evalExpectation.s"
	@echo "... src/twoPhonons/setUpGlobalHamiltonian.o"
	@echo "... src/twoPhonons/setUpGlobalHamiltonian.i"
	@echo "... src/twoPhonons/setUpGlobalHamiltonian.s"
	@echo "... src/twoPhonons/setupBasicOpertors.o"
	@echo "... src/twoPhonons/setupBasicOpertors.i"
	@echo "... src/twoPhonons/setupBasicOpertors.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

