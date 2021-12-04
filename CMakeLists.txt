#specefiy required version - this is the newest one
cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#using intel compiler
set(CMAKE_CXX_COMPILER mpiicpc)
set(CMAKE_C_COMPILER mpiicc)

enable_language(CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

enable_testing()

#set(BLA_VENDOR Intel10_64lp_seq)
#find_package(LAPACK REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(MKL CONFIG REQUIRED)
#find_package(MPI REQUIRED)

set(project_sources
        src/mainSqueeze.cpp
        src/twoPhonons/setupBasicOpertors.cpp
        src/twoPhonons/setUpGlobalHamiltonian.cpp
        src/matrixOperations.cpp
        src/twoPhonons/evalExpectation.cpp
        src/onePhonon/setupBasicOpertorsOnePh.cpp
        src/setupElectronicOperatorsSmall.cpp
        src/onePhonon/setUpGlobalHamiltonianOnePh.cpp
        src/onePhonon/evalExpectationOnePh.cpp
        src/twoPhonons/evalDoccInGS.cpp
        src/onePhonon/evalDoccInGSOnePh.cpp
        src/onePhonon/timeStep.cpp
        src/onePhonon/timeEvolution.cpp
        src/onePhonon/calcGSOnePh.cpp
        src/checkSomeStuff.cpp
        )

set(header_files
        include/globals.h
        include/utils.h
        include/twoPhonons/setupBasicOperators.h
        include/twoPhonons/setUpGlobalHamiltonian.cpp.h
        include/twoPhonons/evalExpectation.h
        include/setupElectronicOperatorsSmall.h
        include/onePhonon/setupBasicOperatorsOnePh.h
        include/onePhonon/setUpGlobalHamiltonianOnePh.h
        include/onePhonon/evalExpectationOnePh.h
        include/twoPhonons/evalDoccInGS.h
        include/onePhonon/evalDoccInGSOnePh.h
        include/onePhonon/timeStep.h
        include/onePhonon/timeEvolution.h
        include/onePhonon/calcGSOnePh.h
        include/checkSomeStuff.h
        )

include_directories("./include")
include_directories("./include/onePhonon")
include_directories("./include/twoPhonon")
#include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(${HDF5_INCLUDE_DIRS})
message(STATUS "---------------------------------------------------------------------")

#add compiler flags
set(C_flags "-O0 -fvar-tracking -ggdb -Wall -Wextra -pedantic -mkl -qopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${C_flags}")

#set project name
project(PhononSqueezing)
#add an executable - using variables
add_executable(${PROJECT_NAME}.out ${project_sources} ${header_files})

target_link_libraries(${PROJECT_NAME}.out ${HDF5_LIBRARIES})
#target_link_libraries(${PROJECT_NAME}.out ${MPI_CXX_LIBRARIES})


