#specefiy required version - this is the newest one
cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#using intel compiler
set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_C_COMPILER icc)

enable_language(CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

enable_testing()

#set(BLA_VENDOR Intel10_64lp_seq)
#find_package(LAPACK REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX)
#find_package(MKL CONFIG REQUIRED)
#find_package(MPI REQUIRED)

set(project_sources
        src/mainSqueeze.cpp
        src/twoPhonons/setupBasicOpertors.cpp
        src/twoPhonons/setUpGlobalHamiltonian.cpp
        src/matrixOperations.cpp
        src/evalExpectation.cpp
        src/onePhonon/setupBasicOpertorsOnePh.cpp
        src/setupElectronicOperatorsSmall.cpp
        src/onePhonon/setUpGlobalHamiltonianOnePh.cpp
        src/twoPhonons/evalDoccInGS.cpp
        src/onePhonon/evalDoccInGSOnePh.cpp
        src/timeStep.cpp
        src/timeEvolution.cpp
        src/calcGS.cpp
        src/checkSomeStuff.cpp
        src/utils/writeStuffToHdf5.cpp
        src/utils/setupOps.cpp
        src/evalGSProps.cpp
        src/utils/makeFilenames.cpp
        src/onlyPhot/setupBasicOpertorsOnlyPhot.cpp
        src/onlyPhot/setUpGlobalHamiltonianOnlyPhot.cpp
        src/onlyPhot/evalGSPropsOnlyPhot.cpp
        src/setupElectronicOperators2Bands.cpp
        src/onePhonon/setHGlobal2Bands.cpp
        src/onePhonon/setupBasicOpertors2Bands.cpp
        src/onePhonon/evalDoccInGS2Bands.cpp
        src/onePhonon/evalGSProps2Bands.cpp
        src/timeEvolution2Bands.cpp
        )

set(header_files
        include/globals.h
        include/utils.h
        include/twoPhonons/setupBasicOperators.h
        include/twoPhonons/setUpGlobalHamiltonian.h
        include/evalExpectation.h
        include/setupElectronicOperatorsSmall.h
        include/onePhonon/setupBasicOperatorsOnePh.h
        include/onePhonon/setUpGlobalHamiltonianOnePh.h
        include/twoPhonons/evalDoccInGS.h
        include/onePhonon/evalDoccInGSOnePh.h
        include/timeStep.h
        include/timeEvolution.h
        include/calcGS.h
        include/checkSomeStuff.h
        include/utils/writeStuffToHdf5.h
        include/utils/setupOps.h
        include/evalGSProps.h
        include/utils/makeFilenames.h
        include/onlyPhot/setupBasicOperatorsOnlyPhot.h
        include/onlyPhot/setUpGlobalHamiltonianOnlyPhot.h
        include/onlyPhot/evalGSPropsOnlyPhot.h
        include/setupElectronicOperators2Bands.h
        include/onePhonon/setupHGlobal2Bands.h
        include/onePhonon/setupBasicOperators2Bands.h
        include/onePhonon/evalDoccInGs2Bands.h
        include/onePhonon/evalGSProps2Bands.h
        include/timeEvolution2Bands.h
        )

include_directories("./include")
include_directories("./include/utils")
include_directories("./include/onePhonon")
include_directories("./include/twoPhonons")
include_directories("./include/onlyPhot")
#include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(${HDF5_INCLUDE_DIRS})
message(STATUS "---------------------------------------------------------------------")

#add compiler flags
set(C_flags "-O3 -fvar-tracking -ggdb -Wall -Wextra -pedantic -qmkl -qopenmp -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${C_flags}")

#set project name
project(PhononSqueezing)
#add an executable - using variables
set(EXNAME "PolaritonSqueezing")
add_executable(${EXNAME}.out ${project_sources} ${header_files})

target_link_libraries(${EXNAME}.out ${HDF5_LIBRARIES})
#target_link_libraries(${PROJECT_NAME}.out ${MPI_CXX_LIBRARIES})


