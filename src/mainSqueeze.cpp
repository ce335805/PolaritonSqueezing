#include <iostream>
#include <iomanip>
#include "vector"
#include "complex"
#include "include/twoPhonons/evalDoccInGS.h"
#include "include/onePhonon/evalDoccInGSOnePh.h"
#include "matrixOperations.h"
#include "timeEvolution.h"
#include "globals.h"
#include "checkSomeStuff.h"
#include <chrono>
#include <algorithm>
#include "evalGSProps.h"
#include "setUpGlobalHamiltonian.h"
#include "evalGSPropsOnlyPhot.h"
#include "calcGS.h"
#include "evalExpectation.h"
#include "evalDoccInGs2Bands.h"
#include "timeEvolution2Bands.h"

#define MKL_Complex16 std::complex<double>

#include "mkl.h"
#include "setUpGlobalHamiltonianOnlyPhot.h"
#include "setupBasicOperatorsOnlyPhot.h"
#include "evalGSProps2Bands.h"


int main() {
  
  std::cout << "Hello World! - let's squeeze some phonons" << std::endl;
  
  std::cout << std::setprecision(12);
  //std::cout << "Bare phonon frequency is: wPh = " << wPh << '\n';
  
  //std::vector<std::complex<double>> HTest;
  //setupGlobalH(HTest);
//
  //std::vector<double> realData(dimHTwoPh * dimHTwoPh, 0.0);
//
  //std::transform(HTest.begin(),
  //               HTest.end(),
  //               realData.begin(),
  //               [](const std::complex<double> entry) -> double {
  //                   return entry.real();
  //               });
//
  //H5::H5File file("data/HTest2.hdf5", H5F_ACC_TRUNC);
//
  //const hsize_t dataShape[2] = {dimHTwoPh, dimHTwoPh};
  //H5::DataSpace dataSpace(2, dataShape);
  //H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  //datatype.setOrder(H5T_ORDER_LE);
//
  //H5::DataSet datasetHTest = file.createDataSet("HTest", datatype, dataSpace);
  //datasetHTest.write(&realData[0], datatype);
  
  
  //evalDoccInGSTwoPh();
  //evalDoccInGSOnePh();
  //evalDoccInGSOnlyPhot();
  //evalDoccInGS2Bands();
  
  //evalGSPropsAsOfWPh();
  
  //calcTimeEvolution2Bands();
  calcTimeEvolutionAsOfWD();
  
  
  //std::cout << "Starting Time Evolution ..." << '\n';
  //auto start = std::chrono::high_resolution_clock::now();
//
  //const bool twoPhonons = true;
  //calcTimeEvolution(twoPhonons);
//
  //auto stop = std::chrono::high_resolution_clock::now();
  //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  //std::cout << "Time Evolution took " << duration.count() << "ms" << '\n';
  
  //const bool twoPhonons = true;
  //evalGSProps(twoPhonons);
  //evalGSPropsTemp();
  
  
  return 0;
}

