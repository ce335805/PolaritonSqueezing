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
#include "H5Cpp.h"
#include "setUpGlobalHamiltonian.h"

#define MKL_Complex16 std::complex<double>
#include "mkl.h"

int main() {

  std::cout << "Hello World! - let's squeeze some phonons" << std::endl;

  std::cout << std::setprecision(12);

  std::cout << "Bare phonon frequency is: wPh = " << wPh << '\n';

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

  std::cout << "Starting Time Evolution ..." << '\n';
  auto start = std::chrono::high_resolution_clock::now();

  const bool twoPhonons = false;
  calcTimeEvolution(twoPhonons);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  std::cout << "Time Evolution took " << duration.count() << "ms" << '\n';

  //const bool twoPhonons = true;
  //evalGSProps(twoPhonons);


  return 0;
}

