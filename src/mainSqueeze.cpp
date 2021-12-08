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

#define MKL_Complex16 std::complex<double>
#include "mkl.h"

int main() {

  std::cout << "Hello World! - let's squeeze some phonons" << std::endl;

  std::cout << std::setprecision(12);

  //evalDoccInGSTwoPh();
  evalDoccInGSOnePh();

  //const bool twoPhonons = false;
  //calcTimeEvolution(twoPhonons);

  return 0;



}

