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


  evalGSProps(1);
  evalGSPropsTemp();
  
  //evalGSPropsAsOfWPh();
  //eigenenergiesAsOfG();
  
  //calcTimeEvolution2Bands();
  //calcTimeEvolutionAsOfWD();
  
  
  return 0;
}

