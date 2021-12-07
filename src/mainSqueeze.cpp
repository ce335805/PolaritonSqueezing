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

  //std::cout << std::fixed;
  //std::cout << std::setprecision(5);

  std::cout << "Hello World! - let's squeeze some phonons" << std::endl;

  //evalDoccInGSTwoPh();
  //evalDoccInGSOnePh();

  const bool twoPhonons = false;
  calcTimeEvolution(twoPhonons);

  //std::vector<std::complex<double>> I({1., 0., 0., 1.});
  //std::vector<std::complex<double>> vec({1., 2.});
  //std::vector<std::complex<double>> vecRes({0., 0.});
  //
          //std::complex<double> alpha(1., 0.);
  //std::complex<double> beta(0., 0.);
  //
          //cblas_zgemv(CblasRowMajor, CblasNoTrans,
                       //            2ul, 2ul, &alpha, I.data(),
          //            2ul, vec.data(), 1,
          //            &beta, vecRes.data(), 1);
  //
          //std::cout << "vecRes[0] = " << vecRes[0] << '\n';
  //std::cout << "vecRes[1] = " << vecRes[1] << '\n';

  //check how to diagonalize matrix via matrix multiplication


  //checkMatrixDiag();



  return 0;
}

