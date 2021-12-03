#include <vector>
#include <complex>

#include "globals.h"

#define MKL_Complex16 std::complex<double>
#include "mkl.h"

#include "include/twoPhonons/evalExpectation.h"

double evalExpectationOnePh(const std::vector<std::complex<double>> &op, const std::vector<std::complex<double>> &vec){

  std::complex<double> alpha (1., 0.);
  std::complex<double> beta (0., 0.);

  std::vector<std::complex<double>> matVec(vec.size(), std::complex<double> (0., 0.));

  cblas_zgemv(CblasRowMajor, CblasNoTrans,
              dimHOnePh, dimHOnePh, &alpha, op.data(),
              dimHOnePh, vec.data(), 1,
              &beta, matVec.data(), 1);

  std::complex<double> expectation (0., 0.);
  for(ulong ind = 0ul; ind < dimHOnePh; ++ind) {
    expectation += std::conj(vec[ind]) * matVec[ind];
  }
  return expectation.real();
}

