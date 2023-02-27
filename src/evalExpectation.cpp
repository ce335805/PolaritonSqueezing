#include <vector>
#include <complex>

#include "globals.h"

#define MKL_Complex16 std::complex<double>
#include "mkl.h"

#include "include/evalExpectation.h"

double evalExpectation(const std::vector<std::complex<double>> &op, const std::vector<std::complex<double>> &vec, const ulong dimH){

  std::complex<double> alpha (1., 0.);
  std::complex<double> beta (0., 0.);

  std::vector<std::complex<double>> matVec(vec.size(), std::complex<double> (0., 0.));

  cblas_zgemv(CblasRowMajor, CblasNoTrans,
              dimH, dimH, &alpha, op.data(),
              dimH, vec.data(), 1,
              &beta, matVec.data(), 1);

  std::complex<double> expectation (0., 0.);
  for(ulong ind = 0ul; ind < dimH; ++ind) {
    expectation += std::conj(vec[ind]) * matVec[ind];
  }
  return expectation.real();
}


void evalAllExpectations(std::vector<double> &expecValues,
                         const std::vector<std::complex<double>> &op,
                         const std::vector<std::complex<double>> &eVecs,
                         const ulong dimH){
  
  std::complex<double> alpha (1., 0.);
  std::complex<double> beta (0., 0.);
  
  std::vector<std::complex<double>> matVec(eVecs.size(), std::complex<double> (0., 0.));
  std::vector<std::complex<double>> matVec2(eVecs.size(), std::complex<double> (0., 0.));
  
  int dimHInt = int(dimH);
  
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHInt, dimHInt, dimHInt, &alpha,
              op.data(), dimHInt,
              eVecs.data(), dimHInt,
              &beta, matVec.data(), dimHInt);
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHInt, dimHInt, dimHInt, &alpha,
              eVecs.data(), dimHInt,
              matVec.data(), dimHInt,
              &beta, matVec2.data(), dimHInt);
  
  for(ulong ind = 0ul; ind < dimH; ++ind) {
    expecValues[ind] = std::real(matVec2[ind * dimH + ind]);
  }
}

