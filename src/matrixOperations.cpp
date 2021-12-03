#include <vector>
#include <complex>
#include <cassert>

#include "matrixOperations.h"
#include "globals.h"

#define MKL_Complex16 std::complex<double>
#include "mkl.h"


void dagger(std::vector<std::complex<double>> &Mat, const ulong dim) {

  mkl_zimatcopy(LAPACK_ROW_MAJOR, 'C', dim, dim, std::complex<double>(1., 0.), Mat.data(), dim, dim);

}

std::vector<double> diagonalize(std::vector<std::complex<double>> &Mat, ulong dim, const char eVecs = 'V') {

  MKL_INT info;
  std::vector<double> eVals(dim, 0.);

  info = LAPACKE_zheev(LAPACK_ROW_MAJOR, eVecs, 'U', int(dim), Mat.data(), int(dim), eVals.data());

  return eVals;

}

void addMatricies(const std::vector<std::complex<double>> &A,
                  const std::vector<std::complex<double>> &B,
                  std::vector<std::complex<double>> &C){
  assert(A.size() == B.size());
  assert(A.size() == C.size());

  for(ulong ind = 0ul; ind < A.size(); ++ind){
    C[ind] = A[ind] + B[ind];
  }
}

void matrixAMinusB(const std::vector<std::complex<double>> &A,
                  const std::vector<std::complex<double>> &B,
                  std::vector<std::complex<double>> &C){
  assert(A.size() == B.size());
  assert(A.size() == C.size());

  for(ulong ind = 0ul; ind < A.size(); ++ind){
    C[ind] = A[ind] - B[ind];
  }
}
