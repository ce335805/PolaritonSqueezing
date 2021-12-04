#include <vector>
#include <complex>
#include <iostream>

#include "matrixOperations.h"

#define MKL_Complex16 std::complex<double>
#include "mkl.h"


void checkMatrixDiag(){

  std::vector<std::complex<double>> A(
          {std::complex<double> (1., 0.), std::complex<double> (3., 0.),
           std::complex<double> (3., 0.), std::complex<double> (4., 0.)}
  );
  std::vector<std::complex<double>> ACopy(A);
  std::vector<std::complex<double>> ADiag(A.size(), std::complex<double> (0., 0.));

  std::vector<double> eValsA(2ul, 0.);
  eValsA = diagonalize(A, 2ul, 'V');

  for (ulong ind = 0ul; ind < eValsA.size(); ++ind) {
    std::cout << eValsA[ind] << '\n';
  }

  std::complex<double> alpha = std::complex<double> (1., 0.);
  std::complex<double> beta = std::complex<double> (0., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              2ul, 2ul, 2ul, &alpha,
              ACopy.data(), 2ul,
              A.data(), 2ul,
              &beta, ADiag.data(), 2ul);


  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              2ul, 2ul, 2ul, &alpha,
              A.data(), 2ul,
              ADiag.data(), 2ul,
              &beta, ADiag.data(), 2ul);

  for (ulong ind = 0ul; ind < eValsA.size(); ++ind) {
    std::cout << ADiag[ind * 2ul + ind].real() << '\n';
  }


}


