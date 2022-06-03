#include <vector>
#include <complex>
#include <iostream>

#include "globals.h"
#include "setUpGlobalHamiltonianOnlyPhot.h"
#include "setupBasicOperatorsOnlyPhot.h"
#include "matrixOperations.h"

#define MKL_Complex16 std::complex<double>
#include "mkl.h"

void setupGlobalHOnlyPhot(std::vector<std::complex<double>> &H){

  H = std::vector<std::complex<double>> (dimHOnlyPhot * dimHOnlyPhot, std::complex<double> (0., 0.));


  std::vector<std::complex<double>> HUncoupled;
  std::vector<std::complex<double>> HEPt;

  setupUncoupledHamiltonianOnlyPhot(HUncoupled);
  setupEPtCoupling(HEPt);

  addMatricies(HUncoupled, HEPt, H);
}

void setupEPtCoupling(std::vector<std::complex<double>> &HEPt){
  
  HEPt = std::vector<std::complex<double>> (dimHOnlyPhot * dimHOnlyPhot, std::complex<double> (0., 0.));

  std::vector<std::complex<double>> A;
  setupAOnlyPhot(A);
  
  std::vector<std::complex<double>> HCoupling;
  setupHCouplingOnlyPhot(HCoupling);

  std::complex<double> alpha (0., gE);
  std::complex<double> beta (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnlyPhot, dimHOnlyPhot, dimHOnlyPhot, &alpha,
              A.data(), dimHOnlyPhot,
              HCoupling.data(), dimHOnlyPhot,
              &beta, HEPt.data(), dimHOnlyPhot);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnlyPhot, dimHOnlyPhot, dimHOnlyPhot, &alpha,
              A.data(), dimHOnlyPhot,
              HCoupling.data(), dimHOnlyPhot,
              &beta, HEPt.data(), dimHOnlyPhot);
  

}
void setupUncoupledHamiltonianOnlyPhot(std::vector<std::complex<double>> &H) {
  
  setupHelectronicOnlyPhot(H);

  std::vector<std::complex<double>> A;
  setupAOnlyPhot(A);

  std::complex<double> alpha (wPt, 0.);
  std::complex<double> beta (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnlyPhot, dimHOnlyPhot, dimHOnlyPhot, &alpha,
              A.data(), dimHOnlyPhot,
              A.data(), dimHOnlyPhot,
              &beta, H.data(), dimHOnlyPhot);
}
