#include <vector>
#include <complex>

#include "matrixOperations.h"
#include "setupHGlobal2Bands.h"
#include "setupBasicOperators2Bands.h"

#include "globals.h"

#define MKL_Complex16 std::complex<double>
#include "mkl.h"

void setupGlobalH2Bands(std::vector<std::complex<double>> &H){
  
  H = std::vector<std::complex<double>> (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));
  
  std::vector<std::complex<double>> HUncoupled;
  std::vector<std::complex<double>> HEBos;
  
  setupHelectronic2Bands(H);
  //setupHEBos(HEBos);
  //addMatricies(HUncoupled, HEBos, H);
}

void setupHEBos(std::vector<std::complex<double>> &HEBos){
  
  HEBos = std::vector<std::complex<double>> (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));
  
  std::vector<std::complex<double>> hopInter0(dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));
  setupInterbandHop0(hopInter0);
  std::vector<std::complex<double>> hopInter1(dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));
  setupInterbandHop0(hopInter1);
  
  std::vector<std::complex<double>> A0;
  setupA0(A0);
  
  std::complex<double> alpha (0., gE);
  std::complex<double> beta (1., 0.);
  
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnlyPhot, dimHOnlyPhot, dimHOnlyPhot, &alpha,
              A0.data(), dimHOnlyPhot,
              hopInter0.data(), dimHOnlyPhot,
              &beta, HEBos.data(), dimHOnlyPhot);
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnlyPhot, dimHOnlyPhot, dimHOnlyPhot, &alpha,
              A0.data(), dimHOnlyPhot,
              hopInter0.data(), dimHOnlyPhot,
              &beta, HEBos.data(), dimHOnlyPhot);
  
  std::vector<std::complex<double>> A1;
  setupA1(A1);
  
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnlyPhot, dimHOnlyPhot, dimHOnlyPhot, &alpha,
              A1.data(), dimHOnlyPhot,
              hopInter0.data(), dimHOnlyPhot,
              &beta, HEBos.data(), dimHOnlyPhot);
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnlyPhot, dimHOnlyPhot, dimHOnlyPhot, &alpha,
              A1.data(), dimHOnlyPhot,
              hopInter0.data(), dimHOnlyPhot,
              &beta, HEBos.data(), dimHOnlyPhot);
  
}

