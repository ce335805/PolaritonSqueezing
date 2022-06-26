#include <vector>
#include <complex>

#include "matrixOperations.h"
#include "setupHGlobal2Bands.h"
#include "setupBasicOperators2Bands.h"

#include "globals.h"

#define MKL_Complex16 std::complex<double>
#include "mkl.h"

void setupGlobalH2Bands(std::vector<std::complex<double>> &H){
  
  static_assert(dimPhoton == dimPhonon, "Phonon and Photon Hilbertspace cutoffs must match!");
  
  H = std::vector<std::complex<double>> (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));
  
  std::vector<std::complex<double>> HUncoupled;
  std::vector<std::complex<double>> HEBos;
  std::vector<std::complex<double>> freePhonon;
  
  setupHelectronic2Bands(HUncoupled);
  setupHEBos(HEBos);
  setupHFreePhon(freePhonon);
  addMatricies(HUncoupled, freePhonon, H);
  addMatricies(H, HEBos, H);
}

void setupHEBos(std::vector<std::complex<double>> &HEBos){
  
  HEBos = std::vector<std::complex<double>> (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));
  
  std::vector<std::complex<double>> hopInter0(dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));
  setupInterbandHop0(hopInter0);
  std::vector<std::complex<double>> hopInter1(dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));
  setupInterbandHop1(hopInter1);
  
  std::vector<std::complex<double>> A0;
  setupA0(A0);
  
  std::complex<double> alpha (gE, 0.);
  std::complex<double> beta (1., 0.);
  
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A0.data(), dimHOnePh,
              hopInter0.data(), dimHOnePh,
              &beta, HEBos.data(), dimHOnePh);
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A0.data(), dimHOnePh,
              hopInter0.data(), dimHOnePh,
              &beta, HEBos.data(), dimHOnePh);
  
  std::vector<std::complex<double>> A1;
  setupA1(A1);
  
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A1.data(), dimHOnePh,
              hopInter1.data(), dimHOnePh,
              &beta, HEBos.data(), dimHOnePh);
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A1.data(), dimHOnePh,
              hopInter1.data(), dimHOnePh,
              &beta, HEBos.data(), dimHOnePh);
  
}

void setupHFreePhon(std::vector<std::complex<double>> &freePhonon){
  
  freePhonon = std::vector<std::complex<double>> (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));
  
  std::vector<std::complex<double>> A0;
  setupA0(A0);
  
  std::complex<double> alpha (wPh, 0.);
  std::complex<double> beta (1., 0.);
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A0.data(), dimHOnePh,
              A0.data(), dimHOnePh,
              &beta, freePhonon.data(), dimHOnePh);
  
  std::vector<std::complex<double>> A1;
  setupA1(A1);
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A1.data(), dimHOnePh,
              A1.data(), dimHOnePh,
              &beta, freePhonon.data(), dimHOnePh);
  
}

