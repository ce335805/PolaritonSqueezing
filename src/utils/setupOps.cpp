#include <vector>
#include <complex>
#include <iostream>
#include <string>

#include "globals.h"
#include "setupBasicOperatorsOnePh.h"
#include "matrixOperations.h"
#include "setupBasicOperators.h"
#include "setupBasicOperatorsOnlyPhot.h"
#include "setupBasicOperators2Bands.h"

#define MKL_Complex16 std::complex<double>

#include "mkl.h"


void setupOpsOnePh(std::vector<std::complex<double>> &dOcc,
                   std::vector<std::complex<double>> &Xpt,
                   std::vector<std::complex<double>> &XptSqr,
                   std::vector<std::complex<double>> &Npt,
                   std::vector<std::complex<double>> &X1ph,
                   std::vector<std::complex<double>> &X1phSqr,
                   std::vector<std::complex<double>> &N1ph,
                   std::vector<std::complex<double>> &X2ph,
                   std::vector<std::complex<double>> &X2phSqr,
                   std::vector<std::complex<double>> &N2ph) {
  
  const ulong dimH = dimHOnePh;
  
  Xpt = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> A(dimH * dimH, std::complex<double>(
      0., 0.));
  setupAOnePh(A);
  std::vector<std::complex<double>> ADag(A);
  dagger(ADag, dimH);
  
  
  addMatricies(ADag, 1. / std::sqrt(wPt), A, 1. / std::sqrt(wPt), Xpt);
  
  Npt = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  std::complex<double> alpha(1., 0.);
  std::complex<double> beta(0., 0.);
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              A.data(), dimH,
              A.data(), dimH,
              &beta, Npt.data(), dimH);
  
  dOcc = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  setupDOccOnePh(dOcc);
  
  XptSqr = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              Xpt.data(), dimH,
              Xpt.data(), dimH,
              &beta, XptSqr.data(), dimH);
  
  X1ph = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> B(dimH * dimH, std::complex<double>(
      0., 0.));
  setupB(B);
  std::vector<std::complex<double>> BDag(B);
  dagger(BDag, dimH);
  
  addMatricies(BDag, 1. / std::sqrt(wPh), B, 1. / std::sqrt(wPh), X1ph);
  
  N1ph = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              B.data(), dimH,
              B.data(), dimH,
              &beta, N1ph.data(), dimH);
  
  X1phSqr = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              X1ph.data(), dimH,
              X1ph.data(), dimH,
              &beta, X1phSqr.data(), dimH);
  
}

void setupOpsTwoPh(std::vector<std::complex<double>> &dOcc,
                   std::vector<std::complex<double>> &Xpt,
                   std::vector<std::complex<double>> &XptSqr,
                   std::vector<std::complex<double>> &Npt,
                   std::vector<std::complex<double>> &X1ph,
                   std::vector<std::complex<double>> &X1phSqr,
                   std::vector<std::complex<double>> &N1ph,
                   std::vector<std::complex<double>> &X2ph,
                   std::vector<std::complex<double>> &X2phSqr,
                   std::vector<std::complex<double>> &N2ph) {
  
  const ulong dimH = dimHTwoPh;
  
  Xpt = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> A(dimH * dimH, std::complex<double>(
      0., 0.));
  setupA(A);
  std::vector<std::complex<double>> ADag(A);
  dagger(ADag, dimH);
  
  
  addMatricies(ADag, 1. / std::sqrt(wPt), A, 1. / std::sqrt(wPt), Xpt);
  
  Npt = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  std::complex<double> alpha(1., 0.);
  std::complex<double> beta(0., 0.);
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              A.data(), dimH,
              A.data(), dimH,
              &beta, Npt.data(), dimH);
  
  dOcc = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  setupDOcc(dOcc);
  
  XptSqr = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              Xpt.data(), dimH,
              Xpt.data(), dimH,
              &beta, XptSqr.data(), dimH);
  
  X1ph = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> B1(dimH * dimH, std::complex<double>(0., 0.));
  setupB1(B1);
  std::vector<std::complex<double>> B1Dag(B1);
  dagger(B1Dag, dimH);
  
  addMatricies(B1Dag, 1. / std::sqrt(wPh), B1, 1. / std::sqrt(wPh), X1ph);
  
  N1ph = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              B1.data(), dimH,
              B1.data(), dimH,
              &beta, N1ph.data(), dimH);
  
  X1phSqr = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              X1ph.data(), dimH,
              X1ph.data(), dimH,
              &beta, X1phSqr.data(), dimH);
  
  X2ph = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> B2(dimH * dimH, std::complex<double>(
      0., 0.));
  setupB2(B2);
  std::vector<std::complex<double>> B2Dag(B2);
  dagger(B2Dag, dimH);
  
  addMatricies(B2Dag, 1. / std::sqrt(wPh), B2, 1. / std::sqrt(wPh), X2ph);
  
  N2ph = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              B2.data(), dimH,
              B2.data(), dimH,
              &beta, N2ph.data(), dimH);
  
  X2phSqr = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              X2ph.data(), dimH,
              X2ph.data(), dimH,
              &beta, X2phSqr.data(), dimH);
  
}


void setupOpsOnlyPhot(std::vector<std::complex<double>> &dOcc,
                      std::vector<std::complex<double>> &Xpt,
                      std::vector<std::complex<double>> &XptSqr,
                      std::vector<std::complex<double>> &Npt) {
  
  const ulong dimH = dimHOnlyPhot;
  
  Xpt = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> A(dimH * dimH, std::complex<double>(
      0., 0.));
  setupAOnlyPhot(A);
  std::vector<std::complex<double>> ADag(A);
  dagger(ADag, dimH);
  
  
  addMatricies(ADag, 1. / std::sqrt(2 * wPt), A, 1. / std::sqrt(2 * wPt), Xpt);
  
  
  XptSqr = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  std::complex<double> alpha(1., 0.);
  std::complex<double> beta(0., 0.);
  
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              Xpt.data(), dimH,
              Xpt.data(), dimH,
              &beta, XptSqr.data(), dimH);
  
  Npt = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              A.data(), dimH,
              A.data(), dimH,
              &beta, Npt.data(), dimH);
  
  dOcc = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  setupDOccOnlyPhot(dOcc);
  
}

void setupOps2Bands(std::vector<std::complex<double>> &dOcc0,
                    std::vector<std::complex<double>> &dOcc1,
                    std::vector<std::complex<double>> &dOccUpDn,
                    std::vector<std::complex<double>> &dOccSigSig,
                    std::vector<std::complex<double>> &n0,
                    std::vector<std::complex<double>> &n1) {
  
  setupDOcc0(dOcc0);
  setupDOcc1(dOcc1);
  setupInterOrbUpDn(dOccUpDn);
  setupInterOrbSigSig(dOccSigSig);
  setupN0(n0);
  setupN1(n1);
}

void setupOps2BandsSpectrum(
    std::vector<std::complex<double>> &nc0,
    std::vector<std::complex<double>> &nd0,
    std::vector<std::complex<double>> &nc1,
    std::vector<std::complex<double>> &nd1,
    std::vector<std::complex<double>> &Nph0,
    std::vector<std::complex<double>> &Nph1
) {
  setupNC0(nc0);
  setupND0(nd0);
  setupNC1(nc1);
  setupND1(nd1);
  
  const ulong dimH = dimHOnePh;
  
  std::complex<double> alpha(1., 0.);
  std::complex<double> beta(0., 0.);
  
  Nph0 = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> A0(dimH * dimH, std::complex<double>(0., 0.));
  setupA0(A0);
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              A0.data(), dimH,
              A0.data(), dimH,
              &beta, Nph0.data(), dimH);
  
  Nph1 = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> A1(dimH * dimH, std::complex<double>(0., 0.));
  setupA1(A1);
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              A1.data(), dimH,
              A1.data(), dimH,
              &beta, Nph1.data(), dimH);
  
}

void setupOps2BandsDrive(
    std::vector<std::complex<double>> &dOcc0,
    std::vector<std::complex<double>> &dOcc1,
    std::vector<std::complex<double>> &dOccUpDn,
    std::vector<std::complex<double>> &dOccSigSig,
    std::vector<std::complex<double>> &n0,
    std::vector<std::complex<double>> &n1,
    std::vector<std::complex<double>> &Xph1,
    std::vector<std::complex<double>> &Xph1Sqr,
    std::vector<std::complex<double>> &Nph1,
    std::vector<std::complex<double>> &Xph2,
    std::vector<std::complex<double>> &Xph2Sqr,
    std::vector<std::complex<double>> &Nph2) {
  
  const ulong dimH = dimHOnePh;
  
  setupOps2Bands(dOcc0, dOcc1, dOccUpDn, dOccSigSig, n0, n1);
  
  Xph1 = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> A(dimH * dimH, std::complex<double>(
      0., 0.));
  setupA0(A);
  std::vector<std::complex<double>> ADag(A);
  dagger(ADag, dimH);
  
  
  addMatricies(ADag, 1. / std::sqrt(wPh), A, 1. / std::sqrt(wPh), Xph1);
  
  Nph1 = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  std::complex<double> alpha(1., 0.);
  std::complex<double> beta(0., 0.);
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              A.data(), dimH,
              A.data(), dimH,
              &beta, Nph1.data(), dimH);
  
  Xph1Sqr = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              Xph1.data(), dimH,
              Xph1.data(), dimH,
              &beta, Xph1Sqr.data(), dimH);
  
  Xph2 = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  std::vector<std::complex<double>> B(dimH * dimH, std::complex<double>(
      0., 0.));
  setupA1(B);
  std::vector<std::complex<double>> BDag(B);
  dagger(BDag, dimH);
  
  addMatricies(BDag, 1. / std::sqrt(wPh), B, 1. / std::sqrt(wPh), Xph2);
  
  Nph2 = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              B.data(), dimH,
              B.data(), dimH,
              &beta, Nph2.data(), dimH);
  
  Xph2Sqr = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              Xph2.data(), dimH,
              Xph2.data(), dimH,
              &beta, Xph2Sqr.data(), dimH);
  
}
