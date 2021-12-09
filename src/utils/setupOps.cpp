#include <vector>
#include <complex>
#include <iostream>
#include <string>

#include "globals.h"
#include "setupBasicOperatorsOnePh.h"
#include "matrixOperations.h"
#include "setUpGlobalHamiltonianOnePh.h"
#include "calcGS.h"
#include "timeStep.h"
#include "setupBasicOperators.h"
#include "setUpGlobalHamiltonian.h"
#include "evalExpectation.h"
#include "writeStuffToHdf5.h"

#include "H5Cpp.h"

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

