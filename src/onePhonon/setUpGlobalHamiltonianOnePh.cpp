#include <vector>
#include <complex>
#include <iostream>

#include "globals.h"
#include "setUpGlobalHamiltonianOnePh.h"
#include "include/onePhonon/setupBasicOperatorsOnePh.h"
#include "matrixOperations.h"

#define MKL_Complex16 std::complex<double>
#include "mkl.h"

void setupGlobalHOnePh(std::vector<std::complex<double>> &H){

  H = std::vector<std::complex<double>> (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));

  std::vector<std::complex<double>> HUncoupled;
  std::vector<std::complex<double>> HEPh;
  std::vector<std::complex<double>> HPhPt;

  setupUncoupledHamiltonianOnePh(HUncoupled);
  setupQuadEPhCouplingOnePh(HEPh);
  setupPhPtCouplingOnePh(HPhPt);

  addMatricies(HUncoupled, HEPh, H);
  addMatricies(H, HPhPt, H);
}

void setupPhPtCouplingOnePh(std::vector<std::complex<double>> &HPhPt){

  HPhPt = std::vector<std::complex<double>> (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));

  std::vector<std::complex<double>> A;
  setupAOnePh(A);
  std::vector<std::complex<double>> ADag(A);
  dagger(ADag, dimHOnePh);

  std::vector<std::complex<double>> APADag(A.size(), std::complex<double> (0., 0.));

  addMatricies(A, ADag, APADag);
  std::complex<double> alpha (wP * wP / (4. * wPt), 0.);
  std::complex<double> beta (0., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              APADag.data(), dimHOnePh,
              APADag.data(), dimHOnePh,
              &beta, HPhPt.data(), dimHOnePh);

  std::vector<std::complex<double>> B;
  setupB(B);
  std::vector<std::complex<double>> BDag(B);
  dagger(BDag, dimHOnePh);

  std::vector<std::complex<double>> BMBDag(B.size(), std::complex<double> (0., 0.));

  matrixAMinusB(B, BDag, BMBDag);

  alpha = std::complex<double> (0., wP * std::sqrt(wPh) / (2. * std::sqrt(wPt)));
  beta = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              APADag.data(), dimHOnePh,
              BMBDag.data(), dimHOnePh,
              &beta, HPhPt.data(), dimHOnePh);
}

void setupQuadEPhCouplingOnePh(std::vector<std::complex<double>> &HEPh){

  HEPh = std::vector<std::complex<double>> (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));

  std::vector<std::complex<double>> dOcc;
  setupDOccOnePh(dOcc);

  std::vector<std::complex<double>> B;
  setupB(B);
  std::vector<std::complex<double>> BDag(B);
  dagger(BDag, dimHOnePh);

  std::vector<std::complex<double>> BPBDag(B.size(), std::complex<double> (0., 0.));

  addMatricies(B, BDag, BPBDag);

  std::complex<double> alpha (gPh / (2. * wPh), 0.);
  std::complex<double> beta (0., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              BPBDag.data(), dimHOnePh,
              BPBDag.data(), dimHOnePh,
              &beta, HEPh.data(), dimHOnePh);

  alpha = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              HEPh.data(), dimHOnePh,
              dOcc.data(), dimHOnePh,
              &beta, HEPh.data(), dimHOnePh);
}

void setupUncoupledHamiltonianOnePh(std::vector<std::complex<double>> &H) {


  std::vector<std::complex<double>> HElectronic;
  setupHelectronicOnePh(HElectronic);

  H = std::vector<std::complex<double>>(HElectronic);
  
  std::vector<std::complex<double>> B;
  setupB(B);

  std::vector<std::complex<double>> A;
  setupAOnePh(A);

  std::complex<double> alpha (wPh, 0.);
  std::complex<double> beta (1., 0.);


  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B.data(), dimHOnePh,
              B.data(), dimHOnePh,
              &beta, H.data(), dimHOnePh);

  alpha = std::complex<double> (wPt, 0.);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A.data(), dimHOnePh,
              A.data(), dimHOnePh,
              &beta, H.data(), dimHOnePh);
}
