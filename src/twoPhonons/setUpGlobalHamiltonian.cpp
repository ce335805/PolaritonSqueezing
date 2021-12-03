#include <vector>
#include <complex>
#include <iostream>

#include "globals.h"
#include "include/twoPhonons/setupBasicOperators.h"
#include "matrixOperations.h"
#include "include/twoPhonons/setUpGlobalHamiltonian.cpp.h"

#define MKL_Complex16 std::complex<double>

#include "mkl.h"

void setupGlobalH(std::vector<std::complex<double>> &H){

  H = std::vector<std::complex<double>> (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));

  std::vector<std::complex<double>> HUncoupled;
  std::vector<std::complex<double>> HEPh;
  std::vector<std::complex<double>> HPhPt;

  setupUncoupledHamiltonian(HUncoupled);
  //setupEPhCoupling(HEPh);
  setupQuadEPhCoupling(HEPh);
  setupPhPtCoupling(HPhPt);

  addMatricies(HUncoupled, HEPh, H);
  addMatricies(H, HPhPt, H);
}

void setupPhPtCoupling(std::vector<std::complex<double>> &HPhPt){

  HPhPt = std::vector<std::complex<double>> (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));

  std::vector<std::complex<double>> A;
  setupA(A);
  std::vector<std::complex<double>> ADag(A);
  dagger(ADag, dimHOnePh);

  std::vector<std::complex<double>> APADag(A.size(), std::complex<double> (0., 0.));

  addMatricies(A, ADag, APADag);
  std::complex<double> alpha (wP * wP / (2. * wPt), 0.);
  std::complex<double> beta (0., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              APADag.data(), dimHOnePh,
              APADag.data(), dimHOnePh,
              &beta, HPhPt.data(), dimHOnePh);

  std::vector<std::complex<double>> B1;
  setupB1(B1);
  std::vector<std::complex<double>> B1Dag(B1);
  dagger(B1Dag, dimHOnePh);

  std::vector<std::complex<double>> B1MB1Dag(B1.size(), std::complex<double> (0., 0.));

  matrixAMinusB(B1, B1Dag, B1MB1Dag);

  alpha = std::complex<double> (0., wP * std::sqrt(wPh) / (2. * std::sqrt(wPt)));
  beta = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              APADag.data(), dimHOnePh,
              B1MB1Dag.data(), dimHOnePh,
              &beta, HPhPt.data(), dimHOnePh);

  std::vector<std::complex<double>> B2;
  setupB2(B2);
  std::vector<std::complex<double>> B2Dag(B2);
  dagger(B2Dag, dimHOnePh);

  std::vector<std::complex<double>> B2MB2Dag(B2.size(), std::complex<double> (0., 0.));

  matrixAMinusB(B2, B2Dag, B2MB2Dag);

  alpha = std::complex<double> (0., wP * std::sqrt(wPh) / (2. * std::sqrt(wPt)));
  beta = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              APADag.data(), dimHOnePh,
              B2MB2Dag.data(), dimHOnePh,
              &beta, HPhPt.data(), dimHOnePh);
}

void setupQuadEPhCoupling(std::vector<std::complex<double>> &HEPh){

  HEPh = std::vector<std::complex<double>> (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));

  std::vector<std::complex<double>> dOcc;
  setupDOcc(dOcc);

  std::vector<std::complex<double>> B1;
  setupB1(B1);
  std::vector<std::complex<double>> B1Dag(B1);
  dagger(B1Dag, dimHOnePh);

  std::vector<std::complex<double>> B1PB1Dag(B1.size(), std::complex<double> (0., 0.));

  addMatricies(B1, B1Dag, B1PB1Dag);

  std::complex<double> alpha (gPh / (2. * wPh), 0.);
  std::complex<double> beta (0., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B1PB1Dag.data(), dimHOnePh,
              B1PB1Dag.data(), dimHOnePh,
              &beta, HEPh.data(), dimHOnePh);

  alpha = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              HEPh.data(), dimHOnePh,
              dOcc.data(), dimHOnePh,
              &beta, HEPh.data(), dimHOnePh);

  std::vector<std::complex<double>> B2;
  setupB2(B2);
  std::vector<std::complex<double>> B2Dag(B2);
  dagger(B2Dag, dimHOnePh);

  std::vector<std::complex<double>> B2PB2Dag(B2.size(), std::complex<double> (0., 0.));

  addMatricies(B2, B2Dag, B2PB2Dag);

  std::vector<std::complex<double>> B2PB2DagSqr(B2.size(), std::complex<double> (0., 0.));

  alpha = std::complex<double> (gPh / (2. * wPh), 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B2PB2Dag.data(), dimHOnePh,
              B2PB2Dag.data(), dimHOnePh,
              &beta, B2PB2DagSqr.data(), dimHOnePh);

  alpha = std::complex<double> (1., 0.);
  beta = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B2PB2DagSqr.data(), dimHOnePh,
              dOcc.data(), dimHOnePh,
              &beta, HEPh.data(), dimHOnePh);

}


void setupEPhCoupling(std::vector<std::complex<double>> &HEPh){

  HEPh = std::vector<std::complex<double>> (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));

  std::vector<std::complex<double>> n1;
  setupOccSiteI(n1, 0ul);

  std::vector<std::complex<double>> B1;
  setupB1(B1);
  std::vector<std::complex<double>> B1Dag(B1);
  dagger(B1Dag, dimHOnePh);

  std::vector<std::complex<double>> B1PB1Dag(B1.size(), std::complex<double> (0., 0.));

  addMatricies(B1, B1Dag, B1PB1Dag);

  std::complex<double> alpha (gPh / (2. * wPh), 0.);
  std::complex<double> beta (0., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B1PB1Dag.data(), dimHOnePh,
              B1PB1Dag.data(), dimHOnePh,
              &beta, HEPh.data(), dimHOnePh);

  alpha = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              HEPh.data(), dimHOnePh,
              n1.data(), dimHOnePh,
              &beta, HEPh.data(), dimHOnePh);

  std::vector<std::complex<double>> n2;
  setupOccSiteI(n2, 1ul);

  std::vector<std::complex<double>> B2;
  setupB2(B2);
  std::vector<std::complex<double>> B2Dag(B2);
  dagger(B2Dag, dimHOnePh);

  std::vector<std::complex<double>> B2PB2Dag(B2.size(), std::complex<double> (0., 0.));

  addMatricies(B2, B2Dag, B2PB2Dag);

  std::vector<std::complex<double>> B2PB2DagSqr(B2.size(), std::complex<double> (0., 0.));

  alpha = std::complex<double> (gPh / (2. * wPh), 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B2PB2Dag.data(), dimHOnePh,
              B2PB2Dag.data(), dimHOnePh,
              &beta, B2PB2DagSqr.data(), dimHOnePh);

  alpha = std::complex<double> (1., 0.);
  beta = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B2PB2DagSqr.data(), dimHOnePh,
              n2.data(), dimHOnePh,
              &beta, HEPh.data(), dimHOnePh);

}


void setupUncoupledHamiltonian(std::vector<std::complex<double>> &H) {


  std::vector<std::complex<double>> HElectronic;
  setupHelectronic(HElectronic);

  H = std::vector<std::complex<double>>(HElectronic);


  std::vector<std::complex<double>> B2;
  setupB2(B2);
  //std::vector<std::complex<double>> B2Dag(B2);
  //dagger(B2Dag, dimHOnePh);

  std::vector<std::complex<double>> B1;
  setupB1(B1);
  //std::vector<std::complex<double>> B1Dag(B1);
  //dagger(B1Dag, dimHOnePh);

  std::vector<std::complex<double>> A;
  setupA(A);
  //std::vector<std::complex<double>> ADag(A);
  //dagger(ADag, dimHOnePh);

  std::complex<double> alpha (wPh, 0.);
  std::complex<double> beta (1., 0.);


  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B2.data(), dimHOnePh,
              B2.data(), dimHOnePh,
              &beta, H.data(), dimHOnePh);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B1.data(), dimHOnePh,
              B1.data(), dimHOnePh,
              &beta, H.data(), dimHOnePh);

  alpha = std::complex<double> (wPt, 0.);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A.data(), dimHOnePh,
              A.data(), dimHOnePh,
              &beta, H.data(), dimHOnePh);

}
