#include <vector>
#include <complex>
#include <iostream>

#include "globals.h"
#include "include/twoPhonons/setupBasicOperators.h"
#include "matrixOperations.h"
#include "include/twoPhonons/setUpGlobalHamiltonian.h"

#define MKL_Complex16 std::complex<double>

#include "mkl.h"

void setupGlobalH(std::vector<std::complex<double>> &H){

  H = std::vector<std::complex<double>> (dimHTwoPh * dimHTwoPh, std::complex<double> (0., 0.));

  std::vector<std::complex<double>> HUncoupled;
  std::vector<std::complex<double>> HEPh;
  std::vector<std::complex<double>> HPhPt;

  setupUncoupledHamiltonian(HUncoupled);
  setupEPhCoupling(HEPh);
  //setupQuadEPhCoupling(HEPh);
  setupPhPtCoupling(HPhPt);

  addMatricies(HUncoupled, HEPh, H);
  addMatricies(H, HPhPt, H);
}

void setupPhPtCoupling(std::vector<std::complex<double>> &HPhPt){

  HPhPt = std::vector<std::complex<double>> (dimHTwoPh * dimHTwoPh, std::complex<double> (0., 0.));

  std::vector<std::complex<double>> A;
  setupA(A);
  std::vector<std::complex<double>> ADag(A);
  dagger(ADag, dimHTwoPh);

  std::vector<std::complex<double>> APADag(A.size(), std::complex<double> (0., 0.));

  addMatricies(A, ADag, APADag);
  std::complex<double> alpha (wP * wP / (2. * wPt), 0.);
  std::complex<double> beta (0., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              APADag.data(), dimHTwoPh,
              APADag.data(), dimHTwoPh,
              &beta, HPhPt.data(), dimHTwoPh);

  std::vector<std::complex<double>> B1;
  setupB1(B1);
  std::vector<std::complex<double>> B1Dag(B1);
  dagger(B1Dag, dimHTwoPh);

  std::vector<std::complex<double>> B1MB1Dag(B1.size(), std::complex<double> (0., 0.));

  matrixAMinusB(B1, B1Dag, B1MB1Dag);

  alpha = std::complex<double> (0., wP * std::sqrt(wPh) / (2. * std::sqrt(wPt)));
  beta = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              APADag.data(), dimHTwoPh,
              B1MB1Dag.data(), dimHTwoPh,
              &beta, HPhPt.data(), dimHTwoPh);

  std::vector<std::complex<double>> B2;
  setupB2(B2);
  std::vector<std::complex<double>> B2Dag(B2);
  dagger(B2Dag, dimHTwoPh);

  std::vector<std::complex<double>> B2MB2Dag(B2.size(), std::complex<double> (0., 0.));

  matrixAMinusB(B2, B2Dag, B2MB2Dag);

  alpha = std::complex<double> (0., wP * std::sqrt(wPh) / (2. * std::sqrt(wPt)));
  beta = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              APADag.data(), dimHTwoPh,
              B2MB2Dag.data(), dimHTwoPh,
              &beta, HPhPt.data(), dimHTwoPh);
}

void setupQuadEPhCoupling(std::vector<std::complex<double>> &HEPh){

  HEPh = std::vector<std::complex<double>> (dimHTwoPh * dimHTwoPh, std::complex<double> (0., 0.));

  std::vector<std::complex<double>> dOcc;
  setupDOcc(dOcc);

  std::vector<std::complex<double>> B1;
  setupB1(B1);
  std::vector<std::complex<double>> B1Dag(B1);
  dagger(B1Dag, dimHTwoPh);

  std::vector<std::complex<double>> B1PB1Dag(B1.size(), std::complex<double> (0., 0.));

  addMatricies(B1, B1Dag, B1PB1Dag);

  std::complex<double> alpha (gPh / (2. * wPh), 0.);
  std::complex<double> beta (0., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              B1PB1Dag.data(), dimHTwoPh,
              B1PB1Dag.data(), dimHTwoPh,
              &beta, HEPh.data(), dimHTwoPh);

  alpha = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              HEPh.data(), dimHTwoPh,
              dOcc.data(), dimHTwoPh,
              &beta, HEPh.data(), dimHTwoPh);

  std::vector<std::complex<double>> B2;
  setupB2(B2);
  std::vector<std::complex<double>> B2Dag(B2);
  dagger(B2Dag, dimHTwoPh);

  std::vector<std::complex<double>> B2PB2Dag(B2.size(), std::complex<double> (0., 0.));

  addMatricies(B2, B2Dag, B2PB2Dag);

  std::vector<std::complex<double>> B2PB2DagSqr(B2.size(), std::complex<double> (0., 0.));

  alpha = std::complex<double> (gPh / (2. * wPh), 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              B2PB2Dag.data(), dimHTwoPh,
              B2PB2Dag.data(), dimHTwoPh,
              &beta, B2PB2DagSqr.data(), dimHTwoPh);

  alpha = std::complex<double> (1., 0.);
  beta = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              B2PB2DagSqr.data(), dimHTwoPh,
              dOcc.data(), dimHTwoPh,
              &beta, HEPh.data(), dimHTwoPh);

}


void setupEPhCoupling(std::vector<std::complex<double>> &HEPh){

  HEPh = std::vector<std::complex<double>> (dimHTwoPh * dimHTwoPh, std::complex<double> (0., 0.));

  std::vector<std::complex<double>> n1;
  setupOccSiteI(n1, 0ul);

  std::vector<std::complex<double>> B1;
  setupB1(B1);
  std::vector<std::complex<double>> B1Dag(B1);
  dagger(B1Dag, dimHTwoPh);

  std::vector<std::complex<double>> B1PB1Dag(B1.size(), std::complex<double> (0., 0.));

  addMatricies(B1, B1Dag, B1PB1Dag);

  std::complex<double> alpha (gPh / (2. * wPh), 0.);
  std::complex<double> beta (0., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              B1PB1Dag.data(), dimHTwoPh,
              B1PB1Dag.data(), dimHTwoPh,
              &beta, HEPh.data(), dimHTwoPh);

  alpha = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              HEPh.data(), dimHTwoPh,
              n1.data(), dimHTwoPh,
              &beta, HEPh.data(), dimHTwoPh);

  std::vector<std::complex<double>> n2;
  setupOccSiteI(n2, 1ul);

  std::vector<std::complex<double>> B2;
  setupB2(B2);
  std::vector<std::complex<double>> B2Dag(B2);
  dagger(B2Dag, dimHTwoPh);

  std::vector<std::complex<double>> B2PB2Dag(B2.size(), std::complex<double> (0., 0.));

  addMatricies(B2, B2Dag, B2PB2Dag);

  std::vector<std::complex<double>> B2PB2DagSqr(B2.size(), std::complex<double> (0., 0.));

  alpha = std::complex<double> (gPh / (2. * wPh), 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              B2PB2Dag.data(), dimHTwoPh,
              B2PB2Dag.data(), dimHTwoPh,
              &beta, B2PB2DagSqr.data(), dimHTwoPh);

  alpha = std::complex<double> (1., 0.);
  beta = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              B2PB2DagSqr.data(), dimHTwoPh,
              n2.data(), dimHTwoPh,
              &beta, HEPh.data(), dimHTwoPh);

}


void setupUncoupledHamiltonian(std::vector<std::complex<double>> &H) {


  std::vector<std::complex<double>> HElectronic;
  setupHelectronic(HElectronic);

  H = std::vector<std::complex<double>>(HElectronic);


  std::vector<std::complex<double>> B2;
  setupB2(B2);
  //std::vector<std::complex<double>> B2Dag(B2);
  //dagger(B2Dag, dimHTwoPh);

  std::vector<std::complex<double>> B1;
  setupB1(B1);
  //std::vector<std::complex<double>> B1Dag(B1);
  //dagger(B1Dag, dimHTwoPh);

  std::vector<std::complex<double>> A;
  setupA(A);
  //std::vector<std::complex<double>> ADag(A);
  //dagger(ADag, dimHTwoPh);

  std::complex<double> alpha (wPh, 0.);
  std::complex<double> beta (1., 0.);


  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              B2.data(), dimHTwoPh,
              B2.data(), dimHTwoPh,
              &beta, H.data(), dimHTwoPh);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              B1.data(), dimHTwoPh,
              B1.data(), dimHTwoPh,
              &beta, H.data(), dimHTwoPh);

  alpha = std::complex<double> (wPt, 0.);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHTwoPh, dimHTwoPh, dimHTwoPh, &alpha,
              A.data(), dimHTwoPh,
              A.data(), dimHTwoPh,
              &beta, H.data(), dimHTwoPh);

}
