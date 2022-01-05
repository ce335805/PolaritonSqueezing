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

  std::complex<double> alpha (wP * wP / (4. * wPt), 0.);
  std::complex<double> beta (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A.data(), dimHOnePh,
              A.data(), dimHOnePh,
              &beta, HPhPt.data(), dimHOnePh);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasConjTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A.data(), dimHOnePh,
              A.data(), dimHOnePh,
              &beta, HPhPt.data(), dimHOnePh);

  alpha = std::complex<double> (wP * wP / (2. * wPt), 0.);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A.data(), dimHOnePh,
              A.data(), dimHOnePh,
              &beta, HPhPt.data(), dimHOnePh);

  std::vector<std::complex<double>> B;
  setupB(B);

  alpha = std::complex<double> (0., wP * std::sqrt(wPh) / (2. * std::sqrt(wPt)));
  beta = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A.data(), dimHOnePh,
              B.data(), dimHOnePh,
              &beta, HPhPt.data(), dimHOnePh);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A.data(), dimHOnePh,
              B.data(), dimHOnePh,
              &beta, HPhPt.data(), dimHOnePh);

  alpha = std::complex<double> (0., - wP * std::sqrt(wPh) / (2. * std::sqrt(wPt)));

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A.data(), dimHOnePh,
              B.data(), dimHOnePh,
              &beta, HPhPt.data(), dimHOnePh);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasConjTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A.data(), dimHOnePh,
              B.data(), dimHOnePh,
              &beta, HPhPt.data(), dimHOnePh);

}

void setupQuadEPhCouplingOnePh(std::vector<std::complex<double>> &HEPh){

  HEPh = std::vector<std::complex<double>> (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));

  std::vector<std::complex<double>> dOcc;
  //setupDOccOnePh(dOcc);
  setupDOccOnePhNoPHS(dOcc);

  std::vector<std::complex<double>> B;
  setupB(B);

  std::vector<std::complex<double>> BBDag(B.size(), std::complex<double> (0., 0.));
  std::vector<std::complex<double>> BDagB(B.size(), std::complex<double> (0., 0.));
  std::vector<std::complex<double>> BB(B.size(), std::complex<double> (0., 0.));
  std::vector<std::complex<double>> BDagBDag(B.size(), std::complex<double> (0., 0.));

  std::complex<double> alpha (1., 0.);
  std::complex<double> beta (0., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B.data(), dimHOnePh,
              B.data(), dimHOnePh,
              &beta, BB.data(), dimHOnePh);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasConjTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B.data(), dimHOnePh,
              B.data(), dimHOnePh,
              &beta, BDagBDag.data(), dimHOnePh);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B.data(), dimHOnePh,
              B.data(), dimHOnePh,
              &beta, BDagB.data(), dimHOnePh);


  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B.data(), dimHOnePh,
              B.data(), dimHOnePh,
              &beta, BBDag.data(), dimHOnePh);

  alpha = std::complex<double>  (gPh, 0.);
  beta = std::complex<double> (1., 0.);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              BB.data(), dimHOnePh,
              dOcc.data(), dimHOnePh,
              &beta, HEPh.data(), dimHOnePh);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              BDagBDag.data(), dimHOnePh,
              dOcc.data(), dimHOnePh,
              &beta, HEPh.data(), dimHOnePh);

  //cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
  //            dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
  //            BBDag.data(), dimHOnePh,
  //            dOcc.data(), dimHOnePh,
  //            &beta, HEPh.data(), dimHOnePh);

  alpha = std::complex<double>  (2. * gPh, 0.);
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              BDagB.data(), dimHOnePh,
              dOcc.data(), dimHOnePh,
              &beta, HEPh.data(), dimHOnePh);

  std::vector<std::complex<double>> unitMat (dimHOnePh * dimHOnePh, std::complex<double> (0., 0.));
  for(ulong ind = 0ul; ind < dimHOnePh; ++ind){
    unitMat[ind * dimHOnePh + ind] = 1.;
  }

  alpha = std::complex<double>  (gPh, 0.);
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              unitMat.data(), dimHOnePh,
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
