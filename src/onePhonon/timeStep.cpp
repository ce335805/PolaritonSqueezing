#include <vector>
#include <complex>
#include <cmath>
#include <iostream>

#include "globals.h"
#include "matrixOperations.h"

#define MKL_Complex16 std::complex<double>

#include "mkl.h"

void calcTimeStep(const double tPoint,
                  const double pumpPrefac,
                  const std::vector<std::complex<double>> &H,
                  const std::vector<std::complex<double>> &Xpt,
                  std::vector<std::complex<double>> &state) {

  std::vector<std::complex<double>> H1(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> H2(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  addMatricies(H, 1., Xpt, pumpPrefac * std::sin(tPoint - (.5 - std::sqrt(3.) / 6.) * dt), H1);
  addMatricies(H, 1., Xpt, pumpPrefac * std::sin(tPoint - (.5 + std::sqrt(3.) / 6.) * dt), H2);

  ////check that GS is really GS/////////
  //std::vector<std::complex<double>> gsAfterApplyingH (dimHOnePh, std::complex<double> (0., 0.));
//
  //std::complex<double> alpha(1., 0.);
  //std::complex<double> beta(0., 0.);
//
  //cblas_zgemv(CblasRowMajor, CblasConjTrans,
  //            dimHOnePh, dimHOnePh, &alpha, H2.data(),
  //            dimHOnePh, state.data(), 1,
  //            &beta, gsAfterApplyingH.data(), 1);
//
  //for(ulong ind = 0ul; ind < dimHOnePh; ++ind){
  //  if(abs(state[ind] - gsAfterApplyingH[ind]) > 1e-10){
  //    std::cout << state[ind] / gsAfterApplyingH[ind] << '\n';
  //  }
  //}
  /////////


  std::vector<std::complex<double>> uMat(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> vMat(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  addMatricies(H1, std::complex<double>((3. - 2. * std::sqrt(2.)) / 12., 0.), H2, std::complex<double>((3. + 2. * std::sqrt(2.)) / 12., 0.), uMat);
  addMatricies(H1, std::complex<double>((3. + 2. * std::sqrt(2.)) / 12., 0.), H2, std::complex<double>((3. - 2. * std::sqrt(2.)) / 12., 0.), vMat);

  std::vector<double> eValsU(dimHOnePh, 0.);
  std::vector<double> eValsV(dimHOnePh, 0.);

  eValsU = diagonalize(uMat, dimHOnePh, 'V');
  eValsV = diagonalize(vMat, dimHOnePh, 'V');

  for (ulong ind1 = 0ul; ind1 < dimHOnePh; ++ind1) {
    for (ulong ind2 = 0ul; ind2 < dimHOnePh; ++ind2) {
      if (ind1 == ind2) {
        //H1[ind1 * dimHOnePh + ind2] = std::complex<double> (1., 0.);
        //H1[ind1 * dimHOnePh + ind2] = eValsU[ind1];
        H1[ind1 * dimHOnePh + ind2] = std::exp(- II * eValsU[ind1]);
        //H2[ind1 * dimHOnePh + ind2] = std::complex<double> (1., 0.);
        //H2[ind1 * dimHOnePh + ind2] = eValsV[ind1];
        H2[ind1 * dimHOnePh + ind2] = std::exp(- II * eValsV[ind1]);
      } else {
        H1[ind1 * dimHOnePh + ind2] = std::complex<double> (0. ,0.);
        H2[ind1 * dimHOnePh + ind2] = std::complex<double> (0., 0.);
      }
    }
  }

  std::complex<double> alpha(1., 0.);
  std::complex<double> beta(0., 0.);

  //alpha = std::complex<double>(1., 0.);
  //beta = std::complex<double>(0., 0.);

  std::vector<std::complex<double>> vecInter(dimHOnePh, std::complex<double> (0., 0.));

  cblas_zgemv(CblasRowMajor, CblasConjTrans,
              dimHOnePh, dimHOnePh, &alpha, vMat.data(),
              dimHOnePh, state.data(), 1,
              &beta, vecInter.data(), 1);

  state = std::vector<std::complex<double>> (vecInter);

  cblas_zgemv(CblasRowMajor, CblasNoTrans,
              dimHOnePh, dimHOnePh, &alpha, H2.data(),
              dimHOnePh, state.data(), 1,
              &beta, vecInter.data(), 1);

  state = std::vector<std::complex<double>> (vecInter);

  cblas_zgemv(CblasRowMajor, CblasNoTrans,
              dimHOnePh, dimHOnePh, &alpha, vMat.data(),
              dimHOnePh, state.data(), 1,
              &beta, vecInter.data(), 1);

  state = std::vector<std::complex<double>> (vecInter);

  cblas_zgemv(CblasRowMajor, CblasConjTrans,
              dimHOnePh, dimHOnePh, &alpha, uMat.data(),
              dimHOnePh, state.data(), 1,
              &beta, vecInter.data(), 1);

  state = std::vector<std::complex<double>> (vecInter);

  cblas_zgemv(CblasRowMajor, CblasNoTrans,
              dimHOnePh, dimHOnePh, &alpha, H1.data(),
              dimHOnePh, state.data(), 1,
              &beta, vecInter.data(), 1);

  state = std::vector<std::complex<double>> (vecInter);

  cblas_zgemv(CblasRowMajor, CblasNoTrans,
              dimHOnePh, dimHOnePh, &alpha, uMat.data(),
              dimHOnePh, state.data(), 1,
              &beta, vecInter.data(), 1);

  state = std::vector<std::complex<double>> (vecInter);

}
