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
                  const std::vector<std::complex<double>> &ODrive,
                  std::vector<std::complex<double>> &state,
                  const ulong dim) {

  std::vector<std::complex<double>> H1(dim * dim, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> H2(dim * dim, std::complex<double>(0., 0.));

  addMatricies(H, 1., ODrive, pumpPrefac * std::sin(wDrive * (tPoint + (.5 - std::sqrt(3.) / 6.) * dt)), H1);
  addMatricies(H, 1., ODrive, pumpPrefac * std::sin(wDrive * (tPoint + (.5 + std::sqrt(3.) / 6.) * dt)), H2);


  std::vector<std::complex<double>> uMat(dim * dim, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> vMat(dim * dim, std::complex<double>(0., 0.));

  addMatricies(H1, std::complex<double>((3. - 2. * std::sqrt(2.)) / 12., 0.), H2, std::complex<double>((3. + 2. * std::sqrt(2.)) / 12., 0.), uMat);
  addMatricies(H1, std::complex<double>((3. + 2. * std::sqrt(2.)) / 12., 0.), H2, std::complex<double>((3. - 2. * std::sqrt(2.)) / 12., 0.), vMat);

  std::vector<double> eValsU(dim, 0.);
  std::vector<double> eValsV(dim, 0.);

  eValsU = diagonalize(uMat, dim, 'V');
  eValsV = diagonalize(vMat, dim, 'V');

  for (ulong ind1 = 0ul; ind1 < dim; ++ind1) {
    for (ulong ind2 = 0ul; ind2 < dim; ++ind2) {
      if (ind1 == ind2) {
        H1[ind1 * dim + ind2] = std::exp(- II * eValsU[ind1] * dt);
        H2[ind1 * dim + ind2] = std::exp(- II * eValsV[ind1] * dt);
      } else {
        H1[ind1 * dim + ind2] = std::complex<double> (0. ,0.);
        H2[ind1 * dim + ind2] = std::complex<double> (0., 0.);
      }
    }
  }

  std::complex<double> alpha(1., 0.);
  std::complex<double> beta(0., 0.);

  std::vector<std::complex<double>> vecInter(dim, std::complex<double> (0., 0.));

  cblas_zgemv(CblasRowMajor, CblasConjTrans,
              dim, dim, &alpha, vMat.data(),
              dim, state.data(), 1,
              &beta, vecInter.data(), 1);

  state = std::vector<std::complex<double>> (vecInter);

  cblas_zgemv(CblasRowMajor, CblasNoTrans,
              dim, dim, &alpha, H2.data(),
              dim, state.data(), 1,
              &beta, vecInter.data(), 1);

  state = std::vector<std::complex<double>> (vecInter);

  cblas_zgemv(CblasRowMajor, CblasNoTrans,
              dim, dim, &alpha, vMat.data(),
              dim, state.data(), 1,
              &beta, vecInter.data(), 1);

  state = std::vector<std::complex<double>> (vecInter);

  cblas_zgemv(CblasRowMajor, CblasConjTrans,
              dim, dim, &alpha, uMat.data(),
              dim, state.data(), 1,
              &beta, vecInter.data(), 1);

  state = std::vector<std::complex<double>> (vecInter);

  cblas_zgemv(CblasRowMajor, CblasNoTrans,
              dim, dim, &alpha, H1.data(),
              dim, state.data(), 1,
              &beta, vecInter.data(), 1);

  state = std::vector<std::complex<double>> (vecInter);

  cblas_zgemv(CblasRowMajor, CblasNoTrans,
              dim, dim, &alpha, uMat.data(),
              dim, state.data(), 1,
              &beta, vecInter.data(), 1);

  state = std::vector<std::complex<double>> (vecInter);

}
