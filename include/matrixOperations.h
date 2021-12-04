#ifndef PHONONSQUEEZING_MATRIXOPERATIONS_H
#define PHONONSQUEEZING_MATRIXOPERATIONS_H

#include <vector>
#include <complex>

void dagger(std::vector<std::complex<double>> &Mat, ulong dim);

std::vector<double> diagonalize(std::vector<std::complex<double>> &Mat, ulong dim, const char eVecs);

void addMatricies(const std::vector<std::complex<double>> &A,
                  const std::vector<std::complex<double>> &B,
                  std::vector<std::complex<double>> &C);

void addMatricies(const std::vector<std::complex<double>> &A,
                  const std::complex<double> alpha,
                  const std::vector<std::complex<double>> &B,
                  const std::complex<double> beta,
                  std::vector<std::complex<double>> &C);

void matrixAMinusB(const std::vector<std::complex<double>> &A,
                   const std::vector<std::complex<double>> &B,
                   std::vector<std::complex<double>> &C);

#endif //PHONONSQUEEZING_MATRIXOPERATIONS_H
