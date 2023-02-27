#ifndef PHONONSQUEEZING_EVALEXPECTATION_H
#define PHONONSQUEEZING_EVALEXPECTATION_H

double evalExpectation(const std::vector<std::complex<double>> &op,
                       const std::vector<std::complex<double>> &vec,
                       const ulong dimH);

void evalAllExpectations(std::vector<double> &expecValues,
                         const std::vector<std::complex<double>> &op,
                         const std::vector<std::complex<double>> &eVecs,
                         const ulong dimH);

#endif //PHONONSQUEEZING_EVALEXPECTATION_H
