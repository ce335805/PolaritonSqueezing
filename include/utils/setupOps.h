#ifndef PHONONSQUEEZING_SETUPOPS_H
#define PHONONSQUEEZING_SETUPOPS_H


void setupOpsOnePh(std::vector<std::complex<double>> &dOcc,
                   std::vector<std::complex<double>> &Xpt,
                   std::vector<std::complex<double>> &XptSqr,
                   std::vector<std::complex<double>> &Npt,
                   std::vector<std::complex<double>> &X1ph,
                   std::vector<std::complex<double>> &X1phSqr,
                   std::vector<std::complex<double>> &N1ph,
                   std::vector<std::complex<double>> &X2ph,
                   std::vector<std::complex<double>> &X2phSqr,
                   std::vector<std::complex<double>> &N2ph
);

void setupOpsTwoPh(std::vector<std::complex<double>> &dOcc,
                   std::vector<std::complex<double>> &Xpt,
                   std::vector<std::complex<double>> &XptSqr,
                   std::vector<std::complex<double>> &Npt,
                   std::vector<std::complex<double>> &X1ph,
                   std::vector<std::complex<double>> &X1phSqr,
                   std::vector<std::complex<double>> &N1ph,
                   std::vector<std::complex<double>> &X2ph,
                   std::vector<std::complex<double>> &X2phSqr,
                   std::vector<std::complex<double>> &N2ph);

void setupOpsOnlyPhot(std::vector<std::complex<double>> &dOcc,
                      std::vector<std::complex<double>> &Xpt,
                      std::vector<std::complex<double>> &XptSqr,
                      std::vector<std::complex<double>> &Npt);

#endif //PHONONSQUEEZING_SETUPOPS_H
