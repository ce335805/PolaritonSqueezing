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

void setupOps2Bands(std::vector<std::complex<double>> &dOcc0,
                    std::vector<std::complex<double>> &dOcc1,
                    std::vector<std::complex<double>> &dOccUpDn,
                    std::vector<std::complex<double>> &dOccSigSig,
                    std::vector<std::complex<double>> &n0,
                    std::vector<std::complex<double>> &n1);

void setupOps2BandsDrive(
    std::vector<std::complex<double>> &dOcc0,
    std::vector<std::complex<double>> &dOcc1,
    std::vector<std::complex<double>> &dOccUpDn,
    std::vector<std::complex<double>> &dOccSigSig,
    std::vector<std::complex<double>> &n0,
    std::vector<std::complex<double>> &n1,
    std::vector<std::complex<double>> &Xph1,
    std::vector<std::complex<double>> &Xph1Sqr,
    std::vector<std::complex<double>> &Nph1,
    std::vector<std::complex<double>> &Xph2,
    std::vector<std::complex<double>> &Xph2Sqr,
    std::vector<std::complex<double>> &Nph2);

void setupOps2BandsSpectrum(
    std::vector<std::complex<double>> &nc0,
    std::vector<std::complex<double>> &nd0,
    std::vector<std::complex<double>> &nc1,
    std::vector<std::complex<double>> &nd1,
    std::vector<std::complex<double>> &Nph0,
    std::vector<std::complex<double>> &Nph1
);

#endif //PHONONSQUEEZING_SETUPOPS_H
