#ifndef PHONONSQUEEZING_WRITESTUFFTOHDF5_H
#define PHONONSQUEEZING_WRITESTUFFTOHDF5_H

void writeStuffToHdf5(
    const std::vector<double> &times,
    const std::vector<double> &pumpFunction,
    const std::vector<double> &dOcc,
    const std::vector<double> &Xpt,
    const std::vector<double> &XptSqr,
    const std::vector<double> &Npt,
    const std::vector<double> &X1ph,
    const std::vector<double> &X1phSqr,
    const std::vector<double> &N1ph,
    const std::vector<double> &X2ph,
    const std::vector<double> &X2phSqr,
    const std::vector<double> &N2ph,
    const std::string &filename,
    const bool twoPhonons);

void writeStuffToHdf5Temps(
    const std::vector<double> &betaArr,
    const std::vector<double> &wPArr,
    const std::vector<double> &dOcc,
    const std::string &filename
);

void writeStuffToHdf5OnlyPhot(
    const std::vector<double> &gArr,
    const std::vector<double> &dOcc,
    const std::vector<double> &Xpt,
    const std::vector<double> &XptSqr,
    const std::vector<double> &Npt,
    const std::vector<double> &eGS,
    const std::string &filename
);

void writeStuffToHdf52Bands(
    const std::vector<double> &wPhArr,
    const std::vector<double> &dOcc0,
    const std::vector<double> &dOcc1,
    const std::vector<double> &dOccUpDn,
    const std::vector<double> &dOccSigSig,
    const std::vector<double> &n0,
    const std::vector<double> &n1,
    const std::string &filename
);

void writeSpectrumToFile(
    const std::vector<double> &gArr,
    const std::vector<double> &spectrum,
    const std::vector<double> &nc0Arr,
    const std::vector<double> &nd0Arr,
    const std::vector<double> &nc1Arr,
    const std::vector<double> &nd1Arr,
    const std::vector<double> &nBos0Arr,
    const std::vector<double> &nBos1Arr,
    const std::string &filename
);


void writeStuffToHdf52BandsTime(
    const std::vector<double> &times,
    const std::vector<double> &pumpFunction,
    const std::vector<double> &dOcc0,
    const std::vector<double> &dOcc1,
    const std::vector<double> &dOccUpDn,
    const std::vector<double> &dOccSigSig,
    const std::vector<double> &Xph1,
    const std::vector<double> &Xph1Sqr,
    const std::vector<double> &Nph1,
    const std::vector<double> &Xph2,
    const std::vector<double> &Xph2Sqr,
    const std::vector<double> &Nph2,
    const std::vector<double> &N0,
    const std::vector<double> &N1,
    const std::string &filename
);

void readInComplex2DArray(std::vector<std::complex<double>> &readInArray, const std::string &fileName);

#endif //PHONONSQUEEZING_WRITESTUFFTOHDF5_H
