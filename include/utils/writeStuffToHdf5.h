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

#endif //PHONONSQUEEZING_WRITESTUFFTOHDF5_H
