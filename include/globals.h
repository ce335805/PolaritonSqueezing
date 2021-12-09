#ifndef PHONONSQUEEZING_GLOBALS_H
#define PHONONSQUEEZING_GLOBALS_H

#include <complex>

constexpr std::complex<double> II = std::complex<double>(0., 1.);
constexpr double PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640;

constexpr double tHop(0.25);
constexpr double wPt(0.1);
constexpr double wPh(0.2);
constexpr double U(1.25);
constexpr double gPh(0.05);
constexpr double wP(0.1);
//extern double wP;

constexpr double wDrive(0.2);
constexpr double fDrive(2.);
constexpr double dt(2. * PI / wDrive / 40.);

constexpr ulong dimPhonon(4ul);
constexpr ulong dimPhoton(4ul);
constexpr ulong dimHTwoPh(4ul * dimPhonon * dimPhonon * dimPhoton);
constexpr ulong dimHOnePh(4ul * dimPhonon * dimPhoton);

#endif //PHONONSQUEEZING_GLOBALS_H
