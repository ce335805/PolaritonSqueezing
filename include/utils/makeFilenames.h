#include <string>

#ifndef PHONONSQUEEZING_MAKEFILENAMES_H
#define PHONONSQUEEZING_MAKEFILENAMES_H

std::string timeEvolName(const bool twoPhonons);

std::string gsPropName(const bool twoPhonons);

std::string gsPropNameTemp(const bool twoPhonons, const bool quadratic);

std::string gsPropNameOnlyPhot();

std::string gsPropName2Bands(const double gSqrOverOmega);

std::string spectrumName();

std::string timeEvolName2Bands();

#endif //PHONONSQUEEZING_MAKEFILENAMES_H
