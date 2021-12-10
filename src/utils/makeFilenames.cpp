#include <vector>
#include <complex>
#include <string>

#include "globals.h"


std::string timeEvolName(const bool twoPhonons){

  std::string fileName ("data/tEvol");

  std::string twoPhonSpeci;
  if(twoPhonons){
    twoPhonSpeci = "2Ph";
  } else {
    twoPhonSpeci = "1Ph";
  }

  std::string gPhStr ("GPH" + std::to_string(int(std::round(100 * std::abs(gPh)))));
  std::string wPStr ("WP" + std::to_string(int(std::round(100 * wP))));
  std::string wDStr ("WD" + std::to_string(int(std::round(100 * wDrive))));
  std::string fDStr ("FD" + std::to_string(int(std::round(100 * fDrive))));
  std::string nBStr ("NB" + std::to_string(int(dimPhonon)));

  fileName = fileName + twoPhonSpeci + gPhStr + wPStr + wDStr + fDStr + nBStr + ".hdf5";

  return fileName;

}




