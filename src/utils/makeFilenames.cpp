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
  std::string tsStr ("TS" + std::to_string(int(timePointsPerDrivingPeriod)));
  std::string wPhStr ("WPh" + std::to_string(int(std::round(1000 * wPh))));

  fileName = fileName + twoPhonSpeci + gPhStr + wPStr + wDStr + fDStr + nBStr + tsStr + wPhStr + ".hdf5";

  return fileName;

}


std::string gsPropName(const bool twoPhonons){

  std::string fileName ("data/gsProp");

  std::string twoPhonSpeci;
  if(twoPhonons){
    twoPhonSpeci = "2Ph";
  } else {
    twoPhonSpeci = "1Ph";
  }

  std::string gPhStr ("GPH" + std::to_string(int(std::round(100 * std::abs(gPh)))));
  std::string nBStr ("NB" + std::to_string(int(dimPhonon)));

  fileName = fileName + twoPhonSpeci + gPhStr + nBStr + ".hdf5";

  return fileName;

}

std::string gsPropNameTemp(const bool twoPhonons, const bool quadratic){

  std::string fileName ("data/gsProp");

  std::string twoPhonSpeci;
  if(twoPhonons){
    twoPhonSpeci = "2Ph";
  } else {
    twoPhonSpeci = "1Ph";
  }

  std::string quadraticSpeci;
  if(quadratic){
    quadraticSpeci = "Quad";
  } else {
    quadraticSpeci = "Lin";
  }

  std::string gPhStr ("GPH" + std::to_string(int(std::round(100 * std::abs(gPh)))));
  std::string nBStr ("NB" + std::to_string(int(dimPhonon)));

  fileName = fileName + twoPhonSpeci + quadraticSpeci + gPhStr + nBStr + "Temp.hdf5";

  return fileName;

}

std::string gsPropNameOnlyPhot(){
  
  std::string fileName ("data/gsPropOnlyPhot");
  
  std::string nBStr ("NB" + std::to_string(int(dimPhoton)));
  
  fileName = fileName + nBStr + ".hdf5";
  
  return fileName;
  
}


