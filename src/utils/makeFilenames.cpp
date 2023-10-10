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

std::string gsPropName2Bands(const double gSqrOverOmega){
  
  std::string fileName ("data/gsProp2Bands");
  
  std::string U0Str ("Ua" + std::to_string(int(10 * U)));
  std::string U1Str ("Ub" + std::to_string(int(10 * U1)));
  std::string UUDStr ("Uud" + std::to_string(int(10 * uUpDn)));
  std::string USSStr ("Uss" + std::to_string(int(10 * uSigSig)));
  std::string eps0Str ("epsA" + std::to_string(int(10 * eps0)));
  std::string eps1Str ("epsB" + std::to_string(int(10 * eps1)));
  std::string gEStr ("gE" + std::to_string(int(std::round(100 * gSqrOverOmega))));
  
  
  fileName = fileName + U0Str + U1Str + UUDStr + USSStr + eps0Str + eps1Str + gEStr + ".hdf5";
  
  return fileName;
  
}

std::string timeEvolName2Bands(){
  
  std::string fileName ("data/gsProp2Bands");
  
  std::string U0Str ("Ua" + std::to_string(int(10 * U)));
  std::string U1Str ("Ub" + std::to_string(int(10 * U1)));
  std::string UUDStr ("Uud" + std::to_string(int(10 * uUpDn)));
  std::string USSStr ("Uss" + std::to_string(int(10 * uSigSig)));
  std::string eps0Str ("epsA" + std::to_string(int(10 * eps0)));
  std::string eps1Str ("epsB" + std::to_string(int(10 * eps1)));
  std::string gEStr ("gE" + std::to_string(int(10 * gE)));
  
  std::string wDStr ("WD" + std::to_string(int(std::round(100 * wDrive))));
  std::string fDStr ("FD" + std::to_string(int(std::round(100 * fDrive))));
  std::string nBStr ("NB" + std::to_string(int(dimPhonon)));
  std::string tsStr ("TS" + std::to_string(int(timePointsPerDrivingPeriod)));
  std::string wPhStr ("WPh" + std::to_string(int(std::round(1000 * wPh))));
  
  
  fileName = fileName + U0Str + U1Str + UUDStr + USSStr + eps0Str + eps1Str + gEStr + wDStr + fDStr + nBStr + tsStr + wPhStr + ".hdf5";
  
  return fileName;
  
}

std::string spectrumName(){
  
  std::string fileName ("data/spectrum2Bands");
  
  std::string tHopStr ("tHop" + std::to_string(int(std::round(100 * tHop))));
  std::string eps0Str ("epsA" + std::to_string(int(10 * eps0)));
  std::string eps1Str ("epsB" + std::to_string(int(10 * eps1)));
  std::string omStr ("Om" + std::to_string(int(std::round(10 * wPt))));
  
  
  fileName = fileName + tHopStr + eps0Str + eps1Str + omStr + ".hdf5";
  
  return fileName;
  
}
