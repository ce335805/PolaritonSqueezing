#include <iostream>

#include "include/twoPhonons/evalDoccInGS.h"
#include "include/onePhonon/evalDoccInGSOnePh.h"

#include <chrono>


int main() {

  std::cout << "Hello World! - let's squeeze some phonons" << std::endl;

  //evalDoccInGSTwoPh();
  evalDoccInGSOnePh();

  return 0;
}

