#include <iostream>

#include "globals.h"
#include "matrixOperations.h"
#include "include/onePhonon/setUpGlobalHamiltonianOnePh.h"
#include "include/onePhonon/setupBasicOperatorsOnePh.h"
#include "include/onePhonon/evalExpectationOnePh.h"

#include <chrono>

void calcGSOnePh(std::vector<std::complex<double>> &gs, std::vector<std::complex<double>> globalH){

  std::cout << "Starting diagonalization ..." << '\n';

  auto start = std::chrono::high_resolution_clock::now();

  std::vector<double> spectrum = diagonalize(globalH, dimHOnePh, 'V');

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  std::cout << "Diagonalization took " << duration.count() << "ms" << '\n';

  std::cout << '\n';

  for(ulong ind = 0ul; ind < dimHOnePh; ++ind){
    gs[ind] = globalH[ind * dimHOnePh];
  }

}