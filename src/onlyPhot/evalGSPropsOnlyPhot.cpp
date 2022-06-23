#include <iostream>

#include "globals.h"
#include "matrixOperations.h"
#include "evalExpectation.h"
#include "setUpGlobalHamiltonianOnlyPhot.h"
#include "setupBasicOperatorsOnlyPhot.h"
#include "setupOps.h"
#include "calcGS.h"
#include "writeStuffToHdf5.h"
#include "makeFilenames.h"

#include <chrono>

//double gE;

void evalDoccInGSOnlyPhot(){
  std::cout << "Using parameters: " << '\n';
  std::cout << "w-Photon = " << wPt << '\n';
  std::cout << "electronic U = " << U << '\n';
  std::cout << "gE = " << gE << '\n';

  std::cout << "\n";

  std::vector<std::complex<double>> H;

  std::cout << "Starting setup of Hamiltonian ..." << '\n';

  auto start = std::chrono::high_resolution_clock::now();

  setupGlobalHOnlyPhot(H);

  std::vector<std::complex<double>> dOcc;
  setupDOccOnlyPhot(dOcc);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  std::cout << "Setup of H took " << duration.count() << "ms" << '\n';

  std::cout << "Starting diagonalization ..." << '\n';

  start = std::chrono::high_resolution_clock::now();


  std::vector<double> spectrum = diagonalize(H, dimHOnlyPhot, 'V');

  std::cout << "GS energy = " << spectrum[0] << '\n';
  std::cout << "GS + 1 energy = " << spectrum[1] << '\n';
  std::cout << "GS + 2 energy = " << spectrum[2] << '\n';
  std::cout << "GS + 3 energy = " << spectrum[3] << '\n';

  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  std::cout << "Diagonalization took " << duration.count() << "ms" << '\n';

  std::cout << '\n';

  std::vector<std::complex<double>> gs (dimHOnlyPhot, std::complex<double> (0., 0.));
  for(ulong ind = 0ul; ind < dimHOnlyPhot; ++ind){
    gs[ind] = H[ind * dimHOnlyPhot];
  }

  double dOccExpec = evalExpectation(dOcc, gs, dimHOnlyPhot);

  std::cout << "<dOcc> = " << dOccExpec + 0.5 << '\n';
}

void evalGSPropsAsOfG() {
  
  const ulong gSteps(101ul);
  std::vector<double> gArr(gSteps, 0.);
  
  for (ulong ind = 0ul; ind < gSteps; ++ind) {
    gArr[ind] = double(ind) / 100.;
  }
  
  std::vector<std::complex<double>> dOcc;
  std::vector<std::complex<double>> Xpt;
  std::vector<std::complex<double>> XptSqr;
  std::vector<std::complex<double>> Npt;
  
  setupOpsOnlyPhot(dOcc, Xpt, XptSqr, Npt);
  
  std::vector<double> dOccExpectation(gSteps, 0.);
  std::vector<double> XptExpectation(gSteps, 0.);
  std::vector<double> XptSqrExpectation(gSteps, 0.);
  std::vector<double> NptExpectation(gSteps, 0.);
  std::vector<double> eGSExpectation(gSteps, 0.);
  
  std::vector<std::complex<double>> gs(dimHOnlyPhot, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> H;
  
  for (ulong gStep = 0ul; gStep < gSteps; ++gStep) {
    
    ////////////////// set g //////////////////////
    //gE = gArr[gStep];
    std::cout << "gE = " << gE << '\n';
    setupGlobalHOnlyPhot(H);
  
    eGSExpectation[gStep] = calcGSWithE(gs, H, dimHOnlyPhot);
    calcGS(gs, H, dimHOnlyPhot);
    dOccExpectation[gStep] = evalExpectation(dOcc, gs, dimHOnlyPhot);
    XptExpectation[gStep] = evalExpectation(Xpt, gs, dimHOnlyPhot);
    XptSqrExpectation[gStep] = evalExpectation(XptSqr, gs, dimHOnlyPhot);
    NptExpectation[gStep] = evalExpectation(Npt, gs, dimHOnlyPhot);
  }
  
  std::string filename;
  filename = gsPropNameOnlyPhot();
  
  writeStuffToHdf5OnlyPhot(gArr,
                   dOccExpectation,
                   XptExpectation,
                   XptSqrExpectation,
                   NptExpectation,
                   eGSExpectation,
                   filename);
}
