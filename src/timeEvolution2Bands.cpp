#include <vector>
#include <complex>
#include <iostream>
#include <string>

#include "globals.h"
#include "matrixOperations.h"
#include "calcGS.h"
#include "timeStep.h"
#include "evalExpectation.h"
#include "writeStuffToHdf5.h"
#include "setupOps.h"
#include "makeFilenames.h"
#include "setupHGlobal2Bands.h"
#include "timeEvolution.h"
#include "timeEvolution2Bands.h"
#include <chrono>

#define MKL_Complex16 std::complex<double>

#include "mkl.h"

double wDrive;
double dt;


void calcTimeEvolutionAsOfWD(){
  
  const ulong wSteps(40ul);
  std::vector<double> wArr(wSteps, 0.);
  
  for (ulong ind = 0ul; ind < wSteps; ++ind) {
    wArr[ind] = double(ind + 1ul) / 4. + 6.;
  }
  
  auto start = std::chrono::high_resolution_clock::now();
  
  
  for (ulong wStep = 0ul; wStep < wSteps; ++wStep) {
  
    auto startLoop = std::chrono::high_resolution_clock::now();
  
    ////////////////// set wPh //////////////////////
    wDrive = wArr[wStep];
    dt = 2. * PI / wDrive / timePointsPerDrivingPeriod;
    std::cout << "wDrive = " << wDrive << '\n';
    calcTimeEvolution2Bands();
  
    auto stopLoop = std::chrono::high_resolution_clock::now();
    auto durationLoop = std::chrono::duration_cast<std::chrono::milliseconds>(stopLoop - startLoop);
  
    std::cout << "This freq. took " << durationLoop.count() << "ms" << '\n';
  }
  
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  
  std::cout << "Entire time-evolution took " << duration.count() << "ms" << '\n';
}

void calcTimeEvolution2Bands() {
  
  const ulong dimH = dimHOnePh;
  
  const ulong timeSteps(timePointsPerDrivingPeriod * 20ul);
  std::vector<double> times(timeSteps, 0.);
  std::vector<double> pumpPreFac(timeSteps, 0.);
  std::vector<double> pumpPreFacOutput(timeSteps, 0.);
  
  const double t0 = 2. * PI / wDrive * 8.;
  const double s = 2. * PI / wDrive * 2.;
  
  for (ulong ind = 0; ind < timeSteps; ++ind) {
    times[ind] = double(ind) * dt;
    pumpPreFac[ind] = fDrive * pumpEnvolope(times[ind], t0, s);
    pumpPreFacOutput[ind] = fDrive * pumpEnvolope(times[ind], t0, s) * std::sin(wDrive * times[ind]);
  }
  std::vector<std::complex<double>> gs(dimH, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> H;
  
  setupGlobalH2Bands(H);
  
  calcGS(gs, H, dimH);
  
  std::vector<std::complex<double>> dOcc0;
  std::vector<std::complex<double>> dOcc1;
  std::vector<std::complex<double>> dOccUpDn;
  std::vector<std::complex<double>> dOccSigSig;
  std::vector<std::complex<double>> n0;
  std::vector<std::complex<double>> n1;
  std::vector<std::complex<double>> Xph1;
  std::vector<std::complex<double>> Xph1Sqr;
  std::vector<std::complex<double>> Npt1;
  std::vector<std::complex<double>> Xph2;
  std::vector<std::complex<double>> Xph2Sqr;
  std::vector<std::complex<double>> Nph2;
  std::vector<std::complex<double>> ODrive(dimH * dimH, std::complex<double>(0., 0.));
  
  std::vector<double> dOcc0Expectation(timeSteps, 0.);
  std::vector<double> dOcc1Expectation(timeSteps, 0.);
  std::vector<double> dOccUpDnExpectation(timeSteps, 0.);
  std::vector<double> dOccSigSigExpectation(timeSteps, 0.);
  std::vector<double> XptExpectation(timeSteps, 0.);
  std::vector<double> XptSqrExpectation(timeSteps, 0.);
  std::vector<double> NptExpectation(timeSteps, 0.);
  std::vector<double> X1phExpectation(timeSteps, 0.);
  std::vector<double> X1phSqrExpectation(timeSteps, 0.);
  std::vector<double> N1phExpectation(timeSteps, 0.);
  std::vector<double> N0Expectation(timeSteps, 0.);
  std::vector<double> N1Expectation(timeSteps, 0.);
  
  
  setupOps2BandsDrive(dOcc0, dOcc1, dOccUpDn, dOccSigSig, n0, n1, Xph1, Xph1Sqr, Npt1, Xph2, Xph2Sqr, Nph2);
  
  
  addMatricies(Xph1, 1. / std::sqrt(2.), Xph2, 1. / std::sqrt(2.), ODrive);
  
  for (ulong timeStep = 0ul; timeStep < timeSteps; ++timeStep) {
    calcTimeStep(times[timeStep], pumpPreFac[timeStep], H, ODrive, gs, dimH);
    
    dOcc0Expectation[timeStep] = evalExpectation(dOcc0, gs, dimH);
    dOcc1Expectation[timeStep] = evalExpectation(dOcc1, gs, dimH);
    dOccUpDnExpectation[timeStep] = evalExpectation(dOccUpDn, gs, dimH);
    dOccSigSigExpectation[timeStep] = evalExpectation(dOccSigSig, gs, dimH);
    XptExpectation[timeStep] = evalExpectation(Xph1, gs, dimH);
    XptSqrExpectation[timeStep] = evalExpectation(Xph1Sqr, gs, dimH);
    NptExpectation[timeStep] = evalExpectation(Npt1, gs, dimH);
    X1phExpectation[timeStep] = evalExpectation(Xph2, gs, dimH);
    X1phSqrExpectation[timeStep] = evalExpectation(Xph2Sqr, gs, dimH);
    N1phExpectation[timeStep] = evalExpectation(Nph2, gs, dimH);
    N0Expectation[timeStep] = evalExpectation(n0, gs, dimH);
    N1Expectation[timeStep] = evalExpectation(n1, gs, dimH);
  }
  
  std::string filename = timeEvolName2Bands();
  
  writeStuffToHdf52BandsTime(times,
                             pumpPreFacOutput,
                             dOcc0Expectation,
                             dOcc1Expectation,
                             dOccUpDnExpectation,
                             dOccSigSigExpectation,
                             XptExpectation,
                             XptSqrExpectation,
                             NptExpectation,
                             X1phExpectation,
                             X1phSqrExpectation,
                             N1phExpectation,
                             N0Expectation,
                             N1Expectation,
                             filename);
}

