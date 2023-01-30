#include <vector>
#include <complex>
#include <iostream>
#include <string>

#include "globals.h"
#include "setupBasicOperatorsOnePh.h"
#include "matrixOperations.h"
#include "setUpGlobalHamiltonianOnePh.h"
#include "calcGS.h"
#include "timeStep.h"
#include "setupBasicOperators.h"
#include "setUpGlobalHamiltonian.h"
#include "evalExpectation.h"
#include "writeStuffToHdf5.h"
#include "setupOps.h"
#include "makeFilenames.h"


#define MKL_Complex16 std::complex<double>

#include "mkl.h"

double pumpEnvolope(const double t, const double t0, const double s) {
  return 1. / std::sqrt(2. * PI * s * s) * std::exp(-0.5 * (t - t0) * (t - t0) / (s * s));
}


void calcTimeEvolution(const bool twoPhonons) {

  const ulong dimH = twoPhonons ? dimHTwoPh : dimHOnePh;

  const ulong timeSteps(timePointsPerDrivingPeriod * 40ul);
  std::vector<double> times(timeSteps, 0.);
  std::vector<double> pumpPreFac(timeSteps, 0.);
  std::vector<double> pumpPreFacOutput(timeSteps, 0.);

  const double t0 = 2. * PI / wPt * 8.;
  const double s = 2. * PI / wPt * 2.;

  for (ulong ind = 0; ind < timeSteps; ++ind) {
    times[ind] = double(ind) * dt;
    pumpPreFac[ind] = fDrive * pumpEnvolope(times[ind], t0, s);
    pumpPreFacOutput[ind] = fDrive * pumpEnvolope(times[ind], t0, s) * std::sin(wDrive * times[ind]);
  }
  std::vector<std::complex<double>> gs(dimH, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> H;

  if (twoPhonons) {
    setupGlobalH(H);
  } else {
    setupGlobalHOnePh(H);
  }

  calcGS(gs, H, dimH);

  std::vector<std::complex<double>> dOcc;
  std::vector<std::complex<double>> Xpt;
  std::vector<std::complex<double>> XptSqr;
  std::vector<std::complex<double>> Npt;
  std::vector<std::complex<double>> X1ph;
  std::vector<std::complex<double>> X1phSqr;
  std::vector<std::complex<double>> N1ph;
  std::vector<std::complex<double>> X2ph;
  std::vector<std::complex<double>> X2phSqr;
  std::vector<std::complex<double>> N2ph;
  std::vector<std::complex<double>> ODrive(dimH * dimH, std::complex<double>(0., 0.));

  std::vector<double> dOccExpectation(timeSteps, 0.);
  std::vector<double> XptExpectation(timeSteps, 0.);
  std::vector<double> XptSqrExpectation(timeSteps, 0.);
  std::vector<double> NptExpectation(timeSteps, 0.);
  std::vector<double> X1phExpectation(timeSteps, 0.);
  std::vector<double> X1phSqrExpectation(timeSteps, 0.);
  std::vector<double> N1phExpectation(timeSteps, 0.);
  std::vector<double> X2phExpectation;
  std::vector<double> X2phSqrExpectation;
  std::vector<double> N2phExpectation;

  if (twoPhonons) {
    X2phExpectation = std::vector<double>(timeSteps, 0.);
    X2phSqrExpectation = std::vector<double>(timeSteps, 0.);
    N2phExpectation = std::vector<double>(timeSteps, 0.);
  }

  if (twoPhonons) {
    setupOpsTwoPh(dOcc, Xpt, XptSqr, Npt, X1ph, X1phSqr, N1ph, X2ph, X2phSqr, N2ph);

    addMatricies(X1ph, 1. / std::sqrt(2.), X2ph, 1. / std::sqrt(2.), ODrive);
    //ODrive = std::vector<std::complex<double>>(Xpt);
  } else {
    setupOpsOnePh(dOcc, Xpt, XptSqr, Npt, X1ph, X1phSqr, N1ph, X2ph, X2phSqr, N2ph);
    //if(std::abs(wP) < 1e-12){
    //ODrive = std::vector<std::complex<double>>(X1ph);
    //} else {
    ODrive = std::vector<std::complex<double>>(Xpt);
    //}
  }
  for (ulong timeStep = 0ul; timeStep < timeSteps; ++timeStep) {
    calcTimeStep(times[timeStep], pumpPreFac[timeStep], H, ODrive, gs, dimH);

    dOccExpectation[timeStep] = evalExpectation(dOcc, gs, dimH);
    XptExpectation[timeStep] = evalExpectation(Xpt, gs, dimH);
    XptSqrExpectation[timeStep] = evalExpectation(XptSqr, gs, dimH);
    NptExpectation[timeStep] = evalExpectation(Npt, gs, dimH);
    X1phExpectation[timeStep] = evalExpectation(X1ph, gs, dimH);
    X1phSqrExpectation[timeStep] = evalExpectation(X1phSqr, gs, dimH);
    N1phExpectation[timeStep] = evalExpectation(N1ph, gs, dimH);
    if (twoPhonons) {
      X2phExpectation[timeStep] = evalExpectation(X2ph, gs, dimH);
      X2phSqrExpectation[timeStep] = evalExpectation(X2phSqr, gs, dimH);
      N2phExpectation[timeStep] = evalExpectation(N2ph, gs, dimH);
    }
  }

  std::string filename = timeEvolName(twoPhonons);

  writeStuffToHdf5(times,
                   pumpPreFacOutput,
                   dOccExpectation,
                   XptExpectation,
                   XptSqrExpectation,
                   NptExpectation,
                   X1phExpectation,
                   X1phSqrExpectation,
                   N1phExpectation,
                   X2phExpectation,
                   X2phSqrExpectation,
                   N2phExpectation,
                   filename,
                   twoPhonons);
}

