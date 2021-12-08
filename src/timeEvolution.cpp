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

#include "H5Cpp.h"

#define MKL_Complex16 std::complex<double>

#include "mkl.h"

double pumpEnvolope(const double t, const double t0, const double s) {
  return 1. / std::sqrt(2. * PI * s) * std::exp(-0.5 * (t - t0) * (t - t0) / (s * s));
}

void setupOpsOnePh(std::vector<std::complex<double>> &dOcc,
                   std::vector<std::complex<double>> &Xpt,
                   std::vector<std::complex<double>> &XptSqr,
                   std::vector<std::complex<double>> &Npt,
                   std::vector<std::complex<double>> &X1ph,
                   std::vector<std::complex<double>> &X1phSqr,
                   std::vector<std::complex<double>> &N1ph,
                   std::vector<std::complex<double>> &X2ph,
                   std::vector<std::complex<double>> &X2phSqr,
                   std::vector<std::complex<double>> &N2ph
);

void setupOpsTwoPh(std::vector<std::complex<double>> &dOcc,
                   std::vector<std::complex<double>> &Xpt,
                   std::vector<std::complex<double>> &XptSqr,
                   std::vector<std::complex<double>> &Npt,
                   std::vector<std::complex<double>> &X1ph,
                   std::vector<std::complex<double>> &X1phSqr,
                   std::vector<std::complex<double>> &N1ph,
                   std::vector<std::complex<double>> &X2ph,
                   std::vector<std::complex<double>> &X2phSqr,
                   std::vector<std::complex<double>> &N2ph);

void writeStuffToHdf5(
        const std::vector<double> &times,
        const std::vector<double> &pumpFunction,
        const std::vector<double> &dOcc,
        const std::vector<double> &Xpt,
        const std::vector<double> &XptSqr,
        const std::vector<double> &Npt,
        const std::vector<double> &X1ph,
        const std::vector<double> &X1phSqr,
        const std::vector<double> &N1ph,
        const std::vector<double> &X2ph,
        const std::vector<double> &X2phSqr,
        const std::vector<double> &N2ph,
        const std::string &filename,
        const bool twoPhonons);

void calcTimeEvolution(const bool twoPhonons) {

  const ulong dimH = twoPhonons ? dimHTwoPh : dimHOnePh;

  const ulong timeSteps(40ul * 30ul);
  std::vector<double> times(timeSteps, 0.);
  std::vector<double> pumpPreFac(timeSteps, 0.);
  std::vector<double> pumpPreFacOutput(timeSteps, 0.);

  const double t0 = 2. * PI / wDrive * 10.;
  const double s = 2. * PI / wDrive * 2.;

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
    addMatricies(X1ph, X2ph, ODrive);
  } else {
    setupOpsOnePh(dOcc, Xpt, XptSqr, Npt, X1ph, X1phSqr, N1ph, X2ph, X2phSqr, N2ph);
    ODrive = std::vector<std::complex<double>>(Xpt);
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

  std::string filename = "data/timeEvolOnePhonWP" + std::to_string(int(100 * wP)) + "N" + std::to_string(dimPhonon) + ".hdf5";

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
void setupOpsOnePh(std::vector<std::complex<double>> &dOcc,
                   std::vector<std::complex<double>> &Xpt,
                   std::vector<std::complex<double>> &XptSqr,
                   std::vector<std::complex<double>> &Npt,
                   std::vector<std::complex<double>> &X1ph,
                   std::vector<std::complex<double>> &X1phSqr,
                   std::vector<std::complex<double>> &N1ph,
                   std::vector<std::complex<double>> &X2ph,
                   std::vector<std::complex<double>> &X2phSqr,
                   std::vector<std::complex<double>> &N2ph) {

  const ulong dimH = dimHOnePh;

  Xpt = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> A(dimH * dimH, std::complex<double>(
          0., 0.));
  setupAOnePh(A);
  std::vector<std::complex<double>> ADag(A);
  dagger(ADag, dimH);


  addMatricies(ADag, 1. / std::sqrt(wPt), A, 1. / std::sqrt(wPt), Xpt);

  Npt = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));

  std::complex<double> alpha(1., 0.);
  std::complex<double> beta(0., 0.);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              A.data(), dimH,
              A.data(), dimH,
              &beta, Npt.data(), dimH);

  dOcc = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));

  setupDOccOnePh(dOcc);

  XptSqr = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));


  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              Xpt.data(), dimH,
              Xpt.data(), dimH,
              &beta, XptSqr.data(), dimH);

  X1ph = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> B(dimH * dimH, std::complex<double>(
          0., 0.));
  setupB(B);
  std::vector<std::complex<double>> BDag(B);
  dagger(BDag, dimH);

  addMatricies(BDag, 1. / std::sqrt(wPh), B, 1. / std::sqrt(wPh), X1ph);

  N1ph = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              B.data(), dimH,
              B.data(), dimH,
              &beta, N1ph.data(), dimH);

  X1phSqr = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              X1ph.data(), dimH,
              X1ph.data(), dimH,
              &beta, X1phSqr.data(), dimH);

}

void setupOpsTwoPh(std::vector<std::complex<double>> &dOcc,
                   std::vector<std::complex<double>> &Xpt,
                   std::vector<std::complex<double>> &XptSqr,
                   std::vector<std::complex<double>> &Npt,
                   std::vector<std::complex<double>> &X1ph,
                   std::vector<std::complex<double>> &X1phSqr,
                   std::vector<std::complex<double>> &N1ph,
                   std::vector<std::complex<double>> &X2ph,
                   std::vector<std::complex<double>> &X2phSqr,
                   std::vector<std::complex<double>> &N2ph) {

  const ulong dimH = dimHTwoPh;

  Xpt = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> A(dimH * dimH, std::complex<double>(
          0., 0.));
  setupA(A);
  std::vector<std::complex<double>> ADag(A);
  dagger(ADag, dimH);


  addMatricies(ADag, 1. / std::sqrt(wPt), A, 1. / std::sqrt(wPt), Xpt);

  Npt = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));

  std::complex<double> alpha(1., 0.);
  std::complex<double> beta(0., 0.);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              A.data(), dimH,
              A.data(), dimH,
              &beta, Npt.data(), dimH);

  dOcc = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));

  setupDOcc(dOcc);

  XptSqr = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));


  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              Xpt.data(), dimH,
              Xpt.data(), dimH,
              &beta, XptSqr.data(), dimH);

  X1ph = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> B1(dimH * dimH, std::complex<double>(
          0., 0.));
  setupB1(B1);
  std::vector<std::complex<double>> B1Dag(B1);
  dagger(B1Dag, dimH);

  addMatricies(B1Dag, 1. / std::sqrt(wPh), B1, 1. / std::sqrt(wPh), X1ph);

  N1ph = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              B1.data(), dimH,
              B1.data(), dimH,
              &beta, N1ph.data(), dimH);

  X1phSqr = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              X1ph.data(), dimH,
              X1ph.data(), dimH,
              &beta, X1phSqr.data(), dimH);

  X2ph = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> B2(dimH * dimH, std::complex<double>(
          0., 0.));
  setupB2(B2);
  std::vector<std::complex<double>> B2Dag(B2);
  dagger(B2Dag, dimH);

  addMatricies(B2Dag, 1. / std::sqrt(wPh), B2, 1. / std::sqrt(wPh), X2ph);

  N2ph = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              B2.data(), dimH,
              B2.data(), dimH,
              &beta, N2ph.data(), dimH);

  X2phSqr = std::vector<std::complex<double>>(dimH * dimH, std::complex<double>(0., 0.));

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimH, dimH, dimH, &alpha,
              X2ph.data(), dimH,
              X2ph.data(), dimH,
              &beta, X2phSqr.data(), dimH);

}

void writeStuffToHdf5(
        const std::vector<double> &times,
        const std::vector<double> &pumpFunction,
        const std::vector<double> &dOcc,
        const std::vector<double> &Xpt,
        const std::vector<double> &XptSqr,
        const std::vector<double> &Npt,
        const std::vector<double> &X1ph,
        const std::vector<double> &X1phSqr,
        const std::vector<double> &N1ph,
        const std::vector<double> &X2ph,
        const std::vector<double> &X2phSqr,
        const std::vector<double> &N2ph,
        const std::string &filename,
        const bool twoPhonons
) {


  H5::H5File file(filename, H5F_ACC_TRUNC);

  const hsize_t dataSize = times.size();
  H5::DataSpace dataSpace(1, &dataSize);
  H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);

  H5::DataSet datasetTimes = file.createDataSet("times", datatype, dataSpace);
  datasetTimes.write(&times[0], datatype);

  H5::DataSet datasetPump = file.createDataSet("pump", datatype, dataSpace);
  datasetPump.write(&pumpFunction[0], datatype);

  H5::DataSet datasetDOcc = file.createDataSet("dOcc", datatype, dataSpace);
  datasetDOcc.write(&dOcc[0], datatype);

  H5::DataSet datasetXpt = file.createDataSet("Xpt", datatype, dataSpace);
  datasetXpt.write(&Xpt[0], datatype);

  H5::DataSet datasetXptSqr = file.createDataSet("XptSqr", datatype, dataSpace);
  datasetXptSqr.write(&XptSqr[0], datatype);

  H5::DataSet datasetNpt = file.createDataSet("Npt", datatype, dataSpace);
  datasetNpt.write(&Npt[0], datatype);

  H5::DataSet datasetXph = file.createDataSet("X1ph", datatype, dataSpace);
  datasetXph.write(&X1ph[0], datatype);

  H5::DataSet datasetXphSqr = file.createDataSet("X1phSqr", datatype, dataSpace);
  datasetXphSqr.write(&X1phSqr[0], datatype);

  H5::DataSet datasetNph = file.createDataSet("N1ph", datatype, dataSpace);
  datasetNph.write(&N1ph[0], datatype);

  if (twoPhonons) {
    H5::DataSet datasetX2ph = file.createDataSet("X2ph", datatype, dataSpace);
    datasetX2ph.write(&X2ph[0], datatype);

    H5::DataSet datasetX2phSqr = file.createDataSet("X2phSqr", datatype, dataSpace);
    datasetX2phSqr.write(&X2phSqr[0], datatype);

    H5::DataSet datasetN2ph = file.createDataSet("N2ph", datatype, dataSpace);
    datasetN2ph.write(&N2ph[0], datatype);
  }

}
