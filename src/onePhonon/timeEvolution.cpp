#include <vector>
#include <complex>
#include <iostream>

#include "globals.h"
#include "setupBasicOperatorsOnePh.h"
#include "matrixOperations.h"
#include "setUpGlobalHamiltonianOnePh.h"
#include "calcGSOnePh.h"
#include "timeStep.h"
#include "evalExpectationOnePh.h"

#include "H5Cpp.h"

#define MKL_Complex16 std::complex<double>

#include "mkl.h"

double pumpEnvolope(const double t, const double t0, const double s) {
  return 1. / std::sqrt(2. * PI * s) * std::exp(-0.5 * (t - t0) * (t - t0) / (s * s));
}

void setupOps(std::vector<std::complex<double>> &dOcc,
              std::vector<std::complex<double>> &Xpt,
              std::vector<std::complex<double>> &XptSqr,
              std::vector<std::complex<double>> &Npt,
              std::vector<std::complex<double>> &Xph,
              std::vector<std::complex<double>> &XphSqr,
              std::vector<std::complex<double>> &Nph);

void writeStuffToHdf5(
        const std::vector<double> &times,
        const std::vector<double> &pumpFunction,
        const std::vector<double> &dOcc,
        const std::vector<double> &Xpt,
        const std::vector<double> &XptSqr,
        const std::vector<double> &Npt,
        const std::vector<double> &Xph,
        const std::vector<double> &XphSqr,
        const std::vector<double> &Nph,
        const std::string &filename);

void calcTimeEvolution() {

  const ulong timeSteps(100ul * 20ul);
  std::vector<double> times(timeSteps, 0.);
  std::vector<double> pumpPreFac(timeSteps, 0.);
  std::vector<double> pumpPreFacOutput(timeSteps, 0.);

  const double t0 = 2. * PI / wDrive * 15.;
  const double s = 2. * PI / wDrive * 3.;

  for (ulong ind = 0; ind < timeSteps; ++ind) {
    times[ind] = double(ind) * dt;
    pumpPreFac[ind] = fDrive * pumpEnvolope(times[ind], t0, s);
    pumpPreFacOutput[ind] = fDrive * pumpEnvolope(times[ind], t0, s) * std::sin(wDrive * times[ind]);
  }
  std::vector<std::complex<double>> gs(dimHOnePh, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> H;
  setupGlobalHOnePh(H);

  calcGSOnePh(gs, H);


  std::vector<std::complex<double>> dOcc;
  std::vector<std::complex<double>> Xpt;
  std::vector<std::complex<double>> XptSqr;
  std::vector<std::complex<double>> Npt;
  std::vector<std::complex<double>> Xph;
  std::vector<std::complex<double>> XphSqr;
  std::vector<std::complex<double>> Nph;

  std::vector<double> dOccExpectation(timeSteps, 0.);
  std::vector<double> XptExpectation(timeSteps, 0.);
  std::vector<double> XptSqrExpectation(timeSteps, 0.);
  std::vector<double> NptExpectation(timeSteps, 0.);
  std::vector<double> XphExpectation(timeSteps, 0.);
  std::vector<double> XphSqrExpectation(timeSteps, 0.);
  std::vector<double> NphExpectation(timeSteps, 0.);

  setupOps(dOcc, Xpt, XptSqr, Npt, Xph, XphSqr, Nph);

  for (ulong timeStep = 0ul; timeStep < timeSteps; ++timeStep) {
    calcTimeStep(times[timeStep], pumpPreFac[timeStep], H, Xph, gs);

    dOccExpectation[timeStep] = evalExpectationOnePh(dOcc, gs);
    XptExpectation[timeStep] = evalExpectationOnePh(Xpt, gs);
    XptSqrExpectation[timeStep] = evalExpectationOnePh(XptSqr, gs);
    NptExpectation[timeStep] = evalExpectationOnePh(Npt, gs);
    XphExpectation[timeStep] = evalExpectationOnePh(Xph, gs);
    XphSqrExpectation[timeStep] = evalExpectationOnePh(XphSqr, gs);
    NphExpectation[timeStep] = evalExpectationOnePh(Nph, gs);
  }

  writeStuffToHdf5(times,
                   pumpPreFacOutput,
                   dOccExpectation,
                   XptExpectation,
                   XptSqrExpectation,
                   NptExpectation,
                   XphExpectation,
                   XphSqrExpectation,
                   NphExpectation,
                   "data/timeEvolutionResultsClassicalDrive.hdf5");
}

void setupOps(std::vector<std::complex<double>> &dOcc,
              std::vector<std::complex<double>> &Xpt,
              std::vector<std::complex<double>> &XptSqr,
              std::vector<std::complex<double>> &Npt,
              std::vector<std::complex<double>> &Xph,
              std::vector<std::complex<double>> &XphSqr,
              std::vector<std::complex<double>> &Nph) {

  Xpt = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  std::vector<std::complex<double>> A(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  setupAOnePh(A);
  std::vector<std::complex<double>> ADag(A);
  dagger(ADag, dimHOnePh);



  addMatricies(ADag, 1. / std::sqrt(wPt), A, 1. / std::sqrt(wPt), Xpt);

  Npt = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  std::complex<double> alpha(1., 0.);
  std::complex<double> beta(0., 0.);

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              A.data(), dimHOnePh,
              A.data(), dimHOnePh,
              &beta, Npt.data(), dimHOnePh);

  dOcc = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  setupDOccOnePh(dOcc);

  XptSqr = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));


  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              Xpt.data(), dimHOnePh,
              Xpt.data(), dimHOnePh,
              &beta, XptSqr.data(), dimHOnePh);

  Xph = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  std::vector<std::complex<double>> B(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));
  setupB(B);
  std::vector<std::complex<double>> BDag(B);
  dagger(BDag, dimHOnePh);

  addMatricies(BDag, 1. / std::sqrt(wPh), B, 1. / std::sqrt(wPh), Xph);

  Nph = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              B.data(), dimHOnePh,
              B.data(), dimHOnePh,
              &beta, Nph.data(), dimHOnePh);

  XphSqr = std::vector<std::complex<double>>(dimHOnePh * dimHOnePh, std::complex<double>(0., 0.));

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dimHOnePh, dimHOnePh, dimHOnePh, &alpha,
              Xph.data(), dimHOnePh,
              Xph.data(), dimHOnePh,
              &beta, XphSqr.data(), dimHOnePh);

}

void writeStuffToHdf5(
        const std::vector<double> &times,
        const std::vector<double> &pumpFunction,
        const std::vector<double> &dOcc,
        const std::vector<double> &Xpt,
        const std::vector<double> &XptSqr,
        const std::vector<double> &Npt,
        const std::vector<double> &Xph,
        const std::vector<double> &XphSqr,
        const std::vector<double> &Nph,
        const std::string &filename
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

  H5::DataSet datasetXph = file.createDataSet("Xph", datatype, dataSpace);
  datasetXph.write(&Xph[0], datatype);

  H5::DataSet datasetXphSqr = file.createDataSet("XphSqr", datatype, dataSpace);
  datasetXphSqr.write(&XphSqr[0], datatype);

  H5::DataSet datasetNph = file.createDataSet("Nph", datatype, dataSpace);
  datasetNph.write(&Nph[0], datatype);

}
