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


void writeParameterToHDF5File(H5::H5File file, const std::string &parameterName, const double parameterValue){

  const hsize_t dataSize = 1ul;
  H5::DataSpace dataSpace(1, &dataSize);
  H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);

  std::vector<double> data (1ul, parameterValue);
  H5::DataSet datasetParameter = file.createDataSet(parameterName, datatype, dataSpace);
  datasetParameter.write(&data[0], datatype);

}

void writeAllPrms(H5::H5File file){

  writeParameterToHDF5File(file, "wPh", wPh);
  writeParameterToHDF5File(file, "wPt", wPt);
  writeParameterToHDF5File(file, "tHop", tHop);
  writeParameterToHDF5File(file, "U", U);
  writeParameterToHDF5File(file, "gPh", gPh);
  writeParameterToHDF5File(file, "wP", wP);

  writeParameterToHDF5File(file, "wDrive", wDrive);
  writeParameterToHDF5File(file, "fDrive", fDrive);
  writeParameterToHDF5File(file, "timePointsPerDrivingPeriod", timePointsPerDrivingPeriod);

  writeParameterToHDF5File(file, "dimPhonon", dimPhonon);
  writeParameterToHDF5File(file, "dimPhoton", dimPhoton);


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


  std::cout << "filename: " << filename << '\n';

  H5::H5File file(filename, H5F_ACC_TRUNC);

  writeAllPrms(file);

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

void writeStuffToHdf5Temps(
        const std::vector<double> &betaArr,
        const std::vector<double> &wPArr,
        const std::vector<double> &dOcc,
        const std::string &filename
) {


  std::cout << "filename: " << filename << '\n';

  H5::H5File file(filename, H5F_ACC_TRUNC);

  writeAllPrms(file);

  H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);

  const hsize_t dataSizeBeta = betaArr.size();
  H5::DataSpace dataSpaceBeta(1, &dataSizeBeta);
  H5::DataSet datasetBetas = file.createDataSet("betas", datatype, dataSpaceBeta);
  datasetBetas.write(&betaArr[0], datatype);

  const hsize_t dataSizeWP = wPArr.size();
  H5::DataSpace dataSpaceWP(1, &dataSizeWP);
  H5::DataSet datasetWP = file.createDataSet("wPs", datatype, dataSpaceWP);
  datasetWP.write(&wPArr[0], datatype);

  const hsize_t dataSizeDOcc = dOcc.size();
  H5::DataSpace dataSpaceDOcc(1, &dataSizeDOcc);
  H5::DataSet datasetDOccc = file.createDataSet("dOccs", datatype, dataSpaceDOcc);
  datasetDOccc.write(&dOcc[0], datatype);
}