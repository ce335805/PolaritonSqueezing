#include <vector>
#include <complex>
#include <iostream>
#include <string>
#include <cassert>

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

void writeStuffToHdf5OnlyPhot(
    const std::vector<double> &gArr,
    const std::vector<double> &dOcc,
    const std::vector<double> &Xpt,
    const std::vector<double> &XptSqr,
    const std::vector<double> &Npt,
    const std::vector<double> &eGS,
    const std::string &filename
) {
  std::cout << "filename: " << filename << '\n';
  
  H5::H5File file(filename, H5F_ACC_TRUNC);
  
  const hsize_t dataSize = gArr.size();
  H5::DataSpace dataSpace(1, &dataSize);
  H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);
  
  H5::DataSet datasetgArr = file.createDataSet("g", datatype, dataSpace);
  datasetgArr.write(&gArr[0], datatype);
  
  H5::DataSet datasetDOcc = file.createDataSet("dOcc", datatype, dataSpace);
  datasetDOcc.write(&dOcc[0], datatype);
  
  H5::DataSet datasetXpt = file.createDataSet("Xpt", datatype, dataSpace);
  datasetXpt.write(&Xpt[0], datatype);
  
  H5::DataSet datasetXptSqr = file.createDataSet("XptSqr", datatype, dataSpace);
  datasetXptSqr.write(&XptSqr[0], datatype);
  
  H5::DataSet datasetNpt = file.createDataSet("Npt", datatype, dataSpace);
  datasetNpt.write(&Npt[0], datatype);
  
  H5::DataSet dataseteGS = file.createDataSet("eGS", datatype, dataSpace);
  dataseteGS.write(&eGS[0], datatype);
  
}

void writeStuffToHdf52Bands(
    const std::vector<double> &wPhArr,
    const std::vector<double> &dOcc0,
    const std::vector<double> &dOcc1,
    const std::vector<double> &dOccUpDn,
    const std::vector<double> &dOccSigSig,
    const std::vector<double> &n0,
    const std::vector<double> &n1,
    const std::string &filename
) {
  std::cout << "filename: " << filename << '\n';
  
  H5::H5File file(filename, H5F_ACC_TRUNC);
  
  const hsize_t dataSize = wPhArr.size();
  H5::DataSpace dataSpace(1, &dataSize);
  H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);
  
  H5::DataSet datasetgArr = file.createDataSet("wPh", datatype, dataSpace);
  datasetgArr.write(&wPhArr[0], datatype);
  
  H5::DataSet datasetDOcc0 = file.createDataSet("dOcc0", datatype, dataSpace);
  datasetDOcc0.write(&dOcc0[0], datatype);
  
  H5::DataSet datasetDOcc1 = file.createDataSet("dOcc1", datatype, dataSpace);
  datasetDOcc1.write(&dOcc1[0], datatype);
  
  H5::DataSet datasetDOccUpDn = file.createDataSet("dOccUpDn", datatype, dataSpace);
  datasetDOccUpDn.write(&dOccUpDn[0], datatype);
  
  H5::DataSet datasetN0 = file.createDataSet("n0", datatype, dataSpace);
  datasetN0.write(&n0[0], datatype);
  
  H5::DataSet datasetN1 = file.createDataSet("n1", datatype, dataSpace);
  datasetN1.write(&n1[0], datatype);
  
  H5::DataSet datasetDOccSigSig = file.createDataSet("dOccSigSig", datatype, dataSpace);
  datasetDOccSigSig.write(&dOccSigSig[0], datatype);
  
}

void writeSpectrumToFile(
    const std::vector<double> &gArr,
    const std::vector<double> &spectrum,
    const std::vector<double> &nc0Arr,
    const std::vector<double> &nd0Arr,
    const std::vector<double> &nc1Arr,
    const std::vector<double> &nd1Arr,
    const std::vector<double> &nBos0Arr,
    const std::vector<double> &nBos1Arr,
    const std::string &filename
)
{
  
  std::cout << "filename: " << filename << '\n';
  
  H5::H5File file(filename, H5F_ACC_TRUNC);
  
  writeAllPrms(file);
  
  const hsize_t dataSize = gArr.size();
  H5::DataSpace dataSpace(1, &dataSize);
  H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);
  
  H5::DataSet datasetGArr = file.createDataSet("gE", datatype, dataSpace);
  datasetGArr.write(&gArr[0], datatype);
  
  const hsize_t dataSize2D = spectrum.size();
  H5::DataSpace dataSpace2D(1, &dataSize2D);
  H5::FloatType datatype2D(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);
  
  H5::DataSet datasetSpectrum = file.createDataSet("spectrum", datatype2D, dataSpace2D);
  datasetSpectrum.write(&spectrum[0], datatype2D);
  
  H5::DataSet datasetnc0 = file.createDataSet("nc0", datatype2D, dataSpace2D);
  datasetnc0.write(nc0Arr.data(), datatype2D);
  
  H5::DataSet datasetnd0 = file.createDataSet("nd0", datatype2D, dataSpace2D);
  datasetnd0.write(nd0Arr.data(), datatype2D);
  
  H5::DataSet datasetnc1 = file.createDataSet("nc1", datatype2D, dataSpace2D);
  datasetnc1.write(nc1Arr.data(), datatype2D);
  
  H5::DataSet datasetnd1 = file.createDataSet("nd1", datatype2D, dataSpace2D);
  datasetnd1.write(nd1Arr.data(), datatype2D);
  
  H5::DataSet datasetnBos0 = file.createDataSet("nBos0", datatype2D, dataSpace2D);
  datasetnBos0.write(nBos0Arr.data(), datatype2D);
  
  H5::DataSet datasetnBos1 = file.createDataSet("nBos1", datatype2D, dataSpace2D);
  datasetnBos1.write(nBos1Arr.data(), datatype2D);
  
  
}

void writeStuffToHdf52BandsTime(
    const std::vector<double> &times,
    const std::vector<double> &pumpFunction,
    const std::vector<double> &dOcc0,
    const std::vector<double> &dOcc1,
    const std::vector<double> &dOccUpDn,
    const std::vector<double> &dOccSigSig,
    const std::vector<double> &Xph1,
    const std::vector<double> &Xph1Sqr,
    const std::vector<double> &Nph1,
    const std::vector<double> &Xph2,
    const std::vector<double> &Xph2Sqr,
    const std::vector<double> &Nph2,
    const std::vector<double> &N0,
    const std::vector<double> &N1,
    const std::string &filename
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
  
  H5::DataSet datasetDOcc0 = file.createDataSet("dOcc0", datatype, dataSpace);
  datasetDOcc0.write(&dOcc0[0], datatype);
  
  H5::DataSet datasetDOcc1 = file.createDataSet("dOcc1", datatype, dataSpace);
  datasetDOcc1.write(&dOcc1[0], datatype);
  
  H5::DataSet datasetDOccUpDn = file.createDataSet("dOccUpDn", datatype, dataSpace);
  datasetDOccUpDn.write(&dOccUpDn[0], datatype);
  
  H5::DataSet datasetDOccSigSig = file.createDataSet("dOccSigSig", datatype, dataSpace);
  datasetDOccSigSig.write(&dOccSigSig[0], datatype);
  
  H5::DataSet datasetN0 = file.createDataSet("n0", datatype, dataSpace);
  datasetN0.write(&N0[0], datatype);
  
  H5::DataSet datasetN1 = file.createDataSet("n1", datatype, dataSpace);
  datasetN1.write(&N1[0], datatype);
  
  H5::DataSet datasetXpt = file.createDataSet("Xph1", datatype, dataSpace);
  datasetXpt.write(&Xph1[0], datatype);
  
  H5::DataSet datasetXptSqr = file.createDataSet("Xph1Sqr", datatype, dataSpace);
  datasetXptSqr.write(&Xph1Sqr[0], datatype);
  
  H5::DataSet datasetNpt = file.createDataSet("Nph1", datatype, dataSpace);
  datasetNpt.write(&Nph1[0], datatype);
  
  H5::DataSet datasetXph = file.createDataSet("Xph2", datatype, dataSpace);
  datasetXph.write(&Xph2[0], datatype);
  
  H5::DataSet datasetXphSqr = file.createDataSet("Xph2Sqr", datatype, dataSpace);
  datasetXphSqr.write(&Xph2Sqr[0], datatype);
  
  H5::DataSet datasetNph = file.createDataSet("Nph2", datatype, dataSpace);
  datasetNph.write(&Nph2[0], datatype);
}




void readInComplex2DArray(std::vector<std::complex<double>> &readInArray, const std::string &fileName) {
  H5::H5File file(fileName, H5F_ACC_RDONLY);

  H5::DataSet realDataset = file.openDataSet("Real");
  H5::DataSet imagDataset = file.openDataSet("Imag");

  H5T_class_t typeClass = realDataset.getTypeClass();
  assert(typeClass == H5T_FLOAT);
  typeClass = imagDataset.getTypeClass();
  assert(typeClass == H5T_FLOAT);

  H5::DataSpace realDataSpace = realDataset.getSpace();
  H5::DataSpace imagDataSpace = imagDataset.getSpace();

  int rank = realDataSpace.getSimpleExtentNdims();
  assert(rank == 2);
  rank = imagDataSpace.getSimpleExtentNdims();
  assert(rank == 2);

  hsize_t dimsOut[2];
  realDataSpace.getSimpleExtentDims(dimsOut, nullptr);
  assert(dimsOut[0] * dimsOut[1] == readInArray.size());
  imagDataSpace.getSimpleExtentDims(dimsOut, nullptr);
  assert(dimsOut[0] * dimsOut[1] == readInArray.size());

  std::vector<double> realInput(readInArray.size(), 0.0);
  std::vector<double> imagInput(readInArray.size(), 0.0);

  const hsize_t inDimension[1] = {readInArray.size()};
  H5::DataSpace memspace(1, inDimension);

  realDataset.read(&realInput[0], H5::PredType::NATIVE_DOUBLE, memspace, realDataSpace);
  imagDataset.read(&imagInput[0], H5::PredType::NATIVE_DOUBLE, memspace, imagDataSpace);

  for (auto ind = 0ul; ind < readInArray.size(); ++ind) {
    readInArray[ind] = std::complex<double>(realInput[ind], imagInput[ind]);
  }
}