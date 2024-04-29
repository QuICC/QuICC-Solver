/**
 * @file WorlandProjector.cpp
 * @brief Source of the interface for a generic FFTW based Worland projector
 */

// External includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "Profiler/Interface.hpp"
#include "QuICC/Transform/Fft/Backend/PfSolve/WorlandProjector.hpp"
#include "Types/Math.hpp"
#include <iostream>
namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace PfSolve {

WorlandProjector::WorlandProjector() {}

WorlandProjector::~WorlandProjector() {
  this->freeBanded(false);
  this->freeBanded(true);
  deleteVkFFT(&mEvenApp);
  deleteVkFFT(&mOddApp);
  for (int i = 0; i < 2; i++) {
#if (VKFFT_BACKEND == 1)
    cudaFree(bufferSolve[i]);
    cudaFree(bufferSolve2[i]);
    cudaFree(bufferTemp0[i]);
    cudaStreamDestroy(pStream[i]);
#elif (VKFFT_BACKEND == 2)
    hipFree(bufferSolve[i]);
    hipFree(bufferSolve2[i]);
    hipFree(bufferTemp0[i]);
    hipStreamDestroy(pStream[i]);
#endif
  }
#if (VKFFT_BACKEND == 1)
  cudaFree(bufferSolve2[2]);
#elif (VKFFT_BACKEND == 2)
  hipFree(bufferSolve2[2]);
#endif
}

void WorlandProjector::init(const SetupType &setup, const int lshift,
                            const int extraN, const bool lshiftOnlyParity,
                            const bool alwaysZeroNegative) const {
  // Get Gpu device
#if (VKFFT_BACKEND == 1)
  cuInit(0);
  cuDeviceGet(&this->device, 0);
#elif (VKFFT_BACKEND == 2)
  hipInit(0);
  hipDeviceGet(&this->device, 0);
#endif
  // Initialize parent
  IWorlandBackend::init(setup, lshift, extraN, lshiftOnlyParity,
                        alwaysZeroNegative);
  this->mpWorkTmp = &this->mInTmp;

  int fwdSize = setup.fwdSize();
  int bwdSize = setup.bwdSize();

  // Set internal spectral resolution
  this->setWSize();

  // Create the plan
  const int *fftSize = &fwdSize;

  // Initialize temporary storage
  int blockSize = std::max(this->mEBlockSize, this->mOBlockSize);
  // Input temporary storage
  this->initStorage(bwdSize, blockSize, 1, this->mInTmp);
  // Output temporary storage
  this->initStorage(fwdSize, blockSize, 1, this->mOutTmp);

  // Set sizes
  this->mPadSize = setup.padSize();

  // Initialise temporary storage
  blockSize = this->mEBlockSize;
  Matrix tmpF = Matrix::Zero(fwdSize, blockSize);
  Matrix tmpB = Matrix::Zero(bwdSize, blockSize);

  // Create the spectral to physical plan
  mEvenConfiguration.FFTdim = 1;
  mEvenConfiguration.size[0] = fwdSize;
  mEvenConfiguration.numberBatches = 2 * blockSize;
  mEvenConfiguration.doublePrecision = 1;
  mEvenConfiguration.performDCT = 3;
  mEvenConfiguration.device = &this->device;
  mEvenConfiguration.num_streams = 1;
  mEvenSize = fwdSize * 2 * blockSize * sizeof(double);
  mEvenConfiguration.bufferSize = &mEvenSize;

  initializeVkFFT(&mEvenApp, mEvenConfiguration);

  // Initialise temporary odd storage
  blockSize = this->mOBlockSize;
  tmpF = Matrix::Zero(fwdSize, blockSize);
  tmpB = Matrix::Zero(bwdSize, blockSize);

  // Create the odd physical to spectral plan
  mOddConfiguration.FFTdim = 1;
  mOddConfiguration.size[0] = bwdSize;
  mOddConfiguration.numberBatches = 2 * blockSize;
  mOddConfiguration.doublePrecision = 1;
  mOddConfiguration.performDCT = 4;
  mOddConfiguration.device = &this->device;
  mOddConfiguration.num_streams = 1;
  mOddSize = bwdSize * 2 * blockSize * sizeof(double);
  mOddConfiguration.bufferSize = &mOddSize;

  initializeVkFFT(&mOddApp, mOddConfiguration);
  // Initialize Jacobi shift matrices
  this->initJ();

  // Initialize Banded matrices
  this->initBanded(false);
  this->initBanded(true);

#if (VKFFT_BACKEND == 1)
  cudaMalloc((void **)&bufferSolve[0], mEvenSize + mOddSize);
  cudaMalloc((void **)&bufferSolve2[0], mEvenSize + mOddSize);
  cudaMalloc((void **)&bufferTemp0[0], mEvenSize + mOddSize);
  cudaStreamCreate(&pStream[0]);

  cudaMalloc((void **)&bufferSolve[1], mEvenSize + mOddSize);
  cudaMalloc((void **)&bufferSolve2[1], mEvenSize + mOddSize);
  cudaMalloc((void **)&bufferTemp0[1], mEvenSize + mOddSize);
  cudaStreamCreate(&pStream[1]);

  cudaMalloc((void **)&bufferSolve2[2], 8 * (mEvenSize + mOddSize));
#elif (VKFFT_BACKEND == 2)
  hipMalloc((void **)&bufferSolve[0], mEvenSize + mOddSize);
  hipMalloc((void **)&bufferSolve2[0], mEvenSize + mOddSize);
  hipMalloc((void **)&bufferTemp0[0], mEvenSize + mOddSize);
  hipStreamCreate(&pStream[0]);

  hipMalloc((void **)&bufferSolve[1], mEvenSize + mOddSize);
  hipMalloc((void **)&bufferSolve2[1], mEvenSize + mOddSize);
  hipMalloc((void **)&bufferTemp0[1], mEvenSize + mOddSize);
  hipStreamCreate(&pStream[1]);

  hipMalloc((void **)&bufferSolve2[2], 8 * (mEvenSize + mOddSize));
#endif
}

int WorlandProjector::lSize(const int l) const { return this->mWSize - l / 2; }

int WorlandProjector::lSize(const int i, const int) const {
  return this->mWSize - i;
}

void WorlandProjector::io(const bool isEven) const {
  this->setPlan(isEven);

  this->io(this->mOutTmp.at(0).data(), this->mInTmp.at(0).data());
}

void WorlandProjector::input(const MatrixZ &in) const {
  Matrix &inTmp = this->mInTmp.at(0);

  PfSolve_JW::PfSolveResult resSolve = PfSolve_JW::PFSOLVE_SUCCESS;
  bufferSolveRes2 = bufferSolve2[2];
  Profiler::RegionStart<2>("applyOperators");

  assert(in.rows() * in.cols() * sizeof(MHDComplex) <= 8 * (mEvenSize + mOddSize));
  IWorlandBackend::transferDataToGPU(&bufferSolveRes2, (void *)in.data(), 0, 0,
                                     in.rows() * in.cols() * sizeof(MHDComplex),
                                     &pStream[currentStream]);
  cudaDeviceSynchronize();
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    if (this->isPhysEven(isEven))
      cudaMemsetAsync(bufferSolveRes, 0, mEvenSize, pStream[currentStream]);
    else
      cudaMemsetAsync(bufferSolveRes, 0, mOddSize, pStream[currentStream]);
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    for (auto &loc : *this->pLoc(isEven)) {
      PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
      PfSolve_JW::PfSolve_MapKey_block mapKey = {};
      mapKey.size[0] = in.rows();
      mapKey.size[1] = std::get<2>(loc);
      mapKey.type = BLOCK_READ_COMPLEX_STRIDED_WRITE_COMPLEX_PACKED;
      mapKey.type += (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION +
                      RUNTIME_INPUTZEROPAD + RUNTIME_OUTPUTZEROPAD +
                      RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE) *
                     1000;

      resSolve =
          PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      if (!tempAppCopy) {
        PfSolve_JW::PfSolveApplication appCopy = {};
        PfSolve_JW::PfSolveConfiguration configurationCopy = {};

        configurationCopy.size[0] = mapKey.size[0];
        configurationCopy.size[1] = mapKey.size[1];
        configurationCopy.block = mapKey.type % 1000;
        configurationCopy.isOutputFormatted = 1;
        configurationCopy.jw_control_bitmask = mapKey.type / 1000;
        configurationCopy.doublePrecision = 1;
        configurationCopy.device = &this->device;
        configurationCopy.num_streams = 1;

        resSolve = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
        resSolve =
            PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
        resSolve =
            PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      }

      PfSolve_JW::PfSolveLaunchParams launchParams = {};
      launchParams.offsetM = std::get<1>(loc) * in.rows();
      launchParams.offsetSolution = start * inTmp.rows();
      launchParams.inputZeropad[0] = 0;
      launchParams.inputZeropad[1] = mapKey.size[0];
      launchParams.outputZeropad[0] = 0;
      launchParams.outputZeropad[1] = mapKey.size[0];
      launchParams.inputBufferStride = mapKey.size[0];
      launchParams.outputBufferStride = inTmp.rows();
      launchParams.buffer = (void **)&this->bufferSolveRes2;
      launchParams.outputBuffer = (void **)&this->bufferSolveRes;
      tempAppCopy->configuration.stream = &pStream[currentStream];

      resSolve = PfSolveAppend(tempAppCopy, -1, &launchParams);

      prevOutputZeropad[currentStream][0][start] =
          launchParams.outputZeropad[0];
      prevOutputZeropad[currentStream][1][start] =
          launchParams.outputZeropad[1];

      start += std::get<2>(loc);
      std::get<3>(loc) = in.rows();
    }
  }
  Profiler::RegionStart<2>("applyImpl");
}

void WorlandProjector::input(const Matrix &in, const bool isEven,
                             const bool needPadding) const {
  if (isPhysEven(isEven)) {
    currentStream = 0;
    bufferSolveRes = bufferSolve[0];
    bufferSolveRes2 = bufferSolve2[0];
    bufferTemp = bufferTemp0[0];
  } else {
    currentStream = 1;
    bufferSolveRes = bufferSolve[1];
    bufferSolveRes2 = bufferSolve2[1];
    bufferTemp = bufferTemp0[1];
  }
  Matrix &inTmp = this->mInTmp.at(0);
  int start = 0;
  for (auto &loc : *this->pLoc(isEven)) {
    inTmp.block(0, start, this->mSpecSize, std::get<2>(loc)) =
        in.block(0, std::get<1>(loc), this->mSpecSize, std::get<2>(loc));
    start += std::get<2>(loc);
    std::get<3>(loc) = this->mSpecSize;
  }

  // Apply padding if required
  if (needPadding) {
    this->applyPadding(inTmp);
  }
  IWorlandBackend::transferDataToGPU(&bufferSolveRes, inTmp.data(), 0, 0,
                                     inTmp.size() * sizeof(MHDFloat),
                                     &pStream[currentStream]);
}

void WorlandProjector::input(const MatrixZ &in, const bool isEven,
                             const bool useReal, const bool needPadding) const {
  WorlandProjector::input(in);
}

void WorlandProjector::output(MatrixZ &rOut) const {
  Profiler::RegionStop<2>("applyImpl");

  Matrix &outTmp = this->mOutTmp.at(0);

  PfSolve_JW::PfSolveResult resSolve = PfSolve_JW::PFSOLVE_SUCCESS;
  bufferSolveRes2 = bufferSolve2[2];

  assert(rOut.rows() * rOut.cols() * sizeof(MHDComplex) <= 8 * (mEvenSize + mOddSize));
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    for (auto &loc : *this->pLoc(isEven)) {
      PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
      PfSolve_JW::PfSolve_MapKey_block mapKey = {};
      mapKey.size[0] = rOut.rows();
      mapKey.size[1] = std::get<2>(loc);
      mapKey.type = BLOCK_READ_COMPLEX_PACKED_WRITE_COMPLEX_STRIDED;
      mapKey.type += (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION +
                      RUNTIME_INPUTZEROPAD + RUNTIME_OUTPUTZEROPAD +
                      RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE) *
                     1000;

      resSolve =
          PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      if (!tempAppCopy) {
        PfSolve_JW::PfSolveApplication appCopy = {};
        PfSolve_JW::PfSolveConfiguration configurationCopy = {};

        configurationCopy.size[0] = mapKey.size[0];
        configurationCopy.size[1] = mapKey.size[1];
        configurationCopy.block = mapKey.type % 1000;
        configurationCopy.isOutputFormatted = 1;
        configurationCopy.jw_control_bitmask = mapKey.type / 1000;
        configurationCopy.doublePrecision = 1;
        configurationCopy.device = &this->device;
        configurationCopy.num_streams = 1;

        resSolve = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
        resSolve =
            PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
        resSolve =
            PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      }

      PfSolve_JW::PfSolveLaunchParams launchParams = {};
      launchParams.offsetM = start * outTmp.rows();
      launchParams.offsetSolution = std::get<1>(loc) * rOut.rows();
      launchParams.inputZeropad[0] = 0;
      launchParams.inputZeropad[1] = mapKey.size[0];
      launchParams.outputZeropad[0] = 0;
      launchParams.outputZeropad[1] = mapKey.size[0];
      launchParams.inputBufferStride = outTmp.rows();
      launchParams.outputBufferStride = rOut.rows();
      launchParams.buffer = (void **)&this->bufferSolveRes;
      launchParams.outputBuffer = (void **)&this->bufferSolveRes2;
      tempAppCopy->configuration.stream = &pStream[currentStream];

      resSolve = PfSolveAppend(tempAppCopy, -1, &launchParams);

      prevOutputZeropad[currentStream][0][start] =
          launchParams.outputZeropad[0];
      prevOutputZeropad[currentStream][1][start] =
          launchParams.outputZeropad[1];
      start += std::get<2>(loc);
    }
  }

  cudaDeviceSynchronize();
  IWorlandBackend::transferDataFromGPU(
      rOut.data(), &bufferSolveRes2, 0, 0,
      rOut.rows() * rOut.cols() * sizeof(MHDComplex), &pStream[currentStream]);
  for (auto loc : this->mZLoc) {
    int s = std::get<1>(loc);
    int cols = std::get<2>(loc);
    rOut.block(0, s, rOut.rows(), cols).setZero();
  }
  Profiler::RegionStop<2>("applyOperators");
}

void WorlandProjector::output(Matrix &rOut, const bool isEven) const {
  Matrix &outTmp = this->mOutTmp.at(0);
  int start = 0;
  IWorlandBackend::transferDataFromGPU(outTmp.data(), &bufferSolveRes, 0, 0,
                                       outTmp.size() * sizeof(MHDFloat),
                                       &pStream[currentStream]);

  for (auto loc : *this->pLoc(isEven)) {
    rOut.block(0, std::get<1>(loc), rOut.rows(), std::get<2>(loc)) =
        outTmp.block(0, start, rOut.rows(), std::get<2>(loc));
    start += std::get<2>(loc);
  }

  // Zero unused values
  for (auto loc : this->mZLoc) {
    int s = std::get<1>(loc);
    int cols = std::get<2>(loc);
    rOut.block(0, s, rOut.rows(), cols).setZero();
  }
}

void WorlandProjector::output(MatrixZ &rOut, const bool isEven,
                              const bool useReal) const {
  WorlandProjector::output(rOut);
}

void WorlandProjector::applyPadding(Matrix &rData, const int extraRows) const {
  if (extraRows >= 0) {
    // Set the padded values to zero
    rData.bottomRows(this->mPadSize - extraRows).setZero();
  } else {
    rData.bottomRows(-extraRows).setZero();
  }
}

void WorlandProjector::applyFft() const {
  for (int isEven = 0; isEven < 2; isEven++) {
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    this->setPlan(isEven);
    VkFFT::VkFFTLaunchParams launchParams = {};
    launchParams.buffer = (void **)&this->bufferSolveRes;
    mpApp->configuration.stream = &pStream[currentStream];
    VkFFT::VkFFTResult rr = VkFFT::VkFFTAppend(mpApp, -1, &launchParams);
  }
}

void WorlandProjector::partialBackwardWorland(
    const int l, const int id, const int i0, const int start, const int cols,
    const std::vector<Matrix> &banded, const std::vector<void *> &bandedGPU,
    Matrix &in) const {
  if (i0 >= l / 2) {
    void **GPUarr = (void **)malloc((i0 - l / 2 + 1) * sizeof(void *));
    for (int i = i0; i >= l / 2; --i) {
      GPUarr[i0 - i] = bandedGPU.at(bandedGPU.size() - 1 - i);
    }
  
    int i_end = l / 2 - 1;
    this->applyPairCombined(in, id, start, cols, this->lSize(i0, l),
                            this->lSize(i_end, l), banded.at(i0),
                            (const void **)GPUarr);
    free(GPUarr);
  }
}

void WorlandProjector::backwardWorland(const bool isEven0,
                                       const unsigned int id) const {
  PfSolve_JW::PfSolveResult resSolve = PfSolve_JW::PFSOLVE_SUCCESS;
  Matrix &inTmp = this->mInTmp.at(id);

  for (int isEven = 0; isEven < 2; isEven++) {
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    void *wTmpGPU = this->workTmpGPU(id);
    const std::vector<Matrix> &banded = this->banded(isEven);
    const std::vector<void *> &bandedGPU = this->bandedGPU(isEven);

    std::vector<LocationType> *pL = this->pLoc(isEven, id);
    if (pL->begin() != pL->end()) {
      int cols = 0;
      int maxL = std::get<4>(*pL->rbegin());
      int start = this->blockSize(isEven);
      int l;
      for (auto loc = pL->rbegin(); loc != pL->rend(); ++loc) {
        l = std::get<4>(*loc);
        this->partialBackwardWorland(l, id, maxL / 2 - 1, start, cols, banded,
                                     bandedGPU, inTmp);
        maxL = l;
        cols += std::get<2>(*loc);
        start -= std::get<2>(*loc);
      }

      // Last iteration if l=0/l=1 not in list
      if (std::get<4>(pL->at(0)) > 1) {
        l = static_cast<int>(!this->isPhysEven(isEven));
        this->partialBackwardWorland(l, id, maxL / 2 - 1, start, cols, banded,
                                     bandedGPU, inTmp);
      }

      // Rescale first mode for l = 0 for FFT
      if (std::get<4>(pL->at(0)) == 0) {

        PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
        PfSolve_JW::PfSolve_MapKey_block mapKey = {};
        mapKey.size[0] = std::get<2>(pL->at(0));
        mapKey.size[1] = 2;
        mapKey.type = (RUNTIME_SCALEC + RUNTIME_INPUTBUFFERSTRIDE +
                       RUNTIME_OUTPUTBUFFERSTRIDE) *
                          1000 +
                      BLOCK_SCALEC + BLOCK_READ_REAL_WRITE_REAL;
        resSolve =
            PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
        // Initialize applications. This function loads shaders, creates
        // pipeline and configures FFT based on configuration file. No buffer
        // allocations inside PfSolve library.
        if (!tempAppCopy) {
          PfSolve_JW::PfSolveApplication appCopy = {};
          PfSolve_JW::PfSolveConfiguration configurationCopy = {};

          // create PfSolve appSolve
          configurationCopy.size[0] = mapKey.size[0];
          configurationCopy.size[1] = mapKey.size[1];
          configurationCopy.block = mapKey.type % 1000;
          configurationCopy.isOutputFormatted = 1;
          configurationCopy.jw_control_bitmask = mapKey.type / 1000;
          configurationCopy.doublePrecision = 1;
          configurationCopy.device = &this->device;
          configurationCopy.num_streams = 1;

          resSolve = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
          resSolve =
              PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
          resSolve =
              PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
        }

        PfSolve_JW::PfSolveLaunchParams launchParams = {};
        launchParams.offsetM = 0;
        launchParams.offsetSolution = 0;
        launchParams.buffer = (void **)&wTmpGPU;
        launchParams.outputBuffer = (void **)&wTmpGPU;
        launchParams.inputBufferStride = inTmp.rows();
        launchParams.outputBufferStride = inTmp.rows();
        launchParams.scaleC = std::sqrt(2.0);
        tempAppCopy->configuration.stream = &pStream[currentStream];
        resSolve = PfSolve_JW::PfSolveAppend(tempAppCopy, -1, &launchParams);
      }

      // reset current l
      this->resetLocations(isEven, id);
    }
  }
}

void WorlandProjector::lowerAlpha(const MHDFloat alpha, const bool isEven0,
                                  const unsigned int id,
                                  const MHDFloat norm) const {
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &wTmp = this->workTmp(id);
    void *wTmpGPU = this->workTmpGPU(id);
    for (auto loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc);
      int cols = std::get<2>(loc);
      Matrix U;
      this->buildShiftU(U, this->lSize(l), alpha, l - 0.5, norm);
      this->applyTriSolve(wTmp, id, start, cols, 1, U);
      start += cols;
    }
  }
}

void WorlandProjector::lowerBeta(const MHDFloat alpha, const bool isEven0,
                                 const unsigned int id,
                                 const MHDFloat norm) const {
  PfSolve_JW::PfSolveResult resSolve = PfSolve_JW::PFSOLVE_SUCCESS;
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &wTmp = this->workTmp(id);
    void *wTmpGPU = this->workTmpGPU(id);
    for (auto &loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc);
      int cols = std::get<2>(loc);
      if (l > 0) {
        Matrix V;
        this->buildShiftV(V, std::get<3>(loc), alpha, l - 0.5, norm);
        MHDFloat scale = (l == 1) ? 1.0 / std::sqrt(2.0) : 1;
        this->applyTriSolve(wTmp, id, start, cols, scale, V);
        std::get<4>(loc)--;
      } else {
        PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
        PfSolve_JW::PfSolve_MapKey_block mapKey = {};
        mapKey.size[0] = this->lSize(l);
        mapKey.size[1] = 2 * cols;
        mapKey.type =
            (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION + RUNTIME_INPUTZEROPAD +
             RUNTIME_OUTPUTZEROPAD + RUNTIME_SCALEC +
             RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE) *
                1000 +
            BLOCK_SCALEC + BLOCK_READ_REAL_WRITE_REAL;
        resSolve =
            PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
        // Initialize applications. This function loads shaders, creates
        // pipeline and configures FFT based on configuration file. No buffer
        // allocations inside PfSolve library.
        if (!tempAppCopy) {
          PfSolve_JW::PfSolveApplication appCopy = {};
          PfSolve_JW::PfSolveConfiguration configurationCopy = {};

          // create PfSolve appSolve
          configurationCopy.size[0] = mapKey.size[0];
          configurationCopy.size[1] = mapKey.size[1];
          configurationCopy.block = mapKey.type % 1000;
          configurationCopy.isOutputFormatted = 1;
          configurationCopy.jw_control_bitmask = mapKey.type / 1000;
          configurationCopy.doublePrecision = 1;
          configurationCopy.device = &this->device;
          configurationCopy.num_streams = 1;

          resSolve = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
          resSolve =
              PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
          resSolve =
              PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
        }

        PfSolve_JW::PfSolveLaunchParams launchParams = {};
        launchParams.offsetM = 2 * start * wTmp.rows();
        launchParams.offsetSolution = 2 * start * wTmp.rows();
        launchParams.buffer = (void **)&this->bufferSolveRes;
        launchParams.outputBuffer = (void **)&this->bufferSolveRes;
        launchParams.scaleC = 1.0 / norm;
        launchParams.inputZeropad[0] =
            prevOutputZeropad[currentStream][0][start];
        launchParams.inputZeropad[1] =
            prevOutputZeropad[currentStream][1][start];
        launchParams.outputZeropad[0] =
            prevOutputZeropad[currentStream][0][start];
        launchParams.outputZeropad[1] =
            prevOutputZeropad[currentStream][1][start];
        launchParams.inputBufferStride = wTmp.rows();
        launchParams.outputBufferStride = wTmp.rows();
        tempAppCopy->configuration.stream = &pStream[currentStream];
        resSolve = PfSolve_JW::PfSolveAppend(tempAppCopy, -1, &launchParams);
        prevOutputZeropad[currentStream][0][start] =
            launchParams.outputZeropad[0];
        prevOutputZeropad[currentStream][1][start] =
            launchParams.outputZeropad[1];

        std::get<4>(loc)++;
      }
      start += cols;
    }
  }
}

void WorlandProjector::lowerR2Beta(const MHDFloat alpha, const bool isEven0,
                                   const unsigned int id,
                                   const MHDFloat norm) const {
  PfSolve_JW::PfSolveResult resSolve = PfSolve_JW::PFSOLVE_SUCCESS;
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &wTmp = this->workTmp(id);
    void *wTmpGPU = this->workTmpGPU(id);
    for (auto &loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc);
      int cols = std::get<2>(loc);
      if (l < 0) {
        PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
        PfSolve_JW::PfSolve_MapKey_block mapKey = {};
        mapKey.size[0] = this->lSize(l);
        mapKey.size[1] = 2 * cols;
        mapKey.type = (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION +
                       RUNTIME_INPUTZEROPAD + RUNTIME_OUTPUTZEROPAD +
                       RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE) *
                          1000 +
                      BLOCK_SCALEC + BLOCK_READ_REAL_WRITE_REAL;
        resSolve =
            PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
        // Initialize applications. This function loads shaders, creates
        // pipeline and configures FFT based on configuration file. No buffer
        // allocations inside PfSolve library.
        if (!tempAppCopy) {
          PfSolve_JW::PfSolveApplication appCopy = {};
          PfSolve_JW::PfSolveConfiguration configurationCopy = {};

          // create PfSolve appSolve
          configurationCopy.size[0] = mapKey.size[0];
          configurationCopy.size[1] = mapKey.size[1];
          configurationCopy.block = mapKey.type % 1000;
          configurationCopy.isOutputFormatted = 1;
          configurationCopy.jw_control_bitmask = mapKey.type / 1000;
          configurationCopy.doublePrecision = 1;
          configurationCopy.scaleC = 0;
          configurationCopy.device = &this->device;
          configurationCopy.num_streams = 1;

          resSolve = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
          resSolve =
              PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
          resSolve =
              PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
        }
        PfSolve_JW::PfSolveLaunchParams launchParams = {};
        launchParams.offsetM = 2 * start * wTmp.rows();
        launchParams.offsetSolution = 2 * start * wTmp.rows();
        launchParams.buffer = (void **)&this->bufferSolveRes;
        launchParams.outputBuffer = (void **)&this->bufferSolveRes;
        launchParams.inputZeropad[0] =
            prevOutputZeropad[currentStream][0][start];
        launchParams.inputZeropad[1] =
            prevOutputZeropad[currentStream][1][start];
        launchParams.outputZeropad[0] =
            prevOutputZeropad[currentStream][0][start];
        launchParams.outputZeropad[1] =
            prevOutputZeropad[currentStream][1][start];
        launchParams.inputBufferStride = wTmp.rows();
        launchParams.outputBufferStride = wTmp.rows();
        tempAppCopy->configuration.stream = &pStream[currentStream];
        resSolve = PfSolve_JW::PfSolveAppend(tempAppCopy, -1, &launchParams);

        prevOutputZeropad[currentStream][0][start] =
            launchParams.outputZeropad[0];
        prevOutputZeropad[currentStream][1][start] =
            launchParams.outputZeropad[1];

      } else if (l > 1) {
        Matrix M;
        this->buildShiftM(M, this->lSize(l), alpha, l - 0.5, norm);
        this->applyTriProduct(wTmp, id, start, cols, M);
        std::get<4>(loc)--;
      } else {
        PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
        PfSolve_JW::PfSolve_MapKey_block mapKey = {};
        mapKey.size[0] = this->lSize(l);
        mapKey.size[1] = 2 * cols;
        mapKey.type = (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION +
                       RUNTIME_INPUTZEROPAD + RUNTIME_OUTPUTZEROPAD +
                       RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE) *
                          1000 +
                      BLOCK_SCALEC + BLOCK_READ_REAL_WRITE_REAL;
        resSolve =
            PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
        // Initialize applications. This function loads shaders, creates
        // pipeline and configures FFT based on configuration file. No buffer
        // allocations inside PfSolve library.
        if (!tempAppCopy) {
          PfSolve_JW::PfSolveApplication appCopy = {};
          PfSolve_JW::PfSolveConfiguration configurationCopy = {};

          // create PfSolve appSolve
          configurationCopy.size[0] = mapKey.size[0];
          configurationCopy.size[1] = mapKey.size[1];
          configurationCopy.block = mapKey.type % 1000;
          configurationCopy.isOutputFormatted = 1;
          configurationCopy.jw_control_bitmask = mapKey.type / 1000;
          configurationCopy.doublePrecision = 1;
          configurationCopy.scaleC = norm;
          configurationCopy.device = &this->device;
          configurationCopy.num_streams = 1;

          resSolve = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
          resSolve =
              PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
          resSolve =
              PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
        }

        PfSolve_JW::PfSolveLaunchParams launchParams = {};
        launchParams.offsetM = 2 * start * wTmp.rows();
        launchParams.offsetSolution = 2 * start * wTmp.rows();
        launchParams.buffer = (void **)&this->bufferSolveRes;
        launchParams.outputBuffer = (void **)&this->bufferSolveRes;
        launchParams.inputZeropad[0] =
            prevOutputZeropad[currentStream][0][start];
        launchParams.inputZeropad[1] =
            prevOutputZeropad[currentStream][1][start];
        launchParams.outputZeropad[0] =
            prevOutputZeropad[currentStream][0][start];
        launchParams.outputZeropad[1] =
            prevOutputZeropad[currentStream][1][start];
        launchParams.inputBufferStride = wTmp.rows();
        launchParams.outputBufferStride = wTmp.rows();
        tempAppCopy->configuration.stream = &pStream[currentStream];
        resSolve = PfSolve_JW::PfSolveAppend(tempAppCopy, -1, &launchParams);
        prevOutputZeropad[currentStream][0][start] =
            launchParams.outputZeropad[0];
        prevOutputZeropad[currentStream][1][start] =
            launchParams.outputZeropad[1];

        std::get<4>(loc)--;
      }
      start += cols;
    }
  }
}

void WorlandProjector::initBanded(const bool isEven) const {
  if (this->pLoc(isEven, 0)->size() > 0) {
    const Matrix &J = this->J(isEven);

    const MHDFloat normV = 1.0;
    const MHDFloat normM = 1.0;
    std::vector<Matrix> *pBanded;
    std::vector<void *> *pBandedGPU;
    if (isPhysEven(isEven)) {
      pBanded = &this->mEBanded;
      pBandedGPU = &this->mEBandedGPU;
    } else {
      pBanded = &this->mOBanded;
      pBandedGPU = &this->mOBandedGPU;
    }

    std::vector<LocationType> *pL = this->pLoc(isEven, 0);
    if (pL->begin() != pL->end()) {
      int maxL = std::get<4>(*pL->rbegin());
      int l;
      for (auto loc = pL->rbegin(); loc != pL->rend(); ++loc) {
        l = std::get<4>(*loc);
#if (VKFFT_BACKEND == 1)
        int warpSize = 32;
#elif (VKFFT_BACKEND == 2)
        int warpSize = 64;
#endif
        int64_t registers_per_thread =
            (uint64_t)ceil(this->lSize(l / 2, l) / (double)warpSize);

        for (int i = maxL / 2 - 1; i >= l / 2; --i) {
          Matrix PS;
          void *tempBuffer;
          this->buildShiftPair(PS, false, i, this->lSize(i, l), J, normV, normM,
                               false);

          // additional GPU operations
          Matrix B2 = PS.transpose();
          double *bb2 = B2.data();
          for (int j = 0; j < B2.rows(); j++) {
            bb2[j + 2 * B2.rows()] = -bb2[j + 2 * B2.rows()];
          }

          Matrix B = B2;
          double *bb = B.data();
          for (int j = 0; j < B.rows(); j++) {
            bb[reorder_indices_PfSolve(j, B.rows(), warpSize,
                                       registers_per_thread) +
               0 * B.rows()] = bb2[j + 0 * B.rows()];
            bb[reorder_indices_PfSolve(j, B.rows(), warpSize,
                                       registers_per_thread) +
               1 * B.rows()] = bb2[j + 1 * B.rows()];
            bb[reorder_indices_PfSolve(j, B.rows(), warpSize,
                                       registers_per_thread) +
               2 * B.rows()] = bb2[j + 2 * B.rows()];
          }

#if (VKFFT_BACKEND == 1)
          cudaMalloc((void **)&tempBuffer, 3 * B.rows() * sizeof(MHDFloat));
#elif (VKFFT_BACKEND == 2)
          hipMalloc((void **)&tempBuffer, 3 * B.rows() * sizeof(MHDFloat));
#endif
          IWorlandBackend::transferDataToGPU(&tempBuffer, bb, 0, 0,
                                             3 * B.rows() * sizeof(MHDFloat),
                                             &pStream[currentStream]);

          pBanded->push_back(PS);
          pBandedGPU->push_back(tempBuffer);
        }
        maxL = l;
      }

      // Last iteration if l=0/l=1 not in list
      if (std::get<4>(pL->at(0)) > 1) {
        l = static_cast<int>(!this->isPhysEven(isEven));
#if (VKFFT_BACKEND == 1)
        int warpSize = 32;
#elif (VKFFT_BACKEND == 2)
        int warpSize = 64;
#endif
        int64_t registers_per_thread =
            (uint64_t)ceil(this->lSize(l / 2, l) / (double)warpSize);
        for (int i = maxL / 2 - 1; i >= l / 2; --i) {
          Matrix PS;
          void *tempBuffer;
          this->buildShiftPair(PS, false, i, this->lSize(i, l), J, normV, normM,
                               false);

          // additional GPU operations
          Matrix B2 = PS.transpose();
          double *bb2 = B2.data();
          for (int j = 0; j < B2.rows(); j++) {
            bb2[j + 2 * B2.rows()] = -bb2[j + 2 * B2.rows()];
          }

          Matrix B = B2;
          double *bb = B.data();
          for (int j = 0; j < B.rows(); j++) {
            bb[reorder_indices_PfSolve(j, B.rows(), warpSize,
                                       registers_per_thread) +
               0 * B.rows()] = bb2[j + 0 * B.rows()];
            bb[reorder_indices_PfSolve(j, B.rows(), warpSize,
                                       registers_per_thread) +
               1 * B.rows()] = bb2[j + 1 * B.rows()];
            bb[reorder_indices_PfSolve(j, B.rows(), warpSize,
                                       registers_per_thread) +
               2 * B.rows()] = bb2[j + 2 * B.rows()];
          }

#if (VKFFT_BACKEND == 1)
          cudaMalloc((void **)&tempBuffer, 3 * B.rows() * sizeof(MHDFloat));
#elif (VKFFT_BACKEND == 2)
          hipMalloc((void **)&tempBuffer, 3 * B.rows() * sizeof(MHDFloat));
#endif
          IWorlandBackend::transferDataToGPU(&tempBuffer, bb, 0, 0,
                                             3 * B.rows() * sizeof(MHDFloat),
                                             &pStream[currentStream]);

          pBanded->push_back(PS);
          pBandedGPU->push_back(tempBuffer);
        }
      }
    }
  }
}

void WorlandProjector::freeBanded(const bool isEven) const {
  if (this->pLoc(isEven, 0)->size() > 0) {
    std::vector<void *> *pBandedGPU;
    if (isPhysEven(isEven)) {
      pBandedGPU = &this->mEBandedGPU;
    } else {
      pBandedGPU = &this->mOBandedGPU;
    }
    for (int i = 0; i < pBandedGPU->size(); ++i) {
#if (VKFFT_BACKEND == 1)
      cudaFree(pBandedGPU->at(i));
#elif (VKFFT_BACKEND == 2)
      hipFree(pBandedGPU->at(i));
#endif
    }
  }
}

} // namespace PfSolve
} // namespace Backend
} // namespace Fft
} // namespace Transform
} // namespace QuICC
