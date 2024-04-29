/**
 * @file WorlandIntegrator.cpp
 * @brief Source of the interface for a generic FFTW based Worland integrator
 */

// External includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "Profiler/Interface.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"
#include "QuICC/Transform/Fft/Backend/PfSolve/WorlandIntegrator.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace PfSolve {

WorlandIntegrator::WorlandIntegrator() {}

WorlandIntegrator::~WorlandIntegrator() {
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

void WorlandIntegrator::init(const SetupType &setup, const int lshift,
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
  this->mpWorkTmp = &this->mOutTmp;

  int fwdSize = setup.fwdSize();
  int bwdSize = setup.bwdSize();

  // Set Internal:: spectral resolution
  this->setWSize();
  assert(this->mWSize <= bwdSize);

  // Set transform scaling
  this->mFftScaling = std::sqrt(Math::PI) / static_cast<MHDFloat>(2 * fwdSize);

  // Initialize main temporary storage
  int blockSize = std::max(this->mEBlockSize, this->mOBlockSize);
  // Input temporary storage
  this->initStorage(fwdSize, blockSize, 1, this->mInTmp);
  // Output temporary storage
  this->initStorage(bwdSize, blockSize, 1, this->mOutTmp);

  // Set sizes
  this->mBwdSize = setup.bwdSize();

  // Create the plan size
  const int *fftSize = &fwdSize;

  // Initialise temporary even storage
  blockSize = this->mEBlockSize;
  Matrix tmpF = Matrix::Zero(fwdSize, blockSize);
  Matrix tmpB = Matrix::Zero(bwdSize, blockSize);

  // Create the even physical to spectral plan
  mEvenConfiguration.FFTdim = 1;
  mEvenConfiguration.size[0] = bwdSize;
  mEvenConfiguration.numberBatches = 2 * blockSize;
  mEvenConfiguration.doublePrecision = 1;
  mEvenConfiguration.performDCT = 2;
  mEvenConfiguration.device = &this->device;
  mEvenConfiguration.num_streams = 1;
  mEvenSize = bwdSize * 2 * blockSize * 8;
  mEvenConfiguration.bufferSize = &mEvenSize;

  initializeVkFFT(&mEvenApp, mEvenConfiguration);

  // Initialise temporary odd storage
  blockSize = this->mOBlockSize;

  tmpF = Matrix::Zero(fwdSize, blockSize);
  tmpB = Matrix::Zero(bwdSize, blockSize);

  // Create the odd physical to spectral plan
  mOddConfiguration.FFTdim = 1;
  mOddConfiguration.size[0] = fwdSize;
  mOddConfiguration.numberBatches = 2 * blockSize;
  mOddConfiguration.doublePrecision = 1;
  mOddConfiguration.performDCT = 4;
  mOddConfiguration.device = &this->device;
  mOddConfiguration.num_streams = 1;
  mOddSize = fwdSize * 2 * blockSize * 8;
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

int WorlandIntegrator::lSize(const int l) const { return this->mWSize - l / 2; }

int WorlandIntegrator::lSize(const int i, const int l) const {
  assert(l >= 0);
  return this->mWSize - i;
}

void WorlandIntegrator::io(const bool isEven) const {
  this->setPlan(isEven);

  this->io(this->mOutTmp.at(0).data(), this->mInTmp.at(0).data());
}

void WorlandIntegrator::input(const MatrixZ &in) const {
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
    }
  }
  Profiler::RegionStart<2>("applyImpl");
}

void WorlandIntegrator::input(const Matrix &in, const bool isEven) const {
  Matrix &inTmp = this->mInTmp.at(0);
  int start = 0;
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
  for (auto loc : *this->pLoc(isEven)) {
    inTmp.block(0, start, inTmp.rows(), std::get<2>(loc)) =
        in.block(0, std::get<1>(loc), in.rows(), std::get<2>(loc));
    start += std::get<2>(loc);
  }
  IWorlandBackend::transferDataToGPU(&bufferSolveRes, inTmp.data(), 0, 0,
                                     inTmp.size() * sizeof(MHDFloat),
                                     &pStream[currentStream]);
}

void WorlandIntegrator::input(const MatrixZ &in, const bool isEven,
                              const bool useReal) const {
  WorlandIntegrator::input(in);
}

void WorlandIntegrator::applyFft() const {
  for (int isEven = 0; isEven < 2; isEven++) {
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    this->setPlan(isEven);
    VkFFT::VkFFTLaunchParams launchParams = {};
    launchParams.buffer = (void **)&this->bufferSolveRes;
    mpApp->configuration.stream = &pStream[currentStream];

    VkFFT::VkFFTResult rr = VkFFTAppend(mpApp, -1, &launchParams);
  }
}

void WorlandIntegrator::partialForwardWorland(
    const int l, const int id, const int i0, const int start, const int cols,
    const std::vector<Matrix> &banded, const std::vector<void *> &bandedGPU,
    Matrix &out) const {
  if (i0 < l / 2) {
    void **GPUarr = (void **)malloc((l / 2 - i0 + 1) * sizeof(void *));
    for (int i = i0; i < l / 2; ++i) {
      GPUarr[i - i0] = bandedGPU.at(i);
    }
  
    int i_end = l / 2;
    this->applyPairCombined(out, id, start, cols, this->lSize(i0, l),
                            this->lSize(i_end, l), banded.at(i0),
                            (const void **)GPUarr);
    free(GPUarr); 
  }
}

void WorlandIntegrator::forwardWorland(const bool isEven0, const int id) const {
  // Reset current l
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    this->resetLocations(isEven, id);

    PfSolve_JW::PfSolveResult resSolve = PfSolve_JW::PFSOLVE_SUCCESS;
    Matrix &outTmp = this->mOutTmp.at(id);
    void *wTmpGPU = this->workTmpGPU(id);
    const std::vector<Matrix> &banded = this->banded(isEven);
    const std::vector<void *> &bandedGPU = this->bandedGPU(isEven);

    int i0 = 0;
    int cols = outTmp.cols();
    for (auto &loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc);
      cols = outTmp.cols() - start;
      if (l < 2) {
        {
          PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
          PfSolve_JW::PfSolve_MapKey_block mapKey = {};
          mapKey.size[0] = this->lSize(l);
          mapKey.size[1] = 2 * std::get<2>(loc);
          mapKey.type =
              (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION + RUNTIME_SCALEC +
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

            resSolve =
                PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
            resSolve =
                PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
            resSolve = PfSolve_JW::checkLibrary_block(&appLibrary, mapKey,
                                                      &tempAppCopy);
          }

          PfSolve_JW::PfSolveLaunchParams launchParams = {};
          launchParams.buffer = (void **)&wTmpGPU;
          launchParams.offsetM = 2 * start * outTmp.rows();
          launchParams.offsetSolution = 2 * start * outTmp.rows();
          launchParams.outputBuffer = (void **)&wTmpGPU;
          launchParams.inputBufferStride = outTmp.rows();
          launchParams.outputBufferStride = outTmp.rows();
          launchParams.scaleC = this->mFftScaling;
          tempAppCopy->configuration.stream = &pStream[currentStream];
          resSolve = PfSolveAppend(tempAppCopy, -1, &launchParams);
        }
        if (l == 0) {
          PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
          PfSolve_JW::PfSolve_MapKey_block mapKey = {};
          mapKey.size[0] = 1;
          mapKey.size[1] = 2 * std::get<2>(loc);
          mapKey.type =
              (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION + RUNTIME_SCALEC +
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

            configurationCopy.size[0] = mapKey.size[0];
            configurationCopy.size[1] = mapKey.size[1];
            configurationCopy.block = mapKey.type % 1000;
            configurationCopy.isOutputFormatted = 1;
            configurationCopy.jw_control_bitmask = mapKey.type / 1000;
            configurationCopy.doublePrecision = 1;
            configurationCopy.device = &this->device;
            configurationCopy.num_streams = 1;

            resSolve =
                PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
            resSolve =
                PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
            resSolve = PfSolve_JW::checkLibrary_block(&appLibrary, mapKey,
                                                      &tempAppCopy);
          }

          PfSolve_JW::PfSolveLaunchParams launchParams = {};
          launchParams.buffer = (void **)&wTmpGPU;
          launchParams.offsetM = 2 * start * outTmp.rows();
          launchParams.offsetSolution = 2 * start * outTmp.rows();
          launchParams.outputBuffer = (void **)&wTmpGPU;
          launchParams.inputBufferStride = outTmp.rows();
          launchParams.outputBufferStride = outTmp.rows();
          launchParams.scaleC = 1.0 / std::sqrt(2.0);
          tempAppCopy->configuration.stream = &pStream[currentStream];
          resSolve = PfSolveAppend(tempAppCopy, -1, &launchParams);
        }
      } else {
        if (i0 == 0) {
          PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
          PfSolve_JW::PfSolve_MapKey_block mapKey = {};
          mapKey.size[0] = this->lSize(0, l);
          mapKey.size[1] = 2 * cols;
          mapKey.type =
              (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION + RUNTIME_SCALEC +
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

            resSolve =
                PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
            resSolve =
                PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
            resSolve = PfSolve_JW::checkLibrary_block(&appLibrary, mapKey,
                                                      &tempAppCopy);
          }

          PfSolve_JW::PfSolveLaunchParams launchParams = {};
          launchParams.buffer = (void **)&wTmpGPU;
          launchParams.offsetM = 2 * start * outTmp.rows();
          launchParams.offsetSolution = 2 * start * outTmp.rows();
          launchParams.outputBuffer = (void **)&wTmpGPU;
          launchParams.inputBufferStride = outTmp.rows();
          launchParams.outputBufferStride = outTmp.rows();
          launchParams.scaleC = this->mFftScaling;
          tempAppCopy->configuration.stream = &pStream[currentStream];
          resSolve = PfSolveAppend(tempAppCopy, -1, &launchParams);
        }
        this->partialForwardWorland(l, id, i0, start, cols, banded, bandedGPU,
                                    outTmp);
        i0 = l / 2;
      }
      std::get<3>(loc) = this->lSize(l);
      start += std::get<2>(loc);
    }
  }
}

void WorlandIntegrator::lowerBeta(const MHDFloat alpha, const bool isEven0,
                                  const int id, const MHDFloat norm) const {
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &outTmp = this->mOutTmp.at(id);
    for (auto &loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc);
      int cols = std::get<2>(loc);
      Matrix V;
      this->buildShiftV(V, this->lSize(l), alpha, l - 0.5, norm);
      MHDFloat scale = (l == 1) ? 1.0 / std::sqrt(2.0) : 1;
      this->applyTriSolve(outTmp, id, start, cols, scale, V);
      std::get<4>(loc)--;
      start += cols;
    }
  }
}

void WorlandIntegrator::raiseBeta(const MHDFloat alpha, const bool isEven0,
                                  const int id, const MHDFloat norm) const {
  PfSolve_JW::PfSolveResult resSolve = PfSolve_JW::PFSOLVE_SUCCESS;
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &outTmp = this->mOutTmp.at(id);
    void *wTmpGPU = this->workTmpGPU(id);
    for (auto &loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc);
      int cols = std::get<2>(loc);
      if (l < 0) {
        PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
        PfSolve_JW::PfSolve_MapKey_block mapKey = {};
        mapKey.size[0] = outTmp.rows();
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
        launchParams.offsetM = 2 * start * mapKey.size[0];
        launchParams.offsetSolution = 2 * start * mapKey.size[0];
        launchParams.buffer = (void **)&wTmpGPU;
        launchParams.outputBuffer = (void **)&wTmpGPU;
        launchParams.inputZeropad[0] =
            prevOutputZeropad[currentStream][0][start];
        launchParams.inputZeropad[1] =
            prevOutputZeropad[currentStream][1][start];
        launchParams.outputZeropad[0] =
            prevOutputZeropad[currentStream][0][start];
        launchParams.outputZeropad[1] =
            prevOutputZeropad[currentStream][1][start];
        launchParams.inputBufferStride = mapKey.size[0];
        launchParams.outputBufferStride = mapKey.size[0];
        tempAppCopy->configuration.stream = &pStream[currentStream];
        resSolve = PfSolveAppend(tempAppCopy, -1, &launchParams);

        prevOutputZeropad[currentStream][0][start] =
            launchParams.outputZeropad[0];
        prevOutputZeropad[currentStream][1][start] =
            launchParams.outputZeropad[1];

      } else {
        Matrix V;
        this->buildShiftV(V, this->lSize(l), alpha, l + 0.5, norm);
        if (l == 0) {
          V(1, 0) *= std::sqrt(2.0);
        }
        this->applyTriProduct(outTmp, id, start, cols, V);
        std::get<4>(loc)++;
      }
      std::get<3>(loc)--;
      start += cols;
    }
  }
}

void WorlandIntegrator::lowerR2Beta(const MHDFloat alpha, const bool isEven0,
                                    const int id, const MHDFloat norm) const {
  PfSolve_JW::PfSolveResult resSolve = PfSolve_JW::PFSOLVE_SUCCESS;
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &outTmp = this->mOutTmp.at(id);
    void *wTmpGPU = this->workTmpGPU(id);
    for (auto &loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc);
      int cols = std::get<2>(loc);
      if (l < 0) {
        PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
        PfSolve_JW::PfSolve_MapKey_block mapKey = {};
        mapKey.size[0] = outTmp.rows();
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
        launchParams.offsetM = 2 * start * mapKey.size[0];
        launchParams.offsetSolution = 2 * start * mapKey.size[0];
        launchParams.buffer = (void **)&wTmpGPU;
        launchParams.outputBuffer = (void **)&wTmpGPU;
        launchParams.inputZeropad[0] =
            prevOutputZeropad[currentStream][0][start];
        launchParams.inputZeropad[1] =
            prevOutputZeropad[currentStream][1][start];
        launchParams.outputZeropad[0] =
            prevOutputZeropad[currentStream][0][start];
        launchParams.outputZeropad[1] =
            prevOutputZeropad[currentStream][1][start];
        launchParams.inputBufferStride = mapKey.size[0];
        launchParams.outputBufferStride = mapKey.size[0];
        tempAppCopy->configuration.stream = &pStream[currentStream];
        resSolve = PfSolveAppend(tempAppCopy, -1, &launchParams);

        prevOutputZeropad[currentStream][0][start] =
            launchParams.outputZeropad[0];
        prevOutputZeropad[currentStream][1][start] =
            launchParams.outputZeropad[1];
      } else {
        Matrix M;
        this->buildShiftM(M, this->lSize(l), alpha, l - 0.5, norm);
        this->applyTriProduct(outTmp, id, start, cols, M);
        std::get<4>(loc)--;
      }
      start += cols;
    }
  }
}

void WorlandIntegrator::raiseR2Beta(const MHDFloat alpha, const bool isEven0,
                                    const int id, const MHDFloat norm,
                                    const bool scaleL0) const {
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &outTmp = this->mOutTmp.at(id);
    for (auto &loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc);
      int cols = std::get<2>(loc);
      Matrix M;
      this->buildShiftM(M, this->lSize(l), alpha, l + 0.5, norm);
      MHDFloat scale = (scaleL0 && l == 0) ? std::sqrt(2.0) : 1;
      this->applyTriSolve(outTmp, id, start, cols, scale, M);
      std::get<4>(loc)++;
      start += cols;
    }
  }
}

void WorlandIntegrator::lowerAlpha(const MHDFloat alpha, const bool isEven0,
                                   const int id, const MHDFloat norm) const {
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &outTmp = this->mOutTmp.at(id);
    for (auto &loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc);
      int cols = std::get<2>(loc);
      Matrix U;
      this->buildShiftU(U, this->lSize(l), alpha, l - 0.5, norm);
      MHDFloat scale = (l == 1) ? 1.0 / std::sqrt(2.0) : 1;
      this->applyTriSolve(outTmp, id, start, cols, scale, U);
      start += cols;
    }
  }
}

void WorlandIntegrator::raiseAlpha(const MHDFloat alpha, const bool isEven0,
                                   const int id, const MHDFloat norm) const {
  throw std::logic_error("Raise alpha operator has not been tested!");
  PfSolve_JW::PfSolveResult resSolve = PfSolve_JW::PFSOLVE_SUCCESS;
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &outTmp = this->mOutTmp.at(id);
    void *wTmpGPU = this->workTmpGPU(id);
    for (auto &loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc);
      int cols = std::get<2>(loc);
      if (l < 0) {
        PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
        PfSolve_JW::PfSolve_MapKey_block mapKey = {};
        mapKey.size[0] = outTmp.rows();
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
        launchParams.offsetM = 2 * start * mapKey.size[0];
        launchParams.offsetSolution = 2 * start * mapKey.size[0];
        launchParams.buffer = (void **)&wTmpGPU;
        launchParams.outputBuffer = (void **)&wTmpGPU;
        launchParams.inputZeropad[0] =
            prevOutputZeropad[currentStream][0][start];
        launchParams.inputZeropad[1] =
            prevOutputZeropad[currentStream][1][start];
        launchParams.outputZeropad[0] =
            prevOutputZeropad[currentStream][0][start];
        launchParams.outputZeropad[1] =
            prevOutputZeropad[currentStream][1][start];
        launchParams.inputBufferStride = mapKey.size[0];
        launchParams.outputBufferStride = mapKey.size[0];
        tempAppCopy->configuration.stream = &pStream[currentStream];
        resSolve = PfSolveAppend(tempAppCopy, -1, &launchParams);

        prevOutputZeropad[currentStream][0][start] =
            launchParams.outputZeropad[0];
        prevOutputZeropad[currentStream][1][start] =
            launchParams.outputZeropad[1];
      } else {
        Matrix U;
        this->buildShiftU(U, this->lSize(l), alpha, l - 0.5, norm);
        this->applyTriProduct(outTmp, id, start, cols, U);
      }
      std::get<3>(loc)--;
      start += cols;
    }
  }
}

void WorlandIntegrator::applyI2(const bool isEven0, const int id) const {
  PfSolve_JW::PfSolveResult resSolve = PfSolve_JW::PFSOLVE_SUCCESS;
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &outTmp = this->mOutTmp.at(id);
    void *wTmpGPU = this->workTmpGPU(id);
    for (auto loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc);
      int cols = std::get<2>(loc);
      int r = this->lSize(l);
      Internal::MHDFloat a = -MHD_MP(0.5);
      Internal::MHDFloat b = -MHD_MP(0.5);
      ::QuICC::SparseSM::Worland::I2 spasm(r, r, a, b, l);

      Matrix bd;
      spasm.buildOp(bd);
      this->applyBandProduct(outTmp, id, start, cols, bd);

      start += cols;
    }
  }
}

void WorlandIntegrator::applyI4(const bool isEven0, const int id) const {
  PfSolve_JW::PfSolveResult resSolve = PfSolve_JW::PFSOLVE_SUCCESS;
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &outTmp = this->mOutTmp.at(id);
    void *wTmpGPU = this->workTmpGPU(id);
    for (auto loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc);
      int cols = std::get<2>(loc);
      int r = this->lSize(l);
      Internal::MHDFloat a = -MHD_MP(0.5);
      Internal::MHDFloat b = -MHD_MP(0.5);
      ::QuICC::SparseSM::Worland::I4 spasm(r, r, a, b, l);

      Matrix bd;
      spasm.buildOp(bd);
      this->applyBandProduct(outTmp, id, start, cols, bd);

      start += cols;
    }
  }
}

void WorlandIntegrator::output(MatrixZ &rOut) const {
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
      mapKey.size[0] = this->mSpecSize;
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
    rOut.block(0, s, this->mSpecSize, cols).setZero();
  }
  Profiler::RegionStop<2>("applyOperators");
}

void WorlandIntegrator::output(Matrix &rOut, const bool isEven) const {
  Matrix &outTmp = this->mOutTmp.at(0);
  int start = 0;
  IWorlandBackend::transferDataFromGPU(outTmp.data(), &bufferSolveRes, 0, 0,
                                       outTmp.size() * sizeof(MHDFloat),
                                       &pStream[currentStream]);

  for (auto loc : *this->pLoc(isEven)) {
    int cols = std::get<2>(loc);
    rOut.block(0, std::get<1>(loc), this->mSpecSize, cols) =
        outTmp.block(0, start, this->mSpecSize, cols);
    start += cols;
  }

  // Zero unused values
  for (auto loc : this->mZLoc) {
    int s = std::get<1>(loc);
    int cols = std::get<2>(loc);
    rOut.block(0, s, this->mSpecSize, cols).setZero();
  }
}

void WorlandIntegrator::output(MatrixZ &rOut, const bool isEven,
                               const bool useReal) const {
  WorlandIntegrator::output(rOut);
}

void WorlandIntegrator::initBanded(const bool isEven) const {
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

    int i0 = 0;
    for (auto &loc : *this->pLoc(isEven, 0)) {
      int l = std::get<4>(loc);
      if (l >= 2) {
#if (VKFFT_BACKEND == 1)
        int warpSize = 32;
#elif (VKFFT_BACKEND == 2)
        int warpSize = 64;
#endif
        int64_t registers_per_thread =
            (uint64_t)ceil(this->lSize(i0, l) / (double)warpSize);
        for (int i = i0; i < l / 2; ++i) {
          Matrix PS;
          void *tempBuffer;
          this->buildShiftPair(PS, true, i, this->lSize(i, l), J, normV, normM);

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
        i0 = l / 2;
      }
    }
  }
}

void WorlandIntegrator::freeBanded(const bool isEven) const {
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
