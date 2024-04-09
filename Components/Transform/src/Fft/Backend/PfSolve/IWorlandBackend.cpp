/**
 * @file IWorlandBackend.cpp
 * @brief Source of the interface for a generic FFTW based Worland integrator
 */

// External includes
//
#include <Eigen/src/misc/blas.h>
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Backend/PfSolve/IWorlandBackend.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace PfSolve {

const MHDFloat IWorlandBackend::UPPER_BANDED = 424242.424242;

const MHDFloat IWorlandBackend::LOWER_BANDED = -424242.424242;

const MHDFloat IWorlandBackend::LOWER_BANDED_PADDED = -424242.434343;

IWorlandBackend::IWorlandBackend() : mFlipped(false) {
  this->mELoc.push_back(LocationVector());
  this->mOLoc.push_back(LocationVector());
}

IWorlandBackend::~IWorlandBackend() {
  for (std::map<PfSolve_JW::PfSolve_MapKey_JonesWorland,
                PfSolve_JW::PfSolveApplication>::iterator it =
           appLibrary.mapJonesWorland.begin();
       it != appLibrary.mapJonesWorland.end(); ++it)
    PfSolve_JW::deletePfSolve(&it->second);
  for (std::map<PfSolve_JW::PfSolve_MapKey_dgbmv,
                PfSolve_JW::PfSolveApplication>::iterator it =
           appLibrary.mapDGBMV.begin();
       it != appLibrary.mapDGBMV.end(); ++it)
    PfSolve_JW::deletePfSolve(&it->second);
  for (std::map<PfSolve_JW::PfSolve_MapKey_block,
                PfSolve_JW::PfSolveApplication>::iterator it =
           appLibrary.mapBlock.begin();
       it != appLibrary.mapBlock.end(); ++it)
    PfSolve_JW::deletePfSolve(&it->second);
}
void IWorlandBackend::transferDataToGPU(void *gpu_buffer, void *cpu_arr,
                                        uint64_t gpuOffset, uint64_t cpuOffset,
                                        uint64_t transferSize, void *stream) {
#if (VKFFT_BACKEND == 1)
  cudaError_t res = cudaSuccess;
  void *buffer = ((void **)gpu_buffer)[0];
  cudaStream_t used_stream = ((cudaStream_t *)stream)[0];
  res = cudaMemcpyAsync((char *)buffer + gpuOffset, (char *)cpu_arr + cpuOffset,
                        transferSize, cudaMemcpyHostToDevice, used_stream);
#elif (VKFFT_BACKEND == 2)
  hipError_t res = hipSuccess;
  void *buffer = ((void **)gpu_buffer)[0];
  hipStream_t used_stream = ((hipStream_t *)stream)[0];
  res = hipMemcpyAsync((char *)buffer + gpuOffset, (char *)cpu_arr + cpuOffset,
                       transferSize, hipMemcpyHostToDevice, used_stream);
#endif
}

void IWorlandBackend::transferDataFromGPU(void *cpu_arr, void *gpu_buffer,
                                          uint64_t cpuOffset,
                                          uint64_t gpuOffset,
                                          uint64_t transferSize, void *stream) {
#if (VKFFT_BACKEND == 1)
  cudaError_t res = cudaSuccess;
  void *buffer = ((void **)gpu_buffer)[0];
  cudaStream_t used_stream = ((cudaStream_t *)stream)[0];
  res = cudaMemcpyAsync((char *)cpu_arr + cpuOffset, (char *)buffer + gpuOffset,
                        transferSize, cudaMemcpyDeviceToHost, used_stream);
#elif (VKFFT_BACKEND == 2)
  hipError_t res = hipSuccess;
  void *buffer = ((void **)gpu_buffer)[0];
  hipStream_t used_stream = ((hipStream_t *)stream)[0];
  res = hipMemcpyAsync((char *)cpu_arr + cpuOffset, (char *)buffer + gpuOffset,
                       transferSize, hipMemcpyDeviceToHost, used_stream);
#endif
}
void IWorlandBackend::initStorage(const int rows, const int cols, const int n,
                                  std::vector<Matrix> &s) const {
  assert(rows > 0);
  assert(cols > 0);
  assert(n > 0);

  // Input temporary storage
  s.reserve(n);
  s.push_back(Matrix(rows, cols));
  s.at(0).setZero(rows, cols);

  // Initialize additional temporary storage
  for (int i = 0; i < n - 1; i++) {
    s.push_back(s.at(0));
  }
}

void IWorlandBackend::addStorage(const int inExtras,
                                 const int outExtras) const {
  assert(inExtras >= 0);
  assert(outExtras >= 0);

  // Initialize additional input temporary storage
  this->mInTmp.reserve(this->mInTmp.size() + inExtras);
  for (int i = 0; i < inExtras; i++) {
    this->mInTmp.push_back(this->mInTmp.at(0));
  }

  // Initialize additional put temporary storage
  this->mOutTmp.reserve(this->mOutTmp.size() + outExtras);
  for (int i = 0; i < outExtras; i++) {
    this->mOutTmp.push_back(this->mOutTmp.at(0));
  }

  // Initialize additional block locations
  this->mELoc.reserve(this->mELoc.size() + std::max(inExtras, outExtras));
  this->mOLoc.reserve(this->mOLoc.size() + std::max(inExtras, outExtras));
  for (int i = 0; i < std::max(inExtras, outExtras); i++) {
    this->mELoc.push_back(LocationVector());
    this->mOLoc.push_back(LocationVector());
  }
}

void IWorlandBackend::init(const SetupType &setup, const int lshift,
                           const int extraN, const bool lshiftOnlyParity,
                           const bool alwaysZeroNegative) const {
  this->mFlipped = (std::abs(lshift) % 2 == 1);
  int lshift_parity = lshift;
  if (lshiftOnlyParity) {
    lshift_parity = 0;
  }
  int lshift_zero = 0;
  if (lshiftOnlyParity && alwaysZeroNegative) {
    lshift_zero = lshift;
  }

  this->mSpecSize = setup.specSize();

  // Compute even and odd block sizes
  this->mEBlockSize = 0;
  this->mOBlockSize = 0;
  int col = 0;
  for (int k = 0; k < setup.slowSize(); ++k) {
    int l = setup.slow(k);
    int loc_l = l + lshift_parity;
    LocationType loc = std::make_tuple(loc_l, col, setup.mult(k), -1, loc_l);
    if (this->inZFilter(l) || std::get<0>(loc) + lshift_zero < 0) {
      this->mZLoc.push_back(loc);
    } else {
      if ((l + lshift) % 2 == 0) {
        this->mEBlockSize += std::get<2>(loc);
        this->mELoc.at(0).push_back(loc);
      } else {
        this->mOBlockSize += std::get<2>(loc);
        this->mOLoc.at(0).push_back(loc);
      }
    }
    col += std::get<2>(loc);
  }

  this->mWExtra = extraN;
}

void IWorlandBackend::setZFilter(const std::set<int> &filter) const {
  this->mZFilter = filter;
}

bool IWorlandBackend::isPhysEven(const bool isSpecEven) const {
  return (isSpecEven ^ this->mFlipped);
}

bool IWorlandBackend::inZFilter(const int l) const {
  return (this->mZFilter.count(l) == 1);
}

int IWorlandBackend::reorder_indices_PfSolve(int i, int M_size, int warpSize,
                                             int used_registers) const {
  int ret = i / used_registers;
  if ((i % used_registers) < (M_size % used_registers)) {
    ret +=
        (i % used_registers) * ((M_size + used_registers - 1) / used_registers);
  } else {
    ret += (M_size % used_registers) *
           ((M_size + used_registers - 1) / used_registers);
    ret += ((i % used_registers) - (M_size % used_registers)) *
           (M_size / used_registers);
  }
  return ret;
}

void IWorlandBackend::setWSize() const {
  int lt = 0;
  int lf = 0;
  if (this->pLoc(true)->rbegin() != this->pLoc(true)->rend()) {
    lt = std::get<0>(*this->pLoc(true)->rbegin());
  }
  if (this->pLoc(false)->rbegin() != this->pLoc(false)->rend()) {
    lf = std::get<0>(*this->pLoc(false)->rbegin());
  }
  // Triangular truncation requirements: WSize = n_lmax_modes + lmax/2
  this->mWSize = this->mSpecSize + this->mWExtra + (std::max(lt, lf)) / 2;
}

IWorlandBackend::LocationVector *
IWorlandBackend::pLoc(const bool isEven, const unsigned int id) const {
  LocationVector *ptr;
  if (this->isPhysEven(isEven)) {
    assert(this->mELoc.size() > id);

    ptr = &this->mELoc.at(id);
  } else {
    assert(this->mOLoc.size() > id);

    ptr = &this->mOLoc.at(id);
  }

  return ptr;
}

void IWorlandBackend::resetLocations(const bool isEven, const int id) const {
  for (auto &loc : *this->pLoc(isEven, id)) {
    std::get<4>(loc) = std::get<0>(loc);
  }
}

Matrix &IWorlandBackend::workTmp(const unsigned int id) const {
  assert(this->mpWorkTmp);
  assert(id >= 0);
  assert(this->mpWorkTmp->size() > id);

  return this->mpWorkTmp->at(id);
}

void *IWorlandBackend::workTmpGPU(const unsigned int id) const {
  assert(id >= 0);
  assert(id < 2);
  if (id == 0)
    return this->bufferSolveRes;
  else
    return this->bufferSolveRes2;
}

void IWorlandBackend::initJ() const {
  MHDFloat alpha = -0.5;

  // Even modes
  if (this->mELoc.at(0).size() > 0) {
    // Max size is L + N + 1
    int lmax = std::get<0>(*this->mELoc.at(0).rbegin());
    int jSize = this->lSize(lmax) + lmax + 1;
    MHDFloat beta = 0.5;
    this->jacobiShiftMatrix(this->mJEven, jSize, alpha, beta);
  }

  // Odd Modes
  if (this->mOLoc.at(0).size() > 0) {
    // Max size is L + N + 1
    int lmax = std::get<0>(*this->mOLoc.at(0).rbegin());
    int jSize = this->lSize(lmax) + lmax + 1;
    MHDFloat beta = 1.5;
    this->jacobiShiftMatrix(this->mJOdd, jSize, alpha, beta);
  }
}

const Matrix &IWorlandBackend::J(const bool isEven) const {
  if (this->isPhysEven(isEven)) {
    return this->mJEven;
  } else {
    return this->mJOdd;
  }
}

const std::vector<Matrix> &IWorlandBackend::banded(const bool isEven) const {
  if (this->isPhysEven(isEven)) {
    return this->mEBanded;
  } else {
    return this->mOBanded;
  }
}

const std::vector<void *> &IWorlandBackend::bandedGPU(const bool isEven) const {
  if (this->isPhysEven(isEven)) {
    return this->mEBandedGPU;
  } else {
    return this->mOBandedGPU;
  }
}

void IWorlandBackend::makeUpper(SparseMatrix &M, const Matrix &Mp) const {
  int n = Mp.cols();
  M.resize(n, n);
  M.reserve(ArrayI::Constant(n, 2));
  M.insert(0, 0) = Mp(1, 0);
  for (int i = 1; i < n; ++i) {
    M.insert(i - 1, i) = Mp(0, i);
    M.insert(i, i) = Mp(1, i);
  }
  M.makeCompressed();
}

void IWorlandBackend::makeLower(SparseMatrix &M, const Matrix &Mp,
                                const bool isSquare) const {
  int n = Mp.cols();
  M.resize(n, n - static_cast<int>(!isSquare));
  M.reserve(ArrayI::Constant(n, 2));
  for (int i = 0; i < n - 1; ++i) {
    M.insert(i, i) = Mp(0, i);
    M.insert(i + 1, i) = Mp(1, i);
  }
  if (isSquare) {
    M.insert(n - 1, n - 1) = Mp(0, n - 1);
  }
  M.makeCompressed();
}

void IWorlandBackend::jacobiShiftMatrix(Matrix &rJ, const int size,
                                        const MHDFloat alpha,
                                        const MHDFloat beta) const {
  MHDFloat a1 = alpha + 1.0;
  MHDFloat ab = alpha + beta;
  MHDFloat ab1 = ab + 1.0;
  MHDFloat ab2 = ab + 2.0;
  Array n = Array::LinSpaced(size, 0.0, static_cast<MHDFloat>(size - 1));
  rJ.resize(size, 4);
  rJ.col(0) = (2.0 * (n.array() + beta) * (n.array() + ab)).sqrt();
  rJ.col(1) = ((2.0 * n.array() + ab) * (2.0 * n.array() + ab1)).sqrt();
  if (ab == 0.0) {
    rJ.col(0)(0) = std::sqrt(beta);
    rJ.col(1)(0) = 1.0;
  }
  rJ.col(2) = (2.0 * (n.array() + 1.0) * (n.array() + a1)).sqrt();
  rJ.col(3) = ((2.0 * n.array() + ab1) * (2.0 * n.array() + ab2)).sqrt();
  if (ab1 == 0.0) {
    rJ.col(0)(1) = 1.0;
    rJ.col(3)(0) = 1.0;
  }
}

void IWorlandBackend::buildShiftU(Matrix &U, const int size,
                                  const MHDFloat alpha, const MHDFloat beta,
                                  const MHDFloat norm) const {
  Matrix J;
  this->jacobiShiftMatrix(J, size, beta, alpha);
  this->buildShiftU(U, 0, size, J, norm);
}

void IWorlandBackend::buildShiftU(SparseMatrix &U, const int size,
                                  const MHDFloat alpha, const MHDFloat beta,
                                  const MHDFloat norm) const {
  Matrix Up;
  this->buildShiftU(Up, size, alpha, beta, norm);
  this->makeUpper(U, Up);
}

void IWorlandBackend::buildShiftU(Matrix &U, const int i, const int size,
                                  const Matrix &J, const MHDFloat norm) const {
  const int n = size;
  U.resize(2, n);

  // Safety asserts
  assert(n - 1 <= J.rows());
  assert(i + n - 1 <= J.rows());
  assert(2 * i + n <= J.rows());
  assert(i + n <= J.rows());

  Eigen::Map<const Array> Voff_num(J.data() + 2 * J.rows(), n - 1);
  Eigen::Map<const Array> Voff_den(J.data() + 3 * J.rows() + i, n - 1);
  U.row(0).rightCols(n - 1) = -Voff_num.array() / Voff_den.array();
  Eigen::Map<const Array> Vdia_num(J.data() + 2 * i, n);
  Eigen::Map<const Array> Vdia_den(J.data() + J.rows() + i, n);
  U.row(1) = Vdia_num.array() / Vdia_den.array();

  // Normalize
  if (norm != 1.0) {
    U *= norm;
  }

  // Set upper triangular flag value
  U(0, 0) = UPPER_BANDED;
}

void IWorlandBackend::buildShiftU(SparseMatrix &U, const int i, const int size,
                                  const Matrix &J, const MHDFloat norm) const {
  Matrix Up;
  this->buildShiftU(Up, i, size, J, norm);
  this->makeUpper(U, Up);
}

void IWorlandBackend::buildShiftV(Matrix &V, const int size,
                                  const MHDFloat alpha, const MHDFloat beta,
                                  const MHDFloat norm) const {
  Matrix J;
  this->jacobiShiftMatrix(J, size, alpha, beta);
  this->buildShiftV(V, 0, size, J, norm);
}

void IWorlandBackend::buildShiftV(SparseMatrix &V, const int size,
                                  const MHDFloat alpha, const MHDFloat beta,
                                  const MHDFloat norm) const {
  Matrix Vp;
  this->buildShiftV(Vp, size, alpha, beta, norm);
  this->makeUpper(V, Vp);
}

void IWorlandBackend::buildShiftV(Matrix &V, const int i, const int size,
                                  const Matrix &J, const MHDFloat norm) const {
  const int n = size;
  V.resize(2, n);

  // Safety asserts
  assert(n - 1 <= J.rows());
  assert(i + n - 1 <= J.rows());
  assert(2 * i + n <= J.rows());
  assert(i + n <= J.rows());

  Eigen::Map<const Array> Voff_num(J.data() + 2 * J.rows(), n - 1);
  Eigen::Map<const Array> Voff_den(J.data() + 3 * J.rows() + i, n - 1);
  V.row(0).rightCols(n - 1) = Voff_num.array() / Voff_den.array();
  Eigen::Map<const Array> Vdia_num(J.data() + 2 * i, n);
  Eigen::Map<const Array> Vdia_den(J.data() + J.rows() + i, n);
  V.row(1) = Vdia_num.array() / Vdia_den.array();

  // Normalize
  if (norm != 1.0) {
    V *= norm;
  }

  // Set upper triangular flag value
  V(0, 0) = UPPER_BANDED;
}

void IWorlandBackend::buildShiftV(SparseMatrix &V, const int i, const int size,
                                  const Matrix &J, const MHDFloat norm) const {
  Matrix Vp;
  this->buildShiftV(Vp, i, size, J, norm);
  this->makeUpper(V, Vp);
}

void IWorlandBackend::buildShiftM(Matrix &M, const int size,
                                  const MHDFloat alpha, const MHDFloat beta,
                                  const MHDFloat norm,
                                  const bool isSquare) const {
  Matrix J;
  this->jacobiShiftMatrix(J, size + 1, alpha, beta - 1);
  this->buildShiftM(M, 0, size, J, norm, isSquare);
}

void IWorlandBackend::buildShiftM(SparseMatrix &M, const int size,
                                  const MHDFloat alpha, const MHDFloat beta,
                                  const MHDFloat norm,
                                  const bool isSquare) const {
  Matrix Mp;
  this->buildShiftM(Mp, size, alpha, beta, norm, isSquare);
  this->makeLower(M, Mp, isSquare);
}

void IWorlandBackend::buildShiftM(Matrix &M, const int i, const int size,
                                  const Matrix &J, const MHDFloat norm,
                                  const bool isSquare) const {
  const int n = size;
  M.resize(2, n);
  M.setConstant(-4242.4242);

  // Safety asserts
  assert(2 * i + 1 + n <= J.rows());
  assert(i + n <= J.rows());
  assert(n - 1 <= J.rows());
  assert(i + 1 + n - 1 <= J.rows());

  Eigen::Map<const Array> Mdia_num(J.data() + 2 * i + 1, n);
  Eigen::Map<const Array> Mdia_den(J.data() + 3 * J.rows() + i, n);
  M.row(0) = Mdia_num.array() / Mdia_den.array();
  Eigen::Map<const Array> Moff_num(J.data() + 2 * J.rows(), n - 1);
  Eigen::Map<const Array> Moff_den(J.data() + 1 * J.rows() + i + 1, n - 1);
  M.row(1).leftCols(n - 1) = Moff_num.array() / Moff_den.array();

  // Normalize
  if (norm != 1.0) {
    M *= norm;
  }

  // Set lower triangular flag value
  if (isSquare) {
    M(1, n - 1) = LOWER_BANDED;
  } else {
    M(1, n - 1) = LOWER_BANDED_PADDED;
  }
}

void IWorlandBackend::buildShiftM(SparseMatrix &M, const int i, const int size,
                                  const Matrix &J, const MHDFloat norm,
                                  const bool isSquare) const {
  Matrix Mp;
  this->buildShiftM(Mp, i, size, J, norm);
  this->makeLower(M, Mp, isSquare);
}

void IWorlandBackend::buildShiftPair(Matrix &PS, const bool isVMOrder,
                                     const int i, const int size,
                                     const Matrix &J, const MHDFloat normV,
                                     const MHDFloat normM,
                                     const bool isSquare) const {
  Matrix V;
  this->buildShiftV(V, i, size, J, normV);

  Matrix M;
  this->buildShiftM(M, i, size, J, normM, isSquare);

  int cols = V.cols();
  PS.resize(3, cols);
  if (isVMOrder) {
    PS.row(0) = V.row(1);
    PS.block(1, 0, 1, cols - 1) = V.block(0, 1, 1, cols - 1);
    PS.block(2, 1, 1, cols - 1) = -M.block(1, 0, 1, cols - 1);
    PS.row(0).array() /= M.row(0).array();
    PS.row(1).array() /= M.row(0).array();
    PS.row(2).array() /= M.row(0).array();
    PS(2, 0) = V(0, 0);
  } else {
    PS.block(0, 1, 1, cols - 1) = M.block(1, 0, 1, cols - 1);
    PS.row(1) = M.row(0);
    PS.block(2, 0, 1, cols - 1) = -V.block(0, 1, 1, cols - 1);
    PS.row(0).array() /= V.row(1).array();
    PS.row(1).array() /= V.row(1).array();
    PS.row(2).array() /= V.row(1).array();
    PS(2, cols - 1) = M(1, cols - 1);
  }
}

void IWorlandBackend::buildShiftPair(SparseMatrix &V, SparseMatrix &M,
                                     const int i, const int size,
                                     const Matrix &J, const MHDFloat normV,
                                     const MHDFloat normM,
                                     const bool isSquare) const {
  this->buildShiftV(V, i, size, J, normV);

  this->buildShiftM(M, i, size, J, normM, isSquare);
}

void IWorlandBackend::setPlan(const bool isEven) const {
  if (this->isPhysEven(isEven)) {
    this->mpApp = &this->mEvenApp;
  } else {
    this->mpApp = &this->mOddApp;
  }
}

int IWorlandBackend::blockSize(const bool isEven) const {
  if (this->isPhysEven(isEven)) {
    return this->mEBlockSize;
  } else {
    return this->mOBlockSize;
  }
}

void IWorlandBackend::applyTriProduct(Matrix &out, const int start,
                                      const int cols,
                                      const SparseMatrix &A) const {
  int r = A.rows();
  int c = A.cols();
  out.block(0, start, r, cols) = A * out.block(0, start, c, cols);
}

void IWorlandBackend::applyTriProduct(Matrix &out, int id, const int start,
                                      const int cols, Matrix &A) const {
  int r = A.cols();
  char UPLO;
  void *wTmpGPU = this->workTmpGPU(id);
  PfSolve_JW::PfSolveResult resFFT = PfSolve_JW::PFSOLVE_SUCCESS;
#if (VKFFT_BACKEND == 1)
  CUresult res = CUDA_SUCCESS;
  cudaError_t res3 = cudaSuccess;
#elif (VKFFT_BACKEND == 2)
  hipError_t res = hipSuccess;
#endif
  if (A(0, 0) == UPPER_BANDED) {
    UPLO = 'U';
  } else {
    assert(A(1, A.cols() - 1) == LOWER_BANDED ||
           A(1, A.cols() - 1) == LOWER_BANDED_PADDED);
    UPLO = 'L';
  }

  Matrix C = A;

  if (A(0, 0) == UPPER_BANDED) {
    C.block(1, 0, 1, A.cols() - 1) = A.block(0, 1, 1, A.cols() - 1);
    C.block(0, 0, 1, A.cols()) = A.block(1, 0, 1, A.cols());
  } else {
    C.block(1, 0, 1, A.cols()) = A.block(0, 0, 1, A.cols());
    C.block(0, 1, 1, A.cols() - 1) = A.block(1, 0, 1, A.cols() - 1);
  }
  Matrix B2 = C.transpose();
  double *bb2 = B2.data();

  Matrix B = B2;
  double *bb = B.data();
#if (VKFFT_BACKEND == 1)
  int warpSize = 32;
#elif (VKFFT_BACKEND == 2)
  int warpSize = 64;
#endif
  int registers_per_thread = (B.rows() + warpSize - 1) / warpSize;
  for (int j = 0; j < B.rows(); j++) {
    bb[reorder_indices_PfSolve(j, B.rows(), warpSize, registers_per_thread) +
       0 * B.rows()] = bb2[j + 0 * B.rows()];
    bb[reorder_indices_PfSolve(j, B.rows(), warpSize, registers_per_thread) +
       1 * B.rows()] = bb2[j + 1 * B.rows()];
  }
  IWorlandBackend::transferDataToGPU(&bufferTemp, B.data(), 0, 0,
                                     2 * B.rows() * sizeof(MHDFloat),
                                     &pStream[currentStream]);

  char TRANS = 'N';
  char DIAG = 'N';
  int K = 1;
  int LDA = 2;
  int INCX = 1;

  PfSolve_JW::PfSolveApplication *tempSolveApp = 0;
  PfSolve_JW::PfSolve_MapKey_JonesWorland mapKey = {};
  mapKey.size[0] = B.rows();
  mapKey.upperBanded = (A(0, 0) == UPPER_BANDED);
  mapKey.type =
      (RUNTIME_OFFSETSOLUTION + RUNTIME_INPUTZEROPAD + RUNTIME_OUTPUTZEROPAD +
       RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE) *
          1000 +
      JWT_USE_PARALLEL_THOMAS + JWT_DISABLE_TRIDIAGONAL_SOLVE;
  PfSolve_JW::checkLibrary_JonesWorland(&appLibrary, mapKey, &tempSolveApp);
  if (!tempSolveApp) {
    PfSolve_JW::PfSolveApplication appSolve = {};
    PfSolve_JW::PfSolveConfiguration configurationSolve = {};
    configurationSolve.size[0] = mapKey.size[0];
    configurationSolve.size[1] = cols;
    configurationSolve.scaleC = 1;

    configurationSolve.jw_type = mapKey.type % 1000;
    configurationSolve.upperBanded = mapKey.upperBanded;
    configurationSolve.jw_control_bitmask = mapKey.type / 1000;
    configurationSolve.doublePrecision = 1;
    configurationSolve.isOutputFormatted = 1;
    configurationSolve.device = &this->device;
    configurationSolve.num_streams = 1;
    resFFT = PfSolve_JW::initializePfSolve(&appSolve, configurationSolve);
    resFFT =
        PfSolve_JW::addToLibrary_JonesWorland(&appLibrary, mapKey, &appSolve);
    PfSolve_JW::checkLibrary_JonesWorland(&appLibrary, mapKey, &tempSolveApp);
  }
  PfSolve_JW::PfSolveLaunchParams launchParams = {};
  launchParams.buffer = (void **)&this->bufferTemp;
  launchParams.outputBuffer = (void **)&wTmpGPU;
  launchParams.inputZeropad[0] = prevOutputZeropad[currentStream][0][start];
  launchParams.inputZeropad[1] = prevOutputZeropad[currentStream][1][start];
  launchParams.outputZeropad[0] = 0;
  launchParams.outputZeropad[1] = mapKey.size[0];
  launchParams.offsetSolution = 2 * start * out.rows();
  launchParams.outputBufferStride = out.rows();
  launchParams.inputBufferStride = mapKey.size[0];
  launchParams.batchSize = 2 * cols;
  tempSolveApp->configuration.stream = &pStream[currentStream];
  resFFT = PfSolve_JW::PfSolveAppend(tempSolveApp, -1, &launchParams);
  prevOutputZeropad[currentStream][0][start] = launchParams.outputZeropad[0];
  prevOutputZeropad[currentStream][1][start] = launchParams.outputZeropad[1];
}

void IWorlandBackend::applyBandProduct(Matrix &out, const int start,
                                       const int cols,
                                       const SparseMatrix &A) const {
  int r = A.rows();
  int c = A.cols();
  out.block(0, start, r, cols) = A * out.block(0, start, c, cols);
}

void IWorlandBackend::applyBandProduct(Matrix &out, int id, const int start,
                                       const int cols, Matrix &A) const {
  int r = A.cols();
  Matrix tmp = out.block(0, start, r, cols);
  int KL = A.bottomRightCorner(1, 1)(0, 0);
  int KU = A.topLeftCorner(1, 1)(0, 0);
  void *wTmpGPU = this->workTmpGPU(id);

  Matrix B = A.transpose();
  double *bb = B.data();
  char TRANS = 'N';
  double ALPHA = 1.0;
  int LDA = KL + KU + 1;
  int INCX = 1;
  double BETA = 0.0;

  PfSolve_JW::PfSolveResult resFFT = PfSolve_JW::PFSOLVE_SUCCESS;
  IWorlandBackend::transferDataToGPU(&bufferTemp, B.data(), 0, 0,
                                     LDA * B.rows() * sizeof(MHDFloat),
                                     &pStream[currentStream]);

  PfSolve_JW::PfSolveApplication *tempSolveApp = 0;
  PfSolve_JW::PfSolve_MapKey_dgbmv mapKey = {};
  mapKey.size[0] = B.rows();
  mapKey.size[1] = 2 * cols;
  mapKey.LDA = LDA;
  mapKey.KU = KU;
  mapKey.KL = KL;

  PfSolve_JW::checkLibrary_dgbmv(&appLibrary, mapKey, &tempSolveApp);
  if (!tempSolveApp) {
    PfSolve_JW::PfSolveApplication appSolve = {};
    PfSolve_JW::PfSolveConfiguration configurationSolve = {};
    configurationSolve.size[0] = mapKey.size[0];
    configurationSolve.size[1] = mapKey.size[1];
    configurationSolve.jw_control_bitmask =
        (RUNTIME_OFFSETSOLUTION + RUNTIME_INPUTZEROPAD + RUNTIME_OUTPUTZEROPAD +
         RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE);
    configurationSolve.scaleC = 1;
    configurationSolve.LDA = mapKey.LDA;
    configurationSolve.KU = mapKey.KU;
    configurationSolve.KL = mapKey.KL;
    configurationSolve.doublePrecision = 1;
    configurationSolve.isOutputFormatted = 1;
    configurationSolve.device = &this->device;
    configurationSolve.num_streams = 1;

    resFFT = PfSolve_JW::initializePfSolve(&appSolve, configurationSolve);
    resFFT = PfSolve_JW::addToLibrary_dgbmv(&appLibrary, mapKey, &appSolve);
    PfSolve_JW::checkLibrary_dgbmv(&appLibrary, mapKey, &tempSolveApp);
  }
  PfSolve_JW::PfSolveLaunchParams launchParams = {};
  launchParams.inputZeropad[0] = prevOutputZeropad[currentStream][0][start];
  launchParams.inputZeropad[1] = prevOutputZeropad[currentStream][1][start];
  launchParams.outputZeropad[0] = 0;
  launchParams.outputZeropad[1] = mapKey.size[0];
  launchParams.outputBufferStride = out.rows();
  launchParams.inputBufferStride = mapKey.size[0];
  launchParams.offsetSolution = 2 * start * out.rows();
  launchParams.buffer = (void **)&this->bufferTemp;
  launchParams.outputBuffer = (void **)&wTmpGPU;
  tempSolveApp->configuration.stream = &pStream[currentStream];
  resFFT = PfSolve_JW::PfSolveAppend(tempSolveApp, -1, &launchParams);
  prevOutputZeropad[currentStream][0][start] = launchParams.outputZeropad[0];
  prevOutputZeropad[currentStream][1][start] = launchParams.outputZeropad[1];
}

void IWorlandBackend::applyTriSolve(Matrix &out, const int start,
                                    const int cols,
                                    const SparseMatrix &A) const {
  int r = A.cols();
  Eigen::SparseLU<SparseMatrix> solver;
  solver.compute(A);
  Matrix tmp = out.block(0, start, r, cols);
  Matrix tmp2 = solver.solve(tmp);
  out.block(0, start, r, cols) = tmp2;
}

void IWorlandBackend::applyTriSolve(Matrix &out, int id, const int start,
                                    const int cols, const MHDFloat scale,
                                    Matrix &A) const {
  char UPLO;
  void *wTmpGPU = this->workTmpGPU(id);
  PfSolve_JW::PfSolveResult resFFT = PfSolve_JW::PFSOLVE_SUCCESS;
#if (VKFFT_BACKEND == 1)
  CUresult res = CUDA_SUCCESS;
  cudaError_t res3 = cudaSuccess;
#elif (VKFFT_BACKEND == 2)
  hipError_t res = hipSuccess;
#endif

  if (A(0, 0) == UPPER_BANDED) {
    UPLO = 'U';
  } else {
    assert(A(1, A.cols() - 1) == LOWER_BANDED);
    UPLO = 'L';
  }

  Matrix C = A;
  double *X2 = out.data() + (start)*out.rows();
  if (A(0, 0) == UPPER_BANDED) {
    C.block(1, 0, 1, A.cols() - 1).array() =
        A.block(0, 1, 1, A.cols() - 1).array() /
        A.block(1, 0, 1, A.cols() - 1).array();
    C.block(0, 0, 1, A.cols()).array() = 1 / A.block(1, 0, 1, A.cols()).array();
  } else {
    C.block(1, 0, 1, A.cols()).array() = 1 / A.block(0, 0, 1, A.cols()).array();
    C.block(0, 1, 1, A.cols() - 1).array() =
        A.block(1, 0, 1, A.cols() - 1).array() /
        A.block(0, 1, 1, A.cols() - 1).array();
  }

  Matrix B2 = C.transpose();
  double *bb2 = B2.data();

  Matrix B = B2;
  double *bb = B.data();
#if (VKFFT_BACKEND == 1)
  int warpSize = 32;
#elif (VKFFT_BACKEND == 2)
  int warpSize = 64;
#endif
  int registers_per_thread = (B.rows() + warpSize - 1) / warpSize;

  for (int j = 0; j < B.rows(); j++) {
    bb[reorder_indices_PfSolve(j, B.rows(), warpSize, registers_per_thread) +
       0 * B.rows()] = bb2[j + 0 * B.rows()];
    bb[reorder_indices_PfSolve(j, B.rows(), warpSize, registers_per_thread) +
       1 * B.rows()] = bb2[j + 1 * B.rows()];
  }
  IWorlandBackend::transferDataToGPU(&bufferTemp, B.data(), 0, 0,
                                     2 * B.rows() * sizeof(MHDFloat),
                                     &pStream[currentStream]);

  char TRANS = 'N';
  char DIAG = 'N';
  int K = 1;
  int LDA = 2;
  int INCX = 1;
  int r = A.cols();

  PfSolve_JW::PfSolveApplication *tempSolveApp = 0;
  PfSolve_JW::PfSolve_MapKey_JonesWorland mapKey = {};
  mapKey.size[0] = B.rows();
  mapKey.upperBanded = (A(0, 0) != UPPER_BANDED);
  mapKey.type = (RUNTIME_OFFSETSOLUTION + RUNTIME_INPUTZEROPAD +
                 RUNTIME_OUTPUTZEROPAD + RUNTIME_SCALEC +
                 RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE) *
                    1000 +
                JWT_USE_PARALLEL_THOMAS + JWT_DIAGONAL_MATVECMUL;
  PfSolve_JW::checkLibrary_JonesWorland(&appLibrary, mapKey, &tempSolveApp);
  if (!tempSolveApp) {
    PfSolve_JW::PfSolveApplication appSolve = {};
    PfSolve_JW::PfSolveConfiguration configurationSolve = {};
    configurationSolve.size[0] = mapKey.size[0];
    configurationSolve.jw_type = mapKey.type % 1000;
    configurationSolve.upperBanded = mapKey.upperBanded;
    configurationSolve.jw_control_bitmask = mapKey.type / 1000;
    configurationSolve.doublePrecision = 1;
    configurationSolve.isOutputFormatted = 1;
    configurationSolve.device = &this->device;
    configurationSolve.num_streams = 1;

    resFFT = PfSolve_JW::initializePfSolve(&appSolve, configurationSolve);
    resFFT =
        PfSolve_JW::addToLibrary_JonesWorland(&appLibrary, mapKey, &appSolve);
    PfSolve_JW::checkLibrary_JonesWorland(&appLibrary, mapKey, &tempSolveApp);
  }

  PfSolve_JW::PfSolveLaunchParams launchParams = {};
  launchParams.buffer = (void **)&this->bufferTemp;
  launchParams.outputBuffer = (void **)&wTmpGPU;
  launchParams.inputZeropad[0] = prevOutputZeropad[currentStream][0][start];
  launchParams.inputZeropad[1] = prevOutputZeropad[currentStream][1][start];
  launchParams.outputZeropad[0] = 0;
  launchParams.outputZeropad[1] = mapKey.size[0];
  launchParams.outputBufferStride = out.rows();
  launchParams.inputBufferStride = mapKey.size[0];
  launchParams.batchSize = 2 * cols;
  launchParams.offsetSolution = 2 * start * out.rows();
  launchParams.scaleC = scale;
  launchParams.offsetM = 0;
  launchParams.offsetV = 0;
  tempSolveApp->configuration.stream = &pStream[currentStream];
  resFFT = PfSolve_JW::PfSolveAppend(tempSolveApp, -1, &launchParams);
  prevOutputZeropad[currentStream][0][start] = launchParams.outputZeropad[0];
  prevOutputZeropad[currentStream][1][start] = launchParams.outputZeropad[1];
}

void IWorlandBackend::applyPair(Matrix &out, const int start, const int cols,
                                const int rows, const SparseMatrix &P,
                                const SparseMatrix &S) const {
  int mr = P.rows();
  int mc = P.cols();

  Matrix tmp = P * out.block(0, start, mc, cols);

  Eigen::SparseLU<SparseMatrix> solver;
  solver.compute(S);
  Matrix tmp2 = solver.solve(tmp);
  out.block(0, start, mr, cols) = tmp2;
}

void IWorlandBackend::applyPair(Matrix &out, const int id, const int start,
                                const int cols, const int rows,
                                const Matrix &PS, const void *PSGPU) const {
  char UPLO;
  void *wTmpGPU = this->workTmpGPU(id);
  PfSolve_JW::PfSolveResult resFFT = PfSolve_JW::PFSOLVE_SUCCESS;
#if (VKFFT_BACKEND == 1)
  CUresult res = CUDA_SUCCESS;
  cudaError_t res3 = cudaSuccess;
#elif (VKFFT_BACKEND == 2)
  hipError_t res = hipSuccess;
#endif
  PfSolve_JW::PfSolveApplication *tempSolveApp = 0;
  PfSolve_JW::PfSolve_MapKey_JonesWorland mapKey = {};
  mapKey.size[0] = PS.cols();
  mapKey.upperBanded = (PS(2, 0) == UPPER_BANDED);
  mapKey.type =
      (RUNTIME_OFFSETSOLUTION + RUNTIME_INPUTZEROPAD + RUNTIME_OUTPUTZEROPAD +
       RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE) *
          1000 + JWT_USE_PARALLEL_THOMAS;
  PfSolve_JW::checkLibrary_JonesWorland(&appLibrary, mapKey, &tempSolveApp);
  if (!tempSolveApp) {
    PfSolve_JW::PfSolveApplication appSolve = {};
    PfSolve_JW::PfSolveConfiguration configurationSolve = {};
    configurationSolve.size[0] = mapKey.size[0];
    configurationSolve.offsetV = 2 * mapKey.size[0];
    configurationSolve.scaleC = 1;
    configurationSolve.jw_type = mapKey.type % 1000;
    configurationSolve.upperBanded = mapKey.upperBanded;
    configurationSolve.jw_control_bitmask = mapKey.type / 1000;
    configurationSolve.doublePrecision = 1;
    configurationSolve.isOutputFormatted = 1;
    configurationSolve.device = &this->device;
    configurationSolve.num_streams = 1;

    resFFT = PfSolve_JW::initializePfSolve(&appSolve, configurationSolve);
    resFFT =
        PfSolve_JW::addToLibrary_JonesWorland(&appLibrary, mapKey, &appSolve);
    PfSolve_JW::checkLibrary_JonesWorland(&appLibrary, mapKey, &tempSolveApp);
  }
  PfSolve_JW::PfSolveLaunchParams launchParams = {};
  launchParams.buffer = (void **)&PSGPU;
  launchParams.outputBuffer = (void **)&wTmpGPU;
  launchParams.inputZeropad[0] = 0;
  launchParams.inputZeropad[1] = mapKey.size[0];
  if (mapKey.upperBanded)
    launchParams.inputZeropad[1]++;
  else
    launchParams.inputZeropad[1]--;
  launchParams.outputZeropad[0] = 0;
  launchParams.outputZeropad[1] = mapKey.size[0];
  launchParams.offsetSolution = 2 * start * out.rows();
  launchParams.offsetM = 0;
  launchParams.offsetV = 0;
  launchParams.outputBufferStride = out.rows();
  launchParams.inputBufferStride = mapKey.size[0];
  launchParams.batchSize = 2 * cols;

  tempSolveApp->configuration.stream = &pStream[currentStream];
  resFFT = PfSolve_JW::PfSolveAppend(tempSolveApp, -1, &launchParams);
  prevOutputZeropad[currentStream][0][start] = launchParams.outputZeropad[0];
  prevOutputZeropad[currentStream][1][start] = launchParams.outputZeropad[1];
}

void IWorlandBackend::applyPairCombined(Matrix &out, const int id,
                                        const int start, const int cols,
                                        const int rows_start,
                                        const int rows_end, const Matrix &PS,
                                        const void **bandedGPU) const {
  char UPLO;
  void *wTmpGPU = this->workTmpGPU(id);
  PfSolve_JW::PfSolveResult resFFT = PfSolve_JW::PFSOLVE_SUCCESS;
#if (VKFFT_BACKEND == 1)
  CUresult res = CUDA_SUCCESS;
  cudaError_t res3 = cudaSuccess;
#elif (VKFFT_BACKEND == 2)
  hipError_t res = hipSuccess;
#endif
  PfSolve_JW::PfSolveApplication *tempSolveApp = 0;
  PfSolve_JW::PfSolve_MapKey_JonesWorland mapKey = {};
  mapKey.size[0] = rows_start * 100000 + rows_end;
  mapKey.upperBanded = (PS(2, 0) == UPPER_BANDED);
  mapKey.type =
      (RUNTIME_OFFSETSOLUTION + RUNTIME_INPUTZEROPAD + RUNTIME_OUTPUTZEROPAD +
       RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE) *
          1000 + JWT_USE_PARALLEL_THOMAS;
  PfSolve_JW::checkLibrary_JonesWorland(&appLibrary, mapKey, &tempSolveApp);
  if (!tempSolveApp) {
    PfSolve_JW::PfSolveApplication appSolve = {};
    PfSolve_JW::PfSolveConfiguration configurationSolve = {};
    configurationSolve.size[0] = mapKey.size[0] / 100000;
    configurationSolve.useMultipleInputBuffers = abs(rows_start - rows_end);
    configurationSolve.numConsecutiveJWIterations = abs(rows_start - rows_end);
    configurationSolve.scaleC = 1;
    configurationSolve.jw_type = mapKey.type % 1000;
    configurationSolve.upperBanded = mapKey.upperBanded;
    configurationSolve.jw_control_bitmask = mapKey.type / 1000;
    configurationSolve.doublePrecision = 1;
    configurationSolve.isOutputFormatted = 1;
    configurationSolve.device = &this->device;
    configurationSolve.num_streams = 1;

    resFFT = PfSolve_JW::initializePfSolve(&appSolve, configurationSolve);
    resFFT =
        PfSolve_JW::addToLibrary_JonesWorland(&appLibrary, mapKey, &appSolve);
    PfSolve_JW::checkLibrary_JonesWorland(&appLibrary, mapKey, &tempSolveApp);
  }
  PfSolve_JW::PfSolveLaunchParams launchParams = {};
  launchParams.buffer = (void **)bandedGPU;
  launchParams.outputBuffer = (void **)&wTmpGPU;
  launchParams.inputZeropad[0] = 0;
  launchParams.inputZeropad[1] = mapKey.size[0] / 100000;
  if (mapKey.upperBanded)
    launchParams.inputZeropad[1]++;
  else
    launchParams.inputZeropad[1]--;
  launchParams.outputZeropad[0] = 0;
  launchParams.outputZeropad[1] = mapKey.size[0] / 100000;
  if (mapKey.upperBanded)
    launchParams.outputZeropad[1] -= (abs(rows_start - rows_end) - 1);
  else
    launchParams.outputZeropad[1] += (abs(rows_start - rows_end) - 1);
  launchParams.offsetSolution = 2 * start * out.rows();
  launchParams.offsetM = 0;
  launchParams.offsetV = 0;
  launchParams.outputBufferStride = out.rows();
  launchParams.inputBufferStride = mapKey.size[0] / 100000;
  launchParams.batchSize = 2 * cols;

  tempSolveApp->configuration.stream = &pStream[currentStream];
  resFFT = PfSolve_JW::PfSolveAppend(tempSolveApp, -1, &launchParams);
  prevOutputZeropad[currentStream][0][start] = launchParams.outputZeropad[0];
  prevOutputZeropad[currentStream][1][start] = launchParams.outputZeropad[1];
}

void IWorlandBackend::io(MHDFloat *out, const MHDFloat *in) const {
  this->mpOut = out;
  this->mpIn = in;
}

void IWorlandBackend::scaleC(const MHDFloat c, const bool isEven0,
                             const unsigned int id) const {
  Matrix &wTmp = this->workTmp(id);

  PfSolve_JW::PfSolveResult resFFT = PfSolve_JW::PFSOLVE_SUCCESS;
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    void *wTmpGPU = this->workTmpGPU(id);
    for (auto loc : *this->pLoc(isEven, id)) {
      int cols = std::get<2>(loc);
      start += cols;
    }
    PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
    PfSolve_JW::PfSolve_MapKey_block mapKey = {};
    mapKey.size[0] = wTmp.rows();
    mapKey.size[1] = 2 * start;
    mapKey.type = (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION + RUNTIME_SCALEC +
                   RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE) *
                      1000 +
                  BLOCK_SCALEC + BLOCK_READ_REAL_WRITE_REAL;
    resFFT = PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
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

      resFFT = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
      resFFT = PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
      resFFT =
          PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
    }

    PfSolve_JW::PfSolveLaunchParams launchParams = {};
    launchParams.offsetM = 0;
    launchParams.offsetSolution = 0;
    launchParams.scaleC = c;
    launchParams.inputBufferStride = wTmp.rows();
    launchParams.outputBufferStride = wTmp.rows();
    launchParams.buffer = (void **)&wTmpGPU;
    launchParams.outputBuffer = (void **)&wTmpGPU;

    tempAppCopy->configuration.stream = &pStream[currentStream];
    resFFT = PfSolve_JW::PfSolveAppend(tempAppCopy, -1, &launchParams);
  }
}

void IWorlandBackend::scaleALPY(const MHDFloat a, const MHDFloat y,
                                const bool isEven0, const int lshift,
                                const unsigned int id) const {
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &wTmp = this->workTmp(id);
    void *wTmpGPU = this->workTmpGPU(id);

    PfSolve_JW::PfSolveResult resFFT = PfSolve_JW::PFSOLVE_SUCCESS;
    for (auto loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc) + lshift;
      int cols = std::get<2>(loc);
      int rows = std::get<3>(loc);
      PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
      PfSolve_JW::PfSolve_MapKey_block mapKey = {};
      mapKey.size[0] = rows;
      mapKey.size[1] = 2 * cols;
      mapKey.type =
          (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION + RUNTIME_INPUTZEROPAD +
           RUNTIME_OUTPUTZEROPAD + RUNTIME_SCALEC + RUNTIME_INPUTBUFFERSTRIDE +
           RUNTIME_OUTPUTBUFFERSTRIDE) *
              1000 +
          BLOCK_SCALEC + BLOCK_READ_REAL_WRITE_REAL;

      resFFT =
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

        resFFT = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
        resFFT = PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
        resFFT =
            PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      }
      PfSolve_JW::PfSolveLaunchParams launchParams = {};
      launchParams.offsetM = 2 * start * wTmp.rows();
      launchParams.offsetSolution = 2 * start * wTmp.rows();
      launchParams.scaleC = a * l + y;
      launchParams.inputZeropad[0] = prevOutputZeropad[currentStream][0][start];
      launchParams.inputZeropad[1] = prevOutputZeropad[currentStream][1][start];
      launchParams.outputZeropad[0] =
          prevOutputZeropad[currentStream][0][start];
      launchParams.outputZeropad[1] =
          prevOutputZeropad[currentStream][1][start];
      launchParams.inputBufferStride = wTmp.rows();
      launchParams.outputBufferStride = wTmp.rows();
      launchParams.buffer = (void **)&wTmpGPU;
      launchParams.outputBuffer = (void **)&wTmpGPU;

      tempAppCopy->configuration.stream = &pStream[currentStream];
      resFFT = PfSolve_JW::PfSolveAppend(tempAppCopy, -1, &launchParams);

      prevOutputZeropad[currentStream][0][start] =
          launchParams.outputZeropad[0];
      prevOutputZeropad[currentStream][1][start] =
          launchParams.outputZeropad[1];
      start += cols;
    }
  }
}

void IWorlandBackend::scaleD(const bool isEven0, const int lshift,
                             const unsigned int id) const {
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    const MHDFloat alpha = -0.5;
    Matrix &wTmp = this->workTmp(id);
    void *wTmpGPU = this->workTmpGPU(id);

    PfSolve_JW::PfSolveResult resFFT = PfSolve_JW::PFSOLVE_SUCCESS;
    for (auto loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc) + lshift;
      int cols = std::get<2>(loc);
      int rows = std::get<3>(loc);
      PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
      PfSolve_JW::PfSolve_MapKey_block mapKey = {};
      mapKey.size[0] = rows;
      mapKey.size[1] = 2 * cols;
      mapKey.type =
          (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION + RUNTIME_INPUTZEROPAD +
           RUNTIME_OUTPUTZEROPAD + RUNTIME_SCALEC + RUNTIME_INPUTBUFFERSTRIDE +
           RUNTIME_OUTPUTBUFFERSTRIDE) *
              1000 +
          BLOCK_SCALED + BLOCK_READ_REAL_WRITE_REAL;
      mapKey.lshift = l;
      resFFT =
          PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      if (!tempAppCopy) {
        PfSolve_JW::PfSolveApplication appCopy = {};
        PfSolve_JW::PfSolveConfiguration configurationCopy = {};

        configurationCopy.size[0] = mapKey.size[0];
        configurationCopy.size[1] = mapKey.size[1];
        configurationCopy.block = mapKey.type % 1000;
        configurationCopy.isOutputFormatted = 1;
        configurationCopy.lshift = l;
        configurationCopy.jw_control_bitmask = mapKey.type / 1000;
        configurationCopy.doublePrecision = 1;
        configurationCopy.device = &this->device;
        configurationCopy.num_streams = 1;

        resFFT = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
        resFFT = PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
        resFFT =
            PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      }
      PfSolve_JW::PfSolveLaunchParams launchParams = {};
      launchParams.offsetM = 2 * start * wTmp.rows();
      launchParams.offsetSolution = 2 * start * wTmp.rows();
      launchParams.scaleC = alpha;
      launchParams.inputZeropad[0] = prevOutputZeropad[currentStream][0][start];
      launchParams.inputZeropad[1] = prevOutputZeropad[currentStream][1][start];
      launchParams.outputZeropad[0] =
          prevOutputZeropad[currentStream][0][start];
      launchParams.outputZeropad[1] =
          prevOutputZeropad[currentStream][1][start];
      launchParams.inputBufferStride = wTmp.rows();
      launchParams.outputBufferStride = wTmp.rows();
      launchParams.buffer = (void **)&wTmpGPU;
      launchParams.outputBuffer = (void **)&wTmpGPU;

      tempAppCopy->configuration.stream = &pStream[currentStream];
      resFFT = PfSolve_JW::PfSolveAppend(tempAppCopy, -1, &launchParams);
      prevOutputZeropad[currentStream][0][start] =
          launchParams.outputZeropad[0];
      prevOutputZeropad[currentStream][1][start] =
          launchParams.outputZeropad[1];
      start += cols;
    }
  }
}

void IWorlandBackend::scaleSphLaplA(const bool isEven0, const int lshift,
                                    const unsigned int id) const {
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    const MHDFloat alpha = -0.5;
    Matrix &wTmp = this->workTmp(id);
    void *wTmpGPU = this->workTmpGPU(id);

    PfSolve_JW::PfSolveResult resFFT = PfSolve_JW::PFSOLVE_SUCCESS;
    for (auto loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc) + lshift;
      int cols = std::get<2>(loc);
      int rows = std::get<3>(loc);
      PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
      PfSolve_JW::PfSolve_MapKey_block mapKey = {};
      mapKey.size[0] = rows;
      mapKey.size[1] = 2 * cols;
      mapKey.type =
          (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION + RUNTIME_INPUTZEROPAD +
           RUNTIME_OUTPUTZEROPAD + RUNTIME_SCALEC + RUNTIME_INPUTBUFFERSTRIDE +
           RUNTIME_OUTPUTBUFFERSTRIDE) *
              1000 +
          BLOCK_SCALESPHLAPLA + BLOCK_READ_REAL_WRITE_REAL;
      mapKey.lshift = l;
      resFFT =
          PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      if (!tempAppCopy) {
        PfSolve_JW::PfSolveApplication appCopy = {};
        PfSolve_JW::PfSolveConfiguration configurationCopy = {};

        configurationCopy.size[0] = mapKey.size[0];
        configurationCopy.size[1] = mapKey.size[1];
        configurationCopy.block = mapKey.type % 1000;
        configurationCopy.isOutputFormatted = 1;
        configurationCopy.lshift = l;
        configurationCopy.jw_control_bitmask = mapKey.type / 1000;
        configurationCopy.doublePrecision = 1;
        configurationCopy.device = &this->device;
        configurationCopy.num_streams = 1;

        resFFT = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
        resFFT = PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
        resFFT =
            PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      }
      PfSolve_JW::PfSolveLaunchParams launchParams = {};
      launchParams.offsetM = 2 * start * wTmp.rows();
      launchParams.offsetSolution = 2 * start * wTmp.rows();
      launchParams.scaleC = alpha;
      launchParams.inputZeropad[0] = prevOutputZeropad[currentStream][0][start];
      launchParams.inputZeropad[1] = prevOutputZeropad[currentStream][1][start];
      launchParams.outputZeropad[0] =
          prevOutputZeropad[currentStream][0][start];
      launchParams.outputZeropad[1] =
          prevOutputZeropad[currentStream][1][start];
      launchParams.inputBufferStride = wTmp.rows();
      launchParams.outputBufferStride = wTmp.rows();
      launchParams.buffer = (void **)&wTmpGPU;
      launchParams.outputBuffer = (void **)&wTmpGPU;

      tempAppCopy->configuration.stream = &pStream[currentStream];
      resFFT = PfSolve_JW::PfSolveAppend(tempAppCopy, -1, &launchParams);
      prevOutputZeropad[currentStream][0][start] =
          launchParams.outputZeropad[0];
      prevOutputZeropad[currentStream][1][start] =
          launchParams.outputZeropad[1];
      start += cols;
    }
  }
}

void IWorlandBackend::scaleSphLaplB(const bool isEven0, const int lshift,
                                    const unsigned int id) const {
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    const MHDFloat alpha = -0.5;
    Matrix &wTmp = this->workTmp(id);
    void *wTmpGPU = this->workTmpGPU(id);

    PfSolve_JW::PfSolveResult resFFT = PfSolve_JW::PFSOLVE_SUCCESS;
    for (auto loc : *this->pLoc(isEven, id)) {
      int l = std::get<4>(loc) + lshift;
      int cols = std::get<2>(loc);
      int rows = std::get<3>(loc);
      PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
      PfSolve_JW::PfSolve_MapKey_block mapKey = {};
      mapKey.size[0] = rows;
      mapKey.size[1] = 2 * cols;
      mapKey.type =
          (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION + RUNTIME_INPUTZEROPAD +
           RUNTIME_OUTPUTZEROPAD + RUNTIME_SCALEC + RUNTIME_INPUTBUFFERSTRIDE +
           RUNTIME_OUTPUTBUFFERSTRIDE) *
              1000 +
          BLOCK_SCALESPHLAPLB + BLOCK_READ_REAL_WRITE_REAL;
      mapKey.lshift = l;
      resFFT =
          PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      if (!tempAppCopy) {
        PfSolve_JW::PfSolveApplication appCopy = {};
        PfSolve_JW::PfSolveConfiguration configurationCopy = {};

        configurationCopy.size[0] = mapKey.size[0];
        configurationCopy.size[1] = mapKey.size[1];
        configurationCopy.block = mapKey.type % 1000;
        configurationCopy.isOutputFormatted = 1;
        configurationCopy.lshift = l;
        configurationCopy.jw_control_bitmask = mapKey.type / 1000;
        configurationCopy.doublePrecision = 1;
        configurationCopy.device = &this->device;
        configurationCopy.num_streams = 1;

        resFFT = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
        resFFT = PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
        resFFT =
            PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      }
      PfSolve_JW::PfSolveLaunchParams launchParams = {};
      launchParams.offsetM = 2 * start * wTmp.rows();
      launchParams.offsetSolution = 2 * start * wTmp.rows();
      launchParams.scaleC = alpha;
      launchParams.inputZeropad[0] = prevOutputZeropad[currentStream][0][start];
      launchParams.inputZeropad[1] = prevOutputZeropad[currentStream][1][start];
      launchParams.outputZeropad[0] =
          prevOutputZeropad[currentStream][0][start];
      launchParams.outputZeropad[1] =
          prevOutputZeropad[currentStream][1][start];
      launchParams.inputBufferStride = wTmp.rows();
      launchParams.outputBufferStride = wTmp.rows();
      launchParams.buffer = (void **)&wTmpGPU;
      launchParams.outputBuffer = (void **)&wTmpGPU;
      tempAppCopy->configuration.stream = &pStream[currentStream];
      resFFT = PfSolve_JW::PfSolveAppend(tempAppCopy, -1, &launchParams);

      prevOutputZeropad[currentStream][0][start] =
          launchParams.outputZeropad[0];
      prevOutputZeropad[currentStream][1][start] =
          launchParams.outputZeropad[1];
      start += cols;
    }
  }
}

void IWorlandBackend::lshift(const unsigned int id, const int lshift,
                             const bool isEven0) const {
  for (int isEven = 0; isEven < 2; isEven++) {
    for (auto &loc : *this->pLoc(isEven, id)) {
      std::get<4>(loc) += lshift;
    }
  }
}

void IWorlandBackend::nshift(const unsigned int id, const int nshift,
                             const bool isEven0) const {
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &wTmp = this->workTmp(id);
    void *wTmpGPU = this->workTmpGPU(id);
    int s = std::abs(nshift);

    PfSolve_JW::PfSolveResult resFFT = PfSolve_JW::PFSOLVE_SUCCESS;
    for (auto loc : *this->pLoc(isEven, id)) {
      int cols = std::get<2>(loc);
      int rows = std::get<3>(loc) - s;
      PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
      PfSolve_JW::PfSolve_MapKey_block mapKey = {};
      mapKey.size[0] = rows;
      mapKey.size[1] = 2 * cols;
      mapKey.type = (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION +
                     RUNTIME_INPUTZEROPAD + RUNTIME_OUTPUTZEROPAD +
                     RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE) *
                        1000 +
                    BLOCK_READ_REAL_WRITE_REAL;
      resFFT =
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

        resFFT = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
        resFFT = PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
        resFFT =
            PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      }
      PfSolve_JW::PfSolveLaunchParams launchParams = {};
      launchParams.inputZeropad[0] = prevOutputZeropad[currentStream][0][start];
      launchParams.inputZeropad[1] = prevOutputZeropad[currentStream][1][start];

      if (nshift > 0) {
        launchParams.offsetM = 2 * start * wTmp.rows();
        launchParams.offsetSolution = 2 * start * wTmp.rows() + s;
        launchParams.outputZeropad[0] = s;
        launchParams.outputZeropad[1] = s + rows;
      } else {
        launchParams.offsetM = 2 * start * wTmp.rows() + s;
        launchParams.offsetSolution = 2 * start * wTmp.rows();
        launchParams.outputZeropad[0] = 0;
        launchParams.outputZeropad[1] = rows;
      }
      launchParams.inputBufferStride = wTmp.rows();
      launchParams.outputBufferStride = wTmp.rows();
      launchParams.buffer = (void **)&wTmpGPU;
      launchParams.outputBuffer = (void **)&wTmpGPU;

      tempAppCopy->configuration.stream = &pStream[currentStream];
      resFFT = PfSolve_JW::PfSolveAppend(tempAppCopy, -1, &launchParams);
      prevOutputZeropad[currentStream][0][start] =
          launchParams.outputZeropad[0];
      prevOutputZeropad[currentStream][1][start] =
          launchParams.outputZeropad[1];
      start += cols;
    }
  }
}

void IWorlandBackend::copy(const int to, const int from, const int nshift,
                           const bool isEven0) const {
  assert(to != from);
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &wTmp = this->workTmp(from);
    Matrix &extraTmp = this->workTmp(to);

    void *wTmpGPU = this->workTmpGPU(from);
    void *extraTmpGPU = this->workTmpGPU(to);
    int s = std::abs(nshift);
    if (to > 0) {
      *this->pLoc(isEven, to) = *this->pLoc(isEven, from);
    }

    PfSolve_JW::PfSolveResult resFFT = PfSolve_JW::PFSOLVE_SUCCESS;
    for (auto loc : *this->pLoc(isEven, to)) {
      int cols = std::get<2>(loc);
      int rows = std::get<3>(loc) - s;
      PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
      PfSolve_JW::PfSolve_MapKey_block mapKey = {};
      mapKey.size[0] = rows;
      mapKey.size[1] = 2 * cols;
      mapKey.type = (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION +
                     RUNTIME_INPUTZEROPAD + RUNTIME_OUTPUTZEROPAD +
                     RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE) *
                        1000 +
                    BLOCK_READ_REAL_WRITE_REAL;
      resFFT =
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

        resFFT = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
        resFFT = PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
        resFFT =
            PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      }

      PfSolve_JW::PfSolveLaunchParams launchParams = {};
      launchParams.inputZeropad[0] = prevOutputZeropad[currentStream][0][start];
      launchParams.inputZeropad[1] = prevOutputZeropad[currentStream][1][start];
      if (nshift > 0) {
        launchParams.offsetM = 2 * start * wTmp.rows();
        launchParams.offsetSolution = 2 * start * extraTmp.rows() + s;
        launchParams.outputZeropad[0] = s;
        launchParams.outputZeropad[1] = s + rows;
      } else {
        launchParams.offsetM = 2 * start * wTmp.rows() + s;
        launchParams.offsetSolution = 2 * start * extraTmp.rows();
        launchParams.outputZeropad[0] = 0;
        launchParams.outputZeropad[1] = rows;
      }

      launchParams.inputBufferStride = wTmp.rows();
      launchParams.outputBufferStride = extraTmp.rows();
      launchParams.buffer = (void **)&wTmpGPU;
      launchParams.outputBuffer = (void **)&extraTmpGPU;
      tempAppCopy->configuration.stream = &pStream[currentStream];
      resFFT = PfSolve_JW::PfSolveAppend(tempAppCopy, -1, &launchParams);

      prevOutputZeropad[currentStream][0][start] =
          launchParams.outputZeropad[0];
      prevOutputZeropad[currentStream][1][start] =
          launchParams.outputZeropad[1];
      start += cols;
    }
  }
}

void IWorlandBackend::add(const int to, const int from, const int nshift,
                          const bool isEven0) const {
  for (int isEven = 0; isEven < 2; isEven++) {
    int start = 0;
    currentStream = isPhysEven(isEven);
    bufferSolveRes = bufferSolve[isPhysEven(isEven)];
    bufferSolveRes2 = bufferSolve2[isPhysEven(isEven)];
    bufferTemp = bufferTemp0[isPhysEven(isEven)];
    Matrix &wTmp = this->workTmp(to);
    Matrix &extraTmp = this->workTmp(from);

    void *wTmpGPU = this->workTmpGPU(to);
    void *extraTmpGPU = this->workTmpGPU(from);

    int s = std::abs(nshift);

    PfSolve_JW::PfSolveResult resFFT = PfSolve_JW::PFSOLVE_SUCCESS;
    for (auto loc : *this->pLoc(isEven, to)) {
      int cols = std::get<2>(loc);
      int rows = std::get<3>(loc) - s;
      PfSolve_JW::PfSolveApplication *tempAppCopy = 0;
      PfSolve_JW::PfSolve_MapKey_block mapKey = {};
      mapKey.size[0] = rows;
      mapKey.size[1] = 2 * cols;
      mapKey.type = (RUNTIME_OFFSETM + RUNTIME_OFFSETSOLUTION +
                     RUNTIME_INPUTBUFFERSTRIDE + RUNTIME_OUTPUTBUFFERSTRIDE) *
                        1000 +
                    BLOCK_ADD_OUTPUTBUFFER_VALUES + BLOCK_READ_REAL_WRITE_REAL;
      resFFT =
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

        resFFT = PfSolve_JW::initializePfSolve(&appCopy, configurationCopy);
        resFFT = PfSolve_JW::addToLibrary_block(&appLibrary, mapKey, &appCopy);
        resFFT =
            PfSolve_JW::checkLibrary_block(&appLibrary, mapKey, &tempAppCopy);
      }
      PfSolve_JW::PfSolveLaunchParams launchParams = {};

      if (nshift > 0) {
        launchParams.offsetM = 2 * start * extraTmp.rows();
        launchParams.offsetSolution = 2 * start * wTmp.rows() + s;
      } else {
        launchParams.offsetM = 2 * start * extraTmp.rows() + s;
        launchParams.offsetSolution = 2 * start * wTmp.rows();
      }
      launchParams.inputBufferStride = extraTmp.rows();
      launchParams.outputBufferStride = wTmp.rows();
      launchParams.buffer = (void **)&extraTmpGPU;
      launchParams.outputBuffer = (void **)&wTmpGPU;
      tempAppCopy->configuration.stream = &pStream[currentStream];
      resFFT = PfSolve_JW::PfSolveAppend(tempAppCopy, -1, &launchParams);
      start += cols;
    }
  }
}

void applyUpperLower(double *X, int q0, const int r, const Matrix &PS) {
  X[q0] = (PS(0, 0) * X[q0] + PS(1, 0) * X[q0 + 1]);
  q0++;

  for (int i = 1; i < r - 1; ++i, ++q0) {
    X[q0] = (PS(0, i) * X[q0] + PS(1, i) * X[q0 + 1] + PS(2, i) * X[q0 - 1]);
  }

  X[q0] = (X[q0] * PS(0, r - 1) + PS(2, r - 1) * X[q0 - 1]);
}

void applyLowerUpper(double *X, int q0, const int r, const Matrix &PS) {
  X[q0] = (PS(0, r - 1) * X[q0 - 1] + PS(1, r - 1) * X[q0]);
  --q0;

  for (int i = r - 2; i > 0; --i, --q0) {
    X[q0] = (PS(0, i) * X[q0 - 1] + PS(1, i) * X[q0] + PS(2, i) * X[q0 + 1]);
  }

  X[q0] = (PS(1, 0) * X[q0] + PS(2, 0) * X[q0 + 1]);
}

void applyUpperLower(Eigen::Ref<Matrix> out, const int r, const Matrix &PS) {
  out.row(0) = (PS(0, 0) * out.row(0) + PS(1, 0) * out.row(1));

  for (int i = 1; i < r - 1; ++i) {
    out.row(i) = (PS(0, i) * out.row(i) + PS(1, i) * out.row(i + 1) +
                  PS(2, i) * out.row(i - 1));
  }
  out.row(r - 1) =
      (PS(0, r - 1) * out.row(r - 1) + PS(2, r - 1) * out.row(r - 2));
}

void applyLowerUpper(Eigen::Ref<Matrix> out, const int r, const Matrix &PS) {
  out.row(r - 1) =
      (PS(0, r - 1) * out.row(r - 2) + PS(1, r - 1) * out.row(r - 1));

  for (int i = r - 2; i > 0; --i) {
    out.row(i) = (PS(0, i) * out.row(i - 1) + PS(1, i) * out.row(i) +
                  PS(2, i) * out.row(i + 1));
  }

  out.row(0) = (PS(1, 0) * out.row(0) + PS(2, 0) * out.row(1));
}

} // namespace PfSolve
} // namespace Backend
} // namespace Fft
} // namespace Transform
} // namespace QuICC
