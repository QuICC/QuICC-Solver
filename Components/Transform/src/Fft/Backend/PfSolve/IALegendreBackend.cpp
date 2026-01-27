/**
 * @file IALegendreBackend.cpp
 * @brief Source of the interface for a generic FFTW based ALegendre integrator
 */

// External includes
//
#include <Eigen/src/misc/blas.h>
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Backend/PfSolve/IALegendreBackend.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace PfSolve_parallALT {

const MHDFloat IALegendreBackend::UPPER_BANDED = 424242.424242;

const MHDFloat IALegendreBackend::LOWER_BANDED = -424242.424242;

const MHDFloat IALegendreBackend::LOWER_BANDED_PADDED = -424242.434343;

IALegendreBackend::IALegendreBackend() : mFlipped(false) {
  this->mELoc.push_back(LocationVector());
  this->mOLoc.push_back(LocationVector());
}

IALegendreBackend::~IALegendreBackend() {

}
void IALegendreBackend::transferDataToGPU(void *gpu_buffer, void *cpu_arr,
                                        uint64_t gpuOffset, uint64_t cpuOffset,
                                        uint64_t transferSize, void *stream) {

}

void IALegendreBackend::transferDataFromGPU(void *cpu_arr, void *gpu_buffer,
                                          uint64_t cpuOffset,
                                          uint64_t gpuOffset,
                                          uint64_t transferSize, void *stream) {

}
void IALegendreBackend::initStorage(const int rows, const int cols, const int n,
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

void IALegendreBackend::addStorage(const int inExtras,
                                 const int outExtras) const {
  assert(inExtras >= 0);
  assert(outExtras >= 0);

}

void IALegendreBackend::init(const SetupType &setup, const int lshift,
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

void IALegendreBackend::setZFilter(const std::set<int> &filter) const {
  this->mZFilter = filter;
}

bool IALegendreBackend::isPhysEven(const bool isSpecEven) const {
  return (isSpecEven ^ this->mFlipped);
}

bool IALegendreBackend::inZFilter(const int l) const {
  return (this->mZFilter.count(l) == 1);
}

int IALegendreBackend::reorder_indices_PfSolve(int i, int M_size, int warpSize,
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

void IALegendreBackend::setWSize() const {
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

IALegendreBackend::LocationVector *
IALegendreBackend::pLoc(const bool isEven, const unsigned int id) const {
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

void IALegendreBackend::resetLocations(const bool isEven, const int id) const {
  for (auto &loc : *this->pLoc(isEven, id)) {
    std::get<4>(loc) = std::get<0>(loc);
  }
}

Matrix &IALegendreBackend::workTmp(const unsigned int id) const {
  assert(this->mpWorkTmp);
  assert(id >= 0);
  assert(this->mpWorkTmp->size() > id);

  return this->mpWorkTmp->at(id);
}

void *IALegendreBackend::workTmpGPU(const unsigned int id) const {
  assert(id >= 0);
  assert(id < 2);
  if (id == 0)
    return this->bufferSolveRes;
  else
    return this->bufferSolveRes2;
}

void IALegendreBackend::initJ() const {
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

const Matrix &IALegendreBackend::J(const bool isEven) const {
  if (this->isPhysEven(isEven)) {
    return this->mJEven;
  } else {
    return this->mJOdd;
  }
}

const std::vector<Matrix> &IALegendreBackend::banded(const bool isEven) const {
  if (this->isPhysEven(isEven)) {
    return this->mEBanded;
  } else {
    return this->mOBanded;
  }
}

const std::vector<void *> &IALegendreBackend::bandedGPU(const bool isEven) const {
  if (this->isPhysEven(isEven)) {
    return this->mEBandedGPU;
  } else {
    return this->mOBandedGPU;
  }
}

void IALegendreBackend::makeUpper(SparseMatrix &M, const Matrix &Mp) const {
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

void IALegendreBackend::makeLower(SparseMatrix &M, const Matrix &Mp,
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

void IALegendreBackend::jacobiShiftMatrix(Matrix &rJ, const int size,
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

void IALegendreBackend::buildShiftU(Matrix &U, const int size,
                                  const MHDFloat alpha, const MHDFloat beta,
                                  const MHDFloat norm) const {
  Matrix J;
  this->jacobiShiftMatrix(J, size, beta, alpha);
  this->buildShiftU(U, 0, size, J, norm);
}

void IALegendreBackend::buildShiftU(SparseMatrix &U, const int size,
                                  const MHDFloat alpha, const MHDFloat beta,
                                  const MHDFloat norm) const {
  Matrix Up;
  this->buildShiftU(Up, size, alpha, beta, norm);
  this->makeUpper(U, Up);
}

void IALegendreBackend::buildShiftU(Matrix &U, const int i, const int size,
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

void IALegendreBackend::buildShiftU(SparseMatrix &U, const int i, const int size,
                                  const Matrix &J, const MHDFloat norm) const {
  Matrix Up;
  this->buildShiftU(Up, i, size, J, norm);
  this->makeUpper(U, Up);
}

void IALegendreBackend::buildShiftV(Matrix &V, const int size,
                                  const MHDFloat alpha, const MHDFloat beta,
                                  const MHDFloat norm) const {
  Matrix J;
  this->jacobiShiftMatrix(J, size, alpha, beta);
  this->buildShiftV(V, 0, size, J, norm);
}

void IALegendreBackend::buildShiftV(SparseMatrix &V, const int size,
                                  const MHDFloat alpha, const MHDFloat beta,
                                  const MHDFloat norm) const {
  Matrix Vp;
  this->buildShiftV(Vp, size, alpha, beta, norm);
  this->makeUpper(V, Vp);
}

void IALegendreBackend::buildShiftV(Matrix &V, const int i, const int size,
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

void IALegendreBackend::buildShiftV(SparseMatrix &V, const int i, const int size,
                                  const Matrix &J, const MHDFloat norm) const {
  Matrix Vp;
  this->buildShiftV(Vp, i, size, J, norm);
  this->makeUpper(V, Vp);
}

void IALegendreBackend::buildShiftM(Matrix &M, const int size,
                                  const MHDFloat alpha, const MHDFloat beta,
                                  const MHDFloat norm,
                                  const bool isSquare) const {
  Matrix J;
  this->jacobiShiftMatrix(J, size + 1, alpha, beta - 1);
  this->buildShiftM(M, 0, size, J, norm, isSquare);
}

void IALegendreBackend::buildShiftM(SparseMatrix &M, const int size,
                                  const MHDFloat alpha, const MHDFloat beta,
                                  const MHDFloat norm,
                                  const bool isSquare) const {
  Matrix Mp;
  this->buildShiftM(Mp, size, alpha, beta, norm, isSquare);
  this->makeLower(M, Mp, isSquare);
}

void IALegendreBackend::buildShiftM(Matrix &M, const int i, const int size,
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

void IALegendreBackend::buildShiftM(SparseMatrix &M, const int i, const int size,
                                  const Matrix &J, const MHDFloat norm,
                                  const bool isSquare) const {
  Matrix Mp;
  this->buildShiftM(Mp, i, size, J, norm);
  this->makeLower(M, Mp, isSquare);
}

void IALegendreBackend::buildShiftPair(Matrix &PS, const bool isVMOrder,
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

void IALegendreBackend::buildShiftPair(SparseMatrix &V, SparseMatrix &M,
                                     const int i, const int size,
                                     const Matrix &J, const MHDFloat normV,
                                     const MHDFloat normM,
                                     const bool isSquare) const {
  this->buildShiftV(V, i, size, J, normV);

  this->buildShiftM(M, i, size, J, normM, isSquare);
}

void IALegendreBackend::setPlan(const bool isEven) const {
  if (this->isPhysEven(isEven)) {
    this->mpApp = &this->mEvenApp;
  } else {
    this->mpApp = &this->mOddApp;
  }
}

int IALegendreBackend::blockSize(const bool isEven) const {
  if (this->isPhysEven(isEven)) {
    return this->mEBlockSize;
  } else {
    return this->mOBlockSize;
  }
}

void IALegendreBackend::applyTriProduct(Matrix &out, const int start,
                                      const int cols,
                                      const SparseMatrix &A) const {
  int r = A.rows();
  int c = A.cols();
  out.block(0, start, r, cols) = A * out.block(0, start, c, cols);
}

void IALegendreBackend::applyTriProduct(Matrix &out, int id, const int start,
                                      const int cols, Matrix &A) const {
  }

void IALegendreBackend::applyBandProduct(Matrix &out, const int start,
                                       const int cols,
                                       const SparseMatrix &A) const {
  int r = A.rows();
  int c = A.cols();
  out.block(0, start, r, cols) = A * out.block(0, start, c, cols);
}

void IALegendreBackend::applyBandProduct(Matrix &out, int id, const int start,
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

 }

void IALegendreBackend::applyTriSolve(Matrix &out, const int start,
                                    const int cols,
                                    const SparseMatrix &A) const {
  int r = A.cols();
  Eigen::SparseLU<SparseMatrix> solver;
  solver.compute(A);
  Matrix tmp = out.block(0, start, r, cols);
  Matrix tmp2 = solver.solve(tmp);
  out.block(0, start, r, cols) = tmp2;
}

void IALegendreBackend::applyTriSolve(Matrix &out, int id, const int start,
                                    const int cols, const MHDFloat scale,
                                    Matrix &A) const {
  }

void IALegendreBackend::applyPair(Matrix &out, const int start, const int cols,
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

void IALegendreBackend::applyPair(Matrix &out, const int id, const int start,
                                const int cols, const int rows,
                                const Matrix &PS, const void *PSGPU) const {
 }

void IALegendreBackend::applyPairCombined(Matrix &out, const int id,
                                        const int start, const int cols,
                                        const int rows_start,
                                        const int rows_end, const Matrix &PS,
                                        const void **bandedGPU) const {
  }

void IALegendreBackend::io(MHDFloat *out, const MHDFloat *in) const {
  this->mpOut = out;
  this->mpIn = in;
}

void IALegendreBackend::scaleC(const MHDFloat c, const bool isEven0,
                             const unsigned int id) const {
  Matrix &wTmp = this->workTmp(id);

  }

void IALegendreBackend::scaleALPY(const MHDFloat a, const MHDFloat y,
                                const bool isEven0, const int lshift,
                                const unsigned int id) const {
  }

void IALegendreBackend::scaleD(const bool isEven0, const int lshift,
                             const unsigned int id) const {
  }

void IALegendreBackend::scaleSphLaplA(const bool isEven0, const int lshift,
                                    const unsigned int id) const {
  }

void IALegendreBackend::scaleSphLaplB(const bool isEven0, const int lshift,
                                    const unsigned int id) const {
  
}

void IALegendreBackend::lshift(const unsigned int id, const int lshift,
                             const bool isEven0) const {
  for (int isEven = 0; isEven < 2; isEven++) {
    for (auto &loc : *this->pLoc(isEven, id)) {
      std::get<4>(loc) += lshift;
    }
  }
}

void IALegendreBackend::nshift(const unsigned int id, const int nshift,
                             const bool isEven0) const {
  
}

void IALegendreBackend::copy(const int to, const int from, const int nshift,
                           const bool isEven0) const {
  assert(to != from);
  
}

void IALegendreBackend::add(const int to, const int from, const int nshift,
                          const bool isEven0) const {
  
}


} // namespace PfSolve
} // namespace Backend
} // namespace Fft
} // namespace Transform
} // namespace QuICC
