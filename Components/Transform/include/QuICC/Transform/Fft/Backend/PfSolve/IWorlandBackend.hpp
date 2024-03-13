/**
 * @file IWorlandBackend.hpp
 * @brief Interface for a generic Worland PFSOLVE based integrator
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_IWORLANDBACKEND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_IWORLANDBACKEND_HPP

// External includes
//
#include <set>

// Project includes
//
#include "QuICC/Transform/Fft/Backend/PfSolve/IPfSolveBackend.hpp"
#include "QuICC/Transform/Fft/Worland/Setup.hpp"
#include "Types/Typedefs.hpp"
#if (VKFFT_BACKEND == 1)
#include <cuComplex.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <nvrtc.h>
#elif (VKFFT_BACKEND == 2)
#ifndef __HIP_PLATFORM_HCC__
#define __HIP_PLATFORM_HCC__
#endif
#include <hip/hip_complex.h>
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <hip/hiprtc.h>
#endif

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace PfSolve {

namespace VkFFT {
#include "vkFFT.h"
}

namespace PfSolve_JW {
#include "pfSolve_AppCollectionManager.hpp"
}
/**
 * @brief Interface for a generic Worland PFSOLVE based integrator
 */
class IWorlandBackend : public IPfSolveBackend {
public:
  /// Typedef for the configuration class
  typedef Worland::Setup SetupType;

  /// Typedef for the configuration class as a shared pointer
  typedef Worland::SharedSetup SharedSetupType;

  /**
   * @brief Typedef for block location
   * <0>: l
   * <1>: start
   * <2>: multiplicity
   * <3>: Accurate resolution
   * <4>: current l
   */
  typedef std::tuple<int, int, int, int, int> LocationType;

  //// Typedef for vector of block locations
  typedef std::vector<LocationType> LocationVector;

  /**
   * @brief Constructor
   */
  IWorlandBackend();

  /**
   * @brief Destructor
   */
  virtual ~IWorlandBackend();

  /**
   * @brief Set zero filter
   */
  void setZFilter(const std::set<int> &filter) const;

  /**
   * @brief Wrapper for Cpu -> Gpu transfers
   * @param gpu_buffer pointer to the output GPU buffer 
   * @param cpu_arr pointer to the input CPU number array's first element
   * @param gpuOffset GPU buffer offset in bytes
   * @param cpuOffset CPU buffer offset in bytes
   * @param transferSize number of bytes transferred
   * @param stream pointer to stream object to perform operation in
   */
  static void transferDataToGPU(void *gpu_buffer, void *cpu_arr,
                                uint64_t gpuOffset, uint64_t cpuOffset,
                                uint64_t transferSize, void *stream);

  /**
   * @brief Wrapper for Gpu -> Cpu transfers
   * @param cpu_arr pointer to the output CPU number array's first element
   * @param gpu_buffer pointer to the input GPU buffer 
   * @param cpuOffset CPU buffer offset in bytes
   * @param gpuOffset GPU buffer offset in bytes
   * @param transferSize number of bytes transferred
   * @param stream pointer to stream object to perform operation in
   */
  static void transferDataFromGPU(void *cpu_arr, void *gpu_buffer,
                                  uint64_t cpuOffset, uint64_t gpuOffset,
                                  uint64_t transferSize, void *stream);

  /**
   * @brief Initialise the PFSOLVE transforms
   */
  virtual void init(const SetupType &setup, const int lshift, const int extraN,
                    const bool lshiftOnlyParity = false,
                    const bool alwaysZeroNegative = false) const;

  /**
   * @brief Initialise the PFSOLVE transforms
   */
  void addStorage(const int inExtras, const int outExtras) const;

  /**
   * @brief Set input and output data pointers for FFT (R2R)
   */
  virtual void io(MHDFloat *out, const MHDFloat *in) const;

  /**
   * @brief Scaling by constant
   */
  void scaleC(const MHDFloat c, const bool isEven0,
              const unsigned int id = 0) const;

  /**
   * @brief Scale coefficients by a*l + y
   */
  void scaleALPY(const MHDFloat a, const MHDFloat y, const bool isEven0,
                 const int lshift = 0, const unsigned int id = 0) const;

  /**
   * @brief Scale coefficients by 2.0(l+n+1)*sqrt((n+1)/(l+n+1))
   */
  void scaleD(const bool isEven0, const int lshift = 0,
              const unsigned int id = 0) const;

  /**
   * @brief Scale coefficients by ??
   */
  void scaleSphLaplA(const bool isEven0, const int lshift = 0,
                     const unsigned int id = 0) const;

  /**
   * @brief Scale coefficients by ??
   */
  void scaleSphLaplB(const bool isEven0, const int lshift = 0,
                     const unsigned int id = 0) const;

  /**
   * @brief Shift l in temporary storage
   */
  void lshift(const unsigned int id, const int lshift,
              const bool isEven0) const;

  /**
   * @brief Shift n expansion in temporary storage
   */
  void nshift(const unsigned int id, const int nshift,
              const bool isEven0) const;

  /**
   * @brief Copy out temporary storage into new extra temporary
   */
  void copy(const int to, const int from, const int nshift,
            const bool isEven0) const;

  /**
   * @brief Add extra temporary to out temporary storage
   */
  void add(const int to, const int from, const int nshift,
           const bool isEven0) const;

protected:
  /**
   * @brief Change internal spectral resolution
   */
  virtual void setWSize() const;

  /**
   * @brief Is physical representation even?
   */
  bool isPhysEven(const bool isSpecEven) const;

  /**
   * @brief Is filtered out by zero filter
   */
  bool inZFilter(const int l) const;

  /**
   * @brief Matrix indices reorder for PfSolve, returns the new index
   * @param i input index
   * @param M_size target matrix size to align the reordering to
   * @param warpSize GPU warp size 
   * @param used_registers how many values of the system each thread will have (num_threads * used_registers >= M_size) 
   */
  int reorder_indices_PfSolve(int i, int M_size, int warpSize,
                              int used_registers) const;
  /**
   * @brief Get blocksize
   */
  int blockSize(const bool isEven) const;

  /**
   * @brief Accurate expansion for given l
   */
  virtual int lSize(const int l) const = 0;

  /**
   * @brief Accurate expansion for given index i and l
   */
  virtual int lSize(const int i, const int l) const = 0;

  /**
   * @brief Get pointer to vector of locations
   */
  LocationVector *pLoc(const bool isEven, const unsigned int id = 0) const;

  /**
   * @brief Reset locations to initial values
   */
  void resetLocations(const bool isEven, const int id = 0) const;

  /**
   * @brief Initialize temporary work storage
   */
  void initStorage(const int rows, const int cols, const int n,
                   std::vector<Matrix> &s) const;

  /**
   * @brief Set even or odd plan
   */
  void setPlan(const bool isEven) const;

  /**
   * @brief Initialize J matrices
   */
  void initJ() const;

  /**
   * @brief Get J matrix based on parity flag
   */
  const Matrix &J(const bool isEven) const;

  /**
   * @brief Get banded matrices based on parity flag
   */
  const std::vector<Matrix> &banded(const bool isEven) const;

  /**
   * @brief Get banded GPU matrices based on parity flag
   */
  const std::vector<void *> &bandedGPU(const bool isEven) const;

  /**
   * @brief Get work temporary
   */
  Matrix &workTmp(const unsigned int id = 0) const;
  
  /**
   * @brief Get work temporary for GPU
   */
  void *workTmpGPU(const unsigned int id = 0) const;

  /**
   * @brief Convert banded storage upper matrix to sparse matrix
   */
  void makeUpper(SparseMatrix &M, const Matrix &Mp) const;

  /**
   * @brief Convert banded storage lower matrix to sparse matrix
   */
  void makeLower(SparseMatrix &M, const Matrix &Mp,
                 const bool isSquare = true) const;

  /**
   * @brief Compute shift matrix for Jacobi polynomials: P^{a,b-1} -> P^{a,b}
   */
  void jacobiShiftMatrix(Matrix &rJ, const int maxN, const MHDFloat alpha,
                         const MHDFloat beta) const;

  /**
   * @brief Build U(upper) for shifting jacobi alpha by +1
   */
  void buildShiftU(Matrix &U, const int maxN, const MHDFloat alpha,
                   const MHDFloat beta, const MHDFloat norm = 1.0) const;

  /**
   * @brief Build U(upper) for shifting jacobi alpha by +1
   */
  void buildShiftU(SparseMatrix &U, const int maxN, const MHDFloat alpha,
                   const MHDFloat beta, const MHDFloat norm = 1.0) const;

  /**
   * @brief Build U(upper) for shifting jacobi alpha by +1
   */
  void buildShiftU(Matrix &U, const int i, const int maxN, const Matrix &J,
                   const MHDFloat norm = 1.0) const;

  /**
   * @brief Build U(upper) for shifting jacobi alpha by +1
   */
  void buildShiftU(SparseMatrix &U, const int i, const int maxN,
                   const Matrix &J, const MHDFloat norm = 1.0) const;

  /**
   * @brief Build V(upper) for shifting jacobi beta by +1
   */
  void buildShiftV(Matrix &V, const int maxN, const MHDFloat alpha,
                   const MHDFloat beta, const MHDFloat norm = 1.0) const;

  /**
   * @brief Build V(upper) for shifting jacobi beta by +1
   */
  void buildShiftV(SparseMatrix &V, const int maxN, const MHDFloat alpha,
                   const MHDFloat beta, const MHDFloat norm = 1.0) const;

  /**
   * @brief Build V(upper) for shifting jacobi beta by +1
   */
  void buildShiftV(Matrix &V, const int i, const int maxN, const Matrix &J,
                   const MHDFloat norm = 1.0) const;

  /**
   * @brief Build V(upper) for shifting jacobi beta by +1
   */
  void buildShiftV(SparseMatrix &V, const int i, const int maxN,
                   const Matrix &J, const MHDFloat norm = 1.0) const;

  /**
   * @brief Build M(lower) for shifting jacobi beta by +1 and multiply by (1+x)
   */
  void buildShiftM(Matrix &M, const int maxN, const MHDFloat alpha,
                   const MHDFloat beta, const MHDFloat norm = 1.0,
                   const bool isSquare = true) const;

  /**
   * @brief Build M(lower) for shifting jacobi beta by +1 and multiply by (1+x)
   */
  void buildShiftM(SparseMatrix &M, const int maxN, const MHDFloat alpha,
                   const MHDFloat beta, const MHDFloat norm = 1.0,
                   const bool isSquare = true) const;

  /**
   * @brief Build M(lower) for shifting jacobi beta by +1 and multiply by (1+x)
   */
  void buildShiftM(Matrix &M, const int i, const int maxN, const Matrix &J,
                   const MHDFloat norm = 1.0, const bool isSquare = true) const;

  /**
   * @brief Build M(lower) for shifting jacobi beta by +1 and multiply by (1+x)
   */
  void buildShiftM(SparseMatrix &M, const int i, const int maxN,
                   const Matrix &J, const MHDFloat norm = 1.0,
                   const bool isSquare = true) const;

  /**
   * @brief Build V(upper),M(lower) pair for shifting jacobi beta by +1
   */
  void buildShiftPair(Matrix &PS, const bool isVMOrder, const int i,
                      const int maxN, const Matrix &J,
                      const MHDFloat normV = 1.0, const MHDFloat normM = 1.0,
                      const bool isSquare = true) const;

  /**
   * @brief Build V(upper),M(lower) pair for shifting jacobi beta by +1
   */
  void buildShiftPair(SparseMatrix &V, SparseMatrix &M, const int i,
                      const int maxN, const Matrix &J,
                      const MHDFloat normV = 1.0, const MHDFloat normM = 1.0,
                      const bool isSquare = true) const;

  /**
   * @brief Apply triangular banded product
   */
  void applyTriProduct(Matrix &out, int id, const int start, const int cols,
                       Matrix &A) const;

  /**
   * @brief Apply triangular sparse product
   */
  void applyTriProduct(Matrix &out, const int start, const int cols,
                       const SparseMatrix &A) const;

  /**
   * @brief Apply triangular banded product
   */
  void applyBandProduct(Matrix &out, int id, const int start, const int cols,
                        Matrix &A) const;

  /**
   * @brief Apply triangular sparse product
   */
  void applyBandProduct(Matrix &out, const int start, const int cols,
                        const SparseMatrix &A) const;

  /**
   * @brief Apply triangular banded linear solve
   */
  void applyTriSolve(Matrix &out, int id, const int start, const int cols,
                     const MHDFloat scale, Matrix &A) const;

  /**
   * @brief Apply triangular sparse linear solve
   */
  void applyTriSolve(Matrix &out, const int start, const int cols,
                     const SparseMatrix &A) const;

  /**
   * @brief Apply pair of banded product and banded linear GPU solve 
   * @param PSGPU pointer to the connection matrix buffer on GPU
   */
  void applyPair(Matrix &out, int id, const int start, const int cols,
                 const int rows, const Matrix &PS, const void *PSGPU) const;

  /**
   * @brief Apply pair of banded product and banded linear GPU solve - multiple connection matrices
   * @param bandedGPU pointer to the first element of the array of connection matrices buffers on GPU
   */
  void applyPairCombined(Matrix &out, int id, const int start, const int cols,
                         const int rows_in, const int rows_out,
                         const Matrix &PS, const void **bandedGPU) const;
  /**
   * @brief Apply pair of sparse product and sparse linear solve
   */
  void applyPair(Matrix &out, const int start, const int cols, const int rows,
                 const SparseMatrix &P, const SparseMatrix &S) const;

  /**
   * @brief Flag for identifying upper triangular matrix
   */
  static const MHDFloat UPPER_BANDED;

  /**
   * @brief Flag for identifying lower triangular matrix
   */
  static const MHDFloat LOWER_BANDED;

  /**
   * @brief Flag for identifying lower rectangular triangular matrix
   */
  static const MHDFloat LOWER_BANDED_PADDED;

  /**
   * @brief Flip parity with respect to l
   */
  mutable bool mFlipped;

  /**
   * @brief VkFFT app pointer that is currently in use
   */
  mutable VkFFT::VkFFTApplication *mpApp;

  /**
   * @brief VkFFT DCT app for even polynomials
   */
  mutable VkFFT::VkFFTConfiguration mEvenConfiguration = {};
  mutable VkFFT::VkFFTApplication mEvenApp = {};
  mutable uint64_t mEvenSize;

  /**
   * @brief VkFFT DCT app for odd polynomials
   */
  mutable VkFFT::VkFFTConfiguration mOddConfiguration = {};
  mutable VkFFT::VkFFTApplication mOddApp = {};
  mutable uint64_t mOddSize;

  /**
   * @brief CUDA stream that id is currently being used for dispatch of kernels
   */
  mutable int currentStream = 0;

  /**
   * @brief Library that collects and manages all kernels related to PfSolve
   */
  mutable PfSolve_JW::PfSolve_AppLibrary appLibrary = {};

  /**
   * @brief Map that keeps track of zeropadding u
   */
  mutable std::map<uint64_t, uint64_t> prevOutputZeropad[2][2];
#if (VKFFT_BACKEND == 1)
  /**
   * @brief CUDA device object
   */
  mutable CUdevice device;

  /**
   * @brief CUDA memory for matrices
   */
  mutable std::vector<void *> bufferSolveOdd;
  mutable std::vector<void *> bufferSolveEven;

  /**
   * @brief CUDA pointer to the memory for the system solved currently
   */
  mutable void *bufferSolveRes = 0;
  mutable void *bufferSolveRes2 = 0;
  mutable void *bufferTemp = 0;

  /**
   * @brief CUDA separate memory for even/odd systems - matches max number of
   * used streams (currently 4, but 2 are used)
   */
  mutable void *bufferSolve[4] = {};
  mutable void *bufferSolve2[4] = {};
  mutable void *bufferTemp0[4] = {};

  /**
   * @brief CUDA streams
   */
  mutable cudaStream_t pStream[4] = {};
#elif (VKFFT_BACKEND == 2)
  mutable hipDevice_t device;

  mutable std::vector<void *> bufferSolveOdd;
  mutable std::vector<void *> bufferSolveEven;
  mutable void *bufferSolveRes = 0;
  mutable void *bufferSolveRes2 = 0;
  mutable void *bufferTemp = 0;

  mutable void *bufferSolve[4] = {};
  mutable void *bufferSolve2[4] = {};
  mutable void *bufferTemp[4] = {};
  mutable hipStream_t pStream[4] = {};
#endif
  /**
   * @brief Spec size
   */
  mutable int mSpecSize;

  /**
   * @brief For some operators, the truncation needs to be extended by a few
   * extra modes i in order to produce mSpecSize accurate modes at the end. For
   * example, the use of the second order quas-inverse I2 requires 3 additional
   * modes due to the 3 super diagonals.
   */
  mutable int mWExtra;

  /**
   * @brief Biggest Worland expansion size (ie. l = 0)
   */
  mutable int mWSize;

  /**
   * @brief Even harmonic degree block size
   */
  mutable int mEBlockSize;

  /**
   * @brief Odd harmonic degree block size
   */
  mutable int mOBlockSize;

  /**
   * @brief Filter to zero modes
   */
  mutable std::set<int> mZFilter;

  /**
   * @brief Location of even harmonic degrees
   */
  mutable std::vector<LocationVector> mELoc;

  /**
   * @brief Location of odd harmonic degrees
   */
  mutable std::vector<LocationVector> mOLoc;

  /**
   * @brief Location of zero-ed harmonic degrees
   */
  mutable LocationVector mZLoc;

  /**
   * @brief Jacobi shift matrix even modes
   */
  mutable Matrix mJEven;

  /**
   * @brief Jacobi shift matrix odd modes
   */
  mutable Matrix mJOdd;

  /**
   * @brief Odd l Banded matrices
   */
  mutable std::vector<Matrix> mOBanded;

  /**
   * @brief Odd l Banded matrices GPU
   */
  mutable std::vector<void *> mOBandedGPU;

  /**
   * @brief Even l Banded matrices
   */
  mutable std::vector<Matrix> mEBanded;

  /**
   * @brief Even l Banded matrices GPU
   */
  mutable std::vector<void *> mEBandedGPU;

  /**
   * @brief Pointer to temporary data
   */
  mutable std::vector<Matrix> *mpWorkTmp;

  /**
   * @brief Temporary data
   */
  mutable std::vector<Matrix> mInTmp;

  /**
   * @brief Temporary data
   */
  mutable std::vector<Matrix> mOutTmp;

  /**
   * @brief Input data pointer
   */
  mutable const MHDFloat *mpIn;

  /**
   * @brief Out data pointer
   */
  mutable MHDFloat *mpOut;

private:
};

/**
 * @brief Apply pair of banded upper product and banded lower linear solve
 */
void applyUpperLower(double *X, const int q0, const int r, const Matrix &PS);

/**
 * @brief Apply pair of banded lower product and banded upper linear solve
 */
void applyLowerUpper(double *X, const int q0, const int r, const Matrix &PS);

/**
 * @brief Apply pair of banded upper product and banded lower linear solve
 */
void applyUpperLower(Eigen::Ref<Matrix> out, const int r, const Matrix &PS);

/**
 * @brief Apply pair of banded lower product and banded upper linear solve
 */
void applyLowerUpper(Eigen::Ref<Matrix> out, const int r, const Matrix &PS);

} // namespace PfSolve
} // namespace Backend
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_IWORLANDBACKEND_HPP
