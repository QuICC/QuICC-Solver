/**
 * @file ALegendreProjector.hpp
 * @brief Interface for a generic ALegendre PFSOLVE based projector
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_ALEGENDREPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_ALEGENDREPROJECTOR_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Backend/PfSolve/IALegendreBackend.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace PfSolve_parallALT {

/**
 * @brief Interface for a generic ALegendre PFSOLVE based projector
 */
class ALegendreProjector : public IALegendreBackend {
public:
  /**
   * @brief Constructor
   */
  ALegendreProjector();

  /**
   * @brief Destructor
   */
  virtual ~ALegendreProjector();

  /**
   * @brief Initialise the PFSOLVE transforms
   */
  virtual void init(const SetupType &setup, const int lshift, const int extraN,
                    const bool lshiftOnlyParity = false,
                    const bool alwaysZeroNegative = false) const override;

  /**
   * @brief Set input and output
   */
  void io(const bool isEven) const;

  /**
   * @brief Set input
   */
  void input(const MatrixZ &in) const;

  /**
   * @brief Set input
   */
  void input(const Matrix &in, const bool isEven,
             const bool needPadding = false) const;

  /**
   * @brief Set input
   */
  void input(const MatrixZ &in, const bool isEven, const bool useReal,
             const bool needPadding = false) const;

  /**
   * @brief Set input and output to internal temporary storage
   */
  void io() const;
  using IALegendreBackend::io;

  /**
   * @brief Set output
   */
  void output(MatrixZ &rOut) const;

  /**
   * @brief Set output
   */
  void output(Matrix &rOut, const bool isEven) const;

  /**
   * @brief Set output
   */
  void output(MatrixZ &rOut, const bool isEven, const bool useReal) const;

  /**
   * @brief Apply FFT
   */
  virtual void applyFft() const override;

  /**
   * @brief Convert ALegendre expansion to  FFT coefficients
   */
  void backwardALegendre(const bool isEven0, const unsigned int id = 0) const;

  /**
   * @brief Lower alpha by 1
   */
  void lowerAlpha(const MHDFloat alpha, const bool isEven0,
                  const unsigned int id = 0,
                  const MHDFloat norm = 1.0 / std::sqrt(2.0)) const;

  /**
   * @brief Lower beta by 1
   */
  void lowerBeta(const MHDFloat alpha, const bool isEven0,
                 const unsigned int id = 0,
                 const MHDFloat norm = 1.0 / std::sqrt(2.0)) const;

  /**
   * @brief Lower beta by 1 and divide by r^2
   */
  void lowerR2Beta(const MHDFloat alpha, const bool isEven0,
                   const unsigned int id = 0,
                   const MHDFloat norm = 1.0 / std::sqrt(2.0)) const;

protected:
  /**
   * @brief Accurate expansion for given l
   */
  virtual int lSize(const int l) const override;

  /**
   * @brief Accurate expansion for given index i and l
   */
  virtual int lSize(const int i, const int l) const override;

private:
  /**
   * @brief Initialize banded matrices
   */
  void initBanded(const bool isEven) const;

  /**
   * @brief Free banded matrices
   */
  void freeBanded(const bool isEven) const;

  /**
   * @brief Apply padding
   */
  void applyPadding(Matrix &rData, const int extraRows = 0) const;

  /**
   * @brief Apply partial backward ALegendre transform
   * @param bandedGPU vector of GPU connection matrices buffers
   */
  void partialBackwardALegendre(const int l, const int id, const int i0,
                              const int start, const int cols,
                              const std::vector<Matrix> &banded,
                              const std::vector<void *> &bandedGPU,
                              Matrix &in) const;

  /**
   * @brief Padding size
   */
  mutable int mPadSize;
};

} // namespace PfSolve
} // namespace Backend
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_WORLANDPROJECTOR_HPP
