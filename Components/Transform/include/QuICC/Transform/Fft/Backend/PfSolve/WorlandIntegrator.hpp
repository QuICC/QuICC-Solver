/**
 * @file WorlandIntegrator.hpp
 * @brief Interface for a generic Worland PFSOLVE based integrator
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_WORLANDINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_WORLANDINTEGRATOR_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Backend/PfSolve/IWorlandBackend.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace PfSolve {

/**
 * @brief Interface for a generic Worland PFSOLVE based integrator
 */
class WorlandIntegrator : public IWorlandBackend {
public:
  /**
   * @brief Constructor
   */
  WorlandIntegrator();

  /**
   * @brief Destructor
   */
  virtual ~WorlandIntegrator();

  /**
   * @brief Initialise the PFSOLVE transforms
   */
  virtual void init(const SetupType &setup, const int lshift, const int extraN,
                    const bool lshiftOnlyParity = false,
                    const bool alwaysZeroNegative = false) const override;

  using IWorlandBackend::io;

  /**
   * @brief Set input and output data pointers for FFT (R2R)
   */
  virtual void io(const bool isEven) const;

  /**
   * @brief Set input and output data pointers for FFT (R2R)
   */
  virtual void input(const MatrixZ &in) const;

  /**
   * @brief Set input and output data pointers for FFT (R2R)
   */
  virtual void input(const Matrix &in, const bool isEven) const;

  /**
   * @brief Set input and output data pointers for FFT (R2R)
   */
  virtual void input(const MatrixZ &out, const bool isEven,
                     const bool useReal) const;

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
   * @brief Convert FFT coefficients to Worland expansion
   */
  void forwardWorland(const bool isEven0, const int id = 0) const;

  /**
   * @brief Lower alpha by 1
   */
  void lowerAlpha(const MHDFloat alpha, const bool isEven0, const int id = 0,
                  const MHDFloat norm = 1.0 / std::sqrt(2.0)) const;

  /**
   * @brief Lower alpha by 1
   */
  void raiseAlpha(const MHDFloat alpha, const bool isEven0, const int id = 0,
                  const MHDFloat norm = 1.0 / std::sqrt(2.0)) const;

  /**
   * @brief Lower beta by 1
   */
  void lowerBeta(const MHDFloat alpha, const bool isEven0, const int id = 0,
                 const MHDFloat norm = 1.0 / std::sqrt(2.0)) const;

  /**
   * @brief Raise beta by 1
   */
  void raiseBeta(const MHDFloat alpha, const bool isEven0, const int id = 0,
                 const MHDFloat norm = 1.0 / std::sqrt(2.0)) const;

  /**
   * @brief Lower beta by 1 and divide by r^2
   */
  void lowerR2Beta(const MHDFloat alpha, const bool isEven0, const int id = 0,
                   const MHDFloat norm = 1.0 / std::sqrt(2.0)) const;

  /**
   * @brief Raise beta by 1 and multiply by r^2
   */
  void raiseR2Beta(const MHDFloat alpha, const bool isEven0, const int id = 0,
                   const MHDFloat norm = 1.0 / std::sqrt(2.0),
                   const bool scaleL0 = false) const;

  /**
   * @brief Apply sparse integration operator I2
   */
  void applyI2(const bool isEven0, const int id = 0) const;

  /**
   * @brief Apply sparse integration operator I4
   */
  void applyI4(const bool isEven0, const int id = 0) const;

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
   * @brief Apply partial backward Worland transform
   * @param bandedGPU vector of GPU connection matrices buffers
   */
  void partialForwardWorland(const int l, const int id, const int i0,
                             const int start, const int cols,
                             const std::vector<Matrix> &banded,
                             const std::vector<void *> &bandedGPU,
                             Matrix &out) const;

  /**
   * @brief Bwd size
   */
  mutable int mBwdSize;

  /**
   * @brief FFT scaling factor
   */
  mutable MHDFloat mFftScaling;
};

} // namespace PfSolve
} // namespace Backend
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_WORLANDINTEGRATOR_HPP
