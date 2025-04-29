/**
 * @file IPfSolveBackend.hpp
 * @brief Interface for a generic PFSOLVE based backend
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_IPFSOLVEBACKEND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_IPFSOLVEBACKEND_HPP

// External includes
//
#include <fftw3.h>

// Project includes
//
#include "Fft/Fftw/Library.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace PfSolve {

/**
 * @brief Interface for a generic PFSOLVE backend
 */
class IPfSolveBackend {
public:
  /**
   * @brief Constructor
   */
  IPfSolveBackend();

  /**
   * @brief Destructor
   */
  virtual ~IPfSolveBackend();

  /**
   * @brief Apply FFT
   *
   * deprecated
   */
  virtual void applyFft() const;

  /**
   * @brief Apply FFT
   *
   * Real to real
   */
  virtual void applyFft(Matrix &, const Matrix &) const;

  /**
   * @brief Apply FFT
   *
   * Complex to real
   */
  virtual void applyFft(Matrix &, const MatrixZ &) const;

  /**
   * @brief Apply FFT
   *
   * Real to complex
   */
  virtual void applyFft(MatrixZ &, const Matrix &) const;

  /**
   * @brief Apply FFT
   *
   * Complex to complex
   */
  virtual void applyFft(MatrixZ &, const MatrixZ &) const;

  /**
   * @brief Get the memory requirements
   */
  virtual MHDFloat requiredStorage() const;

protected:
  /**
   * @brief Plan for the transform
   */
  mutable fftw_plan mPlan;

private:
  /**
   * @brief Initialise the PFSOLVE library
   */
  void initLibrary() const;

  /**
   * @brief Cleanup the FFT
   */
  void cleanupFft();
};

} // namespace PfSolve
} // namespace Backend
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_IPFSOLVEBACKEND_HPP
