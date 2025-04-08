/**
 * @file I2D1_I2.hpp
 * @brief Implementation of the Chebyshev based I^2 of D integrator, but 0 mode
 * is I^2 of P integrator, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I2D1_I2_BASE_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I2D1_I2_BASE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Tags.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/ILinearMapIntegrator.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

template <typename Impl> class I2D1_I2;

/**
 * @brief Implementation of the Chebyshev based I^2 of D integrator, but 0 mode
 * is I^2 of P integrator, with linear map y = ax + b
 */
template <>
class I2D1_I2<base_t> : public ILinearMapIntegrator
{
public:
   /**
    * @brief Constructor
    */
   I2D1_I2() = default;

   /**
    * @brief Destructor
    */
   ~I2D1_I2() = default;

protected:
   /**
    * @brief Sparse matrix operator
    */
   mutable SparseMatrix mOp;

   /**
    * @brief Sparse matrix operator for mean
    */
   mutable SparseMatrix mMeanOp;

private:
   /**
    * @brief Initialize solver operators
    */
   void initOperator() const final;

   /**
    * @brief Apply post FFT operator
    *
    * @param rOut Output values
    */
   void applyPostOperator(Matrix& rOut) const final;

   /**
    * @brief Apply pre FFT operator for component wise openerations
    *
    * @param in   Input values
    * @param useReal Real vs Imag flag
    */
   void applyPreOperator(Matrix& tmp, const MatrixZ& in,
      const bool useReal) const final;

   /**
    * @brief Apply post FFT operator for component wise operations
    *
    * @param rOut Output values
    * @param useReal Real vs Imag flag
    */
   void applyPostOperator(MatrixZ& rOut, const Matrix& tmp,
      const bool useReal) const final;
};

} // namespace Integrator
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I2D1_I2_BASE_HPP
