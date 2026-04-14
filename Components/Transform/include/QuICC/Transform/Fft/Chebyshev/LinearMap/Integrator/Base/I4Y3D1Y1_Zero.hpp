/**
 * @file I4Y3D1Y1_Zero.hpp
 * @brief Implementation of the Chebyshev based I^4 Y^4 of 1/Y D Y integrator,
 * but 0 mode is zeroed, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I4Y3D1Y1_ZERO_BASE_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I4Y3D1Y1_ZERO_BASE_HPP

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

template <typename Impl> class I4Y3D1Y1_Zero;

/**
 * @brief Implementation of the Chebyshev based I^4 Y^4 of 1/Y D Y integrator,
 * but 0 mode is zeroed, with linera map y = ax + b
 */
template <>
class I4Y3D1Y1_Zero<base_t> : public ILinearMapIntegrator
{
public:
   /**
    * @brief Constructor
    */
   I4Y3D1Y1_Zero() = default;

   /**
    * @brief Destructor
    */
   ~I4Y3D1Y1_Zero() = default;

protected:
   /**
    * @brief Sparse matrix operator
    */
   mutable SparseMatrix mOp;

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
   void applyPostOperator(Eigen::Ref<Matrix> rOut) const final;

   /**
    * @brief Apply pre FFT operator for component wise openerations
    *
    * @param in   Input values
    * @param useReal Real vs Imag flag
    */
   void applyPreOperator(Matrix& tmp, const Eigen::Ref<const MatrixZ>& in,
      const bool useReal) const final;

   /**
    * @brief Apply post FFT operator for component wise operations
    *
    * @param rOut Output values
    * @param useReal Real vs Imag flag
    */
   void applyPostOperator(Eigen::Ref<MatrixZ> rOut, const Matrix& tmp,
      const bool useReal) const final;
};

} // namespace Integrator
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I4Y3D1Y1_ZERO_BASE_HPP
