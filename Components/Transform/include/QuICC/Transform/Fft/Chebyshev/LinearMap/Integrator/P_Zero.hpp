/**
 * @file P_Zero.hpp
 * @brief Implementation of the Chebyshev based P integrator, with linear map y
 * = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_P_ZERO_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_P_ZERO_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/IChebyshevIntegrator.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

/**
 * @brief Implementation of the Chebyshev based P integrator, with linear map y
 * = ax + b
 */
class P_Zero : public IChebyshevIntegrator
{
public:
   /**
    * @brief Constructor
    */
   P_Zero() = default;

   /**
    * @brief Destructor
    */
   ~P_Zero() = default;

protected:
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

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_P_ZERO_HPP
