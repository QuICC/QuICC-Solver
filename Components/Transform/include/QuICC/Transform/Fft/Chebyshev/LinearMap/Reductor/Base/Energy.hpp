/**
 * @file Energy.hpp
 * @brief Implementation of the Chebyshev based energy reductor, with linear map
 * y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGY_BASE_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGY_BASE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/ILinearMapEnergy.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Tags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

template <typename Impl> class Energy;

/**
 * @brief Implementation of the Chebyshev based energy reductor, with linear map
 * y = ax + b
 */
template <> class Energy<base_t> : public ILinearMapEnergy
{
public:
   /**
    * @brief Constructor
    */
   Energy() = default;

   /**
    * @brief Destructor
    */
   ~Energy() = default;

protected:
private:
   /**
    * @brief Apply pre FFT operator
    *
    * @param rOut Output values
    * @param in   Input values
    */
   void applyPreOperator(Matrix& tmp, const Matrix& in) const final;

   /**
    * @brief Apply post FFT operator
    *
    * @param rOut Output values
    */
   void applyPostOperator(Matrix& rOut, const Matrix& tmp) const final;

   /**
    * @brief Apply pre FFT operator for component wise openerations
    *
    * @param in   Input values
    * @param useReal Real vs Imag flag
    */
   void applyPreOperator(Matrix& tmp, const MatrixZ& in,
      const bool useReal) const final;
};

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGY_BASE_HPP
