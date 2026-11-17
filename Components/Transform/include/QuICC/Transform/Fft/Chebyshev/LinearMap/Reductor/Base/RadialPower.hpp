/**
 * @file RadialPower.hpp
 * @brief Implementation of the Chebyshev based radial power reductor, with
 * linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_RADIALPOWER_BASE_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_RADIALPOWER_BASE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/ILinearMapRadialPower.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Tags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

template <typename Impl> class RadialPower;

/**
 * @brief Implementation of the Chebyshev based radial power reductor, with
 * linear map y = ax + b
 */
template <> class RadialPower<base_t> : public ILinearMapRadialPower
{
public:
   /**
    * @brief Constructor
    */
   RadialPower() = default;

   /**
    * @brief Destructor
    */
   ~RadialPower() = default;

protected:
private:
   /**
    * @brief Apply pre FFT operator
    *
    * @param rOut Output values
    * @param in   Input values
    */
   void applyPreOperator(Matrix& tmp, const Eigen::Ref<const Matrix>& in) const final;

   /**
    * @brief Apply post FFT operator
    *
    * @param rOut Output values
    */
   void applyPostOperator(Eigen::Ref<Matrix> rOut, const Matrix& tmp) const final;

   /**
    * @brief Apply pre FFT operator for component wise openerations
    *
    * @param in   Input values
    * @param useReal Real vs Imag flag
    */
   void applyPreOperator(Matrix& tmp, const Eigen::Ref<const MatrixZ>& in,
      const bool useReal) const final;
};

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_RADIALPOWER_BASE_HPP
