/**
 * @file EnergyD1.hpp
 * @brief Implementation of the Chebyshev based D^1 energy reductor, with linear
 * map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGYD1_BASE_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGYD1_BASE_HPP

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

template <typename Impl> class EnergyD1;

/**
 * @brief Implementation of the Chebyshev based D^1 energy reductor, with linear
 * map y = ax + b
 */
template <> class EnergyD1<base_t> : public ILinearMapEnergy
{
public:
   /**
    * @brief Constructor
    */
   EnergyD1() = default;

   /**
    * @brief Destructor
    */
   ~EnergyD1() = default;

protected:
   /**
    * @brief Initialise operator
    */
   void initOperator() const final;

   /**
    * @brief Initialize storage
    */
   void initBackend() const final;

private:
   /**
    * @brief Apply pre FFT operator
    *
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

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_REDUCTOR_ENERGYD1_BASE_HPP
