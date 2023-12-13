/**
 * @file I2DivR1_Zero.hpp
 * @brief Implementation of the Worland based I2 1/R1 parallel integrator but 0
 * mode is zeroed
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_KOKKOS_I2DIVR1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_KOKKOS_I2DIVR1_ZERO_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/KokkosIWorlandIntegrator.hpp"
#include "QuICC/Transform/Poly/Worland/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

template <class Impl> class I2DivR1_Zero;
/**
 * @brief Implementation of the associated Worland based I2 1/R1 integrator but
 * 0 mode is zeroed.
 */
template <> class I2DivR1_Zero<kokkos_t> : public KokkosIWorlandIntegrator
{
public:
   /**
    * @brief Constructor
    */
   I2DivR1_Zero();

   /**
    * @brief Destructor
    */
   ~I2DivR1_Zero() = default;

   void applyUnitOperator(const OpMatrixLZ& rOut, const OpMatrixLZL& in,
      const OpVectorI& scan, const int totalOpsCols) const final;

protected:
   /**
    * @brief Make operator
    */
   void makeOperator(Matrix& op, const Internal::Array& igrid,
      const Internal::Array& iweights, const int i) const final;
};

} // namespace Integrator
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_PI2DIVR1_ZERO_HPP
