/**
 * @file I4DivR1D1R1_I2.hpp
 * @brief Implementation of the Worland based parallel I4  1/R1 D R1 but 0 mode
 * is I2 P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_KOKKOS_I4DIVR1D1R1_I2_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_KOKKOS_I4DIVR1D1R1_I2_HPP

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

template <class Impl> class I4DivR1D1R1_I2;
/**
 * @brief Implementation of the Worland based I4 1/R1 D R1 integrator but 0 mode
 * is I2 P integrator
 */
template <> class I4DivR1D1R1_I2<kokkos_t> : public KokkosIWorlandIntegrator
{
public:
   /**
    * @brief Constructor
    */
   I4DivR1D1R1_I2() = default;

   /**
    * @brief Destructor
    */
   ~I4DivR1D1R1_I2() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_PI4DIVR1D1R1_I2_HPP
