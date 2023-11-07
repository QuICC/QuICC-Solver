/**
 * @file I6DivR1D1R1_I4.hpp
 * @brief Implementation of the Worland based parallel I6  1/R1 D R1 but 0 mode
 * is I4 P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_KOKKOS_I6DIVR1D1R1_I4_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_KOKKOS_I6DIVR1D1R1_I4_HPP

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


template <class Impl> class I6DivR1D1R1_I4;

/**
 * @brief Implementation of the associated Worland based I6 1/R1 D R1 integrator
 * but 0 mode is I4 P integrator.
 */
template <> class I6DivR1D1R1_I4<kokkos_t> : public KokkosIWorlandIntegrator
{
public:
   /**
    * @brief Constructor
    */
   I6DivR1D1R1_I4() = default;

   /**
    * @brief Destructor
    */
   ~I6DivR1D1R1_I4() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_PI6DIVR1D1R1_I4_HPP
