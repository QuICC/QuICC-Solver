/**
 * @file P.hpp
 * @brief Implementation of the Worland based parallel P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_KOKKOS_P_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_KOKKOS_P_HPP

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

template <class Impl> class P;

/**
 * @brief Implementation of the associated Worland based P integrator
 */
template <> class P<kokkos_t> : public KokkosIWorlandIntegrator
{
public:
   /**
    * @brief Constructor
    */
   P();

   /**
    * @brief Destructor
    */
   ~P() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_KOKKOS_P_HPP
