/**
 * @file DivR1CylLaplh_Zero.hpp
 * @brief Implementation of the Worland based 1/R cylindrical horizontal
 * laplacian projector but 0 mode is zeroed
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_KOKKOS_DIVR1CYLLAPLH_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_KOKKOS_DIVR1CYLLAPLH_ZERO_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/KokkosIWorlandProjector.hpp"
#include "QuICC/Transform/Poly/Worland/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

template <class Impl> class DivR1CylLaplh_Zero;
/**
 * @brief Implementation of the associated Worland based D projector
 */
template <> class DivR1CylLaplh_Zero<kokkos_t> : public KokkosIWorlandProjector
{
public:
   /**
    * @brief Constructor
    */
   DivR1CylLaplh_Zero() = default;

   /**
    * @brief Destructor
    */
   ~DivR1CylLaplh_Zero() = default;

   void applyUnitOperator(const OpMatrixLZ& rOut, const OpMatrixLZ& in,
      const OpVectorI& scan, const int totalOpsCols) const final;

protected:
   /**
    * @brief Make operator
    */
   void makeOperator(Matrix& op, const Internal::Array& igrid,
      const Internal::Array& iweights, const int i) const final;
};

} // namespace Projector
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_PDIVR1CYLLAPLH_ZERO_HPP
