/**
 * @file D1_P.hpp
 * @brief Implementation of the associated Worland based D parallel projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_KOKKOS_D1_P_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_KOKKOS_D1_P_HPP

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

template <class Impl> class D1_P;
/**
 * @brief Implementation of the associated Worland based D projector
 */
template <> class D1_P<kokkos_t> : public KokkosIWorlandProjector
{
public:
   /**
    * @brief Constructor
    */
   D1_P() = default;

   /**
    * @brief Destructor
    */
   ~D1_P() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_PD1_P_HPP
