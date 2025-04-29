/**
 * @file CylLaplh.hpp
 * @brief Implementation of the Worland based cylindrical horizontal laplacian
 * projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_KOKKOS_CYLLAPLH_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_KOKKOS_CYLLAPLH_HPP

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

template <class Impl> class CylLaplh;
/**
 * @brief Implementation of the associated Worland based D projector
 */
template <> class CylLaplh<kokkos_t> : public KokkosIWorlandProjector
{
public:
   /**
    * @brief Constructor
    */
   CylLaplh() = default;

   /**
    * @brief Destructor
    */
   ~CylLaplh() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_PCYLLAPLH_HPP
