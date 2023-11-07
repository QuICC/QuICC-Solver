/**
 * @file CylLaplh_DivR1D1R1.hpp
 * @brief Implementation of the Worland based cylindrical horizontal laplacian
 * projector but 0 mode 1/R D R projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_KOKKOS_CYLLAPLH_DIVR1D1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_KOKKOS_CYLLAPLH_DIVR1D1R1_HPP

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

template <class Impl> class CylLaplh_DivR1D1R1;
/**
 * @brief Implementation of the associated Worland based D projector
 */
template <> class CylLaplh_DivR1D1R1<kokkos_t> : public KokkosIWorlandProjector
{
public:
   /**
    * @brief Constructor
    */
   CylLaplh_DivR1D1R1() = default;

   /**
    * @brief Destructor
    */
   ~CylLaplh_DivR1D1R1() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_PCYLLAPLH_DIVR1D1R1_HPP
