/**
 * @file D1CylLaplh_D1DivR1D1R1.hpp
 * @brief Implementation of the Worland based D of cylindrical horizontal
 * laplacian projector but 0 mode is D 1/R D R
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_KOKKOS_D1CYLLAPLH_D1DIVR1D1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_KOKKOS_D1CYLLAPLH_D1DIVR1D1R1_HPP


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

template <class Impl> class D1CylLaplh_D1DivR1D1R1;

/**
 * @brief Implementation of the associated Worland based D projector
 */
template <>
class D1CylLaplh_D1DivR1D1R1<kokkos_t> : public KokkosIWorlandProjector
{
public:
   /**
    * @brief Constructor
    */
   D1CylLaplh_D1DivR1D1R1() = default;

   /**
    * @brief Destructor
    */
   ~D1CylLaplh_D1DivR1D1R1() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_PD1CYLLAPLH_D1DIVR1D1R1_HPP
