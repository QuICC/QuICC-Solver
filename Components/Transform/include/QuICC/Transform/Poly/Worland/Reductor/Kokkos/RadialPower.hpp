/**
 * @file RadialPower.hpp
 * @brief Implementation of the Worland based power spectrum operator on radial
 * grid
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_RADIALPOWER_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_RADIALPOWER_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/KokkosIWorlandRadialPower.hpp"
#include "QuICC/Transform/Poly/Worland/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

template <class Impl> class RadialPower;
/**
 * @brief Implementation of the Worland based power spectrum operator on radial
 * grid
 */
template <> class RadialPower<kokkos_t> : public KokkosIWorlandRadialPower
{
public:
   /**
    * @brief Constructor
    */
   RadialPower();

   /**
    * @brief Destructor
    */
   virtual ~RadialPower() = default;

   virtual void applyUnitOperator(const OpMatrixL& rOut, const OpMatrixLZ& in,
      const OpVectorI& scan, const int totalOpsCols) const override;

protected:
   /**
    * @brief Make operator
    */
   virtual void makeOperator(Matrix& op, const Internal::Array& igrid,
      const Internal::Array& iweights, const int i) const override;
};

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_PRADIALPOWER_HPP
