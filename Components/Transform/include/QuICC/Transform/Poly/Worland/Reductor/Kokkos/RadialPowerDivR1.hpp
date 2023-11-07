/**
 * @file RadialPowerDivR1.hpp
 * @brief Implementation of the Worland based 1/R power spectrum operator on
 * radial grid
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_RadialPowerDIVR1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_RadialPowerDIVR1_HPP

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

template <class Impl> class RadialPowerDivR1;
/**
 * @brief Implementation of the Worland based power spectrum operator on radial
 * grid
 */
template <> class RadialPowerDivR1<kokkos_t> : public KokkosIWorlandRadialPower
{
public:
   /**
    * @brief Constructor
    */
   RadialPowerDivR1();

   /**
    * @brief Destructor
    */
   virtual ~RadialPowerDivR1() = default;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_PRadialPowerDivR1_HPP
