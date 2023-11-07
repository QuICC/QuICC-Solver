/**
 * @file PowerSLaplR2.hpp
 * @brief Implementation of the Worland based spherical Laplacian R^2 power
 * spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_POWERSLAPLR2_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_POWERSLAPLR2_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/KokkosIWorlandPower.hpp"
#include "QuICC/Transform/Poly/Worland/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

template <class Impl> class PowerSLaplR2;
/**
 * @brief Implementation of the Worland based power spectrum operator on  grid
 */
template <> class PowerSLaplR2<kokkos_t> : public KokkosIWorlandPower
{
public:
   /**
    * @brief Constructor
    */
   PowerSLaplR2();

   /**
    * @brief Destructor
    */
   virtual ~PowerSLaplR2() = default;

   virtual void applyUnitOperator(const OpMatrixL& rOut, const OpMatrixLZ& in,
      const OpVectorI& scan, const int totalOpsCols) const override;

protected:
   /**
    * @brief Make operator
    */
   virtual void makeOperator(Matrix& op, Matrix& eop,
      const Internal::Array& igrid, const Internal::Array& iweights,
      const int i) const override;
};

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWERSLAPLR2_HPP
