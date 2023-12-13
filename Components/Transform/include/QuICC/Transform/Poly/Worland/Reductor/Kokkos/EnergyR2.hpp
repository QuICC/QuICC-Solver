/**
 * @file Energy.hpp
 * @brief Implementation of the Worland based R^2 energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGYR2_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGYR2_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/KokkosEnergyReductor.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/PowerR2.hpp"
#include "QuICC/Transform/Poly/Worland/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

template <class Impl> class EnergyR2;

/**
 * @brief Implementation of the Worland based R^2 energy operator
 */
template <>
class EnergyR2<kokkos_t> : public KokkosEnergyReductor<PowerR2<kokkos_t>>
{
public:
   /**
    * @brief Constructor
    */
   EnergyR2();

   /**
    * @brief Destructor
    */
   virtual ~EnergyR2() = default;

protected:
private:
};

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGY_HPP
