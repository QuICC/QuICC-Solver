/**
 * @file Energy.hpp
 * @brief Implementation of the Worland based D R energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGYD1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGYD1R1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/KokkosEnergyReductor.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/PowerD1R1.hpp"
#include "QuICC/Transform/Poly/Worland/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

template <class Impl> class EnergyD1R1;

/**
 * @brief Implementation of the Worland based R^2 energy operator
 */
template <>
class EnergyD1R1<kokkos_t> : public KokkosEnergyReductor<PowerD1R1<kokkos_t>>
{
public:
   /**
    * @brief Constructor
    */
   EnergyD1R1();

   /**
    * @brief Destructor
    */
   virtual ~EnergyD1R1() = default;

protected:
private:
};

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGY_HPP
