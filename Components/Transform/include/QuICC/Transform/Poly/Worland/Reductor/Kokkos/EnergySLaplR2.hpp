/**
 * @file EnergySLaplR2.hpp
 * @brief Implementation of the Worland based spherical Laplacian R^2 energy
 * operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGYSLAPLR2_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGYSLAPLR2_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/KokkosEnergyReductor.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/PowerSLaplR2.hpp"
#include "QuICC/Transform/Poly/Worland/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

template <class Impl> class EnergySLaplR2;

/**
 * @brief Implementation of the Worland based R^2 energy operator
 */
template <>
class EnergySLaplR2<kokkos_t>
    : public KokkosEnergyReductor<PowerSLaplR2<kokkos_t>>
{
public:
   /**
    * @brief Constructor
    */
   EnergySLaplR2();

   /**
    * @brief Destructor
    */
   virtual ~EnergySLaplR2() = default;

protected:
private:
};

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGYSLAPLR2_HPP
