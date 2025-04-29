/**
 * @file Energy.hpp
 * @brief Implementation of the Worland based energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGY_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGY_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/KokkosEnergyReductor.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/Power.hpp"
#include "QuICC/Transform/Poly/Worland/Tags.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

template <class Impl> class Energy;

/**
 * @brief Implementation of the Worland based energy operator
 */
template <>
class Energy<kokkos_t> : public KokkosEnergyReductor<Power<kokkos_t>>
{
public:
   /**
    * @brief Constructor
    */
   Energy();

   /**
    * @brief Destructor
    */
   virtual ~Energy() = default;

protected:
private:
};

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_ENERGY_HPP
