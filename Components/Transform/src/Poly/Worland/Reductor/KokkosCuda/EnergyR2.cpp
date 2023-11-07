/**
 * @file Energy.cpp
 * @brief Source of the implementation of the Worland energy operator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/EnergyR2.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

EnergyR2<kokkos_t>::EnergyR2() : KokkosEnergyReductor<PowerR2>()
{
   this->setProfileTag();
}

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
