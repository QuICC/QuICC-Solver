/**
 * @file Energy.cpp
 * @brief Source of the implementation of the Worland D R energy operator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/EnergyD1R1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

EnergyD1R1<kokkos_t>::EnergyD1R1() : KokkosEnergyReductor<PowerD1R1>()
{
   this->setProfileTag();
}

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
