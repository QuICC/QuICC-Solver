/**
 * @file Energy.cpp
 * @brief Source of the implementation of the Worland energy operator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/Energy.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

Energy<kokkos_t>::Energy() : KokkosEnergyReductor<Power>()
{
   this->setProfileTag();
}

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
