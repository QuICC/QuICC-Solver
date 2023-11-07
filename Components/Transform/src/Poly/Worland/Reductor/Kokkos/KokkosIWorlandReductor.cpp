/**
 * @file IWorlandReductor.cpp
 * @brief Source of the interface to a Worland based reduction operator (e.g.
 * energy)
 */

// System includes
//

// Project includes
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/KokkosIWorlandReductor.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

KokkosIWorlandReductor::KokkosIWorlandReductor() : IWorlandReductor()
{
   this->mProfileTag += "-Reductor";
}

KokkosIWorlandReductor::~KokkosIWorlandReductor() {}

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
