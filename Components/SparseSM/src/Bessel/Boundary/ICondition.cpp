/**
 * @file Condition.cpp
 * @brief Source of the interface to a generic Bessel boundary condition
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Bessel/Boundary/ICondition.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Boundary {

ICondition::ICondition(const BesselKind type, const int l) : IDiags(type, l) {}

} // namespace Boundary
} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC
