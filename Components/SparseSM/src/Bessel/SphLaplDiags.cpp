/**
 * @file SphLaplDiags.cpp
 * @brief Source of the implementation of the full sphere Bessel SphLapl sparse
 * operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/SphLaplDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

SphLaplDiags::SphLaplDiags(const BesselKind type, const int l) : IDiags(type, l)
{}

} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC
