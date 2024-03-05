/**
 * @file SphLapl2Diags.cpp
 * @brief Source of the implementation of the full sphere Bessel SphLapl2 sparse
 * operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/SphLapl2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

SphLapl2Diags::SphLapl2Diags(const BesselKind type, const int l) :
    IDiags(type, l)
{}

} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC
