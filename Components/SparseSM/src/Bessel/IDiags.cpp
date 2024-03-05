/**
 * @file IDiags.cpp
 * @brief Source of the interface for Bessel sparse operators
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/IDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

IDiags::IDiags(const BesselKind type, const int l) :
    mType(type), mL(static_cast<Scalar_t>(l))
{}

IDiags::Scalar_t IDiags::l() const
{
   return this->mL;
}

BesselKind IDiags::type() const
{
   return this->mType;
}

} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC
