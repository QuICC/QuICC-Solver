/**
 * @file IdDiags.cpp
 * @brief Source of the implementation of the full sphere Bessel I2 sparse
 * operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/IdDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

IdDiags::IdDiags(const BesselKind type, const int l) : IDiags(type, l) {}

IdDiags::ACoeff_t IdDiags::d0(const ACoeff_t& n) const
{
   ACoeff_t val = ACoeff_t::Ones(n.size());

   return val;
}

} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC
