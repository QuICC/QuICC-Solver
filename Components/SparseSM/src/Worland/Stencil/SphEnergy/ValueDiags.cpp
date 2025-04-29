/**
 * @file ValueDiags.cpp
 * @brief Source of the implementation of the full sphere Worland ValueDiags
 * sparse operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/SphEnergy/ValueDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

namespace SphEnergy {

ValueDiags::ValueDiags(const Scalar_t alpha, const int l) :
    QuICC::SparseSM::Worland::Stencil::ValueDiags(alpha, MHD_MP(0.5), l)
{}

ValueDiags::ACoeff_t ValueDiags::d_1(const ACoeff_t& n) const
{
   ACoeff_t val = -ACoeff_t::Ones(n.size());

   return this->normalizeDiag(n, -1) * val;
}

ValueDiags::ACoeff_t ValueDiags::d0(const ACoeff_t& n) const
{
   ACoeff_t val = ACoeff_t::Ones(n.size());

   return this->normalizeDiag(n, 0) * val;
}

} // namespace SphEnergy
} // namespace Stencil
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
