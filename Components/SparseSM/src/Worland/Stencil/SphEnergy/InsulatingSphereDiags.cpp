/**
 * @file InsulatingSphereDiags.cpp
 * @brief Source of the implementation of the full sphere Worland insulating
 * sphereDiags sparse operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/SphEnergy/InsulatingSphereDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

namespace SphEnergy {

InsulatingSphereDiags::InsulatingSphereDiags(const Scalar_t alpha,
   const int l) :
    QuICC::SparseSM::Worland::Stencil::InsulatingSphereDiags(alpha, MHD_MP(0.5),
       l)
{}

InsulatingSphereDiags::ACoeff_t InsulatingSphereDiags::d_1(
   const ACoeff_t& n) const
{
   auto l1 = this->l();

   ACoeff_t num = -n * (2.0 * l1 + 2.0 * n - 1.0);
   ACoeff_t den = (n + 1.0) * (2.0 * l1 + 2.0 * n + 1.0);

   ACoeff_t val = num / den;

   return this->normalizeDiag(n, -1) * val;
}

InsulatingSphereDiags::ACoeff_t InsulatingSphereDiags::d0(
   const ACoeff_t& n) const
{
   ACoeff_t val = ACoeff_t::Ones(n.size());

   return this->normalizeDiag(n, 0) * val;
}

} // namespace SphEnergy
} // namespace Stencil
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
