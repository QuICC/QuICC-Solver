/**
 * @file D1Diags.cpp
 * @brief Source of the implementation of the full sphere Worland D1Diags sparse
 * operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/SphEnergy/D1Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

namespace SphEnergy {

D1Diags::D1Diags(const Scalar_t alpha, const int l) :
    QuICC::SparseSM::Worland::Stencil::D1Diags(alpha, MHD_MP(0.5), l)
{}

D1Diags::ACoeff_t D1Diags::d_1(const ACoeff_t& n) const
{
   auto l1 = this->l();

   ACoeff_t num = -(2.0 * l1 * n - l1 + 2.0 * n * n - n - 1.0);
   ACoeff_t den = (2.0 * l1 * n + l1 + 2.0 * n * n + 3.0 * n);

   return this->normalizeDiag(n, -1) * (num / den);
}

D1Diags::ACoeff_t D1Diags::d0(const ACoeff_t& n) const
{
   ACoeff_t val = ACoeff_t::Ones(n.size());

   return this->normalizeDiag(n, 0) * val;
}

} // namespace SphEnergy
} // namespace Stencil
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
