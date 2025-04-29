/**
 * @file ValueD2Diags.cpp
 * @brief Source of the implementation of the full sphere Worland ValueD2Diags
 * sparse operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/SphEnergy/ValueD2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

namespace SphEnergy {

ValueD2Diags::ValueD2Diags(const Scalar_t alpha, const int l) :
    QuICC::SparseSM::Worland::Stencil::ValueD2Diags(alpha, MHD_MP(0.5), l)
{}

ValueD2Diags::ACoeff_t ValueD2Diags::d_2(const ACoeff_t& n) const
{
   auto l1 = this->l();

   ACoeff_t num = (2.0 * l1 + 4.0 * n - 3.0) *
                  (2.0 * l1 * n - 2.0 * l1 + 2.0 * n * n - 3.0 * n - 1.0);
   ACoeff_t den =
      ((2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 * n + 2.0 * n * n + n - 2.0));

   return this->normalizeDiag(n, -2) * (num / den);
}

ValueD2Diags::ACoeff_t ValueD2Diags::d_1(const ACoeff_t& n) const
{
   auto l1 = this->l();

   ACoeff_t num = -(2.0 * l1 + 4.0 * n + 3.0) *
                  (4.0 * l1 * n + 2.0 * l1 + 4.0 * n * n + 6.0 * n + 1.0);
   ACoeff_t den = ((2.0 * l1 + 4.0 * n + 5.0) *
                   (2.0 * l1 * n + 2.0 * l1 + 2.0 * n * n + 5.0 * n + 1.0));

   return this->normalizeDiag(n, -1) * (num / den);
}

ValueD2Diags::ACoeff_t ValueD2Diags::d0(const ACoeff_t& n) const
{
   ACoeff_t val = ACoeff_t::Ones(n.size());

   return this->normalizeDiag(n, 0) * val;
}

} // namespace SphEnergy
} // namespace Stencil
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
