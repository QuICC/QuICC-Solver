/**
 * @file I4D3mDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4D3mDiags
 * sparse operator
 *
 * Diagonal bodies generated from the in-tree pyquicc recurrence generator
 * (docs/agent/tools/genradops.py); see
 * docs/plans/stage4-5-generated-radial-ops.md. l-3 coupling -> normalizeDiag
 * third arg = -3. No l=0 patch (velocity, l>=1), no correctQ1/zeroLast
 * (TRUNCATE_QI=OFF).
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I4D3mDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

using namespace Internal::Literals;

I4D3mDiags::I4D3mDiags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I4D3mDiags(alpha, MHD_MP(-0.5), l, q)
{
   if (q > 0)
   {
      throw std::logic_error(
         "I4D3m: truncation (q>0) not supported - build with TRUNCATE_QI=OFF");
   }
}

I4D3mDiags::ACoeff_t I4D3mDiags::d_1(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = 128.0_mp * (l<1>() + n - 4.0_mp) * (l<1>() + n - 3.0_mp) *
         (l<1>() + n - 2.0_mp) * (l<1>() + n - 1.0_mp) /
         ((l<1>() + 2.0_mp * n - 5.0_mp) * (l<1>() + 2.0_mp * n - 4.0_mp) *
            (l<1>() + 2.0_mp * n - 3.0_mp) * (l<1>() + 2.0_mp * n - 2.0_mp) *
            (l<1>() + 2.0_mp * n - 1.0_mp));
   return this->normalizeDiag(n, -1, -3) * val;
}

I4D3mDiags::ACoeff_t I4D3mDiags::d0(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -64.0_mp * (l<1>() + n - 3.0_mp) * (l<1>() + n - 2.0_mp) *
         (l<1>() + n - 1.0_mp) * (2.0_mp * l<1>() - 6.0_mp * n - 3.0_mp) /
         ((l<1>() + 2.0_mp * n - 4.0_mp) * (l<1>() + 2.0_mp * n - 3.0_mp) *
            (l<1>() + 2.0_mp * n - 2.0_mp) * (l<1>() + 2.0_mp * n - 1.0_mp) *
            (l<1>() + 2.0_mp * n + 1.0_mp));
   return this->normalizeDiag(n, 0, -3) * val;
}

I4D3mDiags::ACoeff_t I4D3mDiags::d1(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -64.0_mp * (2.0_mp * n + 1.0_mp) * (l<1>() + n - 2.0_mp) *
         (l<1>() + n - 1.0_mp) * (4.0_mp * l<1>() - 2.0_mp * n - 7.0_mp) /
         ((l<1>() + 2.0_mp * n - 3.0_mp) * (l<1>() + 2.0_mp * n - 2.0_mp) *
            (l<1>() + 2.0_mp * n - 1.0_mp) * (l<1>() + 2.0_mp * n + 1.0_mp) *
            (l<1>() + 2.0_mp * n + 2.0_mp));
   return this->normalizeDiag(n, 1, -3) * val;
}

I4D3mDiags::ACoeff_t I4D3mDiags::d2(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -32.0_mp * (2.0_mp * n + 1.0_mp) * (2.0_mp * n + 3.0_mp) *
         (l<1>() + n - 1.0_mp) * (6.0_mp * l<1>() + 2.0_mp * n - 7.0_mp) /
         ((l<1>() + 2.0_mp * n - 2.0_mp) * (l<1>() + 2.0_mp * n - 1.0_mp) *
            (l<1>() + 2.0_mp * n + 1.0_mp) * (l<1>() + 2.0_mp * n + 2.0_mp) *
            (l<1>() + 2.0_mp * n + 3.0_mp));
   return this->normalizeDiag(n, 2, -3) * val;
}

I4D3mDiags::ACoeff_t I4D3mDiags::d3(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -8.0_mp * (2.0_mp * n + 1.0_mp) * (2.0_mp * n + 3.0_mp) *
         (2.0_mp * n + 5.0_mp) * (8.0_mp * l<1>() + 6.0_mp * n - 3.0_mp) /
         ((l<1>() + 2.0_mp * n - 1.0_mp) * (l<1>() + 2.0_mp * n + 1.0_mp) *
            (l<1>() + 2.0_mp * n + 2.0_mp) * (l<1>() + 2.0_mp * n + 3.0_mp) *
            (l<1>() + 2.0_mp * n + 4.0_mp));
   return this->normalizeDiag(n, 3, -3) * val;
}

I4D3mDiags::ACoeff_t I4D3mDiags::d4(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -4.0_mp * (2.0_mp * n + 1.0_mp) * (2.0_mp * n + 3.0_mp) *
         (2.0_mp * n + 5.0_mp) * (2.0_mp * n + 7.0_mp) *
         (2.0_mp * l<1>() + 2.0_mp * n + 1.0_mp) /
         ((l<1>() + n) * (l<1>() + 2.0_mp * n + 1.0_mp) *
            (l<1>() + 2.0_mp * n + 2.0_mp) * (l<1>() + 2.0_mp * n + 3.0_mp) *
            (l<1>() + 2.0_mp * n + 4.0_mp) * (l<1>() + 2.0_mp * n + 5.0_mp));
   return this->normalizeDiag(n, 4, -3) * val;
}

} // namespace Chebyshev
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
