/**
 * @file I4D2mDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4D2mDiags
 * sparse operator
 *
 * Diagonal bodies generated from the in-tree pyquicc recurrence generator
 * (docs/agent/tools/genradops.py); see
 * docs/plans/stage4-5-generated-radial-ops.md. l-2 coupling -> normalizeDiag
 * third arg = -2. No l=0 patch (velocity, l>=1), no correctQ1/zeroLast
 * (TRUNCATE_QI=OFF).
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I4D2mDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

using namespace Internal::Literals;

I4D2mDiags::I4D2mDiags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I4D2mDiags(alpha, MHD_MP(-0.5), l, q)
{
   if (q > 0)
   {
      throw std::logic_error(
         "I4D2m: truncation (q>0) not supported - build with TRUNCATE_QI=OFF");
   }
}

I4D2mDiags::ACoeff_t I4D2mDiags::d_2(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = 64.0_mp * (l<1>() + n - 4.0_mp) * (l<1>() + n - 3.0_mp) *
         (l<1>() + n - 2.0_mp) * (l<1>() + n - 1.0_mp) /
         ((l<1>() + 2.0_mp * n - 6.0_mp) * (l<1>() + 2.0_mp * n - 5.0_mp) *
            (l<1>() + 2.0_mp * n - 4.0_mp) * (l<1>() + 2.0_mp * n - 3.0_mp) *
            (l<1>() + 2.0_mp * n - 2.0_mp) * (l<1>() + 2.0_mp * n - 1.0_mp));
   return this->normalizeDiag(n, -2, -2) * val;
}

I4D2mDiags::ACoeff_t I4D2mDiags::d_1(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -64.0_mp * (l<1>() + n - 3.0_mp) * (l<1>() + n - 2.0_mp) *
         (l<1>() + n - 1.0_mp) * (2.0_mp * l<1>() - 2.0_mp * n - 1.0_mp) /
         ((l<1>() + 2.0_mp * n - 5.0_mp) * (l<1>() + 2.0_mp * n - 4.0_mp) *
            (l<1>() + 2.0_mp * n - 3.0_mp) * (l<1>() + 2.0_mp * n - 2.0_mp) *
            (l<1>() + 2.0_mp * n - 1.0_mp) * (l<1>() + 2.0_mp * n + 1.0_mp));
   return this->normalizeDiag(n, -1, -2) * val;
}

I4D2mDiags::ACoeff_t I4D2mDiags::d0(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = 16.0_mp * (l<1>() + n - 2.0_mp) * (l<1>() + n - 1.0_mp) *
         (4.0_mp * l<2>() - 24.0_mp * l<1>() * n - 8.0_mp * l<1>() -
            4.0_mp * n * n + 24.0_mp * n + 13.0_mp) /
         ((l<1>() + 2.0_mp * n - 4.0_mp) * (l<1>() + 2.0_mp * n - 3.0_mp) *
            (l<1>() + 2.0_mp * n - 2.0_mp) * (l<1>() + 2.0_mp * n - 1.0_mp) *
            (l<1>() + 2.0_mp * n + 1.0_mp) * (l<1>() + 2.0_mp * n + 2.0_mp));
   return this->normalizeDiag(n, 0, -2) * val;
}

I4D2mDiags::ACoeff_t I4D2mDiags::d1(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = 32.0_mp * (2.0_mp * n + 1.0_mp) * (l<1>() + n - 1.0_mp) *
         (4.0_mp * l<2>() - 4.0_mp * l<1>() * n - 10.0_mp * l<1>() -
            4.0_mp * n * n + 9.0_mp) /
         ((l<1>() + 2.0_mp * n - 3.0_mp) * (l<1>() + 2.0_mp * n - 2.0_mp) *
            (l<1>() + 2.0_mp * n - 1.0_mp) * (l<1>() + 2.0_mp * n + 1.0_mp) *
            (l<1>() + 2.0_mp * n + 2.0_mp) * (l<1>() + 2.0_mp * n + 3.0_mp));
   return this->normalizeDiag(n, 1, -2) * val;
}

I4D2mDiags::ACoeff_t I4D2mDiags::d2(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = 4.0_mp * (2.0_mp * n + 1.0_mp) * (2.0_mp * n + 3.0_mp) *
         (24.0_mp * l<2>() + 16.0_mp * l<1>() * n - 32.0_mp * l<1>() -
            4.0_mp * n * n - 24.0_mp * n + 13.0_mp) /
         ((l<1>() + 2.0_mp * n - 2.0_mp) * (l<1>() + 2.0_mp * n - 1.0_mp) *
            (l<1>() + 2.0_mp * n + 1.0_mp) * (l<1>() + 2.0_mp * n + 2.0_mp) *
            (l<1>() + 2.0_mp * n + 3.0_mp) * (l<1>() + 2.0_mp * n + 4.0_mp));
   return this->normalizeDiag(n, 2, -2) * val;
}

I4D2mDiags::ACoeff_t I4D2mDiags::d3(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = 4.0_mp * (2.0_mp * n + 1.0_mp) * (2.0_mp * n + 3.0_mp) *
         (2.0_mp * n + 5.0_mp) * (2.0_mp * l<1>() + 2.0_mp * n + 1.0_mp) *
         (4.0_mp * l<1>() + 2.0_mp * n - 1.0_mp) /
         ((l<1>() + n) * (l<1>() + 2.0_mp * n - 1.0_mp) *
            (l<1>() + 2.0_mp * n + 1.0_mp) * (l<1>() + 2.0_mp * n + 2.0_mp) *
            (l<1>() + 2.0_mp * n + 3.0_mp) * (l<1>() + 2.0_mp * n + 4.0_mp) *
            (l<1>() + 2.0_mp * n + 5.0_mp));
   return this->normalizeDiag(n, 3, -2) * val;
}

I4D2mDiags::ACoeff_t I4D2mDiags::d4(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = (2.0_mp * n + 1.0_mp) * (2.0_mp * n + 3.0_mp) * (2.0_mp * n + 5.0_mp) *
         (2.0_mp * n + 7.0_mp) * (2.0_mp * l<1>() + 2.0_mp * n + 1.0_mp) *
         (2.0_mp * l<1>() + 2.0_mp * n + 3.0_mp) /
         ((l<1>() + n) * (l<1>() + n + 1.0_mp) *
            (l<1>() + 2.0_mp * n + 1.0_mp) * (l<1>() + 2.0_mp * n + 2.0_mp) *
            (l<1>() + 2.0_mp * n + 3.0_mp) * (l<1>() + 2.0_mp * n + 4.0_mp) *
            (l<1>() + 2.0_mp * n + 5.0_mp) * (l<1>() + 2.0_mp * n + 6.0_mp));
   return this->normalizeDiag(n, 4, -2) * val;
}

} // namespace Chebyshev
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
