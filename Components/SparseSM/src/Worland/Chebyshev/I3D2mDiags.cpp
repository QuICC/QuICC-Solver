/**
 * @file I3D2mDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I3D2mDiags
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
#include "QuICC/SparseSM/Worland/Chebyshev/I3D2mDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

using namespace Internal::Literals;

I3D2mDiags::I3D2mDiags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I3D2mDiags(alpha, MHD_MP(-0.5), l, q)
{
   if (q > 0)
   {
      throw std::logic_error(
         "I3D2m: truncation (q>0) not supported - build with TRUNCATE_QI=OFF");
   }
}

I3D2mDiags::ACoeff_t I3D2mDiags::d_1(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = 32.0_mp * (l<1>() + n - 3.0_mp) * (l<1>() + n - 2.0_mp) *
         (l<1>() + n - 1.0_mp) /
         ((l<1>() + 2.0_mp * n - 4.0_mp) * (l<1>() + 2.0_mp * n - 3.0_mp) *
            (l<1>() + 2.0_mp * n - 2.0_mp) * (l<1>() + 2.0_mp * n - 1.0_mp));
   return this->normalizeDiag(n, -1, -2) * val;
}

I3D2mDiags::ACoeff_t I3D2mDiags::d0(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -32.0_mp * (l<1>() - 2.0_mp * n - 1.0_mp) * (l<1>() + n - 2.0_mp) *
         (l<1>() + n - 1.0_mp) /
         ((l<1>() + 2.0_mp * n - 3.0_mp) * (l<1>() + 2.0_mp * n - 2.0_mp) *
            (l<1>() + 2.0_mp * n - 1.0_mp) * (l<1>() + 2.0_mp * n + 1.0_mp));
   return this->normalizeDiag(n, 0, -2) * val;
}

I3D2mDiags::ACoeff_t I3D2mDiags::d1(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -48.0_mp * (l<1>() - 1.0_mp) * (2.0_mp * n + 1.0_mp) *
         (l<1>() + n - 1.0_mp) /
         ((l<1>() + 2.0_mp * n - 2.0_mp) * (l<1>() + 2.0_mp * n - 1.0_mp) *
            (l<1>() + 2.0_mp * n + 1.0_mp) * (l<1>() + 2.0_mp * n + 2.0_mp));
   return this->normalizeDiag(n, 1, -2) * val;
}

I3D2mDiags::ACoeff_t I3D2mDiags::d2(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -8.0_mp * (2.0_mp * n + 1.0_mp) * (2.0_mp * n + 3.0_mp) *
         (3.0_mp * l<1>() + 2.0_mp * n - 1.0_mp) /
         ((l<1>() + 2.0_mp * n - 1.0_mp) * (l<1>() + 2.0_mp * n + 1.0_mp) *
            (l<1>() + 2.0_mp * n + 2.0_mp) * (l<1>() + 2.0_mp * n + 3.0_mp));
   return this->normalizeDiag(n, 2, -2) * val;
}

I3D2mDiags::ACoeff_t I3D2mDiags::d3(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -2.0_mp * (2.0_mp * n + 1.0_mp) * (2.0_mp * n + 3.0_mp) *
         (2.0_mp * n + 5.0_mp) * (2.0_mp * l<1>() + 2.0_mp * n + 1.0_mp) /
         ((l<1>() + n) * (l<1>() + 2.0_mp * n + 1.0_mp) *
            (l<1>() + 2.0_mp * n + 2.0_mp) * (l<1>() + 2.0_mp * n + 3.0_mp) *
            (l<1>() + 2.0_mp * n + 4.0_mp));
   return this->normalizeDiag(n, 3, -2) * val;
}

} // namespace Chebyshev
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
