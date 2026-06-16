/**
 * @file I4laplDmDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4laplDmDiags
 * sparse operator
 *
 * Diagonal bodies generated from the in-tree pyquicc recurrence generator
 * (docs/agent/tools/genradops.py); see
 * docs/plans/stage4-5-generated-radial-ops.md. l-1 coupling -> normalizeDiag
 * third arg = -1. No l=0 patch (velocity, l>=1), no correctQ1/zeroLast
 * (TRUNCATE_QI=OFF).
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I4laplDmDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

using namespace Internal::Literals;

I4laplDmDiags::I4laplDmDiags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I4laplDmDiags(alpha, MHD_MP(-0.5), l, q)
{
   if (q > 0)
   {
      throw std::logic_error("I4laplDm: truncation (q>0) not supported - build "
                             "with TRUNCATE_QI=OFF");
   }
}

I4laplDmDiags::ACoeff_t I4laplDmDiags::d_2(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = 64.0_mp * (l<1>() + n - 3.0_mp) * (l<1>() + n - 2.0_mp) *
         (l<1>() + n - 1.0_mp) * (2.0_mp * l<1>() + 2.0_mp * n - 5.0_mp) /
         ((l<1>() + 2.0_mp * n - 5.0_mp) * (l<1>() + 2.0_mp * n - 4.0_mp) *
            (l<1>() + 2.0_mp * n - 3.0_mp) * (l<1>() + 2.0_mp * n - 2.0_mp) *
            (l<1>() + 2.0_mp * n - 1.0_mp));
   return this->normalizeDiag(n, -2, -1) * val;
}

I4laplDmDiags::ACoeff_t I4laplDmDiags::d_1(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -32.0_mp * (l<1>() + n - 2.0_mp) * (l<1>() + n - 1.0_mp) *
         (4.0_mp * l<2>() - 8.0_mp * l<1>() * n - 8.0_mp * l<1>() -
            12.0_mp * n * n + 20.0_mp * n + 13.0_mp) /
         ((l<1>() + 2.0_mp * n - 4.0_mp) * (l<1>() + 2.0_mp * n - 3.0_mp) *
            (l<1>() + 2.0_mp * n - 2.0_mp) * (l<1>() + 2.0_mp * n - 1.0_mp) *
            (l<1>() + 2.0_mp * n + 1.0_mp));
   return this->normalizeDiag(n, -1, -1) * val;
}

I4laplDmDiags::ACoeff_t I4laplDmDiags::d0(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -32.0_mp * (l<1>() + n - 1.0_mp) *
         (16.0_mp * l<2>() * n + 4.0_mp * l<2>() + 8.0_mp * l<1>() * n * n -
            32.0_mp * l<1>() * n - 22.0_mp * l<1>() - 8.0_mp * n * n * n -
            12.0_mp * n * n + 26.0_mp * n + 15.0_mp) /
         ((l<1>() + 2.0_mp * n - 3.0_mp) * (l<1>() + 2.0_mp * n - 2.0_mp) *
            (l<1>() + 2.0_mp * n - 1.0_mp) * (l<1>() + 2.0_mp * n + 1.0_mp) *
            (l<1>() + 2.0_mp * n + 2.0_mp));
   return this->normalizeDiag(n, 0, -1) * val;
}

I4laplDmDiags::ACoeff_t I4laplDmDiags::d1(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -16.0_mp * (2.0_mp * n + 1.0_mp) *
         (24.0_mp * l<2>() * n + 24.0_mp * l<2>() + 32.0_mp * l<1>() * n * n +
            8.0_mp * l<1>() * n - 48.0_mp * l<1>() + 8.0_mp * n * n * n -
            12.0_mp * n * n - 26.0_mp * n + 15.0_mp) /
         ((l<1>() + 2.0_mp * n - 2.0_mp) * (l<1>() + 2.0_mp * n - 1.0_mp) *
            (l<1>() + 2.0_mp * n + 1.0_mp) * (l<1>() + 2.0_mp * n + 2.0_mp) *
            (l<1>() + 2.0_mp * n + 3.0_mp));
   return this->normalizeDiag(n, 1, -1) * val;
}

I4laplDmDiags::ACoeff_t I4laplDmDiags::d2(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -4.0_mp * (2.0_mp * n + 1.0_mp) * (2.0_mp * n + 3.0_mp) *
         (2.0_mp * l<1>() + 2.0_mp * n + 1.0_mp) *
         (16.0_mp * l<1>() * n + 28.0_mp * l<1>() + 12.0_mp * n * n +
            20.0_mp * n - 13.0_mp) /
         ((l<1>() + n) * (l<1>() + 2.0_mp * n - 1.0_mp) *
            (l<1>() + 2.0_mp * n + 1.0_mp) * (l<1>() + 2.0_mp * n + 2.0_mp) *
            (l<1>() + 2.0_mp * n + 3.0_mp) * (l<1>() + 2.0_mp * n + 4.0_mp));
   return this->normalizeDiag(n, 2, -1) * val;
}

I4laplDmDiags::ACoeff_t I4laplDmDiags::d3(const ACoeff_t& n) const
{
   ACoeff_t val;
   val =
      -2.0_mp * (2.0_mp * n + 1.0_mp) * (2.0_mp * n + 3.0_mp) *
      (2.0_mp * n + 5.0_mp) * (2.0_mp * n + 5.0_mp) *
      (2.0_mp * l<1>() + 2.0_mp * n + 1.0_mp) *
      (2.0_mp * l<1>() + 2.0_mp * n + 3.0_mp) /
      ((l<1>() + n) * (l<1>() + n + 1.0_mp) * (l<1>() + 2.0_mp * n + 1.0_mp) *
         (l<1>() + 2.0_mp * n + 2.0_mp) * (l<1>() + 2.0_mp * n + 3.0_mp) *
         (l<1>() + 2.0_mp * n + 4.0_mp) * (l<1>() + 2.0_mp * n + 5.0_mp));
   return this->normalizeDiag(n, 3, -1) * val;
}

} // namespace Chebyshev
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
