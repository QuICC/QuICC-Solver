/**
 * @file I4D2pDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4D2pDiags
 * sparse operator
 *
 * Diagonal bodies generated from the in-tree pyquicc recurrence generator
 * (docs/agent/tools/genradops.py); see docs/plans/stage4-5-generated-radial-ops.md.
 * l+2 coupling -> normalizeDiag third arg = 2. No l=0 patch (velocity, l>=1),
 * no correctQ1/zeroLast (TRUNCATE_QI=OFF). Bodies are verbatim generator output
 * (single line) -- run clang-format to wrap if desired.
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I4D2pDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

using namespace Internal::Literals;

I4D2pDiags::I4D2pDiags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I4D2pDiags(alpha, MHD_MP(-0.5), l, q)
{
   if (q > 0)
   {
      throw std::logic_error(
         "I4D2p: truncation (q>0) not supported - build with TRUNCATE_QI=OFF");
   }
}

I4D2pDiags::ACoeff_t I4D2pDiags::d_4(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = 16.0_mp*(l<1>() + n - 2.0_mp)*(l<1>() + n - 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n - 5.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n - 3.0_mp)/((l<1>() + 2.0_mp*n - 6.0_mp)*(l<1>() + 2.0_mp*n - 5.0_mp)*(l<1>() + 2.0_mp*n - 4.0_mp)*(l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp));
   return this->normalizeDiag(n, -4, 2) * val;
}

I4D2pDiags::ACoeff_t I4D2pDiags::d_3(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = -16.0_mp*(l<1>() + n - 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n - 3.0_mp)*(4.0_mp*l<2>() - 4.0_mp*n*n + 8.0_mp*n + 5.0_mp)/((l<1>() + 2.0_mp*n - 5.0_mp)*(l<1>() + 2.0_mp*n - 4.0_mp)*(l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp));
   return this->normalizeDiag(n, -3, 2) * val;
}

I4D2pDiags::ACoeff_t I4D2pDiags::d_2(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = 4.0_mp*(2.0_mp*l<1>() + 2.0_mp*n - 1.0_mp)*(8.0_mp*l<3>() - 40.0_mp*l<2>()*n + 20.0_mp*l<2>() - 56.0_mp*l<1>()*n*n + 56.0_mp*l<1>()*n + 74.0_mp*l<1>() - 8.0_mp*n*n*n + 12.0_mp*n*n + 2.0_mp*n - 3.0_mp)/((l<1>() + 2.0_mp*n - 4.0_mp)*(l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp));
   return this->normalizeDiag(n, -2, 2) * val;
}

I4D2pDiags::ACoeff_t I4D2pDiags::d_1(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = 8.0_mp*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)*(16.0_mp*l<3>()*n - 24.0_mp*l<2>() - 32.0_mp*l<1>()*n*n*n + 40.0_mp*l<1>()*n - 16.0_mp*n*n*n*n + 40.0_mp*n*n - 9.0_mp)/((l<1>() + n)*(l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp));
   return this->normalizeDiag(n, -1, 2) * val;
}

I4D2pDiags::ACoeff_t I4D2pDiags::d0(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = (2.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 3.0_mp)*(48.0_mp*l<2>()*n + 24.0_mp*l<2>() + 32.0_mp*l<1>()*n*n + 32.0_mp*l<1>()*n - 72.0_mp*l<1>() - 8.0_mp*n*n*n - 12.0_mp*n*n + 2.0_mp*n + 3.0_mp)/((l<1>() + n)*(l<1>() + n + 1.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp)*(l<1>() + 2.0_mp*n + 4.0_mp));
   return this->normalizeDiag(n, 0, 2) * val;
}

I4D2pDiags::ACoeff_t I4D2pDiags::d1(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = (2.0_mp*n + 1.0_mp)*(2.0_mp*n + 3.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 3.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 5.0_mp)*(8.0_mp*l<1>()*n + 8.0_mp*l<1>() + 4.0_mp*n*n + 8.0_mp*n - 5.0_mp)/((l<1>() + n)*(l<1>() + n + 1.0_mp)*(l<1>() + n + 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp)*(l<1>() + 2.0_mp*n + 4.0_mp)*(l<1>() + 2.0_mp*n + 5.0_mp));
   return this->normalizeDiag(n, 1, 2) * val;
}

I4D2pDiags::ACoeff_t I4D2pDiags::d2(const ACoeff_t& n) const
{
   ACoeff_t val;
   val = (2.0_mp*n + 1.0_mp)*(2.0_mp*n + 3.0_mp)*(2.0_mp*n + 3.0_mp)*(2.0_mp*n + 5.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 3.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 5.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 7.0_mp)/(4.0_mp*(l<1>() + n)*(l<1>() + n + 1.0_mp)*(l<1>() + n + 2.0_mp)*(l<1>() + n + 3.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp)*(l<1>() + 2.0_mp*n + 4.0_mp)*(l<1>() + 2.0_mp*n + 5.0_mp)*(l<1>() + 2.0_mp*n + 6.0_mp));
   return this->normalizeDiag(n, 2, 2) * val;
}

} // namespace Chebyshev
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
