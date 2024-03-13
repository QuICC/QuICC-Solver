/**
 * @file I4Lapl2Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I4Lapl2Diags
 * sparse operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/SphEnergy/I4Lapl2Diags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace SphEnergy {

   using namespace Internal::Literals;

I4Lapl2Diags::I4Lapl2Diags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I4Lapl2Diags(alpha, 0.5_mp, l, q)
{
   // q <= 2 is equivalent to no truncation (already zero rows)

   if (q > 2)
   {
      throw std::logic_error("I4Lapl2: Truncation for q>2 is not implemented");
   }
}

I4Lapl2Diags::ACoeff_t I4Lapl2Diags::d_2(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 256.0 * (2.0 * l1 + 2.0 * n - 5.0) * (2.0 * l1 + 2.0 * n - 3.0) *
         (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 5.0) * (2.0 * l1 + 4.0 * n - 3.0) *
            (2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0));

   return this->normalizeDiag(n, -2) * val;
}

I4Lapl2Diags::ACoeff_t I4Lapl2Diags::d_1(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 2048.0 * (n + 1.0) * (2.0 * l1 + 2.0 * n - 3.0) *
         (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 3.0) * (2.0 * l1 + 4.0 * n - 1.0) *
            (2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0));

   return this->normalizeDiag(n, -1) * val;
}

I4Lapl2Diags::ACoeff_t I4Lapl2Diags::d0(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 6144.0 * (n + 1.0) * (n + 2.0) * (2.0 * l1 + 2.0 * n - 1.0) *
         (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0) *
            (2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0));

   return this->normalizeDiag(n, 0) * val;
}

I4Lapl2Diags::ACoeff_t I4Lapl2Diags::d1(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 8192.0 * (n + 1.0) * (n + 2.0) * (n + 3.0) *
         (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0) *
            (2.0 * l1 + 4.0 * n + 7.0) * (2.0 * l1 + 4.0 * n + 9.0));

   return this->normalizeDiag(n, 1) * val;
}

I4Lapl2Diags::ACoeff_t I4Lapl2Diags::d2(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 4096.0 * (n + 1.0) * (n + 2.0) * (n + 3.0) * (n + 4.0) /
         ((2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0) *
            (2.0 * l1 + 4.0 * n + 9.0) * (2.0 * l1 + 4.0 * n + 11.0));

   return this->normalizeDiag(n, 2) * val;
}

} // namespace SphEnergy
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
