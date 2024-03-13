/**
 * @file I4Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I4Diags sparse
 * operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/SphEnergy/I4Diags.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace SphEnergy {

   using namespace Internal::Literals;

I4Diags::I4Diags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I4Diags(alpha, 0.5_mp, l, q)
{}

I4Diags::ACoeff_t I4Diags::d_4(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 256.0 * (2.0 * l1 + 2.0 * n - 5.0) * (2.0 * l1 + 2.0 * n - 3.0) *
         (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 13.0) * (2.0 * l1 + 4.0 * n - 11.0) *
            (2.0 * l1 + 4.0 * n - 9.0) * (2.0 * l1 + 4.0 * n - 7.0) *
            (2.0 * l1 + 4.0 * n - 5.0) * (2.0 * l1 + 4.0 * n - 3.0) *
            (2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0));

   // Truncate operator
   this->zeroLast(val, this->mQ - 2);

   return this->normalizeDiag(n, -4) * val;
}

I4Diags::ACoeff_t I4Diags::d_3(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = -1024.0 * (2.0 * l1 + 1.0) * (2.0 * l1 + 2.0 * n - 3.0) *
         (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 11.0) * (2.0 * l1 + 4.0 * n - 9.0) *
            (2.0 * l1 + 4.0 * n - 7.0) * (2.0 * l1 + 4.0 * n - 5.0) *
            (2.0 * l1 + 4.0 * n - 3.0) * (2.0 * l1 + 4.0 * n - 1.0) *
            (2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0));

   // Truncate operator
   this->zeroLast(val, this->mQ - 1);

   return this->normalizeDiag(n, -3) * val;
}

I4Diags::ACoeff_t I4Diags::d_2(const ACoeff_t& n) const
{
   auto l1 = this->l();
   auto l2 = Internal::Math::pow(l1, 2);
   ACoeff_t val;

   val =
      512.0 * (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) *
      (12.0 * l2 - 8.0 * l1 * n + 16.0 * l1 - 8.0 * n.pow(2) + 4.0 * n + 21.0) /
      ((2.0 * l1 + 4.0 * n - 9.0) * (2.0 * l1 + 4.0 * n - 7.0) *
         (2.0 * l1 + 4.0 * n - 5.0) * (2.0 * l1 + 4.0 * n - 3.0) *
         (2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0) *
         (2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0));

   // Truncate operator
   this->zeroLast(val, this->mQ);

   return this->normalizeDiag(n, -2) * val;
}

I4Diags::ACoeff_t I4Diags::d_1(const ACoeff_t& n) const
{
   auto l1 = this->l();
   auto l2 = Internal::Math::pow(l1, 2);
   ACoeff_t val;

   val =
      -1024.0 * (2.0 * l1 + 1.0) * (2.0 * l1 + 2.0 * n + 1.0) *
      (4.0 * l2 - 12.0 * l1 * n + 4.0 * l1 - 12.0 * n.pow(2) - 6.0 * n + 21.0) /
      ((2.0 * l1 + 4.0 * n - 7.0) * (2.0 * l1 + 4.0 * n - 5.0) *
         (2.0 * l1 + 4.0 * n - 3.0) * (2.0 * l1 + 4.0 * n - 1.0) *
         (2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0) *
         (2.0 * l1 + 4.0 * n + 7.0) * (2.0 * l1 + 4.0 * n + 9.0));

   // Truncate operator
   this->zeroLast(val, this->mQ + 1);

   return this->normalizeDiag(n, -1) * val;
}

I4Diags::ACoeff_t I4Diags::d0(const ACoeff_t& n) const
{
   auto l1 = this->l();
   auto l2 = Internal::Math::pow(l1, 2);
   auto l3 = Internal::Math::pow(l1, 3);
   auto l4 = Internal::Math::pow(l1, 4);
   ACoeff_t val;

   val =
      256.0 *
      (16.0 * l4 - 192.0 * l3 * n - 64.0 * l3 - 96.0 * l2 * n.pow(2) -
         384.0 * l2 * n + 56.0 * l2 + 192.0 * l1 * n.pow(3) +
         192.0 * l1 * n.pow(2) - 336.0 * l1 * n + 16.0 * l1 + 96.0 * n.pow(4) +
         288.0 * n.pow(3) + 24.0 * n.pow(2) - 288.0 * n - 15.0) /
      ((2.0 * l1 + 4.0 * n - 5.0) * (2.0 * l1 + 4.0 * n - 3.0) *
         (2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0) *
         (2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0) *
         (2.0 * l1 + 4.0 * n + 9.0) * (2.0 * l1 + 4.0 * n + 11.0));

   // Truncate operator
   this->zeroLast(val, this->mQ + 2);

   return this->normalizeDiag(n, 0) * val;
}

I4Diags::ACoeff_t I4Diags::d1(const ACoeff_t& n) const
{
   auto l1 = this->l();
   auto l2 = Internal::Math::pow(l1, 2);
   ACoeff_t val;

   val =
      2048.0 * (2.0 * l1 + 1.0) * (n + 1.0) *
      (4.0 * l2 - 12.0 * l1 * n - 8.0 * l1 - 12.0 * n.pow(2) - 30.0 * n + 3.0) /
      ((2.0 * l1 + 4.0 * n - 3.0) * (2.0 * l1 + 4.0 * n - 1.0) *
         (2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0) *
         (2.0 * l1 + 4.0 * n + 7.0) * (2.0 * l1 + 4.0 * n + 9.0) *
         (2.0 * l1 + 4.0 * n + 11.0) * (2.0 * l1 + 4.0 * n + 13.0));

   // Truncate operator
   this->zeroLast(val, this->mQ + 3);

   return this->normalizeDiag(n, 1) * val;
}

I4Diags::ACoeff_t I4Diags::d2(const ACoeff_t& n) const
{
   auto l1 = this->l();
   auto l2 = Internal::Math::pow(l1, 2);
   ACoeff_t val;

   val = 2048.0 * (n + 1.0) * (n + 2.0) *
         (12.0 * l2 - 8.0 * l1 * n - 8.0 * n.pow(2) - 28.0 * n - 3.0) /
         ((2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0) *
            (2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0) *
            (2.0 * l1 + 4.0 * n + 9.0) * (2.0 * l1 + 4.0 * n + 11.0) *
            (2.0 * l1 + 4.0 * n + 13.0) * (2.0 * l1 + 4.0 * n + 15.0));

   // Truncate operator
   this->zeroLast(val, this->mQ + 4);

   return this->normalizeDiag(n, 2) * val;
}

I4Diags::ACoeff_t I4Diags::d3(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 8192.0 * (2.0 * l1 + 1.0) * (n + 1.0) * (n + 2.0) * (n + 3.0) /
         ((2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0) *
            (2.0 * l1 + 4.0 * n + 7.0) * (2.0 * l1 + 4.0 * n + 9.0) *
            (2.0 * l1 + 4.0 * n + 11.0) * (2.0 * l1 + 4.0 * n + 13.0) *
            (2.0 * l1 + 4.0 * n + 15.0) * (2.0 * l1 + 4.0 * n + 17.0));

   // Truncate operator
   this->zeroLast(val, this->mQ + 5);

   return this->normalizeDiag(n, 3) * val;
}

I4Diags::ACoeff_t I4Diags::d4(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 4096.0 * (n + 1.0) * (n + 2.0) * (n + 3.0) * (n + 4.0) /
         ((2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0) *
            (2.0 * l1 + 4.0 * n + 9.0) * (2.0 * l1 + 4.0 * n + 11.0) *
            (2.0 * l1 + 4.0 * n + 13.0) * (2.0 * l1 + 4.0 * n + 15.0) *
            (2.0 * l1 + 4.0 * n + 17.0) * (2.0 * l1 + 4.0 * n + 19.0));

   // Truncate operator
   this->zeroLast(val, this->mQ + 6);

   return this->normalizeDiag(n, 4) * val;
}

} // namespace SphEnergy
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
