/**
 * @file I4QmDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4QmDiags
 * sparse operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/SphEnergy/I4QmDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace SphEnergy {

   using namespace Internal::Literals;

I4QmDiags::I4QmDiags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I4QmDiags(alpha, 0.5_mp, l, q),
    mI4(alpha, l, 0)
{
   if (q > 2)
   {
      throw std::logic_error("I4Qm: Truncation for q != 2 is not implemented");
   }
}

I4QmDiags::ACoeff_t I4QmDiags::d_3(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = -256.0 * (2.0 * l1 + 2.0 * n - 5.0) * (2.0 * l1 + 2.0 * n - 3.0) *
         (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 11.0) * (2.0 * l1 + 4.0 * n - 9.0) *
            (2.0 * l1 + 4.0 * n - 7.0) * (2.0 * l1 + 4.0 * n - 5.0) *
            (2.0 * l1 + 4.0 * n - 3.0) * (2.0 * l1 + 4.0 * n - 1.0) *
            (2.0 * l1 + 4.0 * n + 1.0));

   return this->normalizeDiag(n, -3, -1) * val;
}

I4QmDiags::ACoeff_t I4QmDiags::d_2(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 256.0 * (2.0 * l1 + 2.0 * n - 3.0) * (2.0 * l1 + 2.0 * n - 1.0) *
         (2.0 * l1 + 2.0 * n + 1.0) * (6.0 * l1 - 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 9.0) * (2.0 * l1 + 4.0 * n - 7.0) *
            (2.0 * l1 + 4.0 * n - 5.0) * (2.0 * l1 + 4.0 * n - 3.0) *
            (2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0) *
            (2.0 * l1 + 4.0 * n + 5.0));

   // Correct if q == 2
   this->correctQ2(val, n, -2);

   return this->normalizeDiag(n, -2, -1) * val;
}

I4QmDiags::ACoeff_t I4QmDiags::d_1(const ACoeff_t& n) const
{
   const auto l1 = this->l();
   const auto l2 = l1 * l1;
   ACoeff_t val;

   val = -768.0 * (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) *
         (4.0 * l2 - 8.0 * l1 * n - 4.0 * n.pow(2) + 7.0) /
         ((2.0 * l1 + 4.0 * n - 7.0) * (2.0 * l1 + 4.0 * n - 5.0) *
            (2.0 * l1 + 4.0 * n - 3.0) * (2.0 * l1 + 4.0 * n - 1.0) *
            (2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0) *
            (2.0 * l1 + 4.0 * n + 7.0));

   // Correct if q == 2
   this->correctQ2(val, n, -1);

   return this->normalizeDiag(n, -1, -1) * val;
}

I4QmDiags::ACoeff_t I4QmDiags::d0(const ACoeff_t& n) const
{
   const auto l1 = this->l();
   const auto l2 = l1 * l1;
   const auto l3 = l2 * l1;
   ACoeff_t val;

   val = 256.0 * (2.0 * l1 + 2.0 * n + 1.0) *
         (8.0 * l3 - 72.0 * l2 * n - 36.0 * l2 - 24.0 * l1 * n.pow(2) -
            24.0 * l1 * n + 46.0 * l1 + 24.0 * n.pow(3) + 36.0 * n.pow(2) -
            18.0 * n - 15.0) /
         ((2.0 * l1 + 4.0 * n - 5.0) * (2.0 * l1 + 4.0 * n - 3.0) *
            (2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0) *
            (2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0) *
            (2.0 * l1 + 4.0 * n + 9.0));

   // Correct if q == 2
   this->correctQ2(val, n, 0);

   return this->normalizeDiag(n, 0, -1) * val;
}

I4QmDiags::ACoeff_t I4QmDiags::d1(const ACoeff_t& n) const
{
   const auto l1 = this->l();
   const auto l2 = l1 * l1;
   const auto l3 = l2 * l1;
   ACoeff_t val;

   val = 2048.0 * (n + 1.0) *
         (8.0 * l3 - 12.0 * l2 * n - 12.0 * l2 - 24.0 * l1 * n.pow(2) -
            48.0 * l1 * n - 2.0 * l1 - 6.0 * n.pow(3) - 18.0 * n.pow(2) -
            9.0 * n + 3.0) /
         ((2.0 * l1 + 4.0 * n - 3.0) * (2.0 * l1 + 4.0 * n - 1.0) *
            (2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0) *
            (2.0 * l1 + 4.0 * n + 7.0) * (2.0 * l1 + 4.0 * n + 9.0) *
            (2.0 * l1 + 4.0 * n + 11.0));

   // Correct if q == 2
   this->correctQ2(val, n, 1);

   return this->normalizeDiag(n, 1, -1) * val;
}

I4QmDiags::ACoeff_t I4QmDiags::d2(const ACoeff_t& n) const
{
   const auto l1 = this->l();
   const auto l2 = l1 * l1;
   ACoeff_t val;

   val = 6144.0 * (n + 1.0) * (n + 2.0) *
         (4.0 * l2 - 2.0 * n.pow(2) - 6.0 * n - 1.0) /
         ((2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0) *
            (2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0) *
            (2.0 * l1 + 4.0 * n + 9.0) * (2.0 * l1 + 4.0 * n + 11.0) *
            (2.0 * l1 + 4.0 * n + 13.0));

   // Correct if q == 2
   this->correctQ2(val, n, 2);

   return this->normalizeDiag(n, 2, -1) * val;
}

I4QmDiags::ACoeff_t I4QmDiags::d3(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 4096.0 * (n + 1.0) * (n + 2.0) * (n + 3.0) * (4.0 * l1 + n + 2.0) /
         ((2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0) *
            (2.0 * l1 + 4.0 * n + 7.0) * (2.0 * l1 + 4.0 * n + 9.0) *
            (2.0 * l1 + 4.0 * n + 11.0) * (2.0 * l1 + 4.0 * n + 13.0) *
            (2.0 * l1 + 4.0 * n + 15.0));

   // Correct if q == 2
   this->correctQ2(val, n, 3);

   return this->normalizeDiag(n, 3, -1) * val;
}

I4QmDiags::ACoeff_t I4QmDiags::d4(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 4096.0 * (n + 1.0) * (n + 2.0) * (n + 3.0) * (n + 4.0) /
         ((2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0) *
            (2.0 * l1 + 4.0 * n + 9.0) * (2.0 * l1 + 4.0 * n + 11.0) *
            (2.0 * l1 + 4.0 * n + 13.0) * (2.0 * l1 + 4.0 * n + 15.0) *
            (2.0 * l1 + 4.0 * n + 17.0));

   // Correct if q == 2
   this->correctQ2(val, n, 4);

   return this->normalizeDiag(n, 4, -1) * val;
}

I4QmDiags::ACoeff_t I4QmDiags::d5(const ACoeff_t& n) const
{
   ACoeff_t val = ACoeff_t::Zero(n.size());

   // Correct if q == 2
   this->correctQ2(val, n, 5);

   return this->normalizeDiag(n, 5, -1) * val;
}

void I4QmDiags::correctQ2(ACoeff_t& val, const ACoeff_t& n, const int k) const
{
   // Index where to apply correction in val
   auto i_ = val.size() - (k + 3);

   // Only correct if truncation q == 2
   if (this->mQ == 2 && i_ >= 0)
   {
      auto l1 = this->l();
      ACoeff_t m = n.bottomRows(1) - 1;
      // Tau matrix entry: I2Qm(-1,-1)/I2(-1,-2)
      ACoeff_t f = -(2.0 * l1 + 4.0 * m - 5.0);
      ACoeff_t nf =
         (this->normalizeDiag(m, -1, -1) / this->normalizeDiag(m, -2)) * f;

      m = n.bottomRows(1) - static_cast<Scalar_t>(k + 2);
      ACoeff_t g;
      switch (k)
      {
      case -2:
         g = this->mI4.d_3(m);
         break;
      case -1:
         g = this->mI4.d_2(m);
         break;
      case 0:
         g = this->mI4.d_1(m);
         break;
      case 1:
         g = this->mI4.d0(m);
         break;
      case 2:
         g = this->mI4.d1(m);
         break;
      case 3:
         g = this->mI4.d2(m);
         break;
      case 4:
         g = this->mI4.d3(m);
         break;
      case 5:
         g = this->mI4.d4(m);
         break;
      default:
         throw std::logic_error("Unknown diagonal for computing correction");
         break;
      }
      ACoeff_t ng = g / this->normalizeDiag(m, k, -1);

      val(i_) -= (nf * ng)(0);
   }
}

} // namespace SphEnergy
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
