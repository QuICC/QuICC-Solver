/**
 * @file I2QpDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I2QpDiags
 * sparse operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/SphEnergy/I2QpDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace SphEnergy {

   using namespace Internal::Literals;

I2QpDiags::I2QpDiags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I2QpDiags(alpha, 0.5_mp, l, q),
    mI2(alpha, l, 0)
{
   if (q > 1)
   {
      throw std::logic_error("I2Qp: Truncation for q>1 is not implemented");
   }
}

I2QpDiags::ACoeff_t I2QpDiags::d_2(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 16.0 * (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 3.0) * (2.0 * l1 + 4.0 * n - 1.0) *
            (2.0 * l1 + 4.0 * n + 1.0));

   return this->normalizeDiag(n, -2, 1) * val;
}

I2QpDiags::ACoeff_t I2QpDiags::d_1(const ACoeff_t& n) const
{
   const auto l1 = this->l();
   ACoeff_t val;

   val = -16.0 * (2.0 * l1 - 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0) *
            (2.0 * l1 + 4.0 * n + 5.0));

   // Correct if truncation q == 1
   this->correctQ1(val, n, -1);

   return this->normalizeDiag(n, -1, 1) * val;
}

I2QpDiags::ACoeff_t I2QpDiags::d0(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = -64.0 * (n + 1.0) * (2.0 * l1 + n + 1.0) /
         ((2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0) *
            (2.0 * l1 + 4.0 * n + 7.0));

   // Correct if truncation q == 1
   this->correctQ1(val, n, 0);

   return this->normalizeDiag(n, 0, 1) * val;
}

I2QpDiags::ACoeff_t I2QpDiags::d1(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = -64.0 * (n + 1.0) * (n + 2.0) /
         ((2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0) *
            (2.0 * l1 + 4.0 * n + 9.0));

   // Correct if truncation q == 1
   this->correctQ1(val, n, 1);

   return this->normalizeDiag(n, 1, 1) * val;
}

I2QpDiags::ACoeff_t I2QpDiags::d2(const ACoeff_t& n) const
{
   ACoeff_t val = ACoeff_t::Zero(n.size());

   // Correct if truncation q == 1
   this->correctQ1(val, n, 2);

   return this->normalizeDiag(n, 2, 1) * val;
}

void I2QpDiags::correctQ1(ACoeff_t& val, const ACoeff_t& n, const int k) const
{
   // Index where to apply correction in val
   auto i_ = val.size() - (k + 2);

   // Only correct if truncation q == 1
   if (this->mQ == 1 && i_ >= 0)
   {
      auto l1 = this->l();
      ACoeff_t m = n.bottomRows(1) + 1.0;
      // Tau coefficient obtained as the ratio of I2QP 2nd subdiagonal/ I2 2nd
      // subdiagonal
      ACoeff_t f = (2.0 * l1 + 4.0 * m - 5.0);
      ACoeff_t nf =
         (this->normalizeDiag(m, -2, 1) / this->normalizeDiag(m, -2)) * f;

      m = n.bottomRows(1) - static_cast<Scalar_t>(k + 1);
      ACoeff_t g;
      switch (k)
      {
      case -1:
         g = this->mI2.d_1(m);
         break;
      case 0:
         g = this->mI2.d0(m);
         break;
      case 1:
         g = this->mI2.d1(m);
         break;
      case 2:
         g = this->mI2.d2(m);
         break;
      default:
         throw std::logic_error("Unknown diagonal for computing correction");
         break;
      }
      ACoeff_t ng = g / this->normalizeDiag(m, k, 1);

      val(i_) -= (nf * ng)(0);
   }
}

} // namespace SphEnergy
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
