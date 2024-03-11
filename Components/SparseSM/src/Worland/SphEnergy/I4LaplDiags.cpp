/**
 * @file I4LaplDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4LaplDiags
 * sparse operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/SphEnergy/I4LaplDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace SphEnergy {

   using namespace Internal::Literals;

I4LaplDiags::I4LaplDiags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I4LaplDiags(alpha, 0.5_mp, l, q),
    mI4(alpha, l, 0)
{
   if (q > 2)
   {
      throw std::logic_error("I4Lapl: truncation for q>2 is not implemented");
   }
}

I4LaplDiags::ACoeff_t I4LaplDiags::d_3(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 256.0 * (2.0 * l1 + 2.0 * n - 5.0) * (2.0 * l1 + 2.0 * n - 3.0) *
         (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 9.0) * (2.0 * l1 + 4.0 * n - 7.0) *
            (2.0 * l1 + 4.0 * n - 5.0) * (2.0 * l1 + 4.0 * n - 3.0) *
            (2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0));

   return this->normalizeDiag(n, -3) * val;
}

I4LaplDiags::ACoeff_t I4LaplDiags::d_2(const ACoeff_t& n) const
{
   auto l1 = this->l();
   auto l2 = l1 * l1;
   ACoeff_t val;

   val = -512.0 * (2.0 * l1 - 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n - 3.0) *
         (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 7.0) * (2.0 * l1 + 4.0 * n - 5.0) *
            (2.0 * l1 + 4.0 * n - 3.0) * (2.0 * l1 + 4.0 * n - 1.0) *
            (2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0));

   // Correct if q == 2
   this->correctQ2(val, n, -2);

   return this->normalizeDiag(n, -2) * val;
}

I4LaplDiags::ACoeff_t I4LaplDiags::d_1(const ACoeff_t& n) const
{
   auto l1 = this->l();
   auto l2 = l1 * l1;
   auto l3 = l2 * l1;
   ACoeff_t val;

   val =
      256.0 * (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) *
      (4.0 * l2 - 24.0 * l1 * n - 16.0 * l1 - 4.0 * n.pow(2) + 8.0 * n + 15.0) /
      ((2.0 * l1 + 4.0 * n - 5.0) * (2.0 * l1 + 4.0 * n - 3.0) *
         (2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0) *
         (2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0));

   // Correct if q == 2
   this->correctQ2(val, n, -1);

   return this->normalizeDiag(n, -1) * val;
}

I4LaplDiags::ACoeff_t I4LaplDiags::d0(const ACoeff_t& n) const
{
   auto l1 = this->l();
   auto l2 = l1 * l1;
   auto l3 = l2 * l1;
   ACoeff_t val;

   val = 2048.0 * (n + 1.0) * (2.0 * l1 + 2.0 * n + 1.0) *
         (4.0 * l2 - 4.0 * l1 * n - 8.0 * l1 - 4.0 * n.pow(2) - 6.0 * n + 3.0) /
         ((2.0 * l1 + 4.0 * n - 3.0) * (2.0 * l1 + 4.0 * n - 1.0) *
            (2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0) *
            (2.0 * l1 + 4.0 * n + 7.0) * (2.0 * l1 + 4.0 * n + 9.0));

   // Correct if q == 2
   this->correctQ2(val, n, 0);

   return this->normalizeDiag(n, 0) * val;
}

I4LaplDiags::ACoeff_t I4LaplDiags::d1(const ACoeff_t& n) const
{
   auto l1 = this->l();
   auto l2 = l1 * l1;
   ACoeff_t val;

   val = 2048.0 * (n + 1.0) * (n + 2.0) *
         (12.0 * l2 + 8.0 * l1 * n - 2.0 * n.pow(2) - 10.0 * n - 3.0) /
         ((2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0) *
            (2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0) *
            (2.0 * l1 + 4.0 * n + 9.0) * (2.0 * l1 + 4.0 * n + 11.0));

   // Correct if q == 2
   this->correctQ2(val, n, 1);

   return this->normalizeDiag(n, 1) * val;
}

I4LaplDiags::ACoeff_t I4LaplDiags::d2(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 8192.0 * (n + 1.0) * (n + 2.0) * (n + 3.0) * (2.0 * l1 + n + 1.0) /
         ((2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0) *
            (2.0 * l1 + 4.0 * n + 7.0) * (2.0 * l1 + 4.0 * n + 9.0) *
            (2.0 * l1 + 4.0 * n + 11.0) * (2.0 * l1 + 4.0 * n + 13.0));

   // Correct if q == 2
   this->correctQ2(val, n, 2);

   return this->normalizeDiag(n, 2) * val;
}

I4LaplDiags::ACoeff_t I4LaplDiags::d3(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 4096.0 * (n + 1.0) * (n + 2.0) * (n + 3.0) * (n + 4.0) /
         ((2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0) *
            (2.0 * l1 + 4.0 * n + 9.0) * (2.0 * l1 + 4.0 * n + 11.0) *
            (2.0 * l1 + 4.0 * n + 13.0) * (2.0 * l1 + 4.0 * n + 15.0));

   // Correct if q == 2
   this->correctQ2(val, n, 3);

   return this->normalizeDiag(n, 3) * val;
}

I4LaplDiags::ACoeff_t I4LaplDiags::d4(const ACoeff_t& n) const
{
   ACoeff_t val = ACoeff_t::Zero(n.size());

   // Correct if q == 2
   this->correctQ2(val, n, 4);

   return this->normalizeDiag(n, 4) * val;
}

I4LaplDiags::ACoeff_t I4LaplDiags::d5(const ACoeff_t& n) const
{
   ACoeff_t val = ACoeff_t::Zero(n.size());

   // Correct if q == 2
   this->correctQ2(val, n, 5);

   return this->normalizeDiag(n, 5) * val;
}

void I4LaplDiags::correctQ2(ACoeff_t& val, const ACoeff_t& n, const int k) const
{
   // Index where to apply correction in val
   auto i_ = val.size() - (k + 3);

   // Only correct if truncation q == 2
   if (this->mQ == 2 && i_ >= 0)
   {
      auto l1 = this->l();
      ACoeff_t m = n.bottomRows(1) + 1.0;
      // Tau coefficient obtained as the ratio of I4Lapl 3rd subdiagonal/ I4 4th
      // subdiagonal
      ACoeff_t f = (-13.0 + 2.0 * l1 + 4.0 * m) * (-11.0 + 2.0 * l1 + 4.0 * m);
      ACoeff_t nf =
         (this->normalizeDiag(m, -3) / this->normalizeDiag(m, -4)) * f;

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
      ACoeff_t ng = g / this->normalizeDiag(m, k);

      val(i_) -= (nf * ng)(0);
   }
}

} // namespace SphEnergy
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
