/**
 * @file I4QpDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4QpDiags
 * sparse operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/SphEnergy/I4QpDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace SphEnergy {

   using namespace Internal::Literals;

I4QpDiags::I4QpDiags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I4QpDiags(alpha, 0.5_mp, l, q),
    mI4(alpha, l, 0)
{
   if (q > 2)
   {
      throw std::logic_error("I4Qp: Truncation for q != 2 is not implemented");
   }
}

I4QpDiags::ACoeff_t I4QpDiags::d_4(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 256.0 * (2.0 * l1 + 2.0 * n - 5.0) * (2.0 * l1 + 2.0 * n - 3.0) *
         (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 11.0) * (2.0 * l1 + 4.0 * n - 9.0) *
            (2.0 * l1 + 4.0 * n - 7.0) * (2.0 * l1 + 4.0 * n - 5.0) *
            (2.0 * l1 + 4.0 * n - 3.0) * (2.0 * l1 + 4.0 * n - 1.0) *
            (2.0 * l1 + 4.0 * n + 1.0));

   return this->normalizeDiag(n, -4, 1) * val;
}

I4QpDiags::ACoeff_t I4QpDiags::d_3(const ACoeff_t& n) const
{
   const auto l1 = this->l();
   const auto l2 = l1 * l1;
   ACoeff_t val;

   val = -256.0 * (2.0 * l1 + 2.0 * n - 3.0) * (2.0 * l1 + 2.0 * n - 1.0) *
         (2.0 * l1 + 2.0 * n + 1.0) * (6.0 * l1 - 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 9.0) * (2.0 * l1 + 4.0 * n - 7.0) *
            (2.0 * l1 + 4.0 * n - 5.0) * (2.0 * l1 + 4.0 * n - 3.0) *
            (2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0) *
            (2.0 * l1 + 4.0 * n + 5.0));

   // Correct if q == 2
   this->correctQ2(val, n, -3);

   return this->normalizeDiag(n, -3, 1) * val;
}

I4QpDiags::ACoeff_t I4QpDiags::d_2(const ACoeff_t& n) const
{
   const auto l1 = this->l();
   const auto l2 = l1 * l1;
   const auto l3 = l2 * l1;
   ACoeff_t val;

   val = 768.0 * (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) *
         (4.0 * l2 - 8.0 * l1 * n - 4.0 * n.pow(2) + 7.0) /
         ((2.0 * l1 + 4.0 * n - 7.0) * (2.0 * l1 + 4.0 * n - 5.0) *
            (2.0 * l1 + 4.0 * n - 3.0) * (2.0 * l1 + 4.0 * n - 1.0) *
            (2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0) *
            (2.0 * l1 + 4.0 * n + 7.0));

   // Correct if q == 2
   this->correctQ2(val, n, -2);

   return this->normalizeDiag(n, -2, 1) * val;
}

I4QpDiags::ACoeff_t I4QpDiags::d_1(const ACoeff_t& n) const
{
   const auto l1 = this->l();
   const auto l2 = l1 * l1;
   const auto l3 = l2 * l1;
   const auto l4 = l2 * l2;
   ACoeff_t val;

   val = -256.0 * (2.0 * l1 + 2.0 * n + 1.0) *
         (8.0 * l3 - 72.0 * l2 * n - 36.0 * l2 - 24.0 * l1 * n.pow(2) -
            24.0 * l1 * n + 46.0 * l1 + 24.0 * n.pow(3) + 36.0 * n.pow(2) -
            18.0 * n - 15.0) /
         ((2.0 * l1 + 4.0 * n - 5.0) * (2.0 * l1 + 4.0 * n - 3.0) *
            (2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0) *
            (2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0) *
            (2.0 * l1 + 4.0 * n + 9.0));

   // Correct if q == 2
   this->correctQ2(val, n, -1);

   return this->normalizeDiag(n, -1, 1) * val;
}

I4QpDiags::ACoeff_t I4QpDiags::d0(const ACoeff_t& n) const
{
   const auto l1 = this->l();
   const auto l2 = l1 * l1;
   const auto l3 = l2 * l1;
   ACoeff_t val;

   val = -2048.0 * (n + 1.0) *
         (8.0 * l3 - 12.0 * l2 * n - 12.0 * l2 - 24.0 * l1 * n.pow(2) -
            48.0 * l1 * n - 2.0 * l1 - 6.0 * n.pow(3) - 18.0 * n.pow(2) -
            9.0 * n + 3.0) /
         ((2.0 * l1 + 4.0 * n - 3.0) * (2.0 * l1 + 4.0 * n - 1.0) *
            (2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0) *
            (2.0 * l1 + 4.0 * n + 7.0) * (2.0 * l1 + 4.0 * n + 9.0) *
            (2.0 * l1 + 4.0 * n + 11.0));

   // Correct if q == 2
   this->correctQ2(val, n, 0);

   return this->normalizeDiag(n, 0, 1) * val;
}

I4QpDiags::ACoeff_t I4QpDiags::d1(const ACoeff_t& n) const
{
   const auto l1 = this->l();
   const auto l2 = l1 * l1;
   ACoeff_t val;

   val = -6144.0 * (n + 1.0) * (n + 2.0) *
         (4.0 * l2 - 2.0 * n.pow(2) - 6.0 * n - 1.0) /
         ((2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0) *
            (2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0) *
            (2.0 * l1 + 4.0 * n + 9.0) * (2.0 * l1 + 4.0 * n + 11.0) *
            (2.0 * l1 + 4.0 * n + 13.0));

   // Correct if q == 2
   this->correctQ2(val, n, 1);

   return this->normalizeDiag(n, 1, 1) * val;
}

I4QpDiags::ACoeff_t I4QpDiags::d2(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = -4096.0 * (n + 1.0) * (n + 2.0) * (n + 3.0) * (4.0 * l1 + n + 2.0) /
         ((2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0) *
            (2.0 * l1 + 4.0 * n + 7.0) * (2.0 * l1 + 4.0 * n + 9.0) *
            (2.0 * l1 + 4.0 * n + 11.0) * (2.0 * l1 + 4.0 * n + 13.0) *
            (2.0 * l1 + 4.0 * n + 15.0));

   // Correct if q == 2
   this->correctQ2(val, n, 2);

   return this->normalizeDiag(n, 2, 1) * val;
}

I4QpDiags::ACoeff_t I4QpDiags::d3(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = -4096.0 * (n + 1.0) * (n + 2.0) * (n + 3.0) * (n + 4.0) /
         ((2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0) *
            (2.0 * l1 + 4.0 * n + 9.0) * (2.0 * l1 + 4.0 * n + 11.0) *
            (2.0 * l1 + 4.0 * n + 13.0) * (2.0 * l1 + 4.0 * n + 15.0) *
            (2.0 * l1 + 4.0 * n + 17.0));

   // Correct if q == 2
   this->correctQ2(val, n, 3);

   return this->normalizeDiag(n, 3, 1) * val;
}

I4QpDiags::ACoeff_t I4QpDiags::d4(const ACoeff_t& n) const
{
   ACoeff_t val = ACoeff_t::Zero(n.size());

   // Correct if q == 2
   this->correctQ2(val, n, 4);

   return this->normalizeDiag(n, 4, 1) * val;
}

I4QpDiags::ACoeff_t I4QpDiags::d5(const ACoeff_t& n) const
{
   ACoeff_t val = ACoeff_t::Zero(n.size());

   // Correct if q == 2
   this->correctQ2(val, n, 5);

   return this->normalizeDiag(n, 5, 1) * val;
}

void I4QpDiags::correctQ2(ACoeff_t& val, const ACoeff_t& n, const int k) const
{
   // Only correct if truncation q == 2
   if (this->mQ == 2)
   {
      auto l1 = this->l();
      const auto l2 = l1 * l1;

      // Truncation requires 3 coefficients from tau operator T:
      // a = T[-2,-2], b = T[-1,-1], c = T[-2,-1]
      ACoeff_t m = n.bottomRows(1) - 1;
      ACoeff_t f = (2.0 * l1 + 4.0 * m - 5.0);
      ACoeff_t ncA =
         (this->normalizeDiag(m, -2, 1) / this->normalizeDiag(m, -2)) * f;

      m = n.bottomRows(1);
      f = (2.0 * l1 + 4.0 * m - 5.0);
      ACoeff_t ncB =
         (this->normalizeDiag(m, -2, 1) / this->normalizeDiag(m, -2)) * f;

      // compute T[-2,-1] term, and simplify common factors
      m = n.bottomRows(1) - 1;
      ACoeff_t num = -(2.0 * l1 - 2.0 * m - 1.0) / ((2.0 * l1 + 4.0 * m + 5.0));
      num *= this->normalizeDiag(m, -1, 1);
      ACoeff_t prod = -2.0 * (2.0 * l1 + 1.0) /
                      ((2.0 * l1 + 4.0 * m - 3.0) * (2.0 * l1 + 4.0 * m + 5.0));
      prod *= this->normalizeDiag(m, -1);
      ACoeff_t den = (2.0 * l1 + 2.0 * m - 1.0) /
                     ((2.0 * l1 + 4.0 * m - 5.0) * (2.0 * l1 + 4.0 * m - 3.0));
      den *= this->normalizeDiag(m, -2);

      ACoeff_t ncC = (num - ncB * prod) / den;
      ;

      m = n.bottomRows(1) - static_cast<Scalar_t>(k + 3);
      ACoeff_t m1 = n.bottomRows(1) - static_cast<Scalar_t>(k + 2);
      ACoeff_t g, h1, h2, ng;
      switch (k)
      {
      case -3:
         g = this->mI4.d_3(m);
         break;
      case -2:
         g = this->mI4.d_2(m);
         h1 = this->mI4.d_3(m1);
         h2 = this->mI4.d_2(m1);
         break;
      case -1:
         g = this->mI4.d_1(m);
         h1 = this->mI4.d_2(m1);
         h2 = this->mI4.d_1(m1);
         break;
      case 0:
         g = this->mI4.d0(m);
         h1 = this->mI4.d_1(m1);
         h2 = this->mI4.d0(m1);
         break;
      case 1:
         g = this->mI4.d1(m);
         h1 = this->mI4.d0(m1);
         h2 = this->mI4.d1(m1);
         break;
      case 2:
         g = this->mI4.d2(m);
         h1 = this->mI4.d1(m1);
         h2 = this->mI4.d2(m1);
         break;
      case 3:
         g = this->mI4.d3(m);
         h1 = this->mI4.d2(m1);
         h2 = this->mI4.d3(m1);
         break;
      case 4:
         g = this->mI4.d4(m);
         h1 = this->mI4.d3(m1);
         h2 = this->mI4.d4(m1);
         break;
      case 5:
         h1 = this->mI4.d4(m1);
         h2 = 0 * m1;
         break;
      default:
         throw std::logic_error("Unknown diagonal for computing correction");
         break;
      }
      // First part of correction
      auto i_ = val.size() - (k + 4);
      if (k < 5 && i_ >= 0)
      {
         ng = (ncA * g) / this->normalizeDiag(m, k, 1);
         val(i_) -= ng(0);
      }
      // Second part of correction
      auto j_ = val.size() - (k + 3);
      if (k > -3 && j_ >= 0)
      {
         ng = (ncB * h2 + ncC * h1) / this->normalizeDiag(m1, k, 1);
         val(j_) -= ng(0);
      }
   }
}

} // namespace SphEnergy
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
