/**
 * @file I3LaplDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I3LaplDiags
 * sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I3LaplDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

using namespace Internal::Literals;

I3LaplDiags::I3LaplDiags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I3LaplDiags(alpha, MHD_MP(-0.5), l, q),
    mI3(alpha, l, 0)
{
   if (q > 2)
   {
      throw std::logic_error("I3Lapl: truncation for q>2 is not implemented");
   }
}

I3LaplDiags::ACoeff_t I3LaplDiags::d_2(const ACoeff_t& n) const
{

   ACoeff_t val;

   val = 16.0_mp * (l<1>() + n - 2.0_mp) * (l<1>() + n - 1.0_mp) *
         (2.0_mp * l<1>() + 2.0_mp * n - 3.0_mp) /
         ((l<1>() + 2.0_mp * n - 4.0_mp) * (l<1>() + 2.0_mp * n - 3.0_mp) *
            (l<1>() + 2.0_mp * n - 2.0_mp) * (l<1>() + 2.0_mp * n - 1.0_mp));

   // Correct if q == 2
   this->correctQ2(val, n, -2);

   return this->normalizeDiag(n, -2) * val;
}

I3LaplDiags::ACoeff_t I3LaplDiags::d_1(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -16.0_mp * (l<1>() + n - 1.0_mp) *
         (2.0_mp * l<2>() - 2.0_mp * l<1>() * n - l<1>() - 4.0_mp * n * n +
            4.0_mp * n + 3.0_mp) /
         ((l<1>() + 2.0_mp * n - 3.0_mp) * (l<1>() + 2.0_mp * n - 2.0_mp) *
            (l<1>() + 2.0_mp * n - 1.0_mp) * (l<1>() + 2.0_mp * n + 1.0_mp));

   // Correct if q == 2
   this->correctQ2(val, n, -1);

   return this->normalizeDiag(n, -1) * val;
}

I3LaplDiags::ACoeff_t I3LaplDiags::d0(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -8.0_mp *
         (12.0_mp * l<2>() * n + 2.0_mp * l<2>() + 12.0_mp * l<1>() * n * n -
            4.0_mp * l<1>() * n - 9.0_mp * l<1>() - 4.0_mp * n * n + 1.0_mp) /
         ((l<1>() + 2.0_mp * n - 2.0_mp) * (l<1>() + 2.0_mp * n - 1.0_mp) *
            (l<1>() + 2.0_mp * n + 1.0_mp) * (l<1>() + 2.0_mp * n + 2.0_mp));

   // Correct if q == 2
   this->correctQ2(val, n, 0);

   return this->normalizeDiag(n, 0) * val;
}

I3LaplDiags::ACoeff_t I3LaplDiags::d1(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -4.0_mp * (2.0_mp * n + 1.0_mp) *
         (2.0_mp * l<1>() + 2.0_mp * n + 1.0_mp) *
         (6.0_mp * l<1>() * n + 5.0_mp * l<1>() + 4.0_mp * n * n + 4.0_mp * n -
            3.0_mp) /
         ((l<1>() + n) * (l<1>() + 2.0_mp * n - 1.0_mp) *
            (l<1>() + 2.0_mp * n + 1.0_mp) * (l<1>() + 2.0_mp * n + 2.0_mp) *
            (l<1>() + 2.0_mp * n + 3.0_mp));

   // Correct if q == 2
   this->correctQ2(val, n, 1);

   return this->normalizeDiag(n, 1) * val;
}

I3LaplDiags::ACoeff_t I3LaplDiags::d2(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -(2.0_mp * n + 1.0_mp) * (2.0_mp * n + 3.0_mp).pow(2) *
         (2.0_mp * l<1>() + 2.0_mp * n + 1.0_mp) *
         (2.0_mp * l<1>() + 2.0_mp * n + 3.0_mp) /
         ((l<1>() + n) * (l<1>() + n + 1.0_mp) *
            (l<1>() + 2.0_mp * n + 1.0_mp) * (l<1>() + 2.0_mp * n + 2.0_mp) *
            (l<1>() + 2.0_mp * n + 3.0_mp) * (l<1>() + 2.0_mp * n + 4.0_mp));

   // Correct if q == 2
   this->correctQ2(val, n, 2);

   return this->normalizeDiag(n, 2) * val;
}

I3LaplDiags::ACoeff_t I3LaplDiags::d3(const ACoeff_t& n) const
{
   ACoeff_t val = ACoeff_t::Zero(n.size());

   // Correct if q == 2
   this->correctQ2(val, n, 3);

   return this->normalizeDiag(n, 3) * val;
}

I3LaplDiags::ACoeff_t I3LaplDiags::d4(const ACoeff_t& n) const
{
   ACoeff_t val = ACoeff_t::Zero(n.size());

   // Correct if q == 2
   this->correctQ2(val, n, 4);

   return this->normalizeDiag(n, 4) * val;
}

void I3LaplDiags::correctQ2(ACoeff_t& val, const ACoeff_t& n, const int k) const
{
   // Index where to apply correction in val
   auto i_ = val.size() - (k + 3);

   // Only correct if truncation q == 2
   if (this->mQ == 2 && i_ >= 0)
   {
      auto l1 = this->l();
      ACoeff_t m = n.bottomRows(1) + 1.0;
      ACoeff_t f = 2.0 * (-8.0 + l1 + 2.0 * m) * (-7.0 + l1 + 2.0 * m) *
                   (-5.0 + 2.0 * l1 + 2.0 * m) / (-4.0 + l1 + m);
      ACoeff_t nf =
         (this->normalizeDiag(m, -3) / this->normalizeDiag(m, -4)) * f;

      m = n.bottomRows(1) - static_cast<Scalar_t>(k + 2);
      ACoeff_t g;
      switch (k)
      {
      case -1:
         g = this->mI3.d_2(m);
         break;
      case 0:
         g = this->mI3.d_1(m);
         break;
      case 1:
         g = this->mI3.d0(m);
         break;
      case 2:
         g = this->mI3.d1(m);
         break;
      case 3:
         g = this->mI3.d2(m);
         break;
      case 4:
         g = this->mI3.d3(m);
         break;
      default:
         throw std::logic_error("Unknown diagonal for computing correction");
         break;
      }
      ACoeff_t ng = g / this->normalizeDiag(m, k);

      val(i_) -= (nf * ng)(0);
   }
}

} // namespace Chebyshev
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
