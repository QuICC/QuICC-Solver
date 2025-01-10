/**
 * @file I3LaplDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I3LaplDiags
 * sparse operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/SphEnergy/I3LaplDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace SphEnergy {

   using namespace Internal::Literals;

I3LaplDiags::I3LaplDiags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I3LaplDiags(alpha, 0.5_mp, l, q),
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

   val = 64.0_mp*(2.0_mp*l<1>() + 2.0_mp*n - 3.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n - 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)/((2.0_mp*l<1>() + 4.0_mp*n - 5.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n - 3.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n - 1.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n + 1.0_mp));

   // Correct if q == 2
   this->correctQ2(val, n, -2);

   return this->normalizeDiag(n, -2) * val;
}

I3LaplDiags::ACoeff_t I3LaplDiags::d_1(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -64.0_mp*(2.0_mp*l<1>() - 4.0_mp*n - 3.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n - 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)/((2.0_mp*l<1>() + 4.0_mp*n - 3.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n - 1.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n + 5.0_mp));

   // Correct if q == 2
   this->correctQ2(val, n, -1);

   return this->normalizeDiag(n, -1) * val;
}

I3LaplDiags::ACoeff_t I3LaplDiags::d0(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -384.0_mp*(2.0_mp*l<1>() - 1.0_mp)*(n + 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)/((2.0_mp*l<1>() + 4.0_mp*n - 1.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n + 5.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n + 7.0_mp));

   // Correct if q == 2
   this->correctQ2(val, n, 0);

   return this->normalizeDiag(n, 0) * val;
}

I3LaplDiags::ACoeff_t I3LaplDiags::d1(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -256.0_mp*(n + 1.0_mp)*(n + 2.0_mp)*(6.0_mp*l<1>() + 4.0_mp*n + 3.0_mp)/((2.0_mp*l<1>() + 4.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n + 5.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n + 7.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n + 9.0_mp));

   // Correct if q == 2
   this->correctQ2(val, n, 1);

   return this->normalizeDiag(n, 1) * val;
}

I3LaplDiags::ACoeff_t I3LaplDiags::d2(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -512.0_mp*(n + 1.0_mp)*(n + 2.0_mp)*(n + 3.0_mp)/((2.0_mp*l<1>() + 4.0_mp*n + 5.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n + 7.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n + 9.0_mp)*(2.0_mp*l<1>() + 4.0_mp*n + 11.0_mp));

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
      // Tau coefficient obtained as the ratio of I3Lapl 3rd subdiagonal/ I3 4th
      // subdiagonal
      ACoeff_t f = (-13.0 + 2.0 * l1 + 4.0 * m) * (-11.0 + 2.0 * l1 + 4.0 * m);
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

} // namespace SphEnergy
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
