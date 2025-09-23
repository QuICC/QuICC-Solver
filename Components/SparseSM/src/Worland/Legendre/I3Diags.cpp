/**
 * @file I3Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I3Diags sparse
 * operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Legendre/I3Diags.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Legendre {

using namespace Internal::Literals;

I3Diags::I3Diags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I3Diags(alpha, MHD_MP(-0.5), l, q)
{
   if (q > 0)
   {
      throw std::logic_error("Truncation for q>0 is not implemented");
   }
}

I3Diags::ACoeff_t I3Diags::d_3(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = 64.0_mp * (2.0_mp * l<1>() + 2.0_mp * n - 5.0_mp) *
         (2.0_mp * l<1>() + 2.0_mp * n - 3.0_mp) *
         (2.0_mp * l<1>() + 2.0_mp * n - 1.0_mp) /
         ((2.0_mp * l<1>() + 4.0_mp * n - 11.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 9.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 7.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 5.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 3.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 1.0_mp));

   return this->normalizeDiag(n, -3) * val;
}

I3Diags::ACoeff_t I3Diags::d_2(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -192.0_mp * (2.0_mp * l<1>() - 1.0_mp) *
         (2.0_mp * l<1>() + 2.0_mp * n - 3.0_mp) *
         (2.0_mp * l<1>() + 2.0_mp * n - 1.0_mp) /
         ((2.0_mp * l<1>() + 4.0_mp * n - 9.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 7.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 5.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 3.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 1.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 3.0_mp));

   return this->normalizeDiag(n, -2) * val;
}

I3Diags::ACoeff_t I3Diags::d_1(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = 192.0_mp * (2.0_mp * l<1>() + 2.0_mp * n - 1.0_mp) *
         (4.0_mp * l<2>() - 4.0_mp * l<1>() * n - 4.0_mp * l<1>() -
            4.0_mp * n * n + 2.0_mp * n + 5.0_mp) /
         ((2.0_mp * l<1>() + 4.0_mp * n - 7.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 5.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 3.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 1.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 3.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 5.0_mp));

   return this->normalizeDiag(n, -1) * val;
}

I3Diags::ACoeff_t I3Diags::d0(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -64.0_mp * (2.0_mp * l<1>() - 1.0_mp) *
         (4.0_mp * l<2>() - 24.0_mp * l<1>() * n - 16.0_mp * l<1>() -
            24.0_mp * n * n - 12.0_mp * n + 15.0_mp) /
         ((2.0_mp * l<1>() + 4.0_mp * n - 5.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 3.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 1.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 3.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 5.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 7.0_mp));

   return this->normalizeDiag(n, 0) * val;
}

I3Diags::ACoeff_t I3Diags::d1(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -384.0_mp * (n + 1.0_mp) *
         (4.0_mp * l<2>() - 4.0_mp * l<1>() * n - 8.0_mp * l<1>() -
            4.0_mp * n * n - 6.0_mp * n + 3.0_mp) /
         ((2.0_mp * l<1>() + 4.0_mp * n - 3.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n - 1.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 3.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 5.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 7.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 9.0_mp));

   return this->normalizeDiag(n, 1) * val;
}

I3Diags::ACoeff_t I3Diags::d2(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -768.0_mp * (2.0_mp * l<1>() - 1.0_mp) * (n + 1.0_mp) * (n + 2.0_mp) /
         ((2.0_mp * l<1>() + 4.0_mp * n - 1.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 3.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 5.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 7.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 9.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 11.0_mp));

   return this->normalizeDiag(n, 2) * val;
}

I3Diags::ACoeff_t I3Diags::d3(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -512.0_mp * (n + 1.0_mp) * (n + 2.0_mp) * (n + 3.0_mp) /
         ((2.0_mp * l<1>() + 4.0_mp * n + 3.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 5.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 7.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 9.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 11.0_mp) *
            (2.0_mp * l<1>() + 4.0_mp * n + 13.0_mp));

   return this->normalizeDiag(n, 3) * val;
}

} // namespace Legendre
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
