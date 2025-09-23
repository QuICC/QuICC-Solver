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
#include "QuICC/SparseSM/Worland/CylEnergy/I3Diags.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace CylEnergy {

using namespace Internal::Literals;

I3Diags::I3Diags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I3Diags(alpha, MHD_MP(0.0), l, q)
{
   if (q > 0)
   {
      throw std::logic_error("Truncation for q>0 is not implemented");
   }
}

I3Diags::ACoeff_t I3Diags::d_3(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = 8.0_mp * (l<1>() + n) * (l<1>() + n - 2.0_mp) * (l<1>() + n - 1.0_mp) /
         ((l<1>() + 2.0_mp * n) * (l<1>() + 2.0_mp * n - 5.0_mp) *
            (l<1>() + 2.0_mp * n - 4.0_mp) * (l<1>() + 2.0_mp * n - 3.0_mp) *
            (l<1>() + 2.0_mp * n - 2.0_mp) * (l<1>() + 2.0_mp * n - 1.0_mp));

   return this->normalizeDiag(n, -3) * val;
}

I3Diags::ACoeff_t I3Diags::d_2(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -24.0_mp * l<1>() * (l<1>() + n) * (l<1>() + n - 1.0_mp) /
         ((l<1>() + 2.0_mp * n) * (l<1>() + 2.0_mp * n - 4.0_mp) *
            (l<1>() + 2.0_mp * n - 3.0_mp) * (l<1>() + 2.0_mp * n - 2.0_mp) *
            (l<1>() + 2.0_mp * n - 1.0_mp) * (l<1>() + 2.0_mp * n + 2.0_mp));

   return this->normalizeDiag(n, -2) * val;
}

I3Diags::ACoeff_t I3Diags::d_1(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = 24.0_mp * (l<1>() + n) * (l<2>() - l<1>() * n - n * n + 1.0_mp) /
         ((l<1>() + 2.0_mp * n) * (l<1>() + 2.0_mp * n - 3.0_mp) *
            (l<1>() + 2.0_mp * n - 2.0_mp) * (l<1>() + 2.0_mp * n - 1.0_mp) *
            (l<1>() + 2.0_mp * n + 2.0_mp) * (l<1>() + 2.0_mp * n + 3.0_mp));

   return this->normalizeDiag(n, -1) * val;
}

I3Diags::ACoeff_t I3Diags::d0(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -8.0_mp * l<1>() *
         (l<2>() - 6.0_mp * l<1>() * n - 3.0_mp * l<1>() - 6.0_mp * n * n -
            6.0_mp * n + 2.0_mp) /
         ((l<1>() + 2.0_mp * n) * (l<1>() + 2.0_mp * n - 2.0_mp) *
            (l<1>() + 2.0_mp * n - 1.0_mp) * (l<1>() + 2.0_mp * n + 2.0_mp) *
            (l<1>() + 2.0_mp * n + 3.0_mp) * (l<1>() + 2.0_mp * n + 4.0_mp));

   return this->normalizeDiag(n, 0) * val;
}

I3Diags::ACoeff_t I3Diags::d1(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -24.0_mp * (n + 1.0_mp) *
         (l<2>() - l<1>() * n - l<1>() - n * n - 2.0_mp * n) /
         ((l<1>() + 2.0_mp * n) * (l<1>() + 2.0_mp * n - 1.0_mp) *
            (l<1>() + 2.0_mp * n + 2.0_mp) * (l<1>() + 2.0_mp * n + 3.0_mp) *
            (l<1>() + 2.0_mp * n + 4.0_mp) * (l<1>() + 2.0_mp * n + 5.0_mp));

   return this->normalizeDiag(n, 1) * val;
}

I3Diags::ACoeff_t I3Diags::d2(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -24.0_mp * l<1>() * (n + 1.0_mp) * (n + 2.0_mp) /
         ((l<1>() + 2.0_mp * n) * (l<1>() + 2.0_mp * n + 2.0_mp) *
            (l<1>() + 2.0_mp * n + 3.0_mp) * (l<1>() + 2.0_mp * n + 4.0_mp) *
            (l<1>() + 2.0_mp * n + 5.0_mp) * (l<1>() + 2.0_mp * n + 6.0_mp));

   return this->normalizeDiag(n, 2) * val;
}

I3Diags::ACoeff_t I3Diags::d3(const ACoeff_t& n) const
{
   ACoeff_t val;

   val = -8.0_mp * (n + 1.0_mp) * (n + 2.0_mp) * (n + 3.0_mp) /
         ((l<1>() + 2.0_mp * n + 2.0_mp) * (l<1>() + 2.0_mp * n + 3.0_mp) *
            (l<1>() + 2.0_mp * n + 4.0_mp) * (l<1>() + 2.0_mp * n + 5.0_mp) *
            (l<1>() + 2.0_mp * n + 6.0_mp) * (l<1>() + 2.0_mp * n + 7.0_mp));

   return this->normalizeDiag(n, 3) * val;
}

} // namespace CylEnergy
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
