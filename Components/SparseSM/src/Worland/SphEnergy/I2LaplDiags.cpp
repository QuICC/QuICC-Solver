/**
 * @file I2LaplDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I2LaplDiags
 * sparse operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/SphEnergy/I2LaplDiags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace SphEnergy {

   using namespace Internal::Literals;

I2LaplDiags::I2LaplDiags(const Scalar_t alpha, const int l, const int q) :
    QuICC::SparseSM::Worland::I2LaplDiags(alpha, 0.5_mp, l, q)
{
   // q <= 1 is equivalent to no truncation (already zero rows)

   if (q > 1)
   {
      throw std::logic_error("I2Lapl: Truncation for q>1 is not implemented");
   }
}

I2LaplDiags::ACoeff_t I2LaplDiags::d_1(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 16.0 * (2.0 * l1 + 2.0 * n - 1.0) * (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n - 1.0) * (2.0 * l1 + 4.0 * n + 1.0));

   return this->normalizeDiag(n, -1) * val;
}

I2LaplDiags::ACoeff_t I2LaplDiags::d0(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 64.0 * (n + 1.0) * (2.0 * l1 + 2.0 * n + 1.0) /
         ((2.0 * l1 + 4.0 * n + 1.0) * (2.0 * l1 + 4.0 * n + 5.0));

   return this->normalizeDiag(n, 0) * val;
}

I2LaplDiags::ACoeff_t I2LaplDiags::d1(const ACoeff_t& n) const
{
   auto l1 = this->l();
   ACoeff_t val;

   val = 64.0 * (n + 1.0) * (n + 2.0) /
         ((2.0 * l1 + 4.0 * n + 5.0) * (2.0 * l1 + 4.0 * n + 7.0));

   return this->normalizeDiag(n, 1) * val;
}

} // namespace SphEnergy
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
