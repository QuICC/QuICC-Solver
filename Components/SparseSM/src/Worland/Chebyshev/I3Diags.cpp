/**
 * @file I3Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I3Diags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I3Diags.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   using namespace Internal::Literals;

   I3Diags::I3Diags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::I3Diags(alpha, MHD_MP(-0.5), l, q)
   {
   }

   I3Diags::ACoeff_t I3Diags::d_3(const ACoeff_t& n) const
   {
      ACoeff_t val;

      if(l<1>() == 0)
      {
         val = 1.0_mp/((2.0_mp*n - 5.0_mp)*(2.0_mp*n - 3.0_mp)*(2.0_mp*n - 1.0_mp));
         val.head(1) = 8.0_mp/((l<1>() + 3.0_mp)*(l<1>() + 4.0_mp)*(l<1>() + 5.0_mp));
      }
      else
      {
         val = 8.0_mp*(l<1>() + n - 3.0_mp)*(l<1>() + n - 2.0_mp)*(l<1>() + n - 1.0_mp)/((l<1>() + 2.0_mp*n - 6.0_mp)*(l<1>() + 2.0_mp*n - 5.0_mp)*(l<1>() + 2.0_mp*n - 4.0_mp)*(l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp));
      }

      // Truncate operator
      this->zeroLast(val, this->mQ-1);

      return this->normalizeDiag(n, -3)*val;
   }

   I3Diags::ACoeff_t I3Diags::d_2(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = -24.0_mp*l<1>()*(l<1>() + n - 2.0_mp)*(l<1>() + n - 1.0_mp)/((l<1>() + 2.0_mp*n - 5.0_mp)*(l<1>() + 2.0_mp*n - 4.0_mp)*(l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp));

      // Truncate operator
      this->zeroLast(val, this->mQ);

      return this->normalizeDiag(n, -2)*val;
   }

   I3Diags::ACoeff_t I3Diags::d_1(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = 6.0_mp*(l<1>() + n - 1.0_mp)*(4.0_mp*l<2>() - 4.0_mp*l<1>()*n + 2.0_mp*l<1>() - 4.0_mp*n*n + 4.0_mp*n + 3.0_mp)/((l<1>() + 2.0_mp*n - 4.0_mp)*(l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp));

      // Truncate operator
      this->zeroLast(val, this->mQ+1);

      return this->normalizeDiag(n, -1)*val;
   }

   I3Diags::ACoeff_t I3Diags::d0(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = -4.0_mp*l<1>()*(2.0_mp*l<2>() - 12.0_mp*l<1>()*n - 12.0_mp*n*n + 7.0_mp)/((l<1>() + 2.0_mp*n - 3.0_mp)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp));

      // Truncate operator
      this->zeroLast(val, this->mQ+2);

      return this->normalizeDiag(n, 0)*val;
   }

   I3Diags::ACoeff_t I3Diags::d1(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = -3.0_mp*(2.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)*(4.0_mp*l<2>() - 4.0_mp*l<1>()*n - 2.0_mp*l<1>() - 4.0_mp*n*n - 4.0_mp*n + 3.0_mp)/(2.0_mp*(l<1>() + n)*(l<1>() + 2.0_mp*n - 2.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp)*(l<1>() + 2.0_mp*n + 4.0_mp));

      // Truncate operator
      this->zeroLast(val, this->mQ+3);

      return this->normalizeDiag(n, 1)*val;
   }

   I3Diags::ACoeff_t I3Diags::d2(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = -3.0_mp*l<1>()*(2.0_mp*n + 1.0_mp)*(2.0_mp*n + 3.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 3.0_mp)/(2.0_mp*(l<1>() + n)*(l<1>() + n + 1.0_mp)*(l<1>() + 2.0_mp*n - 1.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp)*(l<1>() + 2.0_mp*n + 4.0_mp)*(l<1>() + 2.0_mp*n + 5.0_mp));

      // Truncate operator
      this->zeroLast(val, this->mQ+4);

      return this->normalizeDiag(n, 2)*val;
   }

   I3Diags::ACoeff_t I3Diags::d3(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = -(2.0_mp*n + 1.0_mp)*(2.0_mp*n + 3.0_mp)*(2.0_mp*n + 5.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 1.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 3.0_mp)*(2.0_mp*l<1>() + 2.0_mp*n + 5.0_mp)/(8.0_mp*(l<1>() + n)*(l<1>() + n + 1.0_mp)*(l<1>() + n + 2.0_mp)*(l<1>() + 2.0_mp*n + 1.0_mp)*(l<1>() + 2.0_mp*n + 2.0_mp)*(l<1>() + 2.0_mp*n + 3.0_mp)*(l<1>() + 2.0_mp*n + 4.0_mp)*(l<1>() + 2.0_mp*n + 5.0_mp)*(l<1>() + 2.0_mp*n + 6.0_mp));

      // Truncate operator
      this->zeroLast(val, this->mQ+5);

      return this->normalizeDiag(n, 3)*val;
   }

}
}
}
}
