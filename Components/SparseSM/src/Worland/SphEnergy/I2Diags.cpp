/** 
 * @file I2.cpp
 * @brief Source of the implementation of the full sphere Worland I2 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/SphEnergy/I2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace SphEnergy {

   I2Diags::I2Diags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::I2Diags(alpha, MHD_MP(0.5), l, q)
   {
   }

   I2Diags::ACoeff_t I2Diags::d_2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 16.0*(2.0*l1 + 2.0*n - 1.0)*(2.0*l1 + 2.0*n + 1.0)/((2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0));

      // Truncate operator
      this->zeroLast(val, this->mQ-1);

      return this->normalizeDiag(n, -2)*val;
   }

   I2Diags::ACoeff_t I2Diags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = -32.0*(2.0*l1 + 1.0)*(2.0*l1 + 2.0*n + 1.0)/((2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0));

      // Truncate operator
      this->zeroLast(val, this->mQ);

      return this->normalizeDiag(n, -1)*val;
   }

   I2Diags::ACoeff_t I2Diags::d0(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = 16.0*(4.0*l2 - 8.0*l1*n - 8.0*n.pow(2) - 12.0*n - 1.0)/((2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0));

      // Truncate operator
      this->zeroLast(val, this->mQ+1);

      return this->normalizeDiag(n, 0)*val;
   }

   I2Diags::ACoeff_t I2Diags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 64.0*(2.0*l1 + 1.0)*(n + 1.0)/((2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0));

      // Truncate operator
      this->zeroLast(val, this->mQ+2);

      return this->normalizeDiag(n, 1)*val;
   }

   I2Diags::ACoeff_t I2Diags::d2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 64.0*(n + 1.0)*(n + 2.0)/((2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0));

      // Truncate operator
      this->zeroLast(val, this->mQ+3);

      return this->normalizeDiag(n, 2)*val;
   }

}
}
}
}
