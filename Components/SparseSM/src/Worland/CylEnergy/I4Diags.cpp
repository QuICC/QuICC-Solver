/** 
 * @file I4Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I4Diags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/CylEnergy/I4Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace CylEnergy {

   I4Diags::I4Diags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::I4Diags(alpha, MHD_MP(0.0), l, q)
   {
      if(q > 0)
      {
         throw std::logic_error("Truncation for q>0 is not implemented");
      }
   }

   I4Diags::ACoeff_t I4Diags::d_4(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 16.0*(l1 + n)*(l1 + n - 3.0)*(l1 + n - 2.0)*(l1 + n - 1.0)/((l1 + 2.0*n)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0));

      return this->normalizeDiag(n, -4)*val;
   }

   I4Diags::ACoeff_t I4Diags::d_3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = -64.0*l1*(l1 + n)*(l1 + n - 2.0)*(l1 + n - 1.0)/((l1 + 2.0*n)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0));

      return this->normalizeDiag(n, -3)*val;
   }

   I4Diags::ACoeff_t I4Diags::d_2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = 32.0*(l1 + n)*(l1 + n - 1.0)*(3.0*l2 - 2.0*l1*n + l1 - 2.0*n.pow(2) + 2.0*n + 4.0)/((l1 + 2.0*n)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0));

      return this->normalizeDiag(n, -2)*val;
   }

   I4Diags::ACoeff_t I4Diags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = -64.0*l1*(l1 + n)*(l2 - 3.0*l1*n - 3.0*n.pow(2) + 5.0)/((l1 + 2.0*n)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0));

      return this->normalizeDiag(n, -1)*val;
   }

   I4Diags::ACoeff_t I4Diags::d0(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      auto l3 = precision::pow(l1, 3);
      auto l4 = precision::pow(l1, 4);
      ACoeff_t val;

      val = 16.0*(l4 - 12.0*l3*n - 6.0*l3 - 6.0*l2*n.pow(2) - 6.0*l2*n + 11.0*l2 + 12.0*l1*n.pow(3) + 18.0*l1*n.pow(2) - 6.0*l1*n - 6.0*l1 + 6.0*n.pow(4) + 12.0*n.pow(3) - 6.0*n.pow(2) - 12.0*n)/((l1 + 2.0*n)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0));

      return this->normalizeDiag(n, 0)*val;
   }

   I4Diags::ACoeff_t I4Diags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = 64.0*l1*(n + 1.0)*(l2 - 3.0*l1*n - 3.0*l1 - 3.0*n.pow(2) - 6.0*n + 2.0)/((l1 + 2.0*n)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0));

      return this->normalizeDiag(n, 1)*val;
   }

   I4Diags::ACoeff_t I4Diags::d2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = 32.0*(n + 1.0)*(n + 2.0)*(3.0*l2 - 2.0*l1*n - 3.0*l1 - 2.0*n.pow(2) - 6.0*n)/((l1 + 2.0*n)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0));

      return this->normalizeDiag(n, 2)*val;
   }

   I4Diags::ACoeff_t I4Diags::d3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 64.0*l1*(n + 1.0)*(n + 2.0)*(n + 3.0)/((l1 + 2.0*n)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0));

      return this->normalizeDiag(n, 3)*val;
   }

   I4Diags::ACoeff_t I4Diags::d4(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 16.0*(n + 1.0)*(n + 2.0)*(n + 3.0)*(n + 4.0)/((l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0)*(l1 + 2.0*n + 9.0));

      return this->normalizeDiag(n, 4)*val;
   }

}
}
}
}
