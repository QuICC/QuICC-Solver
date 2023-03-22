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
#include "QuICC/SparseSM/Worland/Legendre/I4Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Legendre {

   I4Diags::I4Diags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::I4Diags(alpha, MHD_MP(-0.5), l, q)
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

      val = 256.0*(2.0*l1 + 2.0*n - 7.0)*(2.0*l1 + 2.0*n - 5.0)*(2.0*l1 + 2.0*n - 3.0)*(2.0*l1 + 2.0*n - 1.0)/((2.0*l1 + 4.0*n - 15.0)*(2.0*l1 + 4.0*n - 13.0)*(2.0*l1 + 4.0*n - 11.0)*(2.0*l1 + 4.0*n - 9.0)*(2.0*l1 + 4.0*n - 7.0)*(2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0));

      return this->normalizeDiag(n, -4)*val;
   }

   I4Diags::ACoeff_t I4Diags::d_3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = -1024.0*(2.0*l1 - 1.0)*(2.0*l1 + 2.0*n - 5.0)*(2.0*l1 + 2.0*n - 3.0)*(2.0*l1 + 2.0*n - 1.0)/((2.0*l1 + 4.0*n - 13.0)*(2.0*l1 + 4.0*n - 11.0)*(2.0*l1 + 4.0*n - 9.0)*(2.0*l1 + 4.0*n - 7.0)*(2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 3.0));

      return this->normalizeDiag(n, -3)*val;
   }

   I4Diags::ACoeff_t I4Diags::d_2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = 512.0*(2.0*l1 + 2.0*n - 3.0)*(2.0*l1 + 2.0*n - 1.0)*(12.0*l2 - 8.0*l1*n - 8.0*l1 - 8.0*n.pow(2) + 12.0*n + 17.0)/((2.0*l1 + 4.0*n - 11.0)*(2.0*l1 + 4.0*n - 9.0)*(2.0*l1 + 4.0*n - 7.0)*(2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 3.0)*(2.0*l1 + 4.0*n + 5.0));

      return this->normalizeDiag(n, -2)*val;
   }

   I4Diags::ACoeff_t I4Diags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = -1024.0*(2.0*l1 - 1.0)*(2.0*l1 + 2.0*n - 1.0)*(4.0*l2 - 12.0*l1*n - 4.0*l1 - 12.0*n.pow(2) + 6.0*n + 21.0)/((2.0*l1 + 4.0*n - 9.0)*(2.0*l1 + 4.0*n - 7.0)*(2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 3.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0));

      return this->normalizeDiag(n, -1)*val;
   }

   I4Diags::ACoeff_t I4Diags::d0(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      auto l3 = precision::pow(l1, 3);
      auto l4 = precision::pow(l1, 4);
      ACoeff_t val;

      val = 256.0*(16.0*l4 - 192.0*l3*n - 128.0*l3 - 96.0*l2*n.pow(2) + 192.0*l2*n + 344.0*l2 + 192.0*l1*n.pow(3) + 384.0*l1*n.pow(2) - 144.0*l1*n - 352.0*l1 + 96.0*n.pow(4) + 96.0*n.pow(3) - 264.0*n.pow(2) - 144.0*n + 105.0)/((2.0*l1 + 4.0*n - 7.0)*(2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 3.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0));

      return this->normalizeDiag(n, 0)*val;
   }

   I4Diags::ACoeff_t I4Diags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = 2048.0*(2.0*l1 - 1.0)*(n + 1.0)*(4.0*l2 - 12.0*l1*n - 16.0*l1 - 12.0*n.pow(2) - 18.0*n + 15.0)/((2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 3.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0));

      return this->normalizeDiag(n, 1)*val;
   }

   I4Diags::ACoeff_t I4Diags::d2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = 2048.0*(n + 1.0)*(n + 2.0)*(12.0*l2 - 8.0*l1*n - 24.0*l1 - 8.0*n.pow(2) - 20.0*n + 9.0)/((2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 3.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0)*(2.0*l1 + 4.0*n + 13.0));

      return this->normalizeDiag(n, 2)*val;
   }

   I4Diags::ACoeff_t I4Diags::d3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 8192.0*(2.0*l1 - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0)/((2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 3.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0)*(2.0*l1 + 4.0*n + 13.0)*(2.0*l1 + 4.0*n + 15.0));

      return this->normalizeDiag(n, 3)*val;
   }

   I4Diags::ACoeff_t I4Diags::d4(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 4096.0*(n + 1.0)*(n + 2.0)*(n + 3.0)*(n + 4.0)/((2.0*l1 + 4.0*n + 3.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0)*(2.0*l1 + 4.0*n + 13.0)*(2.0*l1 + 4.0*n + 15.0)*(2.0*l1 + 4.0*n + 17.0));

      return this->normalizeDiag(n, 4)*val;
   }

}
}
}
}
