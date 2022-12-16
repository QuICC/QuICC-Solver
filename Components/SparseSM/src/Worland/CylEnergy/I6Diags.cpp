/** 
 * @file I6Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I6Diags sparse operator
 */

// Systethis->l() includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/CylEnergy/I6Diags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace CylEnergy {

   I6Diags::I6Diags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::I6Diags(alpha, MHD_MP(0.0), l, q)
   {
      if(q > 0)
      {
         throw std::logic_error("Truncation for q>0 is not implemented");
      }
   }

   I6Diags::ACoeff_t I6Diags::d_6(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 64.0*(l1 + n)*(l1 + n - 5.0)*(l1 + n - 4.0)*(l1 + n - 3.0)*(l1 + n - 2.0)*(l1 + n - 1.0)/((l1 + 2.0*n)*(l1 + 2.0*n - 11.0)*(l1 + 2.0*n - 10.0)*(l1 + 2.0*n - 9.0)*(l1 + 2.0*n - 8.0)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0));

      return this->normalizeDiag(n, -6)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_5(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = -384.0*l1*(l1 + n)*(l1 + n - 4.0)*(l1 + n - 3.0)*(l1 + n - 2.0)*(l1 + n - 1.0)/((l1 + 2.0*n)*(l1 + 2.0*n - 10.0)*(l1 + 2.0*n - 9.0)*(l1 + 2.0*n - 8.0)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0));

      return this->normalizeDiag(n, -5)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_4(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = 192.0*(l1 + n)*(l1 + n - 3.0)*(l1 + n - 2.0)*(l1 + n - 1.0)*(5.0*l2 - 2.0*l1*n + 3.0*l1 - 2.0*n.pow(2) + 6.0*n + 8.0)/((l1 + 2.0*n)*(l1 + 2.0*n - 9.0)*(l1 + 2.0*n - 8.0)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0));

      return this->normalizeDiag(n, -4)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = -640.0*l1*(l1 + n)*(l1 + n - 2.0)*(l1 + n - 1.0)*(2.0*l2 - 3.0*l1*n + 3.0*l1 - 3.0*n.pow(2) + 6.0*n + 13.0)/((l1 + 2.0*n)*(l1 + 2.0*n - 8.0)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0));

      return this->normalizeDiag(n, -3)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      auto l3 = precision::pow(l1, 3);
      auto l4 = precision::pow(l1, 4);
      ACoeff_t val;

      val = 960.0*(l1 + n)*(l1 + n - 1.0)*(l4 - 4.0*l3*n + 2.0*l3 - 3.0*l2*n.pow(2) + 3.0*l2*n + 17.0*l2 + 2.0*l1*n.pow(3) - 3.0*l1*n.pow(2) - 7.0*l1*n + 4.0*l1 + n.pow(4) - 2.0*n.pow(3) - 7.0*n.pow(2) + 8.0*n + 12.0)/((l1 + 2.0*n)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0));

      return this->normalizeDiag(n, -2)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      auto l3 = precision::pow(l1, 3);
      auto l4 = precision::pow(l1, 4);
      ACoeff_t val;

      val = -384.0*l1*(l1 + n)*(l4 - 10.0*l3*n + 35.0*l2 + 20.0*l1*n.pow(3) - 70.0*l1*n + 10.0*n.pow(4) - 70.0*n.pow(2) + 84.0)/((l1 + 2.0*n)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0));

      return this->normalizeDiag(n, -1)*val;
   }

   I6Diags::ACoeff_t I6Diags::d0(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      auto l3 = precision::pow(l1, 3);
      auto l4 = precision::pow(l1, 4);
      auto l5 = precision::pow(l1, 5);
      auto l6 = precision::pow(l1, 6);
      ACoeff_t val;

      val = 64.0*(l6 - 30.0*l5*n - 15.0*l5 + 60.0*l4*n.pow(2) + 60.0*l4*n + 85.0*l4 + 160.0*l3*n.pow(3) + 240.0*l3*n.pow(2) - 370.0*l3*n - 225.0*l3 + 30.0*l2*n.pow(4) + 60.0*l2*n.pow(3) - 270.0*l2*n.pow(2) - 300.0*l2*n + 274.0*l2 - 60.0*l1*n.pow(5) - 150.0*l1*n.pow(4) + 200.0*l1*n.pow(3) + 450.0*l1*n.pow(2) - 80.0*l1*n - 120.0*l1 - 20.0*n.pow(6) - 60.0*n.pow(5) + 100.0*n.pow(4) + 300.0*n.pow(3) - 80.0*n.pow(2) - 240.0*n)/((l1 + 2.0*n)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0));

      return this->normalizeDiag(n, 0)*val;
   }

   I6Diags::ACoeff_t I6Diags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      auto l3 = precision::pow(l1, 3);
      auto l4 = precision::pow(l1, 4);
      ACoeff_t val;

      val = 384.0*l1*(n + 1.0)*(l4 - 10.0*l3*n - 10.0*l3 + 35.0*l2 + 20.0*l1*n.pow(3) + 60.0*l1*n.pow(2) - 10.0*l1*n - 50.0*l1 + 10.0*n.pow(4) + 40.0*n.pow(3) - 10.0*n.pow(2) - 100.0*n + 24.0)/((l1 + 2.0*n)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0));

      return this->normalizeDiag(n, 1)*val;
   }

   I6Diags::ACoeff_t I6Diags::d2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      auto l3 = precision::pow(l1, 3);
      auto l4 = precision::pow(l1, 4);
      ACoeff_t val;

      val = 960.0*(n + 1.0)*(n + 2.0)*(l4 - 4.0*l3*n - 6.0*l3 - 3.0*l2*n.pow(2) - 9.0*l2*n + 11.0*l2 + 2.0*l1*n.pow(3) + 9.0*l1*n.pow(2) + 5.0*l1*n - 6.0*l1 + n.pow(4) + 6.0*n.pow(3) + 5.0*n.pow(2) - 12.0*n)/((l1 + 2.0*n)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0)*(l1 + 2.0*n + 9.0));

      return this->normalizeDiag(n, 2)*val;
   }

   I6Diags::ACoeff_t I6Diags::d3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = 640.0*l1*(n + 1.0)*(n + 2.0)*(n + 3.0)*(2.0*l2 - 3.0*l1*n - 6.0*l1 - 3.0*n.pow(2) - 12.0*n + 4.0)/((l1 + 2.0*n)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0)*(l1 + 2.0*n + 9.0)*(l1 + 2.0*n + 10.0));

      return this->normalizeDiag(n, 3)*val;
   }

   I6Diags::ACoeff_t I6Diags::d4(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = 192.0*(n + 1.0)*(n + 2.0)*(n + 3.0)*(n + 4.0)*(5.0*l2 - 2.0*l1*n - 5.0*l1 - 2.0*n.pow(2) - 10.0*n)/((l1 + 2.0*n)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0)*(l1 + 2.0*n + 9.0)*(l1 + 2.0*n + 10.0)*(l1 + 2.0*n + 11.0));

      return this->normalizeDiag(n, 4)*val;
   }

   I6Diags::ACoeff_t I6Diags::d5(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 384.0*l1*(n + 1.0)*(n + 2.0)*(n + 3.0)*(n + 4.0)*(n + 5.0)/((l1 + 2.0*n)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0)*(l1 + 2.0*n + 9.0)*(l1 + 2.0*n + 10.0)*(l1 + 2.0*n + 11.0)*(l1 + 2.0*n + 12.0));

      return this->normalizeDiag(n, 5)*val;
   }

   I6Diags::ACoeff_t I6Diags::d6(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 64.0*(n + 1.0)*(n + 2.0)*(n + 3.0)*(n + 4.0)*(n + 5.0)*(n + 6.0)/((l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0)*(l1 + 2.0*n + 9.0)*(l1 + 2.0*n + 10.0)*(l1 + 2.0*n + 11.0)*(l1 + 2.0*n + 12.0)*(l1 + 2.0*n + 13.0));

      return this->normalizeDiag(n, 6)*val;
   }

}
}
}
}
