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
#include "QuICC/SparseSM/Worland/SphEnergy/I6Diags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace SphEnergy {

   I6Diags::I6Diags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::I6Diags(alpha, MHD_MP(0.5), l, q)
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

      val = 4096.0*(2.0*l1 + 2.0*n - 9.0)*(2.0*l1 + 2.0*n - 7.0)*(2.0*l1 + 2.0*n - 5.0)*(2.0*l1 + 2.0*n - 3.0)*(2.0*l1 + 2.0*n - 1.0)*(2.0*l1 + 2.0*n + 1.0)/((2.0*l1 + 4.0*n - 21.0)*(2.0*l1 + 4.0*n - 19.0)*(2.0*l1 + 4.0*n - 17.0)*(2.0*l1 + 4.0*n - 15.0)*(2.0*l1 + 4.0*n - 13.0)*(2.0*l1 + 4.0*n - 11.0)*(2.0*l1 + 4.0*n - 9.0)*(2.0*l1 + 4.0*n - 7.0)*(2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0));

      return this->normalizeDiag(n, -6)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_5(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = -24576.0*(2.0*l1 + 1.0)*(2.0*l1 + 2.0*n - 7.0)*(2.0*l1 + 2.0*n - 5.0)*(2.0*l1 + 2.0*n - 3.0)*(2.0*l1 + 2.0*n - 1.0)*(2.0*l1 + 2.0*n + 1.0)/((2.0*l1 + 4.0*n - 19.0)*(2.0*l1 + 4.0*n - 17.0)*(2.0*l1 + 4.0*n - 15.0)*(2.0*l1 + 4.0*n - 13.0)*(2.0*l1 + 4.0*n - 11.0)*(2.0*l1 + 4.0*n - 9.0)*(2.0*l1 + 4.0*n - 7.0)*(2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0));

      return this->normalizeDiag(n, -5)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_4(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = 12288.0*(2.0*l1 + 2.0*n - 5.0)*(2.0*l1 + 2.0*n - 3.0)*(2.0*l1 + 2.0*n - 1.0)*(2.0*l1 + 2.0*n + 1.0)*(20.0*l2 - 8.0*l1*n + 32.0*l1 - 8.0*n.pow(2) + 20.0*n + 43.0)/((2.0*l1 + 4.0*n - 17.0)*(2.0*l1 + 4.0*n - 15.0)*(2.0*l1 + 4.0*n - 13.0)*(2.0*l1 + 4.0*n - 11.0)*(2.0*l1 + 4.0*n - 9.0)*(2.0*l1 + 4.0*n - 7.0)*(2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0));

      return this->normalizeDiag(n, -4)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = -81920.0*(2.0*l1 + 1.0)*(2.0*l1 + 2.0*n - 3.0)*(2.0*l1 + 2.0*n - 1.0)*(2.0*l1 + 2.0*n + 1.0)*(4.0*l2 - 6.0*l1*n + 10.0*l1 - 6.0*n.pow(2) + 9.0*n + 30.0)/((2.0*l1 + 4.0*n - 15.0)*(2.0*l1 + 4.0*n - 13.0)*(2.0*l1 + 4.0*n - 11.0)*(2.0*l1 + 4.0*n - 9.0)*(2.0*l1 + 4.0*n - 7.0)*(2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0));

      return this->normalizeDiag(n, -3)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      auto l3 = precision::pow(l1, 3);
      auto l4 = precision::pow(l1, 4);
      ACoeff_t val;

      val = 61440.0*(2.0*l1 + 2.0*n - 1.0)*(2.0*l1 + 2.0*n + 1.0)*(16.0*l4 - 64.0*l3*n + 64.0*l3 - 48.0*l2*n.pow(2) - 48.0*l2*n + 344.0*l2 + 32.0*l1*n.pow(3) - 96.0*l1*n.pow(2) - 112.0*l1*n + 368.0*l1 + 16.0*n.pow(4) - 16.0*n.pow(3) - 148.0*n.pow(2) + 76.0*n + 297.0)/((2.0*l1 + 4.0*n - 13.0)*(2.0*l1 + 4.0*n - 11.0)*(2.0*l1 + 4.0*n - 9.0)*(2.0*l1 + 4.0*n - 7.0)*(2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0));

      return this->normalizeDiag(n, -2)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      auto l3 = precision::pow(l1, 3);
      auto l4 = precision::pow(l1, 4);
      ACoeff_t val;

      val = -24576.0*(2.0*l1 + 1.0)*(2.0*l1 + 2.0*n + 1.0)*(16.0*l4 - 160.0*l3*n + 32.0*l3 - 240.0*l2*n + 584.0*l2 + 320.0*l1*n.pow(3) - 1240.0*l1*n + 568.0*l1 + 160.0*n.pow(4) + 160.0*n.pow(3) - 1120.0*n.pow(2) - 580.0*n + 1485.0)/((2.0*l1 + 4.0*n - 11.0)*(2.0*l1 + 4.0*n - 9.0)*(2.0*l1 + 4.0*n - 7.0)*(2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0)*(2.0*l1 + 4.0*n + 13.0));

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

      val = 4096.0*(64.0*l6 - 1920.0*l5*n - 768.0*l5 + 3840.0*l4*n.pow(2) - 960.0*l4*n + 3280.0*l4 + 10240.0*l3*n.pow(3) + 23040.0*l3*n.pow(2) - 20800.0*l3*n - 5760.0*l3 + 1920.0*l2*n.pow(4) + 19200.0*l2*n.pow(3) + 11520.0*l2*n.pow(2) - 51360.0*l2*n + 2956.0*l2 - 3840.0*l1*n.pow(5) - 7680.0*l1*n.pow(4) + 24320.0*l1*n.pow(3) + 24960.0*l1*n.pow(2) - 40760.0*l1*n + 1488.0*l1 - 1280.0*n.pow(6) - 5760.0*n.pow(5) + 2080.0*n.pow(4) + 27840.0*n.pow(3) + 7120.0*n.pow(2) - 25500.0*n - 945.0)/((2.0*l1 + 4.0*n - 9.0)*(2.0*l1 + 4.0*n - 7.0)*(2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0)*(2.0*l1 + 4.0*n + 13.0)*(2.0*l1 + 4.0*n + 15.0));

      return this->normalizeDiag(n, 0)*val;
   }

   I6Diags::ACoeff_t I6Diags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      auto l3 = precision::pow(l1, 3);
      auto l4 = precision::pow(l1, 4);
      ACoeff_t val;

      val = 49152.0*(2.0*l1 + 1.0)*(n + 1.0)*(16.0*l4 - 160.0*l3*n - 128.0*l3 - 240.0*l2*n + 344.0*l2 + 320.0*l1*n.pow(3) + 960.0*l1*n.pow(2) - 280.0*l1*n - 352.0*l1 + 160.0*n.pow(4) + 800.0*n.pow(3) + 320.0*n.pow(2) - 1700.0*n + 105.0)/((2.0*l1 + 4.0*n - 7.0)*(2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0)*(2.0*l1 + 4.0*n + 13.0)*(2.0*l1 + 4.0*n + 15.0)*(2.0*l1 + 4.0*n + 17.0));

      return this->normalizeDiag(n, 1)*val;
   }

   I6Diags::ACoeff_t I6Diags::d2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      auto l3 = precision::pow(l1, 3);
      auto l4 = precision::pow(l1, 4);
      ACoeff_t val;

      val = 245760.0*(n + 1.0)*(n + 2.0)*(16.0*l4 - 64.0*l3*n - 64.0*l3 - 48.0*l2*n.pow(2) - 240.0*l2*n + 56.0*l2 + 32.0*l1*n.pow(3) + 96.0*l1*n.pow(2) - 112.0*l1*n + 16.0*l1 + 16.0*n.pow(4) + 112.0*n.pow(3) + 140.0*n.pow(2) - 196.0*n - 15.0)/((2.0*l1 + 4.0*n - 5.0)*(2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0)*(2.0*l1 + 4.0*n + 13.0)*(2.0*l1 + 4.0*n + 15.0)*(2.0*l1 + 4.0*n + 17.0)*(2.0*l1 + 4.0*n + 19.0));

      return this->normalizeDiag(n, 2)*val;
   }

   I6Diags::ACoeff_t I6Diags::d3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = 655360.0*(2.0*l1 + 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0)*(4.0*l2 - 6.0*l1*n - 8.0*l1 - 6.0*n.pow(2) - 27.0*n + 3.0)/((2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0)*(2.0*l1 + 4.0*n + 13.0)*(2.0*l1 + 4.0*n + 15.0)*(2.0*l1 + 4.0*n + 17.0)*(2.0*l1 + 4.0*n + 19.0)*(2.0*l1 + 4.0*n + 21.0));

      return this->normalizeDiag(n, 3)*val;
   }

   I6Diags::ACoeff_t I6Diags::d4(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = 196608.0*(n + 1.0)*(n + 2.0)*(n + 3.0)*(n + 4.0)*(20.0*l2 - 8.0*l1*n - 8.0*n.pow(2) - 44.0*n - 5.0)/((2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0)*(2.0*l1 + 4.0*n + 13.0)*(2.0*l1 + 4.0*n + 15.0)*(2.0*l1 + 4.0*n + 17.0)*(2.0*l1 + 4.0*n + 19.0)*(2.0*l1 + 4.0*n + 21.0)*(2.0*l1 + 4.0*n + 23.0));

      return this->normalizeDiag(n, 4)*val;
   }

   I6Diags::ACoeff_t I6Diags::d5(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 786432.0*(2.0*l1 + 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0)*(n + 4.0)*(n + 5.0)/((2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0)*(2.0*l1 + 4.0*n + 13.0)*(2.0*l1 + 4.0*n + 15.0)*(2.0*l1 + 4.0*n + 17.0)*(2.0*l1 + 4.0*n + 19.0)*(2.0*l1 + 4.0*n + 21.0)*(2.0*l1 + 4.0*n + 23.0)*(2.0*l1 + 4.0*n + 25.0));

      return this->normalizeDiag(n, 5)*val;
   }

   I6Diags::ACoeff_t I6Diags::d6(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 262144.0*(n + 1.0)*(n + 2.0)*(n + 3.0)*(n + 4.0)*(n + 5.0)*(n + 6.0)/((2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0)*(2.0*l1 + 4.0*n + 9.0)*(2.0*l1 + 4.0*n + 11.0)*(2.0*l1 + 4.0*n + 13.0)*(2.0*l1 + 4.0*n + 15.0)*(2.0*l1 + 4.0*n + 17.0)*(2.0*l1 + 4.0*n + 19.0)*(2.0*l1 + 4.0*n + 21.0)*(2.0*l1 + 4.0*n + 23.0)*(2.0*l1 + 4.0*n + 25.0)*(2.0*l1 + 4.0*n + 27.0));

      return this->normalizeDiag(n, 6)*val;
   }

}
}
}
}
