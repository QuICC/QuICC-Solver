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
#include "QuICC/SparseSM/Worland/Chebyshev/I6Diags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   I6Diags::I6Diags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::I6Diags(alpha, MHD_MP(-0.5), l, q)
   {
   }

   I6Diags::ACoeff_t I6Diags::d_6(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      if(l1 == 0)
      {
         val = 1.0/((l1 + 2.0*n - 11.0)*(l1 + 2.0*n - 9.0)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 1.0));
         val.head(1) = 2.0/10395.0;
      } else
      {
         val = 64.0*(l1 + n - 6.0)*(l1 + n - 5.0)*(l1 + n - 4.0)*(l1 + n - 3.0)*(l1 + n - 2.0)*(l1 + n - 1.0)/((l1 + 2.0*n - 12.0)*(l1 + 2.0*n - 11.0)*(l1 + 2.0*n - 10.0)*(l1 + 2.0*n - 9.0)*(l1 + 2.0*n - 8.0)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0));
      }

      // Truncate operator
      this->zeroLast(val, this->mQ-3);

      return this->normalizeDiag(n, -6)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_5(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = -384.0*l1*(l1 + n - 5.0)*(l1 + n - 4.0)*(l1 + n - 3.0)*(l1 + n - 2.0)*(l1 + n - 1.0)/((l1 + 2.0*n - 11.0)*(l1 + 2.0*n - 10.0)*(l1 + 2.0*n - 9.0)*(l1 + 2.0*n - 8.0)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0));

      // Truncate operator
      this->zeroLast(val, this->mQ-2);

      return this->normalizeDiag(n, -5)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_4(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      ACoeff_t ln = l1 + n;
      ACoeff_t l2n = l1 + 2.0*n;
      ACoeff_t val;

      val = 96.0*(ln - 4.0)*(ln - 3.0)*(ln - 2.0)*(ln - 1.0)*(10.0*l2 - 4.0*l1*n + 8.0*l1 - 4.0*n.pow(2) + 16.0*n + 9.0)/((l2n - 10.0)*(l2n - 9.0)*(l2n - 8.0)*(l2n - 7.0)*(l2n - 6.0)*(l2n - 5.0)*(l2n - 4.0)*(l2n - 3.0)*(l2n - 2.0)*(l2n - 1.0)*(l2n + 1.0)*(l2n + 2.0));

      // Truncate operator
      this->zeroLast(val, this->mQ-1);

      return this->normalizeDiag(n, -4)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      ACoeff_t val;

      val = -160.0*l1*(l1 + n - 3.0)*(l1 + n - 2.0)*(l1 + n - 1.0)*(8.0*l2 - 12.0*l1*n + 18.0*l1 - 12.0*n.pow(2) + 36.0*n + 37.0)/((l1 + 2.0*n - 9.0)*(l1 + 2.0*n - 8.0)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0));

      // Truncate operator
      this->zeroLast(val, this->mQ);

      return this->normalizeDiag(n, -3)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l2*l1;
      auto l4 = l2*l2;
      ACoeff_t val;

      val = 60.0*(l1 + n - 2.0)*(l1 + n - 1.0)*(16.0*l4 - 64.0*l3*n + 64.0*l3 - 48.0*l2*n.pow(2) + 96.0*l2*n + 236.0*l2 + 32.0*l1*n.pow(3) - 96.0*l1*n.pow(2) - 40.0*l1*n + 104.0*l1 + 16.0*n.pow(4) - 64.0*n.pow(3) - 40.0*n.pow(2) + 208.0*n + 105.0)/((l1 + 2.0*n - 8.0)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0));

      // Truncate operator
      this->zeroLast(val, this->mQ+1);

      return this->normalizeDiag(n, -2)*val;
   }

   I6Diags::ACoeff_t I6Diags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l2*l1;
      auto l4 = l2*l2;
      ACoeff_t val;

      val = -48.0*l1*(l1 + n - 1.0)*(8.0*l4 - 80.0*l3*n + 40.0*l3 + 280.0*l2 + 160.0*l1*n.pow(3) - 240.0*l1*n.pow(2) - 440.0*l1*n + 260.0*l1 + 80.0*n.pow(4) - 160.0*n.pow(3) - 440.0*n.pow(2) + 520.0*n + 537.0)/((l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0));

      // Truncate operator
      this->zeroLast(val, this->mQ+2);

      return this->normalizeDiag(n, -1)*val;
   }

   I6Diags::ACoeff_t I6Diags::d0(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l2*l1;
      auto l4 = l2*l2;
      auto l5 = l3*l2;
      auto l6 = l3*l3;
      ACoeff_t val;

      val = 4.0*(16.0*l6 - 480.0*l5*n + 960.0*l4*n.pow(2) + 1120.0*l4 + 2560.0*l3*n.pow(3) - 7840.0*l3*n + 480.0*l2*n.pow(4) - 5040.0*l2*n.pow(2) + 5614.0*l2 - 960.0*l1*n.pow(5) + 5600.0*l1*n.pow(3) - 5180.0*l1*n - 320.0*n.pow(6) + 2800.0*n.pow(4) - 5180.0*n.pow(2) + 1125.0)/((l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0));

      // Truncate operator
      this->zeroLast(val, this->mQ+3);

      return this->normalizeDiag(n, 0)*val;
   }

   I6Diags::ACoeff_t I6Diags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l2*l1;
      auto l4 = l2*l2;
      ACoeff_t val;

      val = 12.0*l1*(2.0*n + 1.0)*(2.0*l1 + 2.0*n + 1.0)*(8.0*l4 - 80.0*l3*n - 40.0*l3 + 280.0*l2 + 160.0*l1*n.pow(3) + 240.0*l1*n.pow(2) - 440.0*l1*n - 260.0*l1 + 80.0*n.pow(4) + 160.0*n.pow(3) - 440.0*n.pow(2) - 520.0*n + 537.0)/((l1 + n)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0));

      // Truncate operator
      this->zeroLast(val, this->mQ+4);

      return this->normalizeDiag(n, 1)*val;
   }

   I6Diags::ACoeff_t I6Diags::d2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l2*l1;
      auto l4 = l2*l2;
      ACoeff_t val;

      val = 15.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(16.0*l4 - 64.0*l3*n - 64.0*l3 - 48.0*l2*n.pow(2) - 96.0*l2*n + 236.0*l2 + 32.0*l1*n.pow(3) + 96.0*l1*n.pow(2) - 40.0*l1*n - 104.0*l1 + 16.0*n.pow(4) + 64.0*n.pow(3) - 40.0*n.pow(2) - 208.0*n + 105.0)/(4.0*(l1 + n)*(l1 + n + 1.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0));

      // Truncate operator
      this->zeroLast(val, this->mQ+5);

      return this->normalizeDiag(n, 2)*val;
   }

   I6Diags::ACoeff_t I6Diags::d3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      ACoeff_t val;

      val = 5.0*l1*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(2.0*l1 + 2.0*n + 5.0)*(8.0*l2 - 12.0*l1*n - 18.0*l1 - 12.0*n.pow(2) - 36.0*n + 37.0)/(2.0*(l1 + n)*(l1 + n + 1.0)*(l1 + n + 2.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0)*(l1 + 2.0*n + 9.0));

      // Truncate operator
      this->zeroLast(val, this->mQ+6);

      return this->normalizeDiag(n, 3)*val;
   }

   I6Diags::ACoeff_t I6Diags::d4(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      ACoeff_t val;

      val = 3.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(2.0*l1 + 2.0*n + 5.0)*(2.0*l1 + 2.0*n + 7.0)*(10.0*l2 - 4.0*l1*n - 8.0*l1 - 4.0*n.pow(2) - 16.0*n + 9.0)/(8.0*(l1 + n)*(l1 + n + 1.0)*(l1 + n + 2.0)*(l1 + n + 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0)*(l1 + 2.0*n + 9.0)*(l1 + 2.0*n + 10.0));

      // Truncate operator
      this->zeroLast(val, this->mQ+7);

      return this->normalizeDiag(n, 4)*val;
   }

   I6Diags::ACoeff_t I6Diags::d5(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 3.0*l1*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*n + 9.0)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(2.0*l1 + 2.0*n + 5.0)*(2.0*l1 + 2.0*n + 7.0)*(2.0*l1 + 2.0*n + 9.0)/(8.0*(l1 + n)*(l1 + n + 1.0)*(l1 + n + 2.0)*(l1 + n + 3.0)*(l1 + n + 4.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0)*(l1 + 2.0*n + 9.0)*(l1 + 2.0*n + 10.0)*(l1 + 2.0*n + 11.0));

      // Truncate operator
      this->zeroLast(val, this->mQ+8);

      return this->normalizeDiag(n, 5)*val;
   }

   I6Diags::ACoeff_t I6Diags::d6(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = (2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*n + 9.0)*(2.0*n + 11.0)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(2.0*l1 + 2.0*n + 5.0)*(2.0*l1 + 2.0*n + 7.0)*(2.0*l1 + 2.0*n + 9.0)*(2.0*l1 + 2.0*n + 11.0)/(64.0*(l1 + n)*(l1 + n + 1.0)*(l1 + n + 2.0)*(l1 + n + 3.0)*(l1 + n + 4.0)*(l1 + n + 5.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0)*(l1 + 2.0*n + 9.0)*(l1 + 2.0*n + 10.0)*(l1 + 2.0*n + 11.0)*(l1 + 2.0*n + 12.0));

      // Truncate operator
      this->zeroLast(val, this->mQ+9);

      return this->normalizeDiag(n, 6)*val;
   }

}
}
}
}
