/** 
 * @file I4Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I4Diags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/Chebyshev/I4Diags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   I4Diags::I4Diags(const Scalar_t alpha, const int l)
      : QuICC::SparseSM::Worland::I4Diags(alpha, MHD_MP(-0.5), l)
   {
   }

   I4Diags::ACoeff_t I4Diags::d_4(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      if(l1 == 0)
      {
         val = 1.0/((l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 1.0));
         val.head(1) = 16.0/((l1 + 4.0)*(l1 + 5.0)*(l1 + 6.0)*(l1 + 7.0));
      } else
      {
         val = 16.0*(l1 + n - 4.0)*(l1 + n - 3.0)*(l1 + n - 2.0)*(l1 + n - 1.0)/((l1 + 2.0*n - 8.0)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0));
      }

      return this->normalizeDiag(n, -4)*val;
   }

   I4Diags::ACoeff_t I4Diags::d_3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = -64.0*l1*(l1 + n - 3.0)*(l1 + n - 2.0)*(l1 + n - 1.0)/((l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0));

      return this->normalizeDiag(n, -3)*val;
   }

   I4Diags::ACoeff_t I4Diags::d_2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      ACoeff_t val;

      val = 16.0*(l1 + n - 2.0)*(l1 + n - 1.0)*(6.0*l2 - 4.0*l1*n + 4.0*l1 - 4.0*n.pow(2) + 8.0*n + 5.0)/((l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0));

      return this->normalizeDiag(n, -2)*val;
   }

   I4Diags::ACoeff_t I4Diags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      ACoeff_t val;

      val = -16.0*l1*(l1 + n - 1.0)*(4.0*l2 - 12.0*l1*n + 6.0*l1 - 12.0*n.pow(2) + 12.0*n + 17.0)/((l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0));

      return this->normalizeDiag(n, -1)*val;
   }

   I4Diags::ACoeff_t I4Diags::d0(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l1*l2;
      auto l4 = l2*l2;
      ACoeff_t val;

      val = 2.0*(8.0*l4 - 96.0*l3*n - 48.0*l2*n.pow(2) + 100.0*l2 + 96.0*l1*n.pow(3) - 120.0*l1*n + 48.0*n.pow(4) - 120.0*n.pow(2) + 27.0)/((l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0));

      return this->normalizeDiag(n, 0)*val;
   }

   I4Diags::ACoeff_t I4Diags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      ACoeff_t val;

      val = 4.0*l1*(2.0*n + 1.0)*(2.0*l1 + 2.0*n + 1.0)*(4.0*l2 - 12.0*l1*n - 6.0*l1 - 12.0*n.pow(2) - 12.0*n + 17.0)/((l1 + n)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0));

      return this->normalizeDiag(n, 1)*val;
   }

   I4Diags::ACoeff_t I4Diags::d2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      ACoeff_t val;

      val = (2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(6.0*l2 - 4.0*l1*n - 4.0*l1 - 4.0*n.pow(2) - 8.0*n + 5.0)/((l1 + n)*(l1 + n + 1.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0));

      return this->normalizeDiag(n, 2)*val;
   }

   I4Diags::ACoeff_t I4Diags::d3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = l1*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(2.0*l1 + 2.0*n + 5.0)/((l1 + n)*(l1 + n + 1.0)*(l1 + n + 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0));

      return this->normalizeDiag(n, 3)*val;
   }

   I4Diags::ACoeff_t I4Diags::d4(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = (2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(2.0*l1 + 2.0*n + 5.0)*(2.0*l1 + 2.0*n + 7.0)/(16.0*(l1 + n)*(l1 + n + 1.0)*(l1 + n + 2.0)*(l1 + n + 3.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0));

      return this->normalizeDiag(n, 4)*val;
   }

}
}
}
}
