/** 
 * @file I4LaplDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4LaplDiags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/Chebyshev/I4LaplDiags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   I4LaplDiags::I4LaplDiags(const Scalar_t alpha, const int l)
      : QuICC::SparseSM::Worland::I4LaplDiags(alpha, MHD_MP(-0.5), l)
   {
   }

   I4LaplDiags::ACoeff_t I4LaplDiags::d_3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 32.0*(l1 + n - 3.0)*(l1 + n - 2.0)*(l1 + n - 1.0)*(2.0*l1 + 2.0*n - 5.0)/((l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0));

      return this->normalizeDiag(n, -3)*val;
   }

   I4LaplDiags::ACoeff_t I4LaplDiags::d_2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      ACoeff_t val;

      val = -32.0*(l1 + n - 2.0)*(l1 + n - 1.0)*(4.0*l2 - 6.0*l1 - 4.0*n.pow(2) + 8.0*n + 5.0)/((l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0));

      return this->normalizeDiag(n, -2)*val;
   }

   I4LaplDiags::ACoeff_t I4LaplDiags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l2*l1;
      ACoeff_t val;

      val = 8.0*(l1 + n - 1.0)*(8.0*l3 - 40.0*l2*n - 4.0*l2 - 56.0*l1*n.pow(2) + 80.0*l1*n + 62.0*l1 - 8.0*n.pow(3) + 36.0*n.pow(2) - 22.0*n - 21.0)/((l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0));

      return this->normalizeDiag(n, -1)*val;
   }

   I4LaplDiags::ACoeff_t I4LaplDiags::d0(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l2*l1;
      ACoeff_t val;

      val = 16.0*(16.0*l3*n + 4.0*l3 - 24.0*l2*n - 24.0*l2 - 32.0*l1*n.pow(3) - 24.0*l1*n.pow(2) + 40.0*l1*n + 14.0*l1 - 16.0*n.pow(4) + 40.0*n.pow(2) - 9.0)/((l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0));

      return this->normalizeDiag(n, 0)*val;
   }

   I4LaplDiags::ACoeff_t I4LaplDiags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      ACoeff_t val;

      val = 2.0*(2.0*n + 1.0)*(2.0*l1 + 2.0*n + 1.0)*(48.0*l2*n + 48.0*l2 + 32.0*l1*n.pow(2) + 8.0*l1*n - 84.0*l1 - 8.0*n.pow(3) - 36.0*n.pow(2) - 22.0*n + 21.0)/((l1 + n)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0));

      return this->normalizeDiag(n, 1)*val;
   }

   I4LaplDiags::ACoeff_t I4LaplDiags::d2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 2.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(8.0*l1*n + 14.0*l1 + 4.0*n.pow(2) + 8.0*n - 5.0)/((l1 + n)*(l1 + n + 1.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0));

      return this->normalizeDiag(n, 2)*val;
   }

   I4LaplDiags::ACoeff_t I4LaplDiags::d3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = (2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0).pow(2)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(2.0*l1 + 2.0*n + 5.0)/(2.0*(l1 + n)*(l1 + n + 1.0)*(l1 + n + 2.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0));

      return this->normalizeDiag(n, 3)*val;
   }

}
}
}
}
