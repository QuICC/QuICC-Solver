/** 
 * @file I4D1R1Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I4D1R1Diags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/Chebyshev/I4D1R1Diags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   I4D1R1Diags::I4D1R1Diags(const Scalar_t alpha, const int l)
      : QuICC::SparseSM::Worland::I4D1R1Diags(alpha, MHD_MP(-0.5), l)
   {
   }

   I4D1R1Diags::ACoeff_t I4D1R1Diags::d_3(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = 2.0/((2.0*n - 5.0)*(2.0*n - 3.0)*(2.0*n - 1.0));

      return this->normalizeDiag(n, -3)*val;
   }

   I4D1R1Diags::ACoeff_t I4D1R1Diags::d_2(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = 1.0/((n + 1.0)*(2.0*n - 3.0)*(2.0*n - 1.0));

      return this->normalizeDiag(n, -2)*val;
   }

   I4D1R1Diags::ACoeff_t I4D1R1Diags::d_1(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = -3.0/(2.0*(n - 2.0)*(n + 1.0)*(2.0*n - 1.0));

      return this->normalizeDiag(n, -1)*val;
   }

   I4D1R1Diags::ACoeff_t I4D1R1Diags::d0(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = -3.0/(4.0*(n - 1.0)*(n + 1.0)*(n + 2.0));

      return this->normalizeDiag(n, 0)*val;
   }

   I4D1R1Diags::ACoeff_t I4D1R1Diags::d1(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = 3.0*(2.0*n + 1.0)/(8.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0));

      return this->normalizeDiag(n, 1)*val;
   }

   I4D1R1Diags::ACoeff_t I4D1R1Diags::d2(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = 3.0*(2.0*n + 1.0)*(2.0*n + 3.0)/(16.0*n*(n + 1.0).pow(2)*(n + 2.0)*(n + 3.0));

      return this->normalizeDiag(n, 2)*val;
   }

   I4D1R1Diags::ACoeff_t I4D1R1Diags::d3(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = -(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)/(32.0*n*(n + 1.0).pow(2)*(n + 2.0).pow(2)*(n + 3.0));

      return this->normalizeDiag(n, 3)*val;
   }

   I4D1R1Diags::ACoeff_t I4D1R1Diags::d4(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = -(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)/(64.0*(n + 1.0).pow(2)*(n + 2.0).pow(2)*(n + 3.0).pow(2)*(n + 4.0));

      return this->normalizeDiag(n, 4)*val;
   }

}
}
}
}
