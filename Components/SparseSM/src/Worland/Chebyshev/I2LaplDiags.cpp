/**
 * @file I2LaplDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I2LaplDiags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/Chebyshev/I2LaplDiags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   I2LaplDiags::I2LaplDiags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::I2LaplDiags(alpha, MHD_MP(-0.5), l, q)
   {
      // q <= 1 is equivalent to no truncation (already zero rows)

      if(q > 1)
      {
         throw std::logic_error("I2Lapl: Truncation for q>1 is not implemented");
      }
   }

   I2LaplDiags::ACoeff_t I2LaplDiags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 8.0*(l1 + n - 1.0)*(2.0*l1 + 2.0*n - 1.0)/((l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0));

      return this->normalizeDiag(n,-1)*val;
   }

   I2LaplDiags::ACoeff_t I2LaplDiags::d0(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 8.0*(4.0*l1*n + 4.0*n.pow(2) - 1.0)/((l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0));

      return this->normalizeDiag(n,0)*val;
   }

   I2LaplDiags::ACoeff_t I2LaplDiags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 2.0*(2.0*n + 1.0).pow(2)*(2.0*l1 + 2.0*n + 1.0)/((l1 + n)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0));

      return this->normalizeDiag(n,1)*val;
   }

}
}
}
}
