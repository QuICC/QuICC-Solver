/**
 * @file I2QmDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I2QmDiags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I2QmDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   I2QmDiags::I2QmDiags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::I2QmDiags(alpha, MHD_MP(-0.5), l, q)
   {
      // q <= 1 is equivalent to no truncation (already zero rows)

      if(q > 1)
      {
         throw std::logic_error("I2Qm: Truncation for q>1 is not implemented");
      }
   }

   I2QmDiags::ACoeff_t I2QmDiags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = -8.0*(l1 + n - 2.0)*(l1 + n - 1.0)/((l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0));

      return this->normalizeDiag(n,-1,-1)*val;
   }

   I2QmDiags::ACoeff_t I2QmDiags::d0(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 4.0*(l1 + n - 1.0)*(2.0*l1 - 2.0*n - 1.0)/((l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0));

      return this->normalizeDiag(n,0,-1)*val;
   }

   I2QmDiags::ACoeff_t I2QmDiags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 2.0*(2.0*n + 1.0)*(4.0*l1 + 2.0*n - 1.0)/((l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0));

      return this->normalizeDiag(n,1,-1)*val;
   }

   I2QmDiags::ACoeff_t I2QmDiags::d2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = (2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l1 + 2.0*n + 1.0)/((l1 + n)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0));

      return this->normalizeDiag(n,2,-1)*val;
   }

} // Chebyshev
} // Worland
} // SparseSM
} // QuICC
