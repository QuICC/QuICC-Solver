/** 
 * @file InsulatingShell.cpp
 * @brief Source of the implementation of insulating boundary in spherical shell
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/InsulatingShell.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

   InsulatingShell::InsulatingShell(const Scalar_t lower, const Scalar_t upper, const Position pos, const int l)
      : ICondition(lower, upper, pos), mL(l)
   {
   }

   InsulatingShell::ACoeff_t InsulatingShell::compute(const int maxN)
   {
      ACoeff_t val = ACoeff_t::Zero(maxN+1);

      const auto a_1 = 2.0/this->a();
      if(this->position() == Position::TOP)
      {
         const auto ab_1 = (this->mL + 1.0)/(this->a()+this->b());
         for(int i = 0; i < val.size(); i++)
         {
            const auto n = static_cast<Scalar_t>(i);
            const auto n2 = n*n;
            val(i) = a_1*n2 + (ab_1)*this->c();
         }
      }
      else
      {
         const auto ab_1 = this->mL/(-this->a()+this->b());
         for(int i = 0; i < val.size(); i++)
         {
            const auto n = static_cast<Scalar_t>(i);
            const auto n2 = n*n;
            val(i) = a_1*n2 + (ab_1)*this->c();
            if(i%2 == 0)
            {
               val(i) *= -1.0;
            }
         }
      }

      // Normalization
      val(0) /= this->c();

      return val;
   }
 
} // Boundary
} // LinearMap
} // Chebyshev
} // Polynomial
} // QuICC
