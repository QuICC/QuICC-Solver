/** 
 * @file StressFreeTorAnelastic.cpp
 * @brief Source of the implementation of boundary r D(* / r) - F   
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/StressFreeTorAnelastic.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

   StressFreeTorAnelastic::StressFreeTorAnelastic(const Scalar_t lower, const Scalar_t upper, const Position pos, const Internal::MHDFloat Fb)
      : ICondition(lower, upper, pos), mFb(Fb)
   {
   }

   StressFreeTorAnelastic::ACoeff_t StressFreeTorAnelastic::compute(const int maxN)
   {
      ACoeff_t val = ACoeff_t::Zero(maxN+1);

      const auto a_1 = 1.0/this->a();
      if(this->position() == Position::TOP)
      {
         const auto ab_1 = 1.0/(this->a()+this->b());
         for(int i = 0; i < val.size(); i++)
         {
            const auto n = static_cast<Scalar_t>(i);
            const auto n2 = n*n;
            val(i) = (a_1*n2 - ab_1 - this->mFb)*this->c();
         }
      }
      else
      {
         const auto ab_1 = 1.0/(-this->a()+this->b());
         for(int i = 0; i < val.size(); i++)
         {
            const auto n = static_cast<Scalar_t>(i);
            const auto n2 = n*n;
            val(i) = (a_1*n2 + ab_1 + this->mFb)*this->c();
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
