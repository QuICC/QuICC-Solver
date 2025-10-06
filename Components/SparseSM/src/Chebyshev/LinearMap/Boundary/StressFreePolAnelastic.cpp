/** 
 * @file StressFreePolAnelastic.cpp
 * @brief Source of the implementation of boundary second derivative
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/StressFreePolAnelastic.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

   StressFreePolAnelastic::StressFreePolAnelastic(const Scalar_t lower, const Scalar_t upper, const Position pos, const Internal::MHDFloat Fb)
         : ICondition(lower, upper, pos), mFb(Fb)
   {
   }

   StressFreePolAnelastic::ACoeff_t StressFreePolAnelastic::compute(const int maxN)
   {
      ACoeff_t val = ACoeff_t::Zero(maxN+1);

      auto cnst = this->c()/(3.0*this->a()*this->a());
      const auto a_1 = 1.0/this->a();

      // at the top (x = 1)
      if(this->position() == Position::TOP)
      {
         for(int i = 1; i < val.size(); i++)
         {
            const auto n = static_cast<Scalar_t>(i);
            const auto n2 = n*n;
            const auto n4 = n2*n2;
            //if(n<2)
            //{
            //   val(i) = -a_1*n2*this->mFb;
            //}
            //else
            //{
               val(i) = cnst*(n4 - n2) -a_1*n2*this->mFb*this->c();
            //}
         }
      }
      // at the bottom (x = -1)
      else
      {
         for(int i = 1; i < val.size(); i++)
         {
            const auto n = static_cast<Scalar_t>(i);
            const auto n2 = n*n;
            const auto n4 = n2*n2;
            if(n<2)
            {
               val(i) = a_1*n2*this->mFb*this->c();
            }
            else
            {
               val(i) = cnst*(n4 - n2) +a_1*n2*this->mFb*this->c();
            }
            if(i%2 == 1)
            {
               val(i) *= -1.0;
            }
         }
      }

      return val;
   }
 
} // Boundary
} // LinearMap
} // Chebyshev
} // Polynomial
} // QuICC
