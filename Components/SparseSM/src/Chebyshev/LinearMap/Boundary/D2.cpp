/** 
 * @file D2.cpp
 * @brief Source of the implementation of boundary second derivative
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/D2.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

   D2::D2(const Scalar_t lower, const Scalar_t upper, const Position pos)
      : ICondition(lower, upper, pos)
   {
   }

   D2::ACoeff_t D2::compute(const int maxN)
   {
      ACoeff_t val = ACoeff_t::Zero(maxN+1);

      auto cnst = this->c()/(3.0*this->a()*this->a());
      for(int i = 2; i < val.size(); i++)
      {
         const auto n = static_cast<Scalar_t>(i);
         const auto n2 = n*n;
         const auto n4 = n2*n2;
         val(i) = cnst*(n4 - n2);
      }

      // at the bottom (x = -1)
      if(this->position() == Position::BOTTOM)
      {
         for(int i = 2; i < val.size(); i++)
         {
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
