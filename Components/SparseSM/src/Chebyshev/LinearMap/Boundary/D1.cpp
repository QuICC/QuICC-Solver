/** 
 * @file D1.cpp
 * @brief Source of the implementation of boundary first derivative
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/D1.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

   D1::D1(const Scalar_t lower, const Scalar_t upper, const Position pos)
      : ICondition(lower, upper, pos)
   {
   }

   D1::ACoeff_t D1::compute(const int maxN)
   {
      ACoeff_t val = ACoeff_t::Zero(maxN+1);

      auto cnst = this->c()/this->a();
      for(int i = 1; i < val.size(); i++)
      {
         const auto n = static_cast<Scalar_t>(i);
         val(i) = cnst*n*n;
      }

      // at the bottom (x = -1)
      if(this->position() == Position::BOTTOM)
      {
         for(int i = 1; i < val.size(); i++)
         {
            if(i%2 == 0)
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
