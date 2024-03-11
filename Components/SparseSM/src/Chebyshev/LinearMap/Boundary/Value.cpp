/** 
 * @file Value.cpp
 * @brief Source of the implementation of boundary value
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/Value.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

   Value::Value(const Scalar_t lower, const Scalar_t upper, const Position pos)
      : ICondition(lower, upper, pos)
   {
   }

   Value::ACoeff_t Value::compute(const int maxN)
   {
      ACoeff_t val = ACoeff_t::Constant(maxN+1, this->c());

      // at the bottom (x = -1)
      if(this->position() == Position::BOTTOM)
      {
         for(int i = 0; i < val.size(); i++)
         {
            if(i%2 == 1)
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
