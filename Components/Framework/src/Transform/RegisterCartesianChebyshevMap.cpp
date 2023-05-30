/**
 * file RegisterCartesianChebyshevMap.cpp
 * @brief Source of the registration of transform operators in a cartesian box
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/RegisterCartesianChebyshevMap.hpp"
#include "QuICC/Transform/DefaultCartesianChebyshevMap.hpp"

namespace QuICC {

namespace Transform {

   RegisterCartesianChebyshevMap::MapVector& RegisterCartesianChebyshevMap::mapper()
   {
      static MapVector v;

      if(v.size() == 0)
      {
         auto sp = std::make_shared<DefaultCartesianChebyshevMap>();
         v.push_back(sp);
      }

      return v;
   }

} // Transform
} // QuICC
