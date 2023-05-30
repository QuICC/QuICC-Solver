/**
 * file RegisterSphereWorlandMap.cpp
 * @brief Source of the registration of transform operators in a sphere
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/RegisterSphereWorlandMap.hpp"
#include "QuICC/Transform/DefaultSphereWorlandMap.hpp"

namespace QuICC {

namespace Transform {

   RegisterSphereWorlandMap::MapVector& RegisterSphereWorlandMap::mapper()
   {
      static MapVector v;

      if(v.size() == 0)
      {
         auto sp = std::make_shared<DefaultSphereWorlandMap>();
         v.push_back(sp);
      }

      return v;
   }

} // Transform
} // QuICC
