/**
 * file RegisterSphereFiniteDiffMap.cpp
 * @brief Source of the registration of transform operators in a sphere
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/RegisterSphereFiniteDiffMap.hpp"
#include "QuICC/Transform/DefaultSphereFiniteDiffMap.hpp"

namespace QuICC {

namespace Transform {

   RegisterSphereFiniteDiffMap::MapVector& RegisterSphereFiniteDiffMap::mapper()
   {
      static MapVector v;

      if(v.size() == 0)
      {
         auto sp = std::make_shared<DefaultSphereFiniteDiffMap>();
         v.push_back(sp);
      }

      return v;
   }

} // Transform
} // QuICC
