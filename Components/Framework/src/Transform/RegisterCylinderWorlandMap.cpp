/**
 * file RegisterCylinderWorlandMap.cpp
 * @brief Source of the registration of transform operators in a cylinder
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/RegisterCylinderWorlandMap.hpp"
#include "QuICC/Transform/DefaultCylinderWorlandMap.hpp"

namespace QuICC {

namespace Transform {

   RegisterCylinderWorlandMap::MapVector& RegisterCylinderWorlandMap::mapper()
   {
      static MapVector v;

      if(v.size() == 0)
      {
         auto sp = std::make_shared<DefaultCylinderWorlandMap>();
         v.push_back(sp);
      }

      return v;
   }

} // Transform
} // QuICC
