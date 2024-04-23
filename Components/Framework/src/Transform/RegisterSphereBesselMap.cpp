/**
 * file RegisterSphereBesselMap.cpp
 * @brief Source of the registration of transform operators for spherical Bessel in a sphere
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/RegisterSphereBesselMap.hpp"
#include "QuICC/Transform/DefaultSphereBesselMap.hpp"

namespace QuICC {

namespace Transform {

   RegisterSphereBesselMap::MapVector& RegisterSphereBesselMap::mapper()
   {
      static MapVector v;

      if(v.size() == 0)
      {
         auto sp = std::make_shared<DefaultSphereBesselMap>();
         v.push_back(sp);
      }

      return v;
   }

} // Transform
} // QuICC
