/**
 * file RegisterAnnulusChebyshevMap.cpp
 * @brief Source of registration of transform operators in a annulus
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/RegisterAnnulusChebyshevMap.hpp"
#include "QuICC/Transform/DefaultAnnulusChebyshevMap.hpp"

namespace QuICC {

namespace Transform {

   RegisterAnnulusChebyshevMap::MapVector& RegisterAnnulusChebyshevMap::mapper()
   {
      static MapVector v;

      if(v.size() == 0)
      {
         auto sp = std::make_shared<DefaultAnnulusChebyshevMap>();
         v.push_back(sp);
      }

      return v;
   }

} // Transform
} // QuICC
