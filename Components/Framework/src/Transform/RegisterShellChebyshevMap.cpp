/**
 * file RegisterShellChebyshevMap.cpp
 * @brief Source of the registration of transform operators in a spherical shell
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/RegisterShellChebyshevMap.hpp"
#include "QuICC/Transform/DefaultShellChebyshevMap.hpp"

namespace QuICC {

namespace Transform {

   RegisterShellChebyshevMap::MapVector& RegisterShellChebyshevMap::mapper()
   {
      static MapVector v;

      if(v.size() == 0)
      {
         auto sp = std::make_shared<DefaultShellChebyshevMap>();
         v.push_back(sp);
      }

      return v;
   }

} // Transform
} // QuICC
