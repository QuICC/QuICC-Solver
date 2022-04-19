/** 
 * @file Coordinator.cpp
 * @brief Source of the spatial scheme coordinator
 */

// System includes
//
#include <string>

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/Coordinator.hpp"

// Project includes
//

namespace QuICC {

namespace SpatialScheme {

   std::map<std::size_t,SharedICreator>& Coordinator::map()
   {
      static std::map<std::size_t, SharedICreator> creatorMap;
      return creatorMap;
   }

   Coordinator::Coordinator()
   {
   }

   Coordinator::~Coordinator()
   {
   }

}
}
