/** 
 * @file Stability.cpp
 * @brief Source of the implementation of the stability node of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Simulation/Stability.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Config {

namespace Simulation {

   const std::string Stability::PARENTTAG = "stability";

   Stability::Stability(const std::vector<std::string>& names)
      : IConfigurationNode(Stability::PARENTTAG)
   {
      this->init(names);
   }

   Stability::~Stability()
   {
   }

   void Stability::init(const std::vector<std::string>& names)
   {
      // Get iterator over vector
      for(auto it = names.cbegin(); it != names.cend(); it++)
      {
         this->sTags().addTag(*it, "invalid");
      }
   }

   void Stability::checkData()
   {
   }

}
}
}
}
