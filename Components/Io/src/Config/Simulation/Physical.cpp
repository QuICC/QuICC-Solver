/**
 * @file Physical.cpp
 * @brief Source of the implementation of the physical node of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Simulation/Physical.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Config {

namespace Simulation {

   const std::string Physical::PARENTTAG = "physical";

   Physical::Physical(const std::vector<std::string>& names)
      : IConfigurationNode(Physical::PARENTTAG)
   {
      this->init(names);
   }

   Physical::~Physical()
   {
   }

   void Physical::init(const std::vector<std::string>& names)
   {
      // Get iterator over vector
      for(auto it = names.cbegin(); it != names.cend(); it++)
      {
         this->fTags().addTag(*it, -1.0);
      }
   }

   void Physical::checkData()
   {
   }

}
}
}
}
