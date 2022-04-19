/** 
 * @file Boundary.cpp
 * @brief Source of the implementation of the boundary node of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Simulation/Boundary.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Config {

namespace Simulation {

   const std::string Boundary::PARENTTAG = "boundary";

   Boundary::Boundary(const std::vector<std::string>& names)
      : IConfigurationNode(Boundary::PARENTTAG)
   {
      this->init(names);
   }

   Boundary::~Boundary()
   {
   }

   void Boundary::init(const std::vector<std::string>& names)
   {
      // Get iterator over vector
      for(auto it = names.cbegin(); it != names.cend(); it++)
      {
         this->iTags().addTag(*it, -1);
      }
   }

   void Boundary::checkData()
   {
   }

}
}
}
}
