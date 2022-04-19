/** 
 * @file Timestepper.cpp
 * @brief Source of the implementation of the setup timestepper node of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Setup/Timestepper.hpp"

// Project includes
//

namespace QuICC {
   
namespace Io {
   
namespace Config {
   
namespace Setup {

   const std::string Timestepper::PARENTTAG = "timestepper";

   Timestepper::Timestepper()
      : IConfigurationNode(Timestepper::PARENTTAG)
   {
      this->init();
   }

   Timestepper::~Timestepper()
   {
   }

   void Timestepper::init()
   {
   }

   void Timestepper::checkData()
   {
   }

}
}
}
}
