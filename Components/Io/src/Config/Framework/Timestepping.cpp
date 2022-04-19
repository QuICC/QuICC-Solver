/** 
 * @file Timestepping.cpp
 * @brief Source of the implementation of the parallel node of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Framework/Timestepping.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Config {

namespace Framework {

   const std::string Timestepping::PARENTTAG = "timestepping";

   Timestepping::Timestepping()
      : IConfigurationNode(Timestepping::PARENTTAG)
   {
      this->init();
   }

   Timestepping::~Timestepping()
   {
   }

   void Timestepping::init()
   {
      this->fTags().addTag("time", -1.0);
      this->fTags().addTag("timestep", -1.0);
      this->fTags().addTag("error", -1.0);
   }

   void Timestepping::checkData()
   {
   }

}
}
}
}
