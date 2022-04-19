/** 
 * @file Run.cpp
 * @brief Source of the implementation of the run node of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Framework/Run.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Config {

namespace Framework {

   const std::string Run::PARENTTAG = "run";

   Run::Run()
      : IConfigurationNode(Run::PARENTTAG)
   {
      this->init();
   }

   Run::~Run()
   {
   }

   void Run::init()
   {
      this->fTags().addTag("sim", 0.0);
      this->fTags().addTag("wall", 0.0);
   }

   void Run::checkData()
   {
   }

}
}
}
}
