/** 
 * @file Parallel.cpp
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
#include "QuICC/Io/Config/Framework/Parallel.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Config {

namespace Framework {

   const std::string Parallel::PARENTTAG = "parallel";

   Parallel::Parallel()
      : IConfigurationNode(Parallel::PARENTTAG)
   {
      this->init();
   }

   Parallel::~Parallel()
   {
   }

   void Parallel::init()
   {
      this->iTags().addTag("cpus", -1);
      this->sTags().addTag("algorithm", "serial");
      this->sTags().addTag("grouper", "transform");
   }

   void Parallel::checkData()
   {
   }

}
}
}
}
