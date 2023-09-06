/** 
 * @file Parallel.cpp
 * @brief Source of the implementation of the parallel node of the configuration
 */

// System includes
//

// Project includes
//
#include "QuICC/Io/Config/Framework/Parallel.hpp"

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

   void Parallel::init()
   {
      this->iTags().addTag("cpus", -1);
      this->sTags().addTag("algorithm", "serial");
      this->sTags().addTag("grouper", "transform");
      this->sTags().addTag("decomposition", "auto");
   }

   void Parallel::checkData()
   {
   }

} // Framework
} // Config
} // Io
} // QuICC
