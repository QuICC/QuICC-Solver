/** 
 * @file Parallel.cpp
 * @brief Source of the implementation of the setup parallel node of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Setup/Parallel.hpp"

// Project includes
//

namespace QuICC {
   
namespace Io {
   
namespace Config {
   
namespace Setup {

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
   }

   void Parallel::checkData()
   {
   }

}
}
}
}
