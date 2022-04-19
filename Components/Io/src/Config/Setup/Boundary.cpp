/** 
 * @file Boundary.cpp
 * @brief Source of the implementation of the setup boundary node of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Setup/Boundary.hpp"

// Project includes
//

namespace QuICC {
   
namespace Io {
   
namespace Config {
   
namespace Setup {

   const std::string Boundary::PARENTTAG = "boundary";

   Boundary::Boundary()
      : IConfigurationNode(Boundary::PARENTTAG)
   {
      this->init();
   }

   Boundary::~Boundary()
   {
   }

   void Boundary::init()
   {
      this->iTags().addTag("scheme", 0);
   }

   void Boundary::checkData()
   {
   }

}
}
}
}
