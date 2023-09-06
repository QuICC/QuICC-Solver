/** 
 * @file Model.cpp
 * @brief Source of the implementation of the setup model node of the configuration
 */

// System includes
//

// Project includes
//
#include "QuICC/Io/Config/Setup/Model.hpp"

namespace QuICC {
   
namespace Io {
   
namespace Config {
   
namespace Setup {

   const std::string Model::PARENTTAG = "model";

   Model::Model()
      : IConfigurationNode(Model::PARENTTAG)
   {
      this->init();
   }

   void Model::init()
   {
      this->sTags().addTag("split_equation", "off");
   }

   void Model::checkData()
   {
   }

}
}
}
}
