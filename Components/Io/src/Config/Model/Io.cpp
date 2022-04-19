/** 
 * @file Io.cpp
 * @brief Source of the implementation of a simple IO node of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Model/Io.hpp"

// Project includes
//

namespace QuICC {
   
namespace Io {
   
namespace Config {
   
namespace Model {

   Io::Io(const std::string tag, const std::map<std::string,int>& options)
      : IConfigurationNode(tag)
   {
      this->init(options);
   }

   Io::~Io()
   {
   }

   void Io::init(const std::map<std::string,int>& options)
   {
      for(auto&& [o,v] : options)
      {
         this->iTags().addTag(o, v);
      }
   }

   void Io::checkData()
   {
   }

}
}
}
}
