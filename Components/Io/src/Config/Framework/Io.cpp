/** 
 * @file Io.cpp
 * @brief Source of the implementation of the IO node of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Framework/Io.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Config {

namespace Framework {

   const std::string Io::PARENTTAG = "io";

   Io::Io()
      : IConfigurationNode(Io::PARENTTAG)
   {
      this->init();
   }

   Io::~Io()
   {
   }

   void Io::init()
   {
      this->fTags().addTag("ascii", 0.0);

      this->fTags().addTag("hdf5", 0.0);
   }

   void Io::checkData()
   {
   }

}
}
}
}
