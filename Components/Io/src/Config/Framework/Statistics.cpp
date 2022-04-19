/** 
 * @file Statistics.cpp
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
#include "QuICC/Io/Config/Framework/Statistics.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Config {

namespace Framework {

   const std::string Statistics::PARENTTAG = "statistics";

   Statistics::Statistics()
      : IConfigurationNode(Statistics::PARENTTAG)
   {
      this->init();
   }

   Statistics::~Statistics()
   {
   }

   void Statistics::init()
   {
      this->fTags().addTag("output_rate", 0.0);

      this->fTags().addTag("time_avg_rate", 0.0);
   }

   void Statistics::checkData()
   {
   }

}
}
}
}
