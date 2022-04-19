/**
 * @file Stop.cpp
 * @brief Source of the Stop RuntimeStatus
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/RuntimeStatus/Stop.hpp"

// Project includes
//

namespace QuICC {

namespace RuntimeStatus {

   std::string Stop::sTag()
   {
      return "stop";
   }

   std::string Stop::sFormatted()
   {
      return "Stop";
   }

   Stop::Stop()
      : IRegisterId<Stop>(Stop::sTag(), Stop::sFormatted())
   {
   }

}
}
