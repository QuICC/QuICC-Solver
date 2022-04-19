/**
 * @file GoOn.cpp
 * @brief Source of the GoOn RuntimeStatus
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/RuntimeStatus/GoOn.hpp"

// Project includes
//

namespace QuICC {

namespace RuntimeStatus {

   std::string GoOn::sTag()
   {
      return "go_on";
   }

   std::string GoOn::sFormatted()
   {
      return "GoOn";
   }

   GoOn::GoOn()
      : IRegisterId<GoOn>(GoOn::sTag(), GoOn::sFormatted())
   {
   }

}
}
