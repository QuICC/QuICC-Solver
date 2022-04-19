/**
 * @file Time.cpp
 * @brief Source of the Time ModelOperator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/ModelOperator/Time.hpp"

// Project includes
//

namespace QuICC {

namespace ModelOperator {

   std::string Time::sTag()
   {
      return "time";
   }

   std::string Time::sFormatted()
   {
      return "Time";
   }

   Time::Time()
      : IRegisterId<Time>(Time::sTag(), Time::sFormatted())
   {
   }

}
}
