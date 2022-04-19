/**
 * @file After.cpp
 * @brief Source of the After SolveTiming
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SolveTiming/After.hpp"

// Project includes
//

namespace QuICC {

namespace SolveTiming {

   std::string After::sTag()
   {
      return "after";
   }

   std::string After::sFormatted()
   {
      return "After";
   }

   After::After()
      : IRegisterId<After>(After::sTag(), After::sFormatted())
   {
   }

}
}
