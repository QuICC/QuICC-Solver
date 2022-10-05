/**
 * @file I6R_1Pm.cpp
 * @brief Source of the forward transform operator Forward::I6R_1Pm
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I6R_1Pm.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I6R_1Pm::sTag()
   {
      return "Fwd::I6R_1Pm";
   }

   std::string I6R_1Pm::sFormatted()
   {
      return "Forward::I6R_1Pm";
   }

   I6R_1Pm::I6R_1Pm()
      : IRegisterId<I6R_1Pm>(I6R_1Pm::sTag(), I6R_1Pm::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
