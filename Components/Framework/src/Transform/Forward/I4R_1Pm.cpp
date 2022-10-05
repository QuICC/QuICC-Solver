/**
 * @file I4R_1Pm.cpp
 * @brief Source of the forward transform operator Forward::I4R_1Pm
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4R_1Pm.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4R_1Pm::sTag()
   {
      return "Fwd::I4R_1Pm";
   }

   std::string I4R_1Pm::sFormatted()
   {
      return "Forward::I4R_1Pm";
   }

   I4R_1Pm::I4R_1Pm()
      : IRegisterId<I4R_1Pm>(I4R_1Pm::sTag(), I4R_1Pm::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
