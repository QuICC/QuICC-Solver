/**
 * @file Q.cpp
 * @brief Source of the forward transform operator Forward::Q
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Q.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Q::sTag()
   {
      return "Fwd::Q";
   }

   std::string Q::sFormatted()
   {
      return "Forward::Q";
   }

   Q::Q()
      : IRegisterId<Q>(Q::sTag(), Q::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
