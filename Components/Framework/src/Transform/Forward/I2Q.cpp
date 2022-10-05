/**
 * @file I2Q.cpp
 * @brief Source of the forward transform operator Forward::I2Q
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I2Q.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I2Q::sTag()
   {
      return "Fwd::I2Q";
   }

   std::string I2Q::sFormatted()
   {
      return "Forward::I2Q";
   }

   I2Q::I2Q()
      : IRegisterId<I2Q>(I2Q::sTag(), I2Q::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
