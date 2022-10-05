/**
 * @file I2S.cpp
 * @brief Source of the forward transform operator Forward::I2S
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I2S.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I2S::sTag()
   {
      return "Fwd::I2S";
   }

   std::string I2S::sFormatted()
   {
      return "Forward::I2S";
   }

   I2S::I2S()
      : IRegisterId<I2S>(I2S::sTag(), I2S::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
