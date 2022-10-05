/**
 * @file I2P.cpp
 * @brief Source of the forward transform operator Forward::I2P
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I2P.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I2P::sTag()
   {
      return "Fwd::I2P";
   }

   std::string I2P::sFormatted()
   {
      return "Forward::I2P";
   }

   I2P::I2P()
      : IRegisterId<I2P>(I2P::sTag(), I2P::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
