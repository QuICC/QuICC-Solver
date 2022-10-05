/**
 * @file I2T.cpp
 * @brief Source of the forward transform operator Forward::I2T
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I2T.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I2T::sTag()
   {
      return "Fwd::I2T";
   }

   std::string I2T::sFormatted()
   {
      return "Forward::I2T";
   }

   I2T::I2T()
      : IRegisterId<I2T>(I2T::sTag(), I2T::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
