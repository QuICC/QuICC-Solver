/**
 * @file I2D1.cpp
 * @brief Source of the forward transform operator Forward::I2D1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I2D1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I2D1::sTag()
   {
      return "Fwd::I2D1";
   }

   std::string I2D1::sFormatted()
   {
      return "Forward::I2D1";
   }

   I2D1::I2D1()
      : IRegisterId<I2D1>(I2D1::sTag(), I2D1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
