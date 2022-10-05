/**
 * @file D1ZP0.cpp
 * @brief Source of the forward transform operator Forward::D1ZP0
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/D1ZP0.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string D1ZP0::sTag()
   {
      return "Fwd::D1ZP0";
   }

   std::string D1ZP0::sFormatted()
   {
      return "Forward::D1ZP0";
   }

   D1ZP0::D1ZP0()
      : IRegisterId<D1ZP0>(D1ZP0::sTag(), D1ZP0::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
