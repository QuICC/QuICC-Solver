/**
 * @file P0.cpp
 * @brief Source of the forward transform operator Forward::P0
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/P0.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string P0::sTag()
   {
      return "Fwd::P0";
   }

   std::string P0::sFormatted()
   {
      return "Forward::P0";
   }

   P0::P0()
      : IRegisterId<P0>(P0::sTag(), P0::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
