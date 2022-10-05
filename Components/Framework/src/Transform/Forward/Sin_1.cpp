/**
 * @file Sin_1.cpp
 * @brief Source of the forward transform operator Forward::Sin_1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Sin_1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Sin_1::sTag()
   {
      return "Fwd::Sin_1";
   }

   std::string Sin_1::sFormatted()
   {
      return "Forward::Sin_1";
   }

   Sin_1::Sin_1()
      : IRegisterId<Sin_1>(Sin_1::sTag(), Sin_1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
