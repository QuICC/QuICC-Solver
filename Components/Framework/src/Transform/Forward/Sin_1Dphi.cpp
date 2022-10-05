/**
 * @file Sin_1Dphi.cpp
 * @brief Source of the forward transform operator Forward::Sin_1Dphi
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Sin_1Dphi.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Sin_1Dphi::sTag()
   {
      return "Fwd::Sin_1Dphi";
   }

   std::string Sin_1Dphi::sFormatted()
   {
      return "Forward::Sin_1Dphi";
   }

   Sin_1Dphi::Sin_1Dphi()
      : IRegisterId<Sin_1Dphi>(Sin_1Dphi::sTag(), Sin_1Dphi::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
